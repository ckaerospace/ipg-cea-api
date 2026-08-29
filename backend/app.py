from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, field_validator
from starlette.responses import Response

from physics.cea_rocket import (
    CEA_VERSION,
    characteristics_sweep,
    co2_enthalpy_20C,
    hinj_for_target_mdot,
    hp_equilibrium_o2,
)
from physics.constants import (
    ADVANCED_SPECIES,
    COMMON_GASES,
    FACILITIES,
    IPG6S,
    TANK7,
    THESIS_PRESETS,
    VENUS_AO_EV,
    NozzleGeometry,
)
from physics.mixture import UnknownSpecies, parse_mixture, validate_species
from physics.pipeline import solve_cached
from public_access import CORS_ORIGIN_REGEX, GRID_MAX, GRID_MIN, ApiKeyMiddleware, api_key

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "frontend" / "dist"

app = FastAPI(
    title="IRS collisionless plume",
    description="NASA CEA rocket stations + Khasawneh–Cai 2-D free-molecular jet.",
    version="1.5.0",
)
# API key first (inner), CORS last so it wraps 401s for grok.me / localhost.
app.add_middleware(ApiKeyMiddleware)
app.add_middleware(
    CORSMiddleware,
    allow_origins=[],
    allow_origin_regex=rf"^(?:{CORS_ORIGIN_REGEX})$",
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["Authorization", "Content-Type", "X-API-Key"],
    expose_headers=["*"],
    max_age=86400,
)


class SolveRequest(BaseModel):
    pinj_Pa: float = Field(100.0, gt=1.0, lt=2e5)
    hinj_MJ_kg: float = Field(23.0, gt=-20.0, lt=80.0)
    mixture: dict[str, float] | None = None
    basis: str = "mole"
    d_c_mm: float = Field(37.0, gt=1.0, lt=500.0)
    d_t_mm: float = Field(20.0, gt=1.0, lt=500.0)
    d_e_mm: float = Field(40.0, gt=1.0, lt=500.0)
    nozzle_name: str = "IPG6-S"
    a_c_mm2: Optional[float] = None
    a_t_mm2: Optional[float] = None
    a_e_mm2: Optional[float] = None
    xmax_m: float = Field(2.0, gt=0.2, lt=6.0)
    ymax_m: float = Field(1.0, gt=0.2, lt=2.0)
    nx: int = Field(65, ge=GRID_MIN, le=GRID_MAX)
    ny: int = Field(65, ge=GRID_MIN, le=GRID_MAX)
    mode: str = Field("enthalpy", description="enthalpy | generator | power | point")
    plume_mode: str = Field("auto", description="auto | collisionless | sudden_freeze")
    p_tank_Pa: float = Field(10.0, gt=0.05, lt=2e5, description="Ambient / tank pressure for the shock overlay")
    # legacy
    gas: Optional[str] = None
    he_mole_frac: float = Field(0.0, ge=0.0, le=0.99)
    power_W: Optional[float] = None
    mdot_mg_s: Optional[float] = None

    @field_validator("nx", "ny", mode="before")
    @classmethod
    def _cap_grid(cls, v):
        n = int(v)
        return max(GRID_MIN, min(n, GRID_MAX))


class MixPreview(BaseModel):
    mixture: dict[str, float]
    basis: str = "mole"


class CharacteristicsRequest(BaseModel):
    pinj_ref_Pa: float = Field(100.0, gt=1.0, lt=2e5)
    mixture: dict[str, float] | None = None
    basis: str = "mole"
    gas: Optional[str] = None
    he_mole_frac: float = Field(0.0, ge=0.0, le=0.99)
    d_c_mm: float = Field(37.0, gt=1.0, lt=500.0)
    d_t_mm: float = Field(20.0, gt=1.0, lt=500.0)
    d_e_mm: float = Field(40.0, gt=1.0, lt=500.0)
    nozzle_name: str = "IPG6-S"
    a_c_mm2: Optional[float] = None
    a_t_mm2: Optional[float] = None
    a_e_mm2: Optional[float] = None
    hinj_min: Optional[float] = Field(None, gt=-20.0, lt=80.0)
    hinj_max: float = Field(40.0, gt=-15.0, lt=80.0)
    n_h: int = Field(29, ge=3, le=80)
    mdot_mg_s_lines: Optional[list[float]] = None
    power_W_lines: Optional[list[float]] = None
    ions: bool = True


def _presets_resolved() -> list[dict[str, Any]]:
    h_co2 = co2_enthalpy_20C() / 1e6
    out = []
    for p in THESIS_PRESETS:
        item = dict(p)
        if "delta_h_MJ_kg" in item and "hinj_MJ_kg" not in item:
            item["hinj_MJ_kg"] = h_co2 + float(item["delta_h_MJ_kg"])
            item["h_ref_MJ_kg"] = h_co2
        fac = FACILITIES[item["facility"]]
        item["geometry"] = fac["geometry"]
        out.append(item)
    return out


def _geometry_from_body(body) -> NozzleGeometry:
    name = getattr(body, "nozzle_name", None) or "custom"
    if getattr(body, "a_c_mm2", None) and body.a_t_mm2 and body.a_e_mm2:
        return NozzleGeometry.from_areas_mm2(name, body.a_c_mm2, body.a_t_mm2, body.a_e_mm2)
    return NozzleGeometry.from_mm(name, body.d_c_mm, body.d_t_mm, body.d_e_mm)


@app.options("/api/{rest:path}")
def api_preflight(rest: str) -> Response:
    """Explicit OPTIONS so Grok Build / browsers always get a preflight on /api/*."""
    return Response(status_code=204)


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {
        "ok": True,
        "cea_version": CEA_VERSION,
        "geometry": IPG6S.as_dict(),
        "tank": TANK7,
        "grid_max": GRID_MAX,
        "api_key_required": bool(api_key()),
    }


@app.get("/api/catalog")
def catalog() -> dict[str, Any]:
    return {
        "facilities": list(FACILITIES.values()),
        "gases": COMMON_GASES,
        "advanced_species": ADVANCED_SPECIES,
        "air_recipe": {"N2": 0.79, "O2": 0.21},
        "presets": _presets_resolved(),
        "venus_ao_eV": VENUS_AO_EV,
        "cea_version": CEA_VERSION,
        "tank": TANK7,
        "units": {
            "pc": "bar (p_inj Pa × 1e-5)",
            "h_inj": "MJ/kg absolute CEA enthalpy",
            "h_ref": "MJ/kg at 298.15 K",
        },
    }


@app.get("/api/presets")
def presets() -> dict[str, Any]:
    return {"presets": _presets_resolved(), "venus_ao_eV": VENUS_AO_EV, "facilities": list(FACILITIES.values())}


@app.post("/api/mixture-preview")
def mixture_preview(body: MixPreview) -> dict[str, Any]:
    try:
        spec = parse_mixture(body.mixture, basis=body.basis)
    except (UnknownSpecies, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return spec.as_dict()


@app.get("/api/species")
def species_check(name: str) -> dict[str, Any]:
    try:
        ok = validate_species(name)
        return {"ok": True, "name": ok}
    except UnknownSpecies as exc:
        return {"ok": False, "name": name, "error": str(exc)}


@app.post("/api/solve")
def solve(req: SolveRequest) -> dict[str, Any]:
    mode = (req.mode or "enthalpy").strip().lower()
    if mode == "point":
        mode = "enthalpy"
    if mode not in ("enthalpy", "generator", "power"):
        raise HTTPException(status_code=400, detail="mode must be enthalpy, generator, power, or point")

    hinj = req.hinj_MJ_kg
    # Generator setup: operators set ṁ (MFC) and measure pinj; invert CEA for hinj.
    # Do NOT use P/ṁ — that is not how IPG enthalpy is recovered.
    if mode == "generator":
        if req.mdot_mg_s is None or req.mdot_mg_s <= 0:
            raise HTTPException(status_code=400, detail="generator mode needs mdot_mg_s > 0 plus pinj_Pa")
    elif mode == "power":
        if req.power_W is not None and req.mdot_mg_s is not None and req.mdot_mg_s > 0:
            hinj = req.power_W / (req.mdot_mg_s * 1e-6) / 1e6
            hinj = float(max(-15.0, min(hinj, 80.0)))

    mixture = req.mixture
    if not mixture and req.gas:
        mixture = None
    else:
        mixture = mixture or {"O2": 1.0}
        mix_key = json.dumps({k: float(v) for k, v in mixture.items() if float(v) > 0}, sort_keys=True)

    if req.a_c_mm2 and req.a_t_mm2 and req.a_e_mm2:
        geom = NozzleGeometry.from_areas_mm2(req.nozzle_name, req.a_c_mm2, req.a_t_mm2, req.a_e_mm2)
        d_c, d_t, d_e = geom.d_c * 1e3, geom.d_t * 1e3, geom.d_e * 1e3
    else:
        d_c, d_t, d_e = req.d_c_mm, req.d_t_mm, req.d_e_mm

    try:
        if mixture is None:
            from physics.mixture import parse_mixture as _pm

            spec = _pm(None, gas=req.gas, he_mole_frac=req.he_mole_frac)
            mix_key = json.dumps(spec.as_dict()["mole_fractions"], sort_keys=True)
            basis = "mole"
        else:
            basis = req.basis

        if mode == "generator":
            geom = NozzleGeometry.from_mm(req.nozzle_name or "custom", d_c, d_t, d_e)
            mix_dict = json.loads(mix_key)
            hinj, _ = hinj_for_target_mdot(
                pinj_Pa=float(req.pinj_Pa),
                mdot_mg_s=float(req.mdot_mg_s),
                geometry=geom,
                mixture=mix_dict,
                basis=basis,
            )

        pmode = (req.plume_mode or "auto").strip().lower()
        if pmode in ("freeze", "collisional"):
            pmode = "sudden_freeze"
        if pmode not in ("auto", "collisionless", "sudden_freeze"):
            raise HTTPException(status_code=400, detail="plume_mode must be auto, collisionless, or sudden_freeze")

        payload = solve_cached(
            float(req.pinj_Pa),
            float(round(hinj, 4)),
            mix_key,
            basis,
            float(d_c),
            float(d_t),
            float(d_e),
            req.nozzle_name or "custom",
            float(req.xmax_m),
            float(req.ymax_m),
            int(req.nx),
            int(req.ny),
            pmode,
            float(req.p_tank_Pa),
        )
        payload["mode"] = mode
        if mode == "generator":
            payload["requested_mdot_mg_s"] = float(req.mdot_mg_s)
        return payload
    except (UnknownSpecies, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"CEA/plume solve failed: {exc}") from exc



@app.post("/api/characteristics")
def characteristics(req: CharacteristicsRequest) -> dict[str, Any]:
    """1-D hinj sweep at pinj_ref → mdot/power isolines and composition kinks."""
    n_h = max(5, min(41, int(req.n_h)))
    mixture = req.mixture
    gas = req.gas
    if not mixture and gas:
        mixture = None
    else:
        mixture = mixture or {"O2": 1.0}
        gas = None
    try:
        geom = _geometry_from_body(req)
        payload = characteristics_sweep(
            pinj_ref_Pa=float(req.pinj_ref_Pa),
            mixture=mixture,
            basis=req.basis,
            gas=gas,
            he_mole_frac=float(req.he_mole_frac),
            geometry=geom,
            hinj_min=req.hinj_min,
            hinj_max=float(req.hinj_max),
            n_h=n_h,
            mdot_mg_s_lines=req.mdot_mg_s_lines,
            power_W_lines=req.power_W_lines,
            ions=bool(req.ions),
        )
        return payload
    except (UnknownSpecies, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"CEA characteristics sweep failed: {exc}") from exc


@app.get("/api/kinks")
def kinks(pinj_Pa: float = 100.0) -> dict[str, Any]:
    hs = [0.8, 1.2, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0, 18.0, 21.0, 23.0, 26.0, 28.0, 30.0, 35.0, 40.0]
    rows = []
    for h in hs:
        r = hp_equilibrium_o2(pinj_Pa, h, ions=True)
        rows.append(
            {
                "hinj_MJ_kg": h,
                "T": r["T"],
                "MW": r["MW"],
                "x_O": r["x_O"],
                "x_O2": r["x_O2"],
                "x_O+": r["x_O+"],
                "converged": r["converged"],
            }
        )
    return {"pinj_Pa": pinj_Pa, "points": rows}


if DIST.is_dir():
    app.mount("/", StaticFiles(directory=str(DIST), html=True), name="ui")
