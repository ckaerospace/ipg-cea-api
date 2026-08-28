"""CEA rocket stations → collisionless plume grid for any nozzle + mixture."""

from __future__ import annotations

import json
from functools import lru_cache
from typing import Any, Mapping

import numpy as np

from .cea_rocket import run_rocket
from .constants import IPG6S, NozzleGeometry
from .derived import probe_quantities, probe_to_dict
from .plume import CollisionlessPlume, marching_squares
from .sudden_freeze import mix_d_hs, sudden_freeze_field

DENSITY_LEVELS = [0.8, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01]
SPEED_LEVELS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.35]
TEMP_LEVELS = [0.5, 0.525, 0.55, 0.575, 0.6, 0.7, 0.8]
HTOT_LEVELS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05]


def solve_operating_point(
    pinj_Pa: float,
    hinj_MJ_kg: float,
    gas: str | None = "O2",
    he_mole_frac: float = 0.0,
    mixture: Mapping[str, float] | None = None,
    basis: str = "mole",
    xmax_m: float = 2.0,
    ymax_m: float = 1.0,
    nx: int = 97,
    ny: int = 81,
    geometry: NozzleGeometry | None = None,
    ions: bool = True,
    include_contours: bool = True,
    plume_mode: str = "collisionless",
) -> dict[str, Any]:
    geom = geometry or IPG6S
    cea_res = run_rocket(
        pinj_Pa=pinj_Pa,
        hinj_MJ_kg=hinj_MJ_kg,
        gas=gas,
        he_mole_frac=he_mole_frac,
        mixture=mixture,
        basis=basis,
        geometry=geom,
        ions=ions,
    )
    ex = cea_res.exit
    plume = CollisionlessPlume.from_exit(
        T0=ex["T0"],
        R_specific=ex["R"],
        U0=ex["U0"],
        n0=ex["n0"],
        H=ex["H"],
    )
    x = np.linspace(0.0, xmax_m, nx)
    y = np.linspace(-ymax_m, ymax_m, ny)
    field = plume.grid(x, y)
    freeze_meta = {
        "mode": "collisionless",
        "kn_gll_exit": None,
        "r_freeze_m": None,
        "theta_jet_deg": None,
        "kn_crit": None,
        "d_hs": None,
    }
    if (plume_mode or "collisionless").strip().lower() in ("sudden_freeze", "freeze", "collisional"):
        d_hs = mix_d_hs(ex.get("mole_fractions") or {})
        field = sudden_freeze_field(
            x, y,
            T0=ex["T0"],
            R_specific=ex["R"],
            U0=ex["U0"],
            n0=ex["n0"],
            H=ex["H"],
            gamma=float(ex.get("gamma") or 1.4),
            Mach_e=float(ex.get("Mach") or 1.2),
            d_hs=d_hs,
            collisionless=plume,
        )
        freeze_meta = {
            "mode": "sudden_freeze",
            "kn_gll_exit": field.get("kn_gll_exit"),
            "r_freeze_m": field.get("r_freeze_m"),
            "theta_jet_deg": field.get("theta_jet_deg"),
            "kn_crit": field.get("kn_crit"),
            "d_hs": field.get("d_hs"),
        }

    # Frozen total enthalpy on the collisionless grid:
    # static h ≈ href + (h_exit − href) * T/T0, then add ½V².
    href_J = float(cea_res.mixture.get("h_ref_MJ_kg", 0.0)) * 1e6
    h_exit_J = float(ex.get("h_kJ_kg") or 0.0) * 1e3
    hinj_J = float(cea_res.hinj_MJ_kg) * 1e6
    t_ratio = np.asarray(field["t_ratio"], dtype=np.float64)
    speed = np.asarray(field["speed"], dtype=np.float64)
    h_static = href_J + (h_exit_J - href_J) * t_ratio
    h_tot = h_static + 0.5 * speed * speed
    h_tot_ratio = h_tot / hinj_J if abs(hinj_J) > 1.0 else np.zeros_like(h_tot)
    field["h_tot_MJ_kg"] = h_tot / 1e6
    field["h_tot_ratio"] = h_tot_ratio

    contours: dict[str, Any] = {}
    if include_contours:
        speed_ratio = field["speed"] / ex["U0"] if ex["U0"] > 0 else field["speed"]
        contours = {
            "n": marching_squares(x, y, field["n_ratio"], DENSITY_LEVELS),
            "speed": marching_squares(x, y, speed_ratio, SPEED_LEVELS),
            "T": marching_squares(x, y, field["t_ratio"], TEMP_LEVELS),
            "h_tot": marching_squares(x, y, h_tot_ratio, HTOT_LEVELS),
        }

    def _f32(a: np.ndarray) -> list[float]:
        return np.asarray(a, dtype=np.float32).ravel(order="C").tolist()

    return {
        "cea": cea_res.as_dict(),
        "plume": {
            "S0": plume.S0,
            "T0": plume.T0,
            "R": plume.R,
            "U0": plume.U0,
            "n0": plume.n0,
            "H": plume.H,
            "thermal": plume.thermal,
            "mode": freeze_meta["mode"],
            "kn_gll_exit": freeze_meta["kn_gll_exit"],
            "r_freeze_m": freeze_meta["r_freeze_m"],
            "theta_jet_deg": freeze_meta["theta_jet_deg"],
            "kn_crit": freeze_meta["kn_crit"],
            "d_hs": freeze_meta["d_hs"],
            "xmax_m": xmax_m,
            "ymax_m": ymax_m,
            "nx": nx,
            "ny": ny,
            "x": x.tolist(),
            "y": y.tolist(),
            "n_ratio": _f32(field["n_ratio"]),
            "u": _f32(field["u"]),
            "v": _f32(field["v"]),
            "t_ratio": _f32(field["t_ratio"]),
            "speed": _f32(field["speed"]),
            "h_tot_MJ_kg": _f32(field["h_tot_MJ_kg"]),
            "h_tot_ratio": _f32(field["h_tot_ratio"]),
            "contours": contours,
        },
        "caveats": [
            "CEA is equilibrium chemistry and can under-predict dissociation/ions relative to the RF plasma.",
            "The default plume is a 2-D planar collisionless jet applied to a round nozzle: H = D_E/2.",
            "plume_mode=sudden_freeze uses isentropic source flow until Kn_GLL=0.05 (Boyd/Bird), then freezes T.",
            "Atomic-oxygen kinetic energy in the tunnel is typically 2–3.8 eV, below Venus aerobraking (8.3 eV).",
            "Catalytic heat flux is a fully-catalytic (γ=1) upper-bound estimate from O recombination when O is present.",
        ],
    }


def interpolate_field(payload: dict[str, Any], xq: float, yq: float) -> dict[str, float]:
    p = payload["plume"]
    x = np.asarray(p["x"])
    y = np.asarray(p["y"])
    nx, ny = int(p["nx"]), int(p["ny"])

    def reshape(key: str) -> np.ndarray:
        return np.asarray(p[key], dtype=np.float64).reshape(ny, nx)

    n = reshape("n_ratio")
    u = reshape("u")
    v = reshape("v")
    t = reshape("t_ratio")
    xq = float(np.clip(xq, x[0], x[-1]))
    yq = float(np.clip(yq, y[0], y[-1]))
    ix = int(np.clip(np.searchsorted(x, xq) - 1, 0, nx - 2))
    iy = int(np.clip(np.searchsorted(y, yq) - 1, 0, ny - 2))
    tx = 0.0 if x[ix + 1] == x[ix] else (xq - x[ix]) / (x[ix + 1] - x[ix])
    ty = 0.0 if y[iy + 1] == y[iy] else (yq - y[iy]) / (y[iy + 1] - y[iy])

    def bl(a: np.ndarray) -> float:
        return float(
            (1 - tx) * (1 - ty) * a[iy, ix]
            + tx * (1 - ty) * a[iy, ix + 1]
            + (1 - tx) * ty * a[iy + 1, ix]
            + tx * ty * a[iy + 1, ix + 1]
        )

    return {"n_ratio": bl(n), "u": bl(u), "v": bl(v), "t_ratio": bl(t)}


def probe_at(payload: dict[str, Any], xq: float, yq: float) -> dict[str, Any]:
    samp = interpolate_field(payload, xq, yq)
    ex = payload["cea"]["exit"]
    reading = probe_quantities(
        x=xq,
        y=yq,
        n_ratio=samp["n_ratio"],
        u=samp["u"],
        v=samp["v"],
        t_ratio=samp["t_ratio"],
        n0=ex["n0"],
        T0=ex["T0"],
        mole_fractions=ex.get("mole_fractions") or {},
        H=ex["H"],
        x_O=ex.get("x_O"),
        x_O2=ex.get("x_O2"),
        x_He=ex.get("x_He"),
    )
    d = probe_to_dict(reading)
    d["n_ratio"] = samp["n_ratio"]
    d["t_ratio"] = samp["t_ratio"]
    d["speed_ratio"] = (d["speed"] / ex["U0"]) if ex["U0"] else 0.0
    return d


@lru_cache(maxsize=48)
def solve_cached(
    pinj_Pa: float,
    hinj_MJ_kg: float,
    mix_key: str,
    basis: str,
    d_c_mm: float,
    d_t_mm: float,
    d_e_mm: float,
    nozzle_name: str,
    xmax_m: float,
    ymax_m: float,
    nx: int,
    ny: int,
    plume_mode: str = "collisionless",
) -> dict[str, Any]:
    mixture = json.loads(mix_key)
    geom = NozzleGeometry.from_mm(nozzle_name, d_c_mm, d_t_mm, d_e_mm)
    return solve_operating_point(
        pinj_Pa=pinj_Pa,
        hinj_MJ_kg=hinj_MJ_kg,
        mixture=mixture,
        basis=basis,
        gas=None,
        xmax_m=xmax_m,
        ymax_m=ymax_m,
        nx=nx,
        ny=ny,
        geometry=geom,
        plume_mode=plume_mode,
    )
