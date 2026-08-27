"""NASA CEA rocket-mode wrapper for assigned-enthalpy injection.

Uses the official modernized NASA CEA Python package (`pip install cea`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

import cea

from .constants import IPG6S, N_A, NozzleGeometry, R_UNIV
from .mixture import MixtureSpec, parse_mixture

CEA_VERSION = getattr(cea, "__version__", "unknown")

_ALWAYS = (
    "O", "O2", "O+", "O2+", "O-", "N", "N2", "N+", "N2+", "NO", "NO+",
    "He", "He+", "Ar", "Ar+", "C", "C+", "CO", "CO2", "C2", "e-", "H", "H2",
)


@dataclass
class Station:
    name: str
    position: str
    T: float
    p_Pa: float
    h_kJ_kg: float
    gamma: float
    Mach: float
    MW: float
    rho: float
    u: float
    sonic: float
    ae_at: float
    mole_fractions: dict[str, float]


@dataclass
class CEAResult:
    version: str
    pinj_Pa: float
    hinj_MJ_kg: float
    hc_over_R: float
    geometry: dict
    mixture: dict
    ions: bool
    converged: bool
    last_error: int
    stations: list[Station]
    exit: dict[str, Any]
    mdot_kg_s: float
    power_W: float
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "pinj_Pa": self.pinj_Pa,
            "hinj_MJ_kg": self.hinj_MJ_kg,
            "delta_h_MJ_kg": self.hinj_MJ_kg - self.mixture.get("h_ref_MJ_kg", 0.0),
            "hc_over_R": self.hc_over_R,
            "geometry": self.geometry,
            "mixture": self.mixture,
            "ions": self.ions,
            "converged": self.converged,
            "last_error": self.last_error,
            "mdot_kg_s": self.mdot_kg_s,
            "mdot_mg_s": self.mdot_kg_s * 1e6,
            "power_W": self.power_W,
            "notes": self.notes,
            "stations": [
                {
                    "name": s.name,
                    "position": s.position,
                    "T": s.T,
                    "p_Pa": s.p_Pa,
                    "h_MJ_kg": s.h_kJ_kg / 1000.0,
                    "h_kJ_kg": s.h_kJ_kg,
                    "gamma": s.gamma,
                    "Mach": s.Mach,
                    "MW": s.MW,
                    "rho": s.rho,
                    "u": s.u,
                    "sonic": s.sonic,
                    "ae_at": s.ae_at,
                    "mole_fractions": s.mole_fractions,
                }
                for s in self.stations
            ],
            "exit": self.exit,
            # legacy keys so older UI/tests keep working
            "gas": "+".join(self.mixture.get("mole_fractions", {})),
            "he_mole_frac": self.mixture.get("mole_fractions", {}).get("He", 0.0),
        }


def _mole_fracs_at(sol: cea.RocketSolution, idx: int, cutoff: float = 1e-8) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, arr in sol.mole_fractions.items():
        val = float(arr[idx])
        if val >= cutoff or name in _ALWAYS:
            out[name] = val
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def _station(sol: cea.RocketSolution, idx: int, name: str, position: str) -> Station:
    mach = float(sol.Mach[idx])
    sonic = float(sol.sonic_velocity[idx])
    u = float(sol.Isp[idx]) if np.isfinite(sol.Isp[idx]) else mach * sonic
    if not np.isfinite(u):
        u = 0.0
    return Station(
        name=name,
        position=position,
        T=float(sol.T[idx]),
        p_Pa=float(sol.P[idx]) * 1e5,
        h_kJ_kg=float(sol.enthalpy[idx]),
        gamma=float(sol.gamma_s[idx]),
        Mach=mach,
        MW=float(sol.MW[idx]),
        rho=float(sol.density[idx]),
        u=u,
        sonic=sonic,
        ae_at=float(sol.ae_at[idx]),
        mole_fractions=_mole_fracs_at(sol, idx),
    )


def run_rocket(
    pinj_Pa: float,
    hinj_MJ_kg: float,
    gas: str | None = "O2",
    he_mole_frac: float = 0.0,
    mixture: Mapping[str, float] | None = None,
    basis: str = "mole",
    geometry: NozzleGeometry | None = None,
    ions: bool = True,
    freeze_at_throat: bool = False,
) -> CEAResult:
    """Assigned-enthalpy rocket problem for an arbitrary IPG nozzle + gas mix."""
    cea.init()
    geom = geometry or IPG6S
    spec: MixtureSpec = parse_mixture(
        mixture, basis=basis, gas=gas, he_mole_frac=he_mole_frac
    )
    pc_bar = pinj_Pa * 1e-5
    h_J_kg = hinj_MJ_kg * 1e6
    hc = h_J_kg / cea.R

    subar = geom.subar
    supar = geom.supar
    subar_arg = None if abs(subar - 1.0) < 1e-6 else subar
    supar_arg = 1.0001 if abs(supar - 1.0) < 1e-6 else supar

    kwargs: dict[str, Any] = dict(
        iac=True,
        hc=hc,
        tc_est=max(1500.0, 300.0 + abs(hinj_MJ_kg) * 200.0),
    )
    if subar_arg is not None:
        kwargs["subar"] = subar_arg
    if supar_arg is not None:
        kwargs["supar"] = supar_arg
    if freeze_at_throat:
        kwargs["n_frz"] = 2

    used_ions = ions
    sol = _rocket_solve(spec, used_ions, pc_bar, kwargs)
    notes = []
    if ions and not _sol_usable(sol):
        notes.append(
            "CEA rocket with ions produced a non-physical state (typical for pure He); "
            "retrying with ions off."
        )
        used_ions = False
        sol = _rocket_solve(spec, False, pc_bar, kwargs)

    notes = [
        f"NASA CEA {CEA_VERSION} rocket application, IAC, ions={used_ions}.",
        "hc = h_inj(J/kg) / cea.R with cea.R = 8314.51 J/(kmol·K).",
        "h_inj is the absolute CEA specific enthalpy (kJ/kg × 1e-3 = MJ/kg).",
        f"Mixture h_ref(298 K) = {spec.h_ref_MJ_kg:.3f} MJ/kg; Δh = h_inj − h_ref.",
        "pc is bar (p_inj Pa × 1e-5). Exit composition is frozen into the plume.",
        f"H = D_E/2 = {geom.H*1e3:.2f} mm for the 2-D collisionless jet.",
        *notes,
    ]

    mach = np.asarray(sol.Mach, dtype=float)
    ae = np.asarray(sol.ae_at, dtype=float)
    n = sol.num_pts
    i_ch = int(np.nanargmin(np.abs(mach)))
    i_th = int(np.nanargmin(np.abs(mach - 1.0)))
    supersonic = np.where(np.nan_to_num(mach, nan=0.0) > 1.05)[0]
    i_ex = int(supersonic[np.argmax(ae[supersonic])]) if len(supersonic) else n - 1

    stations = [
        _station(sol, i_ch, "chamber", "2"),
        _station(sol, i_th, "throat", "3"),
        _station(sol, i_ex, "exit", "4"),
    ]

    u_t = stations[1].u if stations[1].u > 1.0 else stations[1].Mach * stations[1].sonic
    mdot = stations[1].rho * u_t * geom.a_t
    power = mdot * hinj_MJ_kg * 1e6

    st_e = stations[2]
    if not np.isfinite(st_e.MW) or st_e.MW < 0.2:
        raise RuntimeError(
            "CEA rocket exit MW is invalid for this mixture and assigned enthalpy. "
            "Try a different h_inj or add a molecular gas (pure He with ions is a known CEA failure)."
        )
    R_sp = R_UNIV / st_e.MW
    n0 = st_e.rho * N_A / (st_e.MW * 1e-3)
    mf = st_e.mole_fractions

    if not sol.converged:
        notes.append(f"CEA did not fully converge (last_error={sol.last_error}). Treat numbers as qualitative.")

    mix_dict = spec.as_dict()
    mix_dict["delta_h_MJ_kg"] = hinj_MJ_kg - spec.h_ref_MJ_kg

    exit_state = {
        "T0": st_e.T,
        "p_Pa": st_e.p_Pa,
        "U0": st_e.u,
        "R": R_sp,
        "n0": n0,
        "rho": st_e.rho,
        "MW": st_e.MW,
        "gamma": st_e.gamma,
        "Mach": st_e.Mach,
        "H": geom.H,
        "S0": st_e.u / np.sqrt(2.0 * R_sp * st_e.T) if st_e.T > 0 else 0.0,
        "mole_fractions": mf,
        "h_kJ_kg": st_e.h_kJ_kg,
        "x_O": mf.get("O", 0.0),
        "x_O2": mf.get("O2", 0.0),
        "x_He": mf.get("He", 0.0),
        "x_N": mf.get("N", 0.0),
        "x_N2": mf.get("N2", 0.0),
        "x_C": mf.get("C", 0.0),
        "x_CO": mf.get("CO", 0.0),
        "x_CO2": mf.get("CO2", 0.0),
        "x_Ar": mf.get("Ar", 0.0),
        "x_ion": sum(v for k, v in mf.items() if k.endswith("+")),
    }

    return CEAResult(
        version=CEA_VERSION,
        pinj_Pa=pinj_Pa,
        hinj_MJ_kg=hinj_MJ_kg,
        hc_over_R=hc,
        geometry=geom.as_dict(),
        mixture=mix_dict,
        ions=used_ions,
        converged=bool(sol.converged),
        last_error=int(sol.last_error),
        stations=stations,
        exit=exit_state,
        mdot_kg_s=float(mdot),
        power_W=float(power),
        notes=notes,
    )


def _rocket_solve(spec: MixtureSpec, ions: bool, pc_bar: float, kwargs: dict[str, Any]):
    reac = cea.Mixture(spec.names, ions=ions)
    prod = cea.Mixture(spec.names, products_from_reactants=True, ions=ions)
    solver = cea.RocketSolver(prod, reactants=reac, ions=ions, transport=True)
    sol = cea.RocketSolution(solver)
    solver.solve(sol, spec.mass_fracs, pc_bar, **kwargs)
    return sol


def _sol_usable(sol: cea.RocketSolution) -> bool:
    if getattr(sol, "num_pts", 0) < 2:
        return False
    mw = np.asarray(sol.MW, dtype=float)
    T = np.asarray(sol.T, dtype=float)
    ok = (
        np.isfinite(mw)
        & (mw > 0.5)
        & np.isfinite(T)
        & (T > 50.0)
        & (T < 8.0e4)
    )
    return int(np.count_nonzero(ok)) >= 2


def co2_enthalpy_20C() -> float:
    """CEA product enthalpy of CO2 at 20 °C, J/kg (thesis: ≈ −8.9 MJ/kg)."""
    cea.init()
    reac = cea.Mixture(["CO2"])
    return float(reac.calc_property(cea.ENTHALPY, np.array([1.0]), 293.15))


def hp_equilibrium_o2(pinj_Pa: float, hinj_MJ_kg: float, ions: bool = True) -> dict[str, Any]:
    cea.init()
    pc = pinj_Pa * 1e-5
    hc = hinj_MJ_kg * 1e6 / cea.R
    reac = cea.Mixture(["O2"], ions=ions)
    prod = cea.Mixture(["O2"], products_from_reactants=True, ions=ions)
    solver = cea.EqSolver(prod, reactants=reac, ions=ions)
    sol = cea.EqSolution(solver)
    solver.solve(sol, cea.HP, hc, pc, np.array([1.0]))
    mf = {k: float(v) for k, v in sol.mole_fractions.items() if float(v) > 1e-10}
    return {
        "converged": bool(sol.converged),
        "T": float(sol.T) if sol.converged else float("nan"),
        "MW": float(sol.MW) if sol.converged else float("nan"),
        "h_kJ_kg": float(sol.enthalpy) if sol.converged else float("nan"),
        "gamma": float(sol.gamma_s) if sol.converged else float("nan"),
        "mole_fractions": mf,
        "x_O": mf.get("O", 0.0),
        "x_O2": mf.get("O2", 0.0),
        "x_O+": mf.get("O+", 0.0),
        "x_e": mf.get("e-", 0.0),
    }
