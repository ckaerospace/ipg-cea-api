"""NASA CEA rocket-mode wrapper for assigned-enthalpy injection.

Uses the official modernized NASA CEA Python package (`pip install cea`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

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
    notes = []
    sol = None
    try:
        sol = _rocket_solve(spec, used_ions, pc_bar, kwargs)
    except Exception:
        sol = None
    if ions and (sol is None or not _sol_usable(sol)):
        notes.append(
            "CEA rocket with ions produced a non-physical state (typical for pure He); "
            "retrying with ions off."
        )
        used_ions = False
        sol = _rocket_solve(spec, False, pc_bar, kwargs)
    if sol is None:
        raise RuntimeError("CEA rocket produced no solution for this assigned enthalpy.")

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


def hinj_for_target_mdot(
    pinj_Pa: float,
    mdot_mg_s: float,
    geometry: NozzleGeometry,
    mixture: Mapping[str, float] | None = None,
    basis: str = "mole",
    gas: str | None = None,
    he_mole_frac: float = 0.0,
    ions: bool = True,
    h_lo: float | None = None,
    h_hi: float = 75.0,
) -> tuple[float, CEAResult]:
    """Invert CEA rocket: operators set mass flow and measure pinj; hinj is solved.

    At fixed pinj and throat area, mdot falls as assigned enthalpy rises
    (hotter, lower density). Matches IPG operation with mass-flow controllers.
    """
    from scipy.optimize import brentq

    target = float(mdot_mg_s) * 1e-6
    if target <= 0:
        raise ValueError("mdot_mg_s must be positive")

    spec = parse_mixture(mixture, basis=basis, gas=gas, he_mole_frac=he_mole_frac)
    lo = spec.h_ref_MJ_kg + 0.3 if h_lo is None else float(h_lo)
    hi = float(h_hi)
    if hi <= lo + 0.5:
        hi = lo + 20.0

    cache: dict[float, CEAResult] = {}

    def run_h(h: float) -> CEAResult:
        key = round(float(h), 4)
        if key not in cache:
            cache[key] = run_rocket(
                pinj_Pa=pinj_Pa,
                hinj_MJ_kg=key,
                gas=gas,
                he_mole_frac=he_mole_frac,
                mixture=mixture,
                basis=basis,
                geometry=geometry,
                ions=ions,
            )
        return cache[key]

    def f(h: float) -> float:
        return run_h(h).mdot_kg_s - target

    # Expand the bracket if the target sits outside.
    m_lo = f(lo)
    m_hi = f(hi)
    tries = 0
    while m_lo * m_hi > 0 and tries < 8:
        if m_lo > 0 and m_hi > 0:
            # both too much flow → need hotter (lower mdot) → raise hi
            hi = min(80.0, hi + 10.0)
            m_hi = f(hi)
        else:
            # both too little flow → cooler
            lo = max(spec.h_ref_MJ_kg - 5.0, lo - 8.0)
            m_lo = f(lo)
        tries += 1
    if m_lo * m_hi > 0:
        raise ValueError(
            f"No CEA enthalpy matches mdot={mdot_mg_s:g} mg/s at pinj={pinj_Pa:g} Pa "
            f"(bracket mdot { (m_lo+target)*1e6:.3g} … {(m_hi+target)*1e6:.3g} mg/s)."
        )
    h_star = float(brentq(f, lo, hi, xtol=2e-4, maxiter=40))
    res = run_h(h_star)
    return h_star, res

# ---------------------------------------------------------------------------
# Characteristics field (pinj–hinj isolines) from a 1-D hinj sweep at pinj_ref.
# Approximation: at fixed hinj, composition and T are only weakly p-dependent,
# so mdot ≈ k(h) * pinj. ~n_h CEA rocket calls, not a 2-D grid.
# ---------------------------------------------------------------------------

DEFAULT_MDOT_MG_S_LINES: tuple[float, ...] = (2.0, 5.0, 8.0, 13.0, 20.0, 30.0, 50.0)
DEFAULT_POWER_W_LINES: tuple[float, ...] = (50.0, 150.0, 300.0, 450.0, 600.0)
CHAR_AXES_PINJ_PA: tuple[float, float] = (0.0, 250.0)
CHAR_AXES_HINJ_MJ_KG: tuple[float, float] = (0.0, 40.0)

_ALWAYS_X = ("O2", "O", "O+", "e-")
_NO_DISSOCIATION = frozenset({"He", "Ar", "Ne", "Kr", "Xe", "O", "N", "C", "H", "e-"})
_PARENT_PRIORITY = (
    ("O2", "O"),
    ("N2", "N"),
    ("CO2", "O"),
    ("CO", "C"),
    ("H2", "H"),
    ("NO", "N"),
    ("H2O", "H"),
)
_CHAR_NOTES = [
    "mdot isolines use mdot ≈ k(h)·pinj at fixed hinj (composition weakly p-dependent).",
    "Kinks mark composition change: energy into dissociation/ionization, T and pinj rise more slowly.",
]


def _parent_from_spec(spec: MixtureSpec) -> tuple[str, str]:
    """Dominant molecule (and its atom) for kink detection."""
    moles = {n: float(x) for n, x in zip(spec.names, spec.mole_fracs)}
    for mol, atom in _PARENT_PRIORITY:
        if moles.get(mol, 0.0) > 0.05:
            return mol, atom
    if not moles:
        return "O2", "O"
    top = max(moles, key=moles.get)
    atom = {"He": "He", "Ar": "Ar", "Ne": "Ne", "Kr": "Kr", "Xe": "Xe"}.get(top, top.rstrip("234"))
    return top, atom


def _interp_crossing(h: np.ndarray, y: np.ndarray, y_cross: float, *, rising: bool) -> float | None:
    """Linear hinj where y crosses y_cross. None if it never happens."""
    h = np.asarray(h, dtype=float)
    y = np.asarray(y, dtype=float)
    if h.size < 2:
        return None
    if rising:
        mask = y >= y_cross
    else:
        mask = y <= y_cross
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return None
    i = int(idx[0])
    if i == 0:
        # already on the far side at the first sample
        return float(h[0])
    y0, y1 = float(y[i - 1]), float(y[i])
    h0, h1 = float(h[i - 1]), float(h[i])
    if not np.isfinite(y0) or not np.isfinite(y1) or y1 == y0:
        return h1
    t = (float(y_cross) - y0) / (y1 - y0)
    t = float(np.clip(t, 0.0, 1.0))
    return h0 + t * (h1 - h0)


def detect_kinks(
    hinj: Sequence[float],
    x_O2: Sequence[float],
    x_O: Sequence[float],
    x_e: Sequence[float],
    parent: str = "O2",
) -> list[dict[str, Any]]:
    """Locate dissociation/ionization kinks from a 1-D hinj composition sweep.

    Finite differences on chamber (or any consistent station) mole fractions:
    - dissociation_start: |dx_mol/dh| first exceeds a small threshold while x_mol is still high
    - dissociation_end: x_mol has dropped and |dx_mol/dh| falls back after its peak
    - ionization_start: x_e rises through ~1e-3 (else a peak in |dx_e/dh|)

    ``x_O2`` / ``x_O`` are the parent-molecule and atom series (N2/N, CO2/O, …).
    Atomic gases (He, Ar, …) skip dissociation kinks.
    """
    h = np.asarray(hinj, dtype=float)
    x_mol = np.nan_to_num(np.asarray(x_O2, dtype=float), nan=0.0)
    x_atom = np.nan_to_num(np.asarray(x_O, dtype=float), nan=0.0)
    xe = np.nan_to_num(np.asarray(x_e, dtype=float), nan=0.0)
    n = int(h.size)
    if n < 5 or np.nanmax(h) - np.nanmin(h) < 0.5:
        return []

    with np.errstate(invalid="ignore", divide="ignore"):
        dmol = np.gradient(x_mol, h)
        de = np.gradient(xe, h)
    abs_dmol = np.abs(dmol)
    abs_de = np.abs(de)

    kinks: list[dict[str, Any]] = []
    parent = parent or "O2"
    do_dissoc = parent not in _NO_DISSOCIATION and float(np.nanmax(x_mol)) > 0.02

    if do_dissoc:
        x0 = float(x_mol[0])
        peak_dmol = float(np.nanmax(abs_dmol)) if np.isfinite(np.nanmax(abs_dmol)) else 0.0
        # Forward differences avoid np.gradient's noisy endpoint slope.
        dh = np.diff(h)
        dmol_fwd = np.abs(np.diff(x_mol) / np.where(np.abs(dh) < 1e-12, np.nan, dh))
        peak_fwd = float(np.nanmax(dmol_fwd)) if dmol_fwd.size and np.isfinite(np.nanmax(dmol_fwd)) else peak_dmol
        thr_start = max(0.0025, 0.04 * peak_fwd) if peak_fwd > 0 else 0.0025
        high = (0.70 * x0) if x0 > 0 else 0.40

        start_h = None
        start_i = None
        for i, slope in enumerate(dmol_fwd):
            # slope[i] lives on the segment h[i] → h[i+1]
            if not np.isfinite(slope) or slope < thr_start:
                continue
            if x_mol[i] < high and x_mol[i + 1] < high:
                continue
            start_i = i + 1
            if slope == dmol_fwd[i] and i > 0 and np.isfinite(dmol_fwd[i - 1]) and dmol_fwd[i - 1] < thr_start:
                t = (thr_start - dmol_fwd[i - 1]) / (slope - dmol_fwd[i - 1] + 1e-30)
                t = float(np.clip(t, 0.0, 1.0))
                start_h = float(h[i]) + t * float(h[i + 1] - h[i])
            else:
                start_h = float(h[i + 1])
            break
        if x0 > 0.05:
            drop = x0 - 0.02 * max(x0, 0.5)
            drop_h = _interp_crossing(h, x_mol, drop, rising=False)
            if drop_h is not None and drop_h > float(h[0]) + 0.3:
                j = int(np.argmin(np.abs(h - drop_h)))
                if x_mol[j] >= high * 0.9:
                    if start_h is None or drop_h < start_h:
                        start_h = drop_h
                        start_i = max(j, 1)
        if start_h is not None and start_h <= float(h[0]) + 0.15:
            start_h = None
            start_i = None
        if start_h is not None:
            # do not interpolate onset across a skipped-CEA hole (>2 MJ/kg)
            ia = int(np.searchsorted(h, start_h))
            ia = min(max(ia, 1), n - 1)
            if float(h[ia] - h[ia - 1]) > 2.0:
                start_h = float(h[ia])
                start_i = ia

        if start_h is not None:
            kinks.append(
                {
                    "hinj_MJ_kg": float(start_h),
                    "kind": "dissociation_start",
                    "label": f"{parent} dissociation start",
                }
            )

        # End: after the |dx/dh| peak, parent has dropped and the slope falls back.
        i0 = start_i if start_i is not None else 0
        peak_i = i0 + int(np.nanargmax(abs_dmol[i0:])) if i0 < n else i0
        # "dropped" ≈ residual parent after dissociation (O2 ~0.02 near 21 MJ/kg).
        dropped = max(0.003, 0.02 * x0) if x0 > 0 else 0.02
        fall_thr = 0.28 * peak_dmol if peak_dmol > 0 else 0.0

        end_h = None
        for i in range(max(peak_i, i0) + 1, n):
            if x_mol[i] <= dropped and abs_dmol[i] <= max(fall_thr, 0.006):
                # interpolate the x_mol drop through `dropped` (more stable than the slope)
                end_h = _interp_crossing(h[: i + 1], x_mol[: i + 1], dropped, rising=False)
                if end_h is None:
                    end_h = float(h[i])
                break
        if end_h is None:
            end_h = _interp_crossing(h[peak_i:], x_mol[peak_i:], dropped, rising=False)
            if end_h is not None:
                # require that we are past the steep region; otherwise leave as that crossing
                j = int(np.argmin(np.abs(h - end_h)))
                if abs_dmol[j] > max(0.55 * peak_dmol, 0.02) and j + 1 < n:
                    # still steep — walk forward until the slope relaxes
                    for i in range(j + 1, n):
                        if abs_dmol[i] <= max(fall_thr, 0.006):
                            end_h = float(h[i])
                            break

        if end_h is not None and (start_h is None or abs(end_h - start_h) > 1.5):
            kinks.append(
                {
                    "hinj_MJ_kg": float(end_h),
                    "kind": "dissociation_end",
                    "label": f"{parent} dissociation end",
                }
            )

    ion_h = _interp_crossing(h, xe, 1e-3, rising=True)
    ion_h_vis = _interp_crossing(h, xe, 1e-2, rising=True)
    diss_end_h = next((k["hinj_MJ_kg"] for k in kinks if k["kind"] == "dissociation_end"), None)
    if ion_h is not None and ion_h_vis is not None and diss_end_h is not None:
        if (ion_h - float(diss_end_h)) < 6.0:
            ion_h = ion_h_vis
    elif ion_h is None:
        ion_h = ion_h_vis
    if ion_h is None and float(np.nanmax(xe)) >= 1e-4 and float(np.nanmax(abs_de)) > 0:
        ion_i = int(np.nanargmax(abs_de))
        ion_h = float(h[ion_i])
    if ion_h is not None:
        too_close = any(abs(ion_h - k["hinj_MJ_kg"]) < 1.0 for k in kinks)
        if not too_close:
            kinks.append(
                {
                    "hinj_MJ_kg": float(ion_h),
                    "kind": "ionization_start",
                    "label": "ionization start",
                }
            )

    # keep in hinj order; drop NaNs
    kinks = [k for k in kinks if np.isfinite(k["hinj_MJ_kg"])]
    kinks.sort(key=lambda k: k["hinj_MJ_kg"])
    return kinks


def isolines_from_k(
    hinj: Sequence[float],
    k_kg_s_Pa: Sequence[float],
    href_MJ_kg: float,
    mdot_mg_s_lines: Sequence[float] | None = None,
    power_W_lines: Sequence[float] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Build mdot and coupled-power isolines from k(h) = mdot(h, pinj_ref) / pinj_ref.

    pinj(h) for target mdot:  mdot_target / k(h)
    Coupled power P = mdot * (hinj - href) * 1e6  →  mdot = P / ((hinj - href)*1e6)
    then pinj = mdot / k(h).
    """
    h = np.asarray(hinj, dtype=float)
    k = np.asarray(k_kg_s_Pa, dtype=float)
    mdot_targets = (
        list(DEFAULT_MDOT_MG_S_LINES) if mdot_mg_s_lines is None else [float(v) for v in mdot_mg_s_lines]
    )
    power_targets = (
        list(DEFAULT_POWER_W_LINES) if power_W_lines is None else [float(v) for v in power_W_lines]
    )

    mdot_isolines: list[dict[str, Any]] = []
    for m_mg in mdot_targets:
        if not np.isfinite(m_mg) or m_mg <= 0:
            continue
        m_kg = float(m_mg) * 1e-6
        pinj: list[float] = []
        hh: list[float] = []
        for hi, ki in zip(h, k):
            if not np.isfinite(ki) or ki <= 0 or not np.isfinite(hi):
                continue
            p = m_kg / float(ki)
            if np.isfinite(p) and p >= 0.0:
                pinj.append(float(p))
                hh.append(float(hi))
        mdot_isolines.append({"mdot_mg_s": float(m_mg), "pinj_Pa": pinj, "hinj_MJ_kg": hh})

    power_isolines: list[dict[str, Any]] = []
    href = float(href_MJ_kg)
    for P in power_targets:
        if not np.isfinite(P) or P <= 0:
            continue
        pinj = []
        hh = []
        for hi, ki in zip(h, k):
            dh = float(hi) - href
            if dh <= 1e-9 or not np.isfinite(ki) or ki <= 0 or not np.isfinite(hi):
                continue
            m_kg = float(P) / (dh * 1e6)
            p = m_kg / float(ki)
            if np.isfinite(p) and p >= 0.0:
                pinj.append(float(p))
                hh.append(float(hi))
        power_isolines.append({"power_W": float(P), "pinj_Pa": pinj, "hinj_MJ_kg": hh})
    return mdot_isolines, power_isolines


def _collect_x_series(rows: list[dict[str, float]]) -> dict[str, list[float]]:
    keys: set[str] = set(_ALWAYS_X)
    for mf in rows:
        for name, val in mf.items():
            try:
                if float(val) > 1e-3:
                    keys.add(str(name))
            except (TypeError, ValueError):
                continue
    ordered: dict[str, list[float]] = {}
    for name in _ALWAYS_X:
        ordered[name] = [float(mf.get(name, 0.0) or 0.0) for mf in rows]
        keys.discard(name)
    for name in sorted(keys):
        ordered[name] = [float(mf.get(name, 0.0) or 0.0) for mf in rows]
    return ordered


def _cap_line_values(vals: Sequence[float] | None, default: Sequence[float], cap: int = 20) -> list[float]:
    if vals is None:
        return [float(v) for v in default]
    out: list[float] = []
    for v in list(vals)[:cap]:
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if np.isfinite(fv) and fv > 0:
            out.append(fv)
    return out


def characteristics_sweep(
    pinj_ref_Pa: float = 100.0,
    mixture: Mapping[str, float] | None = None,
    basis: str = "mole",
    gas: str | None = None,
    he_mole_frac: float = 0.0,
    geometry: NozzleGeometry | None = None,
    hinj_min: float | None = None,
    hinj_max: float = 40.0,
    n_h: int = 29,
    mdot_mg_s_lines: Sequence[float] | None = None,
    power_W_lines: Sequence[float] | None = None,
    ions: bool = True,
) -> dict[str, Any]:
    """Sweep hinj at one pinj_ref and return k(h), composition, kinks, and isolines."""
    geom = geometry or IPG6S
    spec = parse_mixture(mixture, basis=basis, gas=gas, he_mole_frac=he_mole_frac)
    href = float(spec.h_ref_MJ_kg)
    pinj_ref = float(pinj_ref_Pa)
    if not np.isfinite(pinj_ref) or pinj_ref <= 0:
        raise ValueError("pinj_ref_Pa must be positive")

    h_min = href + 0.2 if hinj_min is None else float(hinj_min)
    h_max = float(hinj_max)
    if not np.isfinite(h_min) or not np.isfinite(h_max) or h_max <= h_min + 0.2:
        raise ValueError("hinj_max must be greater than hinj_min")
    n = int(n_h) if n_h is not None else 29
    n = max(5, min(41, n))
    grid = np.linspace(h_min, h_max, n)

    mdot_lines = _cap_line_values(mdot_mg_s_lines, DEFAULT_MDOT_MG_S_LINES)
    power_lines = _cap_line_values(power_W_lines, DEFAULT_POWER_W_LINES)

    hinj_ok: list[float] = []
    T_ch: list[float] = []
    MW_ch: list[float] = []
    mdot_mg: list[float] = []
    k_arr: list[float] = []
    x_ch_rows: list[dict[str, float]] = []
    T_ex: list[float] = []
    MW_ex: list[float] = []
    x_ex_rows: list[dict[str, float]] = []

    parent, atom = _parent_from_spec(spec)

    for h in grid:
        try:
            res = run_rocket(
                pinj_Pa=pinj_ref,
                hinj_MJ_kg=float(h),
                gas=gas,
                he_mole_frac=he_mole_frac,
                mixture=mixture,
                basis=basis,
                geometry=geom,
                ions=bool(ions),
            )
        except Exception:
            continue
        mdot = float(res.mdot_kg_s)
        if not np.isfinite(mdot) or mdot <= 0:
            continue
        st_ch = next((s for s in res.stations if s.name == "chamber"), res.stations[0] if res.stations else None)
        if st_ch is None or not np.isfinite(st_ch.MW) or st_ch.MW < 0.2:
            continue
        ex = res.exit or {}
        if not np.isfinite(float(ex.get("MW", float("nan")))):
            continue
        hinj_ok.append(float(res.hinj_MJ_kg))
        T_ch.append(float(st_ch.T))
        MW_ch.append(float(st_ch.MW))
        mdot_mg.append(mdot * 1e6)
        k_arr.append(mdot / pinj_ref)
        x_ch_rows.append(dict(st_ch.mole_fractions or {}))
        T_ex.append(float(ex.get("T0", float("nan"))))
        MW_ex.append(float(ex.get("MW", float("nan"))))
        x_ex_rows.append(dict(ex.get("mole_fractions") or {}))

    if len(hinj_ok) < 3:
        raise RuntimeError(
            f"CEA characteristics sweep produced only {len(hinj_ok)} usable hinj points "
            f"at pinj_ref={pinj_ref:g} Pa (n_h={n}). Try a different enthalpy window."
        )

    x_chamber = _collect_x_series(x_ch_rows)
    x_exit = _collect_x_series(x_ex_rows)
    # union of keys so chamber/exit x dicts align for plotting
    all_keys: list[str] = list(_ALWAYS_X)
    extra = sorted(set(x_chamber) | set(x_exit) - set(_ALWAYS_X))
    all_keys.extend(extra)
    x_chamber = {k: x_chamber.get(k, [0.0] * len(hinj_ok)) for k in all_keys}
    x_exit = {k: x_exit.get(k, [0.0] * len(hinj_ok)) for k in all_keys}

    x_mol = x_chamber.get(parent, x_chamber.get("O2", [0.0] * len(hinj_ok)))
    x_atom = x_chamber.get(atom, x_chamber.get("O", [0.0] * len(hinj_ok)))
    x_e = x_chamber.get("e-", [0.0] * len(hinj_ok))
    kinks = detect_kinks(hinj_ok, x_mol, x_atom, x_e, parent=parent)

    mdot_isolines, power_isolines = isolines_from_k(
        hinj_ok,
        k_arr,
        href,
        mdot_mg_s_lines=mdot_lines,
        power_W_lines=power_lines,
    )

    return {
        "pinj_ref_Pa": pinj_ref,
        "href_MJ_kg": href,
        "geometry": geom.as_dict(),
        "hinj_MJ_kg": hinj_ok,
        "k_kg_s_Pa": k_arr,
        "chamber": {
            "T": T_ch,
            "MW": MW_ch,
            "mdot_mg_s": mdot_mg,
            "x": x_chamber,
        },
        "exit": {
            "T0": T_ex,
            "MW": MW_ex,
            "x": x_exit,
        },
        "kinks": kinks,
        "mdot_isolines": mdot_isolines,
        "power_isolines": power_isolines,
        "axes": {
            "pinj_Pa": [CHAR_AXES_PINJ_PA[0], CHAR_AXES_PINJ_PA[1]],
            "hinj_MJ_kg": [CHAR_AXES_HINJ_MJ_KG[0], CHAR_AXES_HINJ_MJ_KG[1]],
        },
        "notes": list(_CHAR_NOTES),
        "ions": bool(ions),
        "parent_molecule": parent,
    }
