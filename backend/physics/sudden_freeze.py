"""Sudden-freeze collisional plume (Bird / Dettleff / Miller).

Near the nozzle the expansion is treated as a planar isentropic source flow
(collisions maintain one T and p = n k T).  When the Boyd gradient-length
Knudsen number Kn_GLL = λ / R reaches KN_CRIT (0.05), translational
temperature freezes and the remainder of the ray is free-molecular
(n ∝ 1/R in this 2-D geometry, T = T_f, bulk speed held at U_f).

Outside the Prandtl–Meyer vacuum turning cone the field falls back to the
Khasawneh–Cai collisionless jet already used in collisionless mode.

References: Bird, Molecular Gas Dynamics (1994); Boyd, Chen & Candler,
Phys. Fluids 7 (1995); Dettleff, Prog. Aerosp. Sci. 28 (1991); Miller in
Scoles, Atomic and Molecular Beam Methods (1988).
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .constants import D_HS, K_B
from .plume import CollisionlessPlume

KN_CRIT = 0.05  # Boyd Kn_GLL / Bird P-order freeze


def resolve_plume_mode(plume_mode: str, kn_exit: float) -> tuple[str, str]:
    """Kn_exit is the only auto trigger. collisionless / sudden_freeze chips override Auto."""
    requested = (plume_mode or "auto").strip().lower()
    if requested in ("freeze", "collisional"):
        requested = "sudden_freeze"
    if requested in ("auto", ""):
        chosen = "sudden_freeze" if kn_exit < KN_CRIT else "collisionless"
        requested = "auto"
    elif requested in ("sudden_freeze", "collisionless"):
        chosen = requested
    else:
        chosen = "sudden_freeze" if kn_exit < KN_CRIT else "collisionless"
        requested = "auto"
    return requested, chosen


def kn_gll_exit(n0: float, H: float, d_hs: float) -> float:
    """λ/H at the nozzle lip. Continuum if this is well below KN_CRIT."""
    n0 = max(float(n0), 1.0)
    H = max(float(H), 1e-6)
    d_hs = max(float(d_hs), 1e-12)
    lam = 1.0 / (np.sqrt(2.0) * np.pi * d_hs * d_hs * n0)
    return float(lam / H)


def mix_d_hs(mole_fractions: Mapping[str, float] | None) -> float:
    num = 0.0
    den = 0.0
    for name, x in (mole_fractions or {}).items():
        xv = float(x)
        if xv <= 0.0:
            continue
        key = name.rstrip("+-")
        d = D_HS.get(name) or D_HS.get(key)
        if d:
            num += xv * d
            den += xv
    return num / den if den > 0 else 3.6e-10


def _A_over_Astar(M: np.ndarray, g: float) -> np.ndarray:
    M = np.maximum(np.asarray(M, dtype=np.float64), 1.001)
    gp = (g + 1.0) / 2.0
    expo = (g + 1.0) / (2.0 * (g - 1.0))
    return (gp ** (-expo)) * (1.0 + (g - 1.0) / 2.0 * M * M) ** expo / M


def _mach_from_area(ar: np.ndarray, g: float, Me: float) -> np.ndarray:
    """Invert A/A_e = ar for supersonic M, given exit Mach Me."""
    target = np.asarray(ar, dtype=np.float64) * _A_over_Astar(Me, g)
    M = np.full(target.shape, max(float(Me), 1.2), dtype=np.float64)
    M = np.where(ar <= 1.0, float(Me), M)
    for _ in range(20):
        f = _A_over_Astar(M, g) - target
        dM = 1e-4 * np.maximum(M, 1.0)
        df = (_A_over_Astar(M + dM, g) - _A_over_Astar(M, g)) / dM
        df = np.where(np.abs(df) < 1e-12, 1e-12, df)
        step = np.clip(f / df, -2.0, 2.0)
        M = np.clip(M - step, 1.001, 80.0)
    return M


def _prandtl_meyer(M: float, g: float) -> float:
    M = max(float(M), 1.0)
    k = np.sqrt((g + 1.0) / (g - 1.0))
    a = np.sqrt(M * M - 1.0)
    return float(k * np.arctan(a / k) - np.arctan(a))


def sudden_freeze_field(
    x: np.ndarray,
    y: np.ndarray,
    *,
    T0: float,
    R_specific: float,
    U0: float,
    n0: float,
    H: float,
    gamma: float,
    Mach_e: float,
    d_hs: float,
    collisionless: CollisionlessPlume,
) -> dict[str, np.ndarray]:
    """Return n_ratio, u, v, t_ratio, speed plus freeze diagnostics."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    X, Y = np.meshgrid(x, y)
    g = float(np.clip(gamma, 1.05, 1.67))
    Me = float(max(Mach_e, 1.05))
    cl = collisionless.evaluate(X, Y)

    R = np.hypot(X, np.maximum(np.abs(Y), 1e-12))
    R_e = max(float(H), 1e-4)
    theta = np.arctan2(Y, np.maximum(X, 1e-12))

    nu_max = 0.5 * np.pi * (np.sqrt((g + 1.0) / (g - 1.0)) - 1.0)
    nu_e = _prandtl_meyer(Me, g)
    theta_jet = float(np.clip(nu_max - nu_e, 0.35, 1.4))  # rad, ~20–80 deg

    ar = np.maximum(R / R_e, 1.0)
    M = _mach_from_area(ar, g, Me)
    T_ratio_isen = (1.0 + (g - 1.0) / 2.0 * Me * Me) / (1.0 + (g - 1.0) / 2.0 * M * M)
    T_isen = T0 * T_ratio_isen
    n_isen = n0 * (T_ratio_isen ** (1.0 / (g - 1.0)))  # n/n0 = (T/T0)^{1/(γ-1)}
    a = np.sqrt(g * R_specific * np.maximum(T_isen, 1.0))
    U_isen = M * a
    # Energy cap: cannot exceed stagnation from exit
    cp = g * R_specific / (g - 1.0)
    e0 = 0.5 * U0 * U0 + cp * T0
    U_max = np.sqrt(np.maximum(2.0 * e0, 0.0))
    U_isen = np.minimum(U_isen, U_max)

    lam = 1.0 / (np.sqrt(2.0) * np.pi * d_hs * d_hs * np.maximum(n_isen, 1.0))
    Kn = lam / np.maximum(R, R_e)
    kn_exit = float(1.0 / (np.sqrt(2.0) * np.pi * d_hs * d_hs * max(n0, 1.0)) / R_e)

    # Freeze radius on each ray: first R where Kn >= KN_CRIT (approx Kn ∝ R^γ for isentropic).
    frozen = Kn >= KN_CRIT
    # Reconstruct T_f, n_f, U_f by evaluating isentropic at R_f ≈ R * (KN_CRIT / Kn)^{something}
    # Kn = λ/R ∝ 1/(n R) and n ∝ R^{-1/(γ-1)} in 2-D source → Kn ∝ R^{1/(γ-1) - 1} wait
    # 2-D: n ∝ 1/R^{1/(γ-1)}? Area ∝ R so A/A_e = R/R_e, n ∝ T^{1/(γ-1)} and T drops with M(R).
    # Practical: along the ray, R_f = R * (KN_CRIT / Kn)^{1} clipped, then recompute isentropic at R_f.
    kn_safe = np.maximum(Kn, 1e-12)
    R_f = np.clip(R * (KN_CRIT / kn_safe), R_e, np.maximum(R, R_e))
    R_f = np.where(frozen, R_f, R)
    ar_f = np.maximum(R_f / R_e, 1.0)
    M_f = _mach_from_area(ar_f, g, Me)
    T_f_ratio = (1.0 + (g - 1.0) / 2.0 * Me * Me) / (1.0 + (g - 1.0) / 2.0 * M_f * M_f)
    T_f = T0 * T_f_ratio
    n_f = n0 * (T_f_ratio ** (1.0 / (g - 1.0)))
    U_f = np.minimum(M_f * np.sqrt(g * R_specific * np.maximum(T_f, 1.0)), U_max)

    n_free = n_f * (R_f / np.maximum(R, R_f))
    T_free = T_f
    U_free = U_f

    n_core = np.where(frozen, n_free, n_isen)
    T_core = np.where(frozen, T_free, T_isen)
    U_core = np.where(frozen, U_free, U_isen)

    cth = X / np.maximum(R, 1e-12)
    sth = Y / np.maximum(R, 1e-12)
    u_core = U_core * cth
    v_core = U_core * sth

    in_cone = np.abs(theta) <= theta_jet
    n_ratio = np.where(in_cone, n_core / max(n0, 1.0), cl["n_ratio"])
    t_ratio = np.where(in_cone, T_core / max(T0, 1.0), cl["t_ratio"])
    u = np.where(in_cone, u_core, cl["u"])
    v = np.where(in_cone, v_core, cl["v"])
    # Exit plane strip: hold CEA exit state across the slit.
    on_lip = (X <= 0.02 * H) & (np.abs(Y) <= H)
    n_ratio = np.where(on_lip, 1.0, n_ratio)
    t_ratio = np.where(on_lip, 1.0, t_ratio)
    u = np.where(on_lip, U0, u)
    v = np.where(on_lip, 0.0, v)

    n_ratio = np.clip(n_ratio, 0.0, 5.0)
    t_ratio = np.clip(t_ratio, 0.0, 5.0)
    speed = np.hypot(u, v)

    r_freeze = float(np.median(R_f[frozen])) if np.any(frozen) else float("inf")
    return {
        "n_ratio": n_ratio,
        "u": u,
        "v": v,
        "t_ratio": t_ratio,
        "speed": speed,
        "kn_gll": Kn,
        "frozen": frozen.astype(np.float64),
        "mode": "sudden_freeze",
        "kn_gll_exit": kn_exit,
        "r_freeze_m": r_freeze if np.isfinite(r_freeze) else None,
        "theta_jet_deg": float(np.degrees(theta_jet)),
        "kn_crit": KN_CRIT,
        "d_hs": float(d_hs),
    }
