"""Engineering shock overlay on a continuum (sudden-freeze) plume.

Not an Euler / Navier–Stokes / DSMC solve.  Underexpanded Mach-disk
station follows the Crist / Addy / Ashkenas–Sherman sonic-orifice family.
Overexpanded lips use a θ–β–M oblique shock.  Rankine–Hugoniot closes
the disk.  A Knudsen / sudden-freeze veto keeps vacuum-like tanks on
today's sudden-freeze field.

Kn_exit = λ/H at the lip is the only mode trigger (see pipeline).  This
module only runs a shock overlay when the chosen mode is continuum.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from .sudden_freeze import KN_CRIT, _mach_from_area, _prandtl_meyer

CRIST_XM_COEF = 0.67
MATCHED_LOG_NPR = 0.05
TRIPLE_X_FRAC = 0.35  # r_tp / x_m in the 0.3–0.4 band
N_BARREL = 25


def npr(p_e_Pa: float, p_tank_Pa: float) -> float:
    pe = max(float(p_e_Pa), 1e-30)
    pt = max(float(p_tank_Pa), 1e-30)
    return pe / pt


def classify_npr_regime(p_e_Pa: float, p_tank_Pa: float) -> str:
    """NPR-only label before freeze / kinetic vetoes."""
    ratio = npr(p_e_Pa, p_tank_Pa)
    if abs(np.log(ratio)) < MATCHED_LOG_NPR:
        return "matched"
    if ratio > 1.0:
        return "underexpanded"
    return "overexpanded"


def crist_mach_disk_x(p_e_Pa: float, p_tank_Pa: float, H: float) -> float:
    """x_m / D_e = 0.67 sqrt(p_e / p_tank), D_e = 2H."""
    de = 2.0 * max(float(H), 1e-9)
    return CRIST_XM_COEF * de * float(np.sqrt(npr(p_e_Pa, p_tank_Pa)))


def kn_gll_at_radius(
    R: float,
    *,
    n0: float,
    H: float,
    d_hs: float,
    gamma: float,
    Mach_e: float,
) -> float:
    """λ/R at radius R on the isentropic source-flow core."""
    R = max(float(R), max(float(H), 1e-6))
    g = float(np.clip(gamma, 1.05, 1.67))
    Me = float(max(Mach_e, 1.05))
    ar = np.array([R / max(float(H), 1e-6)], dtype=np.float64)
    M = float(_mach_from_area(ar, g, Me)[0])
    t_ratio = (1.0 + 0.5 * (g - 1.0) * Me * Me) / (1.0 + 0.5 * (g - 1.0) * M * M)
    n = max(float(n0), 1.0) * (t_ratio ** (1.0 / (g - 1.0)))
    d = max(float(d_hs), 1e-12)
    lam = 1.0 / (np.sqrt(2.0) * np.pi * d * d * max(n, 1.0))
    return float(lam / R)


def isentropic_core_at_R(
    R: float,
    *,
    T0: float,
    n0: float,
    U0: float,
    H: float,
    gamma: float,
    Mach_e: float,
    R_specific: float,
) -> dict[str, float]:
    g = float(np.clip(gamma, 1.05, 1.67))
    Me = float(max(Mach_e, 1.05))
    R = max(float(R), max(float(H), 1e-6))
    ar = np.array([R / max(float(H), 1e-6)], dtype=np.float64)
    M = float(_mach_from_area(ar, g, Me)[0])
    t_ratio = (1.0 + 0.5 * (g - 1.0) * Me * Me) / (1.0 + 0.5 * (g - 1.0) * M * M)
    T = float(T0) * t_ratio
    n = float(n0) * (t_ratio ** (1.0 / (g - 1.0)))
    a = float(np.sqrt(g * max(float(R_specific), 1e-12) * max(T, 1.0)))
    cp = g * float(R_specific) / (g - 1.0)
    e0 = 0.5 * float(U0) * float(U0) + cp * float(T0)
    U_max = float(np.sqrt(max(2.0 * e0, 0.0)))
    U = min(M * a, U_max)
    return {"M": M, "T": T, "n": n, "U": U, "t_ratio": t_ratio}


def rankine_hugoniot_normal(M1: float, gamma: float) -> dict[str, float]:
    """Normal-shock jump. Downstream is subsonic; n and T rise; U drops."""
    g = float(np.clip(gamma, 1.05, 1.67))
    M1 = float(max(M1, 1.01))
    gp, gm = g + 1.0, g - 1.0
    M2sq = (1.0 + 0.5 * gm * M1 * M1) / (g * M1 * M1 - 0.5 * gm)
    M2 = float(np.sqrt(max(M2sq, 1e-8)))
    rho_ratio = (gp * M1 * M1) / (gm * M1 * M1 + 2.0)
    p_ratio = 1.0 + (2.0 * g / gp) * (M1 * M1 - 1.0)
    T_ratio = p_ratio / rho_ratio
    return {
        "M2": M2,
        "rho_ratio": float(rho_ratio),
        "p_ratio": float(p_ratio),
        "T_ratio": float(T_ratio),
        "U_ratio": float(1.0 / rho_ratio),
    }


def _theta_from_beta(M: float, beta: float, g: float) -> float:
    s2 = np.sin(beta) ** 2
    num = 2.0 * (1.0 / np.tan(beta)) * (M * M * s2 - 1.0)
    den = M * M * (g + np.cos(2.0 * beta)) + 2.0
    if abs(den) < 1e-14:
        return 0.0
    return float(np.arctan(num / den))


def _theta_max(M: float, g: float) -> float:
    M = max(float(M), 1.01)
    mu = float(np.arcsin(min(1.0, 1.0 / M)))
    betas = np.linspace(mu + 1e-3, 0.5 * np.pi - 1e-3, 64)
    return float(max(_theta_from_beta(M, float(b), g) for b in betas))


def _beta_from_pressure_ratio(M1: float, p2_over_p1: float, g: float) -> float | None:
    """Wave angle β from the oblique-shock pressure ratio. None if detached."""
    M1 = max(float(M1), 1.01)
    pr = max(float(p2_over_p1), 1.0)
    p_ns = 1.0 + (2.0 * g / (g + 1.0)) * (M1 * M1 - 1.0)
    if pr >= p_ns * 0.999:
        return None
    s2 = (((g + 1.0) / (2.0 * g)) * (pr - 1.0) + 1.0) / (M1 * M1)
    if s2 >= 1.0 or s2 <= 1.0 / (M1 * M1) + 1e-12:
        return None
    beta = float(np.arcsin(np.sqrt(s2)))
    mu = float(np.arcsin(min(1.0, 1.0 / M1)))
    if beta <= mu:
        return None
    return beta


def _oblique_post(M1: float, beta: float, g: float) -> dict[str, float]:
    theta = _theta_from_beta(M1, beta, g)
    Mn1 = M1 * float(np.sin(beta))
    rh = rankine_hugoniot_normal(Mn1, g)
    dbeta = max(beta - theta, 1e-4)
    M2 = rh["M2"] / float(np.sin(dbeta))
    return {"theta": theta, "M2": float(M2), **rh}


def _barrel_underexpanded(H: float, x_m: float, r_tp: float) -> list[list[float]]:
    """Smooth lip → triple (upper half). Quadratic Bézier with a barrel bulge."""
    x_m = max(float(x_m), 1e-9)
    H = max(float(H), 0.0)
    r_tp = max(float(r_tp), 0.0)
    r_mid = max(H * 1.1, 0.50 * x_m)
    t = np.linspace(0.0, 1.0, N_BARREL)
    x = (1.0 - t) ** 2 * 0.0 + 2.0 * (1.0 - t) * t * (0.45 * x_m) + t ** 2 * x_m
    y = (1.0 - t) ** 2 * H + 2.0 * (1.0 - t) * t * r_mid + t ** 2 * r_tp
    y = np.maximum(y, 0.0)
    return [[float(a), float(b)] for a, b in zip(x, y)]


def _barrel_overexpanded(H: float, x_end: float, r_end: float, beta: float) -> list[list[float]]:
    """Lip → triple / axis along the incident wave (upper half, y >= 0)."""
    x_end = max(float(x_end), 1e-9)
    H = max(float(H), 0.0)
    r_end = max(float(r_end), 0.0)
    t = np.linspace(0.0, 1.0, N_BARREL)
    x = t * x_end
    y_line = H - x * float(np.tan(max(beta, 1e-4)))
    # Blend the straight wave into the end radius so Mach-stem cases stay smooth.
    y = (1.0 - t) * np.maximum(y_line, 0.0) + t * r_end
    y = np.maximum(y, 0.0)
    return [[float(a), float(b)] for a, b in zip(x, y)]


def _empty_shock(
    p_e_Pa: float,
    p_tank_Pa: float,
    *,
    regime: str,
    reason: str,
) -> dict[str, Any]:
    return {
        "p_tank_Pa": float(p_tank_Pa),
        "p_e_Pa": float(p_e_Pa),
        "npr": float(npr(p_e_Pa, p_tank_Pa)),
        "regime": regime,
        "x_mach_disk_m": None,
        "r_triple_m": None,
        "shock_applied": False,
        "shock_reason": reason,
        "barrel_xy": [],
        "disk_y0": None,
        "disk_y1": None,
    }


def plan_shock(
    *,
    p_e_Pa: float,
    p_tank_Pa: float,
    H: float,
    kn_exit: float,
    r_freeze_m: float | None,
    n0: float,
    d_hs: float,
    gamma: float,
    Mach_e: float,
    mode: str,
    kn_crit: float = KN_CRIT,
) -> dict[str, Any]:
    """Geometry + apply/veto decision. Does not mutate a field."""
    pe = float(p_e_Pa)
    pt = float(p_tank_Pa)
    H = max(float(H), 1e-9)
    ratio = npr(pe, pt)
    npr_regime = classify_npr_regime(pe, pt)
    g = float(np.clip(gamma, 1.05, 1.67))
    Me = float(max(Mach_e, 1.05))
    r_freeze = float(r_freeze_m) if r_freeze_m is not None and np.isfinite(r_freeze_m) else None

    # Thesis / collisionless chip: hard override. p_tank is diagnostics only —
    # never a barrel, disk, or field mutation, even at huge NPR or low Kn.
    if mode == "collisionless":
        regime = "vacuum" if npr_regime == "underexpanded" else npr_regime
        if npr_regime == "matched":
            regime = "matched"
        return _empty_shock(pe, pt, regime=regime, reason="collisionless")
    if kn_exit >= kn_crit and mode != "sudden_freeze":
        # Auto path already chose collisionless; keep the same empty overlay.
        regime = "vacuum" if npr_regime == "underexpanded" else npr_regime
        if npr_regime == "matched":
            regime = "matched"
        return _empty_shock(pe, pt, regime=regime, reason="collisionless")

    if npr_regime == "matched":
        return _empty_shock(pe, pt, regime="matched", reason="matched")

    def _veto(x_wave: float) -> bool:
        if r_freeze is not None and r_freeze < x_wave:
            return True
        kn_wave = kn_gll_at_radius(
            x_wave, n0=n0, H=H, d_hs=d_hs, gamma=g, Mach_e=Me
        )
        return kn_wave >= kn_crit

    if npr_regime == "underexpanded":
        x_m = crist_mach_disk_x(pe, pt, H)
        r_tp = TRIPLE_X_FRAC * x_m
        if _veto(x_m):
            return _empty_shock(pe, pt, regime="vacuum", reason="freeze_before_disk")
        barrel = _barrel_underexpanded(H, x_m, r_tp)
        return {
            "p_tank_Pa": pt,
            "p_e_Pa": pe,
            "npr": float(ratio),
            "regime": "underexpanded",
            "x_mach_disk_m": float(x_m),
            "r_triple_m": float(r_tp),
            "shock_applied": True,
            "shock_reason": "underexpanded",
            "barrel_xy": barrel,
            "disk_y0": float(-r_tp),
            "disk_y1": float(r_tp),
            "kind": "underexpanded_disk",
            "beta": None,
            "theta": None,
        }

    # Overexpanded: lip oblique shock from the pressure ratio.
    p2p1 = pt / max(pe, 1e-30)
    beta = _beta_from_pressure_ratio(Me, p2p1, g)
    if beta is None:
        # Detached / would-be normal: small near-exit disk (Mach reflection).
        x_m = max(1.5 * H, 0.08)
        r_tp = 0.22 * H
        mach_reflect = True
        theta = 0.0
        beta_use = 0.5 * np.pi
    else:
        post = _oblique_post(Me, beta, g)
        theta = post["theta"]
        M2 = post["M2"]
        th_max = _theta_max(max(M2, 1.01), g)
        mach_reflect = theta > th_max - 1e-3
        beta_use = beta
        if mach_reflect:
            x_m = 0.70 * H / max(float(np.tan(beta)), 1e-3)
            r_tp = 0.22 * H
        else:
            x_m = H / max(float(np.tan(beta)), 1e-3)
            r_tp = 0.0

    x_m = float(max(x_m, 1e-6))
    if _veto(x_m):
        return _empty_shock(pe, pt, regime="overexpanded", reason="freeze_before_disk")

    if mach_reflect:
        barrel = _barrel_overexpanded(H, x_m, r_tp, float(beta_use))
        return {
            "p_tank_Pa": pt,
            "p_e_Pa": pe,
            "npr": float(ratio),
            "regime": "overexpanded",
            "x_mach_disk_m": x_m,
            "r_triple_m": float(r_tp),
            "shock_applied": True,
            "shock_reason": "overexpanded_mach",
            "barrel_xy": barrel,
            "disk_y0": float(-r_tp),
            "disk_y1": float(r_tp),
            "kind": "overexpanded_mach",
            "beta": float(beta_use),
            "theta": float(theta),
        }

    barrel = _barrel_overexpanded(H, x_m, 0.0, float(beta_use))
    return {
        "p_tank_Pa": pt,
        "p_e_Pa": pe,
        "npr": float(ratio),
        "regime": "overexpanded",
        "x_mach_disk_m": None,
        "r_triple_m": None,
        "shock_applied": True,
        "shock_reason": "overexpanded_regular",
        "barrel_xy": barrel,
        "disk_y0": None,
        "disk_y1": None,
        "kind": "overexpanded_regular",
        "beta": float(beta_use),
        "theta": float(theta),
        "x_reflect_m": x_m,
    }


def _public_meta(plan: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "p_tank_Pa": plan["p_tank_Pa"],
        "p_e_Pa": plan["p_e_Pa"],
        "npr": plan["npr"],
        "regime": plan["regime"],
        "x_mach_disk_m": plan["x_mach_disk_m"],
        "r_triple_m": plan["r_triple_m"],
        "shock_applied": bool(plan["shock_applied"]),
        "shock_reason": plan["shock_reason"],
        "barrel_xy": list(plan.get("barrel_xy") or []),
        "disk_y0": plan["disk_y0"],
        "disk_y1": plan["disk_y1"],
    }


def apply_shock_to_field(
    field: dict[str, Any],
    x: np.ndarray,
    y: np.ndarray,
    plan: Mapping[str, Any],
    *,
    T0: float,
    n0: float,
    U0: float,
    H: float,
    gamma: float,
    Mach_e: float,
    R_specific: float,
) -> dict[str, Any]:
    """Rankine–Hugoniot disk and/or inward deflection. No-op if not applied."""
    if not plan.get("shock_applied"):
        return field
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    X, Y = np.meshgrid(x, y)
    n_ratio = np.asarray(field["n_ratio"], dtype=np.float64).copy()
    t_ratio = np.asarray(field["t_ratio"], dtype=np.float64).copy()
    u = np.asarray(field["u"], dtype=np.float64).copy()
    v = np.asarray(field["v"], dtype=np.float64).copy()
    g = float(np.clip(gamma, 1.05, 1.67))
    n0 = max(float(n0), 1.0)
    T0 = max(float(T0), 1.0)

    kind = plan.get("kind") or ""
    if kind == "underexpanded_disk" and plan.get("x_mach_disk_m") is not None:
        x_m = float(plan["x_mach_disk_m"])
        r_tp = float(plan["r_triple_m"] or 0.0)
        core = isentropic_core_at_R(
            x_m, T0=T0, n0=n0, U0=U0, H=H, gamma=g, Mach_e=Mach_e, R_specific=R_specific
        )
        rh = rankine_hugoniot_normal(core["M"], g)
        n2 = core["n"] * rh["rho_ratio"]
        T2 = core["T"] * rh["T_ratio"]
        U2 = core["U"] * rh["U_ratio"]
        wake = (X >= x_m) & (np.abs(Y) <= r_tp)
        n_ratio = np.where(wake, n2 / n0, n_ratio)
        t_ratio = np.where(wake, T2 / T0, t_ratio)
        u = np.where(wake, U2, u)
        v = np.where(wake, 0.0, v)
    elif kind.startswith("overexpanded"):
        theta = float(plan.get("theta") or 0.0)
        beta = float(plan.get("beta") or (0.4))
        x_end = float(plan.get("x_mach_disk_m") or plan.get("x_reflect_m") or 0.0)
        tanb = float(np.tan(max(beta, 1e-4)))
        y_shock = np.maximum(H - np.maximum(X, 0.0) * tanb, 0.0)
        inside = (X > 0.02 * H) & (X <= max(x_end, 0.0) + 1e-9) & (np.abs(Y) <= y_shock + 1e-12)
        # Post-oblique (or post-normal if detached) from the exit state.
        if kind == "overexpanded_mach" and (plan.get("beta") is None or beta > 1.2):
            rh = rankine_hugoniot_normal(max(float(Mach_e), 1.05), g)
            n2 = n0 * rh["rho_ratio"]
            T2 = T0 * rh["T_ratio"]
            U2 = float(U0) * rh["U_ratio"]
            n_ratio = np.where(inside, n2 / n0, n_ratio)
            t_ratio = np.where(inside, T2 / T0, t_ratio)
            u = np.where(inside, U2, u)
            v = np.where(inside, 0.0, v)
        else:
            post = _oblique_post(max(float(Mach_e), 1.05), beta, g)
            n2 = n0 * post["rho_ratio"]
            T2 = T0 * post["T_ratio"]
            U2 = float(U0) * post["U_ratio"]
            cth, sth = float(np.cos(theta)), float(np.sin(theta))
            u_in = U2 * cth
            v_in = -np.sign(Y) * U2 * sth
            n_ratio = np.where(inside, n2 / n0, n_ratio)
            t_ratio = np.where(inside, T2 / T0, t_ratio)
            u = np.where(inside, u_in, u)
            v = np.where(inside, v_in, v)
        if kind == "overexpanded_mach" and plan.get("x_mach_disk_m") is not None:
            x_m = float(plan["x_mach_disk_m"])
            r_tp = float(plan["r_triple_m"] or 0.0)
            core = isentropic_core_at_R(
                max(x_m, H), T0=T0, n0=n0, U0=U0, H=H, gamma=g, Mach_e=Mach_e, R_specific=R_specific
            )
            rh = rankine_hugoniot_normal(core["M"], g)
            n2 = core["n"] * rh["rho_ratio"]
            T2 = core["T"] * rh["T_ratio"]
            U2 = core["U"] * rh["U_ratio"]
            wake = (X >= x_m) & (np.abs(Y) <= r_tp)
            n_ratio = np.where(wake, n2 / n0, n_ratio)
            t_ratio = np.where(wake, T2 / T0, t_ratio)
            u = np.where(wake, U2, u)
            v = np.where(wake, 0.0, v)

    n_ratio = np.clip(n_ratio, 0.0, 8.0)
    t_ratio = np.clip(t_ratio, 0.0, 8.0)
    speed = np.hypot(u, v)
    out = dict(field)
    out["n_ratio"] = n_ratio
    out["t_ratio"] = t_ratio
    out["u"] = u
    out["v"] = v
    out["speed"] = speed
    return out


def attach_shock_overlay(
    field: dict[str, Any],
    x: np.ndarray,
    y: np.ndarray,
    *,
    p_e_Pa: float,
    p_tank_Pa: float,
    H: float,
    T0: float,
    n0: float,
    U0: float,
    R_specific: float,
    gamma: float,
    Mach_e: float,
    d_hs: float,
    kn_exit: float,
    r_freeze_m: float | None,
    mode: str,
    kn_crit: float = KN_CRIT,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Plan + optional field overlay. Always returns public plume extras."""
    plan = plan_shock(
        p_e_Pa=p_e_Pa,
        p_tank_Pa=p_tank_Pa,
        H=H,
        kn_exit=kn_exit,
        r_freeze_m=r_freeze_m,
        n0=n0,
        d_hs=d_hs,
        gamma=gamma,
        Mach_e=Mach_e,
        mode=mode,
        kn_crit=kn_crit,
    )
    if plan.get("shock_applied") and mode != "collisionless":
        field = apply_shock_to_field(
            field, x, y, plan,
            T0=T0, n0=n0, U0=U0, H=H, gamma=gamma,
            Mach_e=Mach_e, R_specific=R_specific,
        )
    return field, _public_meta(plan)
