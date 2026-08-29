"""Optional calorimeter disk on the plume axis.

Incident state is sampled from the existing plume field at (x, y=0)
*before* inserting a body. Two engineering closures, switched on Kn_obj
and plume mode — not Euler/NS/DSMC.

Kinetic (collisionless mode, or Kn_obj >= 0.05):
  Local drifting-Maxwellian incident fluxes on the forward face, with
  fully diffuse re-emission at Tw (accommodation α = 1). Number flux,
  normal pressure and translational energy follow the standard
  gaskinetic forward-face integrals used by Khasawneh, Liu & Cai
  (Phys. Fluids 22, 117101, 2010, doi:10.1063/1.3490409) and the
  collisionless jet of Cai & Boyd (J. Spacecraft 44(3) 2007,
  doi:10.2514/1.25893). Translational energy only; no catalysis.

Continuum (sudden_freeze / auto-continuum and Kn_obj < 0.05):
  Billig (1967) sphere-like bow standoff on a blunt face,
  modified-Newtonian pressure, Sutton–Graves stagnation heating for a
  sphere of radius R. Not a 2-D NS body.

Collisionless geometric shadow downstream of the disk is omitted in v1
(the plume arrays are the undisturbed field).
"""

from __future__ import annotations

from typing import Any, Mapping

from math import erf

import numpy as np

from .constants import K_B
from .sudden_freeze import KN_CRIT, mix_d_hs

SQRT_PI = np.sqrt(np.pi)
ATM_PA = 101325.0


def _finite(x: float, default: float = 0.0) -> float:
    v = float(x)
    return v if np.isfinite(v) else default

# Sutton & Graves, NASA TR R-376 (1971), Table II.
# K in kg s^{-1} m^{-3/2} atm^{-1/2} for q = K (p_s/R_n)^{1/2} (h_s − h_w).
_SG_K = {
    "N2": 0.1112,
    "O2": 0.1201,
    "H2": 0.0395,
    "He": 0.0797,
    "Ne": 0.1474,
    "Ar": 0.1495,
    "CO2": 0.1210,
    "NH3": 0.0990,
    "CH4": 0.0807,
}

# Undissociated parent used to look up K (TR R-376 is for the cold mixture).
_SG_PARENT = {
    "O": "O2",
    "O2": "O2",
    "O+": "O2",
    "O2+": "O2",
    "N": "N2",
    "N2": "N2",
    "N+": "N2",
    "C": "CO2",
    "CO": "CO2",
    "CO2": "CO2",
    "C+": "CO2",
    "He": "He",
    "He+": "He",
    "Ar": "Ar",
    "Ar+": "Ar",
    "H": "H2",
    "H2": "H2",
    "Ne": "Ne",
    "CH4": "CH4",
    "NH3": "NH3",
}


def mean_free_path(n: float, d_hs: float) -> float:
    n = max(float(n), 0.0)
    d_hs = max(float(d_hs), 1e-12)
    if n <= 0.0:
        return float("inf")
    return float(1.0 / (np.sqrt(2.0) * np.pi * d_hs * d_hs * n))


def kn_obj(n: float, r_m: float, d_hs: float) -> float:
    lam = mean_free_path(n, d_hs)
    r_m = max(float(r_m), 1e-9)
    if not np.isfinite(lam):
        return float("inf")
    return float(lam / r_m)


def _chi_psi(s: float) -> tuple[float, float, float]:
    """Forward-face speed-ratio auxiliaries. s = U_n / sqrt(2 R T)."""
    s = float(s)
    e = float(np.exp(-s * s))
    er = float(1.0 + erf(s))
    chi = e + SQRT_PI * s * er
    # Normal-momentum coefficient: p_i / (n k T)
    p_hat = (s / SQRT_PI) * e + (s * s + 0.5) * er
    # Energy auxiliary; ψ(0) = 2 so the 1/(4√π) prefactor recovers 2 kT per particle.
    psi = (s * s + 2.0) * e + SQRT_PI * s * (s * s + 2.5) * er
    return chi, p_hat, psi


def kinetic_face(
    n: float,
    T: float,
    U_n: float,
    R_specific: float,
    Tw: float,
) -> tuple[float, float]:
    """Incident + diffuse re-emitted pressure and translational heat (α = 1).

    Returns (p_w, q_w) with q_w the heat *to* the wall (incident energy minus
    wall re-emission). Empty/undefined states return (0, 0).
    """
    n = max(float(n), 0.0)
    T = max(float(T), 0.0)
    Tw = max(float(Tw), 1.0)
    R_specific = max(float(R_specific), 1.0)
    U_n = max(float(U_n), 0.0)
    if n <= 0.0 or T <= 1e-6:
        return 0.0, 0.0

    m = K_B / R_specific
    c_m = np.sqrt(2.0 * R_specific * T)
    s = U_n / c_m if c_m > 0.0 else 0.0
    chi, p_hat, psi = _chi_psi(s)

    p_i = n * K_B * T * p_hat
    p_r = 0.5 * n * K_B * chi * np.sqrt(T * Tw)
    p_w = float(p_i + p_r)

    # E_i = n m (2RT)^{3/2} ψ / (4 √π)   (3-D translational)
    E_i = n * m * (c_m ** 3) * psi / (4.0 * SQRT_PI)
    N_i = n * np.sqrt(R_specific * T / (2.0 * np.pi)) * chi
    E_r = N_i * 2.0 * K_B * Tw
    q_w = float(E_i - E_r)
    return p_w, q_w


def rayleigh_pitot(p: float, M: float, gamma: float) -> float:
    """Total pressure behind a normal shock (M>1) or isentropic pitot (M<=1)."""
    p = max(float(p), 0.0)
    g = float(np.clip(gamma, 1.05, 1.67))
    M = max(float(M), 0.0)
    if p <= 0.0:
        return 0.0
    if M <= 1.0:
        return float(p * (1.0 + 0.5 * (g - 1.0) * M * M) ** (g / (g - 1.0)))
    # Rayleigh–Pitot: p_t2/p_1
    num = ((g + 1.0) / 2.0) * M * M
    den = (2.0 * g * M * M - (g - 1.0)) / (g + 1.0)
    if den <= 0.0:
        return p
    return float(p * (num ** (g / (g - 1.0))) * (1.0 / den) ** (1.0 / (g - 1.0)))


def billig_standoff_sphere(M: float) -> float:
    """Δ/R for a spherical nose. Billig, J. Spacecraft 4(6) 1967, doi:10.2514/3.28969."""
    M = max(float(M), 1.0)
    return float(0.143 * np.exp(3.24 / (M * M)))


def billig_rc_sphere(M: float) -> float:
    """Shock vertex radius of curvature / R. Billig 1967, sphere-cone fit."""
    M = max(float(M), 1.05)
    return float(1.143 * np.exp(0.54 / (M - 1.0) ** 1.2))


def billig_bow_polyline(
    x_face: float,
    R: float,
    M: float,
    n_pts: int = 33,
    y_max_R: float = 4.0,
) -> list[list[float]]:
    """Upper-half Billig hyperbola in plume (x, y) coordinates.

    Shock vertex on the axis at x_face − Δ; x increases downstream.
    x = (x_face − Δ) + R_c cot²θ (sqrt(1 + y² tan²θ / R_c²) − 1)
    with θ = arcsin(1/M) the asymptotic Mach angle.
    """
    if M <= 1.02 or R <= 0.0:
        return []
    delta = billig_standoff_sphere(M) * R
    rc = min(billig_rc_sphere(M) * R, 80.0 * R)
    theta = float(np.arcsin(1.0 / max(M, 1.02)))
    cot = 1.0 / np.tan(theta) if theta > 1e-6 else 1e6
    tan = np.tan(theta)
    x0 = x_face - delta
    ys = np.linspace(0.0, y_max_R * R, n_pts)
    out: list[list[float]] = []
    for y in ys:
        x = x0 + rc * cot * cot * (np.sqrt(1.0 + (y * tan / rc) ** 2) - 1.0)
        out.append([float(x), float(y)])
    return out


def sutton_graves_K(mass_fractions: Mapping[str, float] | None) -> float:
    """Mixture K from TR R-376 Table II via 1/K² = Σ y_i / K_i² (eq. 44 style)."""
    y_parent: dict[str, float] = {}
    for name, y in (mass_fractions or {}).items():
        yv = float(y)
        if yv <= 0.0:
            continue
        parent = _SG_PARENT.get(name) or _SG_PARENT.get(name.rstrip("+-"))
        if parent is None or parent not in _SG_K:
            continue
        y_parent[parent] = y_parent.get(parent, 0.0) + yv
    s = sum(y_parent.values())
    if s <= 0.0:
        return _SG_K["O2"]
    acc = 0.0
    for parent, y in y_parent.items():
        ki = _SG_K[parent]
        acc += (y / s) / (ki * ki)
    return float(acc ** (-0.5)) if acc > 0.0 else _SG_K["O2"]


def sutton_graves_q(
    p_stag_Pa: float,
    R_n: float,
    h_s: float,
    h_w: float,
    mass_fractions: Mapping[str, float] | None,
) -> float:
    """Stagnation convective heat, Sutton & Graves NASA TR R-376 eq. (33).

    q = K (p_s / R_n)^{1/2} (h_s − h_w), p_s in atm, h in J/kg, q in W/m².
    Cold-wall: clamp (h_s − h_w) at 0. Sphere of radius R_n.
    """
    R_n = max(float(R_n), 1e-6)
    p_atm = max(float(p_stag_Pa), 0.0) / ATM_PA
    dh = max(float(h_s) - float(h_w), 0.0)
    if p_atm <= 0.0 or dh <= 0.0:
        return 0.0
    K = sutton_graves_K(mass_fractions)
    return float(K * np.sqrt(p_atm / R_n) * dh)


def _mass_fracs_from_moles(mole_fractions: Mapping[str, float] | None) -> dict[str, float]:
    """Rough mass fractions from mole fractions using integer atomic masses."""
    amu = {
        "O": 16.0, "O2": 32.0, "O+": 16.0, "O2+": 32.0,
        "N": 14.0, "N2": 28.0, "N+": 14.0,
        "C": 12.0, "CO": 28.0, "CO2": 44.0, "C+": 12.0,
        "He": 4.0, "He+": 4.0, "Ar": 40.0, "Ar+": 40.0,
        "H": 1.0, "H2": 2.0, "Ne": 20.0, "CH4": 16.0, "NH3": 17.0,
        "NO": 30.0, "e-": 5.5e-4,
    }
    num: dict[str, float] = {}
    den = 0.0
    for k, x in (mole_fractions or {}).items():
        xv = float(x)
        if xv <= 0.0:
            continue
        m = amu.get(k)
        if m is None:
            continue
        w = xv * m
        num[k] = num.get(k, 0.0) + w
        den += w
    if den <= 0.0:
        return {}
    return {k: v / den for k, v in num.items()}


def evaluate_sample_disk(
    *,
    x_m: float,
    r_m: float,
    Tw_K: float,
    n_inf: float,
    T_inf: float,
    U_inf: float,
    R_specific: float,
    gamma: float,
    plume_mode: str,
    mole_fractions: Mapping[str, float] | None = None,
    mass_fractions: Mapping[str, float] | None = None,
    h_tot_J_kg: float | None = None,
) -> dict[str, Any]:
    """Face-average probe reading. p_stag/q_stag are continuum-only."""
    x_m = float(x_m)
    r_m = max(float(r_m), 1e-6)
    Tw_K = float(Tw_K)
    n_inf = max(float(n_inf), 0.0)
    T_inf = max(float(T_inf), 0.0)
    U_inf = max(float(U_inf), 0.0)
    R_specific = max(float(R_specific), 1.0)
    gamma = float(np.clip(gamma, 1.05, 1.67))
    d_hs = mix_d_hs(mole_fractions)
    kn = kn_obj(n_inf, r_m, d_hs)
    if not np.isfinite(kn):
        kn = 1e6
    kn = float(min(max(kn, 0.0), 1e6))

    mode = (plume_mode or "collisionless").strip().lower()
    continuum_mode = mode in ("sudden_freeze", "freeze", "collisional")
    use_continuum = continuum_mode and kn < KN_CRIT

    notes = [
        "Translational heat only; no chemistry catalysis in v1.",
        "Incident state is sampled on the undisturbed plume field (no body in the grid).",
        "T0 on the plume is nozzle-exit translational T, not chamber T.",
        "Collisionless geometric shadow downstream of the disk is omitted in v1.",
    ]

    a = np.sqrt(max(gamma * R_specific * max(T_inf, 1.0), 1.0))
    M = U_inf / a if a > 0.0 else 0.0
    p_inf = n_inf * K_B * T_inf if T_inf > 0.0 else 0.0

    if use_continuum:
        p_t2 = rayleigh_pitot(p_inf, M, gamma)
        # Flat forward face: θ = 0, so modified Newtonian is uniform on the face.
        p_stag = p_t2
        p_w = p_inf + (p_t2 - p_inf) * 1.0  # cos²θ = 1
        cp = gamma * R_specific / (gamma - 1.0)
        h_w = cp * Tw_K
        if h_tot_J_kg is not None and np.isfinite(h_tot_J_kg):
            h_s = float(h_tot_J_kg)
        else:
            h_s = cp * T_inf + 0.5 * U_inf * U_inf
        y_mass = dict(mass_fractions or {}) or _mass_fracs_from_moles(mole_fractions)
        q_stag = sutton_graves_q(p_stag, r_m, h_s, h_w, y_mass)
        # v1 face average: Lees-type cosine drop on an equivalent spherical nose,
        # frontal-area mean of q/q_s = cos θ is 2/3. Not a resolved NS face.
        q_w = (2.0 / 3.0) * q_stag
        bow = billig_bow_polyline(x_m, r_m, M) if M > 1.02 else []
        notes.append(
            "Continuum closure: Billig 1967 sphere standoff "
            "(Δ/R = 0.143 exp(3.24/M²), doi:10.2514/3.28969), "
            "modified Newtonian p = p_inf + (p_t2 − p_inf) cos²θ, "
            "Sutton–Graves NASA TR R-376 (1971) stagnation q for a sphere of radius R."
        )
        notes.append("Not Euler/NS/DSMC; stagnation + a simple face average only.")
        model = "newtonian_billig"
        regime = "continuum"
        return {
            "x_m": x_m,
            "r_mm": r_m * 1e3,
            "Tw_K": Tw_K,
            "kn_obj": _finite(kn, 1e6),
            "regime": regime,
            "n_inf": _finite(n_inf),
            "T_inf": _finite(T_inf),
            "U_inf": _finite(U_inf),
            "p_w_Pa": _finite(p_w),
            "q_w_W_m2": _finite(q_w),
            "p_stag_Pa": _finite(p_stag),
            "q_stag_W_m2": _finite(q_stag),
            "bow_xy": bow,
            "model": model,
            "notes": notes,
        }

    p_w, q_w = kinetic_face(n_inf, T_inf, U_inf, R_specific, Tw_K)
    notes.append(
        "Kinetic closure: Khasawneh–Cai / Cai–Boyd drifting-Maxwellian "
        "incident flux on the forward face with diffuse re-emission at Tw, α=1. "
        "No bow shock in this path."
    )
    return {
        "x_m": x_m,
        "r_mm": r_m * 1e3,
        "Tw_K": Tw_K,
        "kn_obj": _finite(kn, 1e6),
        "regime": "kinetic",
        "n_inf": _finite(n_inf),
        "T_inf": _finite(T_inf),
        "U_inf": _finite(U_inf),
        "p_w_Pa": _finite(p_w),
        "q_w_W_m2": _finite(q_w),
        "p_stag_Pa": None,
        "q_stag_W_m2": None,
        "bow_xy": [],
        "model": "khasawneh_diffuse",
        "notes": notes,
    }


