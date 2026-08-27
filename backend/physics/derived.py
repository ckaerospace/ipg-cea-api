"""Probe-derived quantities for whatever mixture is frozen at the nozzle exit."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .constants import (
    ATOM_MASS_U,
    D_HS,
    E_CHARGE,
    ESA_MLI_FLUENCE_CM2,
    K_B,
    M_U,
    N_A,
    O_RECOMBINATION_EV,
    VENUS_AO_EV,
)

# Species we surface on the probe when they are present at the exit.
ATOM_KEYS = ("O", "N", "C", "He", "Ar", "H")
MOLECULE_KEYS = ("O2", "N2", "CO2", "CO", "NO", "H2")
ION_KEYS = ("O+", "N+", "C+", "Ar+", "He+", "NO+", "e-")


@dataclass
class ProbeReading:
    x: float
    y: float
    n: float
    U: float
    V: float
    speed: float
    T: float
    knudsen: float
    q_kin_W_m2: float
    species: dict[str, dict]
    ao_energy_eV: float | None
    ao_flux_m2s: float | None
    mli_minutes: float | None
    q_cat_W_m2: float | None
    notes: list[str] = field(default_factory=list)


def _hs_diameter(mf: dict[str, float]) -> float:
    num = 0.0
    den = 0.0
    for name, x in mf.items():
        if x <= 0:
            continue
        d = D_HS.get(name.rstrip("+-"))
        if d is None:
            continue
        num += x * d
        den += x
    return num / den if den > 0 else 3.3e-10


def _mass_u(mf: dict[str, float]) -> float:
    """Mean particle mass in u from mole fractions (rough, atoms/molecules)."""
    amu = {
        "O": 16.0, "O2": 32.0, "O+": 16.0,
        "N": 14.0, "N2": 28.0, "N+": 14.0,
        "C": 12.0, "CO": 28.0, "CO2": 44.0, "C+": 12.0,
        "He": 4.0, "He+": 4.0, "Ar": 40.0, "Ar+": 40.0,
        "H": 1.0, "H2": 2.0, "NO": 30.0, "e-": 5.5e-4,
    }
    num = den = 0.0
    for k, x in mf.items():
        if x <= 0:
            continue
        m = amu.get(k)
        if m is None:
            continue
        num += x * m
        den += x
    return num / den if den > 0 else 16.0


def probe_quantities(
    x: float,
    y: float,
    n_ratio: float,
    u: float,
    v: float,
    t_ratio: float,
    n0: float,
    T0: float,
    mole_fractions: dict[str, float],
    H: float,
    gamma_cat: float = 1.0,
    x_O: float | None = None,
    x_O2: float | None = None,
    x_He: float | None = None,
) -> ProbeReading:
    mf = dict(mole_fractions or {})
    if x_O is not None and "O" not in mf:
        mf["O"] = x_O
    if x_O2 is not None and "O2" not in mf:
        mf["O2"] = x_O2
    if x_He is not None and "He" not in mf:
        mf["He"] = x_He

    n = max(n_ratio, 0.0) * n0
    T = max(t_ratio, 0.0) * T0
    speed = float(np.hypot(u, v))
    u_in = max(u, 0.0)

    d_hs = _hs_diameter(mf)
    lam = 1.0 / (np.sqrt(2.0) * np.pi * d_hs**2 * n) if n > 0 else float("inf")
    kn = lam / (2.0 * H) if np.isfinite(lam) else float("inf")

    m_avg = _mass_u(mf) * M_U
    q_kin = n * u_in * (0.5 * m_avg * speed * speed + 2.0 * K_B * T)

    species: dict[str, dict] = {}
    for key in ATOM_KEYS + MOLECULE_KEYS + ION_KEYS:
        xf = float(mf.get(key, 0.0))
        if xf < 5e-4 and key not in ("O", "N", "C"):
            continue
        if xf < 1e-6:
            continue
        entry: dict = {
            "x": xf,
            "n": n * xf,
            "flux_m2s": n * xf * u_in,
        }
        if key in ATOM_MASS_U:
            m = ATOM_MASS_U[key] * M_U
            entry["energy_eV"] = 0.5 * m * speed * speed / E_CHARGE if speed > 0 else 0.0
        species[key] = entry

    xO = float(mf.get("O", 0.0))
    ao_eV = ao_flux = mli_min = q_cat = None
    notes: list[str] = []
    if xO > 1e-4:
        ao_eV = species.get("O", {}).get("energy_eV")
        ao_flux = n * xO * u_in
        flux_cm2 = ao_flux / 1e4
        mli_min = (ESA_MLI_FLUENCE_CM2 / flux_cm2 / 60.0) if flux_cm2 > 0 else None
        q_cat = gamma_cat * ao_flux * (O_RECOMBINATION_EV * E_CHARGE)
        if ao_eV is not None and ao_eV < VENUS_AO_EV:
            notes.append(
                f"AO kinetic energy {ao_eV:.2f} eV is below Venus aerobraking ({VENUS_AO_EV} eV)."
            )
    if kn < 1.0:
        notes.append(f"Kn = {kn:.2f} < 1: collisionless assumption is marginal at this point.")

    return ProbeReading(
        x=x,
        y=y,
        n=n,
        U=u,
        V=v,
        speed=speed,
        T=T,
        knudsen=float(kn) if np.isfinite(kn) else 1e6,
        q_kin_W_m2=q_kin,
        species=species,
        ao_energy_eV=ao_eV,
        ao_flux_m2s=ao_flux,
        mli_minutes=mli_min,
        q_cat_W_m2=q_cat,
        notes=notes,
    )


def probe_to_dict(p: ProbeReading) -> dict:
    return {
        "x": p.x,
        "y": p.y,
        "n": p.n,
        "U": p.U,
        "V": p.V,
        "speed": p.speed,
        "T": p.T,
        "knudsen": p.knudsen,
        "q_kin_W_m2": p.q_kin_W_m2,
        "species": p.species,
        "ao_energy_eV": p.ao_energy_eV,
        "ao_flux_m2s": p.ao_flux_m2s,
        "ao_flux_cm2s": None if p.ao_flux_m2s is None else p.ao_flux_m2s / 1e4,
        "mli_minutes": p.mli_minutes,
        "q_cat_W_m2": p.q_cat_W_m2,
        "n_AO": p.species.get("O", {}).get("n", 0.0),
        "notes": p.notes,
    }
