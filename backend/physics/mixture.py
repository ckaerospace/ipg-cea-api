"""Working-gas mixtures for assigned-enthalpy CEA (no .inp files)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

import numpy as np

import cea

from .constants import R_UNIV

COMMON_CEA = ("O2", "N2", "CO2", "He", "Ar", "H2", "CO", "NO", "Ne", "Kr", "Xe", "NH3", "CH4", "H2O", "N2O")


class UnknownSpecies(ValueError):
    pass


def normalize(fracs: Mapping[str, float], drop_zero: bool = True) -> dict[str, float]:
    out = {str(k): float(v) for k, v in fracs.items() if float(v) > 0 or not drop_zero}
    out = {k: v for k, v in out.items() if v > 1e-12}
    s = sum(out.values())
    if s <= 0:
        raise ValueError("mixture is empty — add at least one working gas")
    return {k: v / s for k, v in out.items()}


def expand_air_recipe(fracs: Mapping[str, float]) -> dict[str, float]:
    """Replace a chip named Air with N2/O2 79/21 by mole (then re-normalize later)."""
    out = dict(fracs)
    air = out.pop("Air", None)
    if air and float(air) > 0:
        out["N2"] = out.get("N2", 0.0) + 0.79 * float(air)
        out["O2"] = out.get("O2", 0.0) + 0.21 * float(air)
    return out


def validate_species(name: str) -> str:
    name = name.strip()
    if not name:
        raise UnknownSpecies("empty species name")
    cea.init()
    try:
        cea.Mixture([name])
    except Exception as exc:
        raise UnknownSpecies(f"'{name}' is not a CEA thermo species") from exc
    return name


def mole_to_mass(names: list[str], moles: np.ndarray) -> np.ndarray:
    cea.init()
    reac = cea.Mixture(names)
    w = np.asarray(reac.moles_to_weights(moles), dtype=np.float64)
    s = w.sum()
    return w / s if s > 0 else w


def mass_to_mole(names: list[str], mass: np.ndarray) -> np.ndarray:
    cea.init()
    reac = cea.Mixture(names)
    m = np.asarray(reac.weights_to_moles(mass), dtype=np.float64)
    s = m.sum()
    return m / s if s > 0 else m


@dataclass
class MixtureSpec:
    names: list[str]
    mole_fracs: np.ndarray
    mass_fracs: np.ndarray
    MW: float  # g/mol
    R: float  # J/(kg·K)
    h_ref_J_kg: float  # CEA enthalpy at 298.15 K
    basis: str

    @property
    def h_ref_MJ_kg(self) -> float:
        return self.h_ref_J_kg / 1e6

    def as_dict(self) -> dict:
        return {
            "basis": self.basis,
            "mole_fractions": {n: float(x) for n, x in zip(self.names, self.mole_fracs)},
            "mass_fractions": {n: float(x) for n, x in zip(self.names, self.mass_fracs)},
            "MW": self.MW,
            "R": self.R,
            "h_ref_MJ_kg": self.h_ref_MJ_kg,
            "h_ref_kJ_kg": self.h_ref_J_kg / 1000.0,
        }


def parse_mixture(
    mixture: Mapping[str, float] | Iterable[Mapping[str, float]] | None,
    basis: str = "mole",
    gas: str | None = None,
    he_mole_frac: float = 0.0,
) -> MixtureSpec:
    """Build a MixtureSpec from UI chips, a dict, or the legacy gas= field."""
    basis = (basis or "mole").lower()
    if basis not in ("mole", "mass"):
        raise ValueError("basis must be 'mole' or 'mass'")

    fracs: dict[str, float]
    if mixture:
        if isinstance(mixture, Mapping):
            fracs = {str(k): float(v) for k, v in mixture.items()}
        else:
            fracs = {}
            for item in mixture:
                name = str(item.get("name") or item.get("id") or "")
                fracs[name] = fracs.get(name, 0.0) + float(item.get("fraction", 0.0))
    elif gas:
        fracs = _legacy_gas(gas, he_mole_frac)
        basis = "mole"
    else:
        fracs = {"O2": 1.0}
        basis = "mole"

    fracs = expand_air_recipe(fracs)
    fracs = normalize(fracs)
    names = [validate_species(n) for n in fracs]
    x = np.array([fracs[n] for n in names], dtype=np.float64)

    cea.init()
    reac = cea.Mixture(names)
    if basis == "mole":
        mole = x / x.sum()
        mass = mole_to_mass(names, mole)
    else:
        mass = x / x.sum()
        mole = mass_to_mole(names, mass)

    mw_vec = np.asarray(reac.moles_to_weights(np.ones(len(names))), dtype=np.float64)
    MW = float(np.dot(mole, mw_vec))  # g/mol
    R = R_UNIV / MW
    h_ref = float(reac.calc_property(cea.ENTHALPY, mass, 298.15))
    return MixtureSpec(
        names=names,
        mole_fracs=mole,
        mass_fracs=mass,
        MW=MW,
        R=R,
        h_ref_J_kg=h_ref,
        basis=basis,
    )


def _legacy_gas(gas: str, he_mole_frac: float) -> dict[str, float]:
    g = gas.upper().replace("+", "").replace("/", "").replace("-", "").replace(" ", "")
    if g in ("O2", "OXYGEN"):
        return {"O2": 1.0}
    if g in ("N2", "NITROGEN"):
        return {"N2": 1.0}
    if g in ("CO2",):
        return {"CO2": 1.0}
    if g in ("HE", "HELIUM"):
        return {"He": 1.0}
    if g in ("AR", "ARGON"):
        return {"Ar": 1.0}
    if g in ("AIR",):
        return {"N2": 0.79, "O2": 0.21}
    if g in ("HEO2", "HEOXYGEN"):
        x = float(np.clip(he_mole_frac, 0.0, 0.99))
        return {"He": x, "O2": 1.0 - x}
    raise ValueError(f"Unsupported legacy gas '{gas}'")
