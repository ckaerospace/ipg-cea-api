"""Geometry, working-gas chips, and thesis operating-point presets."""

from __future__ import annotations

from dataclasses import asdict, dataclass

# CODATA / SI
R_UNIV = 8314.51  # J/(kmol·K), matches NASA CEA
N_A = 6.02214076e23
K_B = 1.380649e-23
E_CHARGE = 1.602176634e-19
M_U = 1.66053906660e-27
M_O = 15.999 * M_U
M_N = 14.007 * M_U
M_C = 12.011 * M_U
M_HE = 4.002602 * M_U
M_AR = 39.948 * M_U

O2_DISSOCIATION_EV = 5.115
O_RECOMBINATION_EV = O2_DISSOCIATION_EV / 2.0

D_HS = {
    "O": 3.0e-10,
    "O2": 3.6e-10,
    "N": 3.0e-10,
    "N2": 3.7e-10,
    "He": 2.2e-10,
    "Ar": 3.6e-10,
    "C": 3.1e-10,
    "CO": 3.6e-10,
    "CO2": 3.9e-10,
    "NO": 3.5e-10,
    "H2": 2.9e-10,
    "H": 2.5e-10,
}

ATOM_MASS_U = {"O": 15.999, "N": 14.007, "C": 12.011, "He": 4.0026, "Ar": 39.948, "H": 1.008}

ESA_MLI_FLUENCE_CM2 = 5.0e21
VENUS_AO_EV = 8.3

PI = 3.141592653589793


@dataclass(frozen=True)
class NozzleGeometry:
    """De Laval (or tube) diameters in metres."""

    name: str
    d_c: float
    d_t: float
    d_e: float

    @property
    def a_c(self) -> float:
        return 0.25 * PI * self.d_c**2

    @property
    def a_t(self) -> float:
        return 0.25 * PI * self.d_t**2

    @property
    def a_e(self) -> float:
        return 0.25 * PI * self.d_e**2

    @property
    def subar(self) -> float:
        return (self.d_c / self.d_t) ** 2

    @property
    def supar(self) -> float:
        return (self.d_e / self.d_t) ** 2

    @property
    def H(self) -> float:
        """Semi-slit height / exit radius for the 2-D collisionless model."""
        return 0.5 * self.d_e

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "d_c_mm": self.d_c * 1e3,
            "d_t_mm": self.d_t * 1e3,
            "d_e_mm": self.d_e * 1e3,
            "a_c_mm2": self.a_c * 1e6,
            "a_t_mm2": self.a_t * 1e6,
            "a_e_mm2": self.a_e * 1e6,
            "H_mm": self.H * 1e3,
            "subar": self.subar,
            "supar": self.supar,
        }

    @classmethod
    def from_mm(
        cls,
        name: str,
        d_c_mm: float,
        d_t_mm: float,
        d_e_mm: float,
    ) -> "NozzleGeometry":
        if d_c_mm <= 0 or d_t_mm <= 0 or d_e_mm <= 0:
            raise ValueError("nozzle diameters must be positive")
        return cls(name=name or "custom", d_c=d_c_mm * 1e-3, d_t=d_t_mm * 1e-3, d_e=d_e_mm * 1e-3)

    @classmethod
    def from_areas_mm2(
        cls,
        name: str,
        a_c_mm2: float,
        a_t_mm2: float,
        a_e_mm2: float,
    ) -> "NozzleGeometry":
        def d_mm(a: float) -> float:
            return 2.0 * (a / PI) ** 0.5

        return cls.from_mm(name, d_mm(a_c_mm2), d_mm(a_t_mm2), d_mm(a_e_mm2))


IPG6S = NozzleGeometry(name="IPG6-S", d_c=37e-3, d_t=20e-3, d_e=40e-3)
IPG3 = NozzleGeometry(name="IPG3", d_c=84e-3, d_t=84e-3, d_e=84e-3)
IPG4 = NozzleGeometry(name="IPG4", d_c=84e-3, d_t=50e-3, d_e=50e-3)

FACILITIES = {
    "IPG6-S": {
        "id": "IPG6-S",
        "label": "IPG6-S",
        "blurb": "Inductive generator, de Laval. Default for this app.",
        "geometry": IPG6S.as_dict(),
        "default_gas": "O2",
        "typical_mixture": {"O2": 1.0},
        "editable": False,
    },
    "IPG3": {
        "id": "IPG3",
        "label": "IPG3",
        "blurb": "Thesis Table 4.1 — constant-area tube, typically O2.",
        "geometry": IPG3.as_dict(),
        "default_gas": "O2",
        "typical_mixture": {"O2": 1.0},
        "editable": False,
    },
    "IPG4": {
        "id": "IPG4",
        "label": "IPG4",
        "blurb": "Convergent nozzle, typically CO2 (Burghaus / Venus).",
        "geometry": IPG4.as_dict(),
        "default_gas": "CO2",
        "typical_mixture": {"CO2": 1.0},
        "editable": False,
    },
    "custom": {
        "id": "custom",
        "label": "Custom",
        "blurb": "Edit chamber, throat, and exit diameters (or areas).",
        "geometry": IPG6S.as_dict() | {"name": "custom"},
        "default_gas": "O2",
        "typical_mixture": {"O2": 1.0},
        "editable": True,
    },
}

# Lab chips. Air is a recipe (N2/O2 79/21 mole), not the CEA "Air" reactant.
COMMON_GASES = [
    {"id": "O2", "label": "O2", "cea": "O2", "recipe": None},
    {"id": "N2", "label": "N2", "cea": "N2", "recipe": None},
    {"id": "CO2", "label": "CO2", "cea": "CO2", "recipe": None},
    {"id": "He", "label": "He", "cea": "He", "recipe": None},
    {"id": "Ar", "label": "Ar", "cea": "Ar", "recipe": None},
    {
        "id": "Air",
        "label": "Air",
        "cea": None,
        "recipe": {"N2": 0.79, "O2": 0.21},
        "note": "N2/O2 79/21 by mole — editable after applying.",
    },
]

ADVANCED_SPECIES = ["H2", "CO", "NO", "Ne", "Kr", "Xe", "NH3", "CH4", "H2O", "N2O"]

TANK7 = {"name": "Tank 7", "diameter_m": 2.0, "length_m": 4.8}

THESIS_PRESETS = [
    {
        "id": "ipg6s-1",
        "label": "IPG6-S (1) · 15 MJ/kg",
        "facility": "IPG6-S",
        "pinj_Pa": 100.0,
        "hinj_MJ_kg": 15.0,
        "basis": "mole",
        "mixture": {"O2": 1.0},
    },
    {
        "id": "ipg6s-2",
        "label": "IPG6-S (2) · 23 MJ/kg",
        "facility": "IPG6-S",
        "pinj_Pa": 100.0,
        "hinj_MJ_kg": 23.0,
        "basis": "mole",
        "mixture": {"O2": 1.0},
    },
    {
        "id": "ipg6s-3",
        "label": "IPG6-S (3) · 30 MJ/kg",
        "facility": "IPG6-S",
        "pinj_Pa": 100.0,
        "hinj_MJ_kg": 30.0,
        "basis": "mole",
        "mixture": {"O2": 1.0},
    },
    {
        "id": "ipg6s-heo2",
        "label": "IPG6-S He/O2 70/30",
        "facility": "IPG6-S",
        "pinj_Pa": 100.0,
        "hinj_MJ_kg": 23.0,
        "basis": "mole",
        "mixture": {"He": 0.7, "O2": 0.3},
    },
    {
        "id": "ipg3-o01",
        "label": "IPG3 O#01",
        "facility": "IPG3",
        "pinj_Pa": 1450.0,
        "hinj_MJ_kg": 17.9,
        "basis": "mole",
        "mixture": {"O2": 1.0},
        "note": "Table 4.2 O#01 · ṁ ≈ 3210 mg/s (measured). CEA ṁ is equilibrium-choked.",
    },
    {
        "id": "ipg4-burghaus",
        "label": "IPG4 Burghaus",
        "facility": "IPG4",
        "pinj_Pa": 2900.0,
        "delta_h_MJ_kg": 26.3,
        "basis": "mole",
        "mixture": {"CO2": 1.0},
        "note": "Δh ≈ 26.3 MJ/kg above CO2 at 20 °C (h_ref ≈ −8.9 MJ/kg). ṁ ≈ 2200 mg/s measured.",
    },
]
