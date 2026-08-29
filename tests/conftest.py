"""Shared fixtures. CEA is stubbed when the NASA package is not installed."""

from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "backend"
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

CEA_AVAILABLE = False
try:
    import cea as _cea  # noqa: F401

    CEA_AVAILABLE = hasattr(_cea, "init") or hasattr(_cea, "RocketSolver")
except ImportError:
    _fake = types.ModuleType("cea")
    _fake.__version__ = "0.0.0"
    sys.modules["cea"] = _fake


def exit_state(
    *,
    n0: float,
    p_Pa: float,
    H: float = 0.02,
    T0: float = 2000.0,
    U0: float = 2500.0,
    R: float = 260.0,
    gamma: float = 1.4,
    Mach: float = 2.0,
    MW: float = 32.0,
) -> dict:
    """Fixture CEA exit. n0 sets Kn_exit; p_Pa is p_e for NPR."""
    return {
        "T0": T0,
        "R": R,
        "U0": U0,
        "n0": n0,
        "H": H,
        "p_Pa": p_Pa,
        "gamma": gamma,
        "Mach": Mach,
        "MW": MW,
        "h_kJ_kg": 0.0,
        "mole_fractions": {"O2": 1.0},
        "x_O": 0.0,
        "x_O2": 1.0,
        "x_He": 0.0,
    }


@pytest.fixture
def dense_exit() -> dict:
    """Kn_exit well below 0.05 so Auto picks sudden_freeze."""
    return exit_state(n0=2.0e23, p_Pa=400.0)


@pytest.fixture
def rarefied_exit() -> dict:
    """Kn_exit >= 0.05 so Auto picks collisionless."""
    return exit_state(n0=5.0e20, p_Pa=200.0)


@pytest.fixture
def near_freeze_exit() -> dict:
    """Continuum at the lip but freeze close-in (tiny ambient → freeze_before_disk)."""
    return exit_state(n0=2.0e21, p_Pa=50.0)
