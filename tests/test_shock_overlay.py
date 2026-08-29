"""Unit tests for the continuum shock overlay (no CEA required)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from physics.sudden_freeze import KN_CRIT, resolve_plume_mode  # noqa: E402
from physics.shocks import (  # noqa: E402
    CRIST_XM_COEF,
    attach_shock_overlay,
    crist_mach_disk_x,
    kn_gll_at_radius,
    plan_shock,
)

# Dense continuum exit so Kn at a nearby disk stays below KN_CRIT.
_CORE = dict(n0=2.0e23, d_hs=3.6e-10, gamma=1.4, Mach_e=2.0, H=0.02)
_FIELD = dict(T0=2000.0, U0=2500.0, R_specific=260.0)
_DENSE = {**_CORE, **_FIELD}


def _plan(**kwargs):
    base = dict(
        kn_exit=0.002,
        r_freeze_m=2.0,
        mode="sudden_freeze",
        kn_crit=KN_CRIT,
        **_CORE,
    )
    base.update(kwargs)
    return plan_shock(**base)


def test_tiny_p_tank_freeze_before_disk():
    """Tiny ambient pressure: x_m sits past r_freeze → no disk."""
    p_e, p_tank = 50.0, 0.08
    H = _DENSE["H"]
    r_freeze = 0.20
    x_m = crist_mach_disk_x(p_e, p_tank, H)
    assert x_m > r_freeze

    plan = _plan(p_e_Pa=p_e, p_tank_Pa=p_tank, r_freeze_m=r_freeze)
    assert plan["shock_applied"] is False
    assert plan["shock_reason"] == "freeze_before_disk"
    assert plan["x_mach_disk_m"] is None
    assert plan["barrel_xy"] == []
    assert plan["regime"] == "vacuum"


def test_moderate_npr_low_kn_mach_disk():
    """Moderate NPR + continuum Kn → Crist station and shock_applied."""
    p_e, p_tank = 400.0, 25.0
    H = _DENSE["H"]
    npr = p_e / p_tank
    expected = CRIST_XM_COEF * (2.0 * H) * math.sqrt(npr)

    kn_exit = kn_gll_at_radius(
        H, n0=_DENSE["n0"], H=H, d_hs=_DENSE["d_hs"],
        gamma=_DENSE["gamma"], Mach_e=_DENSE["Mach_e"],
    )
    assert kn_exit < KN_CRIT

    plan = _plan(p_e_Pa=p_e, p_tank_Pa=p_tank, kn_exit=kn_exit, r_freeze_m=2.0)
    assert plan["shock_applied"] is True
    assert plan["regime"] == "underexpanded"
    assert plan["x_mach_disk_m"] is not None
    assert abs(plan["x_mach_disk_m"] - expected) / expected < 0.02
    assert plan["barrel_xy"]
    assert plan["barrel_xy"][0][1] >= 0.0
    assert plan["disk_y0"] < 0.0 < plan["disk_y1"]


def test_overexpanded_regime():
    """p_tank > p_e → overexpanded."""
    plan = _plan(p_e_Pa=10.0, p_tank_Pa=50.0, kn_exit=0.001, r_freeze_m=2.0)
    assert plan["regime"] == "overexpanded"
    assert plan["npr"] == pytest.approx(10.0 / 50.0)


def test_high_kn_collisionless_no_barrel_even_if_npr_huge():
    """Kn_exit >= 0.05 → collisionless; no barrel even at huge NPR."""
    kn_exit = 0.20
    requested, chosen = resolve_plume_mode("auto", kn_exit)
    assert requested == "auto"
    assert chosen == "collisionless"
    assert kn_exit >= KN_CRIT

    plan = _plan(
        p_e_Pa=200.0,
        p_tank_Pa=0.08,
        kn_exit=kn_exit,
        r_freeze_m=5.0,
        mode=chosen,
    )
    assert plan["shock_applied"] is False
    assert plan["barrel_xy"] == []
    assert plan["x_mach_disk_m"] is None
    assert plan["shock_reason"] == "collisionless"


def test_chips_override_auto():
    assert resolve_plume_mode("collisionless", 0.001) == ("collisionless", "collisionless")
    assert resolve_plume_mode("sudden_freeze", 0.20) == ("sudden_freeze", "sudden_freeze")
    assert resolve_plume_mode("auto", 0.001)[1] == "sudden_freeze"
    assert resolve_plume_mode("", 0.20)[0] == "auto"


def test_collisionless_hard_override_ignores_p_tank_for_field():
    """Thesis sends collisionless: no barrel/disk even at low Kn and huge NPR."""
    p_e, p_tank = 400.0, 0.08
    plan = _plan(
        p_e_Pa=p_e,
        p_tank_Pa=p_tank,
        kn_exit=0.001,
        r_freeze_m=5.0,
        mode="collisionless",
    )
    assert plan["shock_applied"] is False
    assert plan["shock_reason"] == "collisionless"
    assert plan["barrel_xy"] == []
    assert plan["x_mach_disk_m"] is None
    assert plan["disk_y0"] is None and plan["disk_y1"] is None
    # Diagnostics still echo NPR if a tank pressure was sent.
    assert plan["p_e_Pa"] == p_e
    assert plan["p_tank_Pa"] == p_tank
    assert plan["npr"] == pytest.approx(p_e / p_tank)


def test_matched_no_shock():
    plan = _plan(p_e_Pa=10.0, p_tank_Pa=10.0, kn_exit=0.001)
    assert plan["regime"] == "matched"
    assert plan["shock_applied"] is False
    assert plan["shock_reason"] == "matched"


def test_attach_overlay_collisionless_leaves_field():
    import numpy as np

    x = np.linspace(0.0, 0.4, 9)
    y = np.linspace(-0.1, 0.1, 7)
    ny, nx = 7, 9
    field = {
        "n_ratio": np.ones((ny, nx)),
        "t_ratio": np.ones((ny, nx)),
        "u": np.full((ny, nx), 1000.0),
        "v": np.zeros((ny, nx)),
        "speed": np.full((ny, nx), 1000.0),
    }
    out, meta = attach_shock_overlay(
        field, x, y,
        p_e_Pa=200.0,
        p_tank_Pa=0.08,
        kn_exit=0.2,
        r_freeze_m=5.0,
        mode="collisionless",
        **_DENSE,
    )
    assert meta["shock_applied"] is False
    assert meta["barrel_xy"] == []
    assert np.allclose(out["n_ratio"], field["n_ratio"])


def test_solve_request_p_tank_default_and_bounds():
    import types

    from pydantic import ValidationError

    if "cea" not in sys.modules:
        fake = types.ModuleType("cea")
        fake.__version__ = "0.0.0"
        sys.modules["cea"] = fake

    from app import SolveRequest

    req = SolveRequest()
    assert req.p_tank_Pa == 10.0
    with pytest.raises(ValidationError):
        SolveRequest(p_tank_Pa=0.05)
    with pytest.raises(ValidationError):
        SolveRequest(p_tank_Pa=2e5)
    ok = SolveRequest(p_tank_Pa=0.051)
    assert ok.p_tank_Pa == pytest.approx(0.051)
