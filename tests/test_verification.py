"""Verification suite: Thesis/Advanced modes, Kn trigger, shocks, HTTP default.

Exit state is fixtured (no live CEA). One smoke test runs if CEA is installed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from physics.pipeline import assemble_plume
from physics.plume import CollisionlessPlume
from physics.shocks import CRIST_XM_COEF
from physics.sudden_freeze import KN_CRIT, kn_gll_exit, mix_d_hs, resolve_plume_mode

from conftest import CEA_AVAILABLE

_GRID = dict(xmax_m=0.8, ymax_m=0.3, nx=17, ny=17, include_contours=False)


def _kn(ex: dict) -> float:
    return kn_gll_exit(ex["n0"], ex["H"], mix_d_hs(ex.get("mole_fractions") or {}))


def test_thesis_collisionless_never_applies_shock(dense_exit):
    """Thesis path: collisionless hard override — no shock, no barrel, kernel field."""
    p_tank = 0.08
    payload = assemble_plume(
        dense_exit, plume_mode="collisionless", p_tank_Pa=p_tank, **_GRID
    )
    pl = payload["plume"]
    assert pl["mode"] == "collisionless"
    assert pl["shock_applied"] is False
    assert pl["barrel_xy"] == []
    assert pl["x_mach_disk_m"] is None
    assert pl["p_e_Pa"] == dense_exit["p_Pa"]
    assert pl["p_tank_Pa"] == p_tank
    assert pl["npr"] == pytest.approx(dense_exit["p_Pa"] / p_tank)

    kernel = CollisionlessPlume.from_exit(
        T0=dense_exit["T0"],
        R_specific=dense_exit["R"],
        U0=dense_exit["U0"],
        n0=dense_exit["n0"],
        H=dense_exit["H"],
    )
    x = np.asarray(pl["x"])
    y = np.asarray(pl["y"])
    ref = kernel.grid(x, y)
    ny, nx = pl["ny"], pl["nx"]
    n = np.asarray(pl["n_ratio"], dtype=np.float64).reshape(ny, nx)
    u = np.asarray(pl["u"], dtype=np.float64).reshape(ny, nx)
    np.testing.assert_allclose(n, ref["n_ratio"], rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(u, ref["u"], rtol=1e-5, atol=1e-6)


def test_auto_high_kn_is_collisionless(rarefied_exit):
    kn = _kn(rarefied_exit)
    assert kn >= KN_CRIT
    requested, chosen = resolve_plume_mode("auto", kn)
    assert requested == "auto"
    assert chosen == "collisionless"
    payload = assemble_plume(rarefied_exit, plume_mode="auto", p_tank_Pa=0.08, **_GRID)
    assert payload["plume"]["mode"] == "collisionless"
    assert payload["plume"]["shock_applied"] is False
    assert payload["plume"]["barrel_xy"] == []


def test_auto_low_kn_is_sudden_freeze(dense_exit):
    kn = _kn(dense_exit)
    assert kn < KN_CRIT
    _, chosen = resolve_plume_mode("auto", kn)
    assert chosen == "sudden_freeze"
    payload = assemble_plume(dense_exit, plume_mode="auto", p_tank_Pa=25.0, **_GRID)
    assert payload["plume"]["mode"] == "sudden_freeze"


def test_freeze_before_disk_tiny_p_tank(near_freeze_exit):
    kn = _kn(near_freeze_exit)
    assert kn < KN_CRIT
    payload = assemble_plume(
        near_freeze_exit, plume_mode="auto", p_tank_Pa=0.08, **_GRID
    )
    pl = payload["plume"]
    assert pl["shock_applied"] is False
    assert pl["shock_reason"] == "freeze_before_disk"
    assert pl["x_mach_disk_m"] is None
    assert pl["barrel_xy"] == []


def test_disk_moderate_npr_low_kn(dense_exit):
    p_e = float(dense_exit["p_Pa"])
    p_tank = 25.0
    de = 2.0 * float(dense_exit["H"])
    expected = CRIST_XM_COEF * de * math.sqrt(p_e / p_tank)
    payload = assemble_plume(
        dense_exit, plume_mode="auto", p_tank_Pa=p_tank, **_GRID
    )
    pl = payload["plume"]
    assert pl["mode"] == "sudden_freeze"
    assert pl["shock_applied"] is True
    assert pl["regime"] == "underexpanded"
    assert pl["x_mach_disk_m"] is not None
    assert abs(pl["x_mach_disk_m"] - expected) / expected <= 0.15


def test_overexpanded_when_p_tank_gt_p_e(dense_exit):
    ex = dict(dense_exit)
    ex["p_Pa"] = 10.0
    payload = assemble_plume(ex, plume_mode="auto", p_tank_Pa=50.0, **_GRID)
    assert payload["plume"]["regime"] == "overexpanded"
    assert payload["plume"]["npr"] == pytest.approx(10.0 / 50.0)


def test_omit_p_tank_defaults_ten_no_422():
    from pydantic import ValidationError

    from app import SolveRequest

    body = SolveRequest.model_validate(
        {"pinj_Pa": 100, "hinj_MJ_kg": 23, "mixture": {"O2": 1.0}}
    )
    assert body.p_tank_Pa == 10.0
    assert "p_tank_Pa" not in body.model_fields_set
    with pytest.raises(ValidationError):
        SolveRequest.model_validate({"p_tank_Pa": 0.05})


def test_omit_p_tank_http_no_422(monkeypatch):
    captured: dict = {}

    def _fake(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        return {"cea": {"exit": {"p_Pa": 5.0}}, "plume": {"shock_applied": False}}

    monkeypatch.setattr("app.solve_cached", _fake)
    from fastapi.testclient import TestClient

    from app import app

    client = TestClient(app)
    res = client.post(
        "/api/solve",
        json={"pinj_Pa": 100, "hinj_MJ_kg": 23, "mixture": {"O2": 1.0}, "nx": 17, "ny": 17},
    )
    assert res.status_code == 200
    assert captured["args"][-1] == 10.0
    assert captured["args"][-2] == "auto"


@pytest.mark.skipif(not CEA_AVAILABLE, reason="NASA CEA not installed")
def test_live_cea_smoke():
    from physics.pipeline import solve_operating_point

    out = solve_operating_point(
        pinj_Pa=100.0,
        hinj_MJ_kg=23.0,
        mixture={"O2": 1.0},
        xmax_m=0.5,
        ymax_m=0.3,
        nx=17,
        ny=17,
        include_contours=False,
        plume_mode="auto",
        p_tank_Pa=10.0,
    )
    pl = out["plume"]
    assert pl["kn_gll_exit"] > 0.0
    assert pl["mode"] in ("collisionless", "sudden_freeze")
    assert "n_ratio" in pl
    assert pl["p_tank_Pa"] == 10.0
    assert "p_e_Pa" in pl
