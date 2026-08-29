"""Sample-disk calorimeter: unit physics plus POST /api/solve wiring."""

from __future__ import annotations

import math

import pytest
from fastapi.testclient import TestClient

from physics.sample_disk import (
    KN_CRIT,
    evaluate_sample_disk,
    kinetic_face,
    rayleigh_pitot,
    sutton_graves_q,
)
from physics.pipeline import attach_sample_disk


# Dense O2-like state so Kn_obj << 0.05 for R = 80 mm.
_CONT_N = 2.0e22  # m^-3
_CONT_T = 800.0
_CONT_U = 2500.0
_CONT_R = 260.0  # J/(kg K), ~O2
_CONT_G = 1.4


def test_kinetic_equilibrium_static_gas_q_near_zero_p_equals_nkT():
    n, T, Tw = 1.0e20, 300.0, 300.0
    p_w, q_w = kinetic_face(n, T, 0.0, 260.0, Tw)
    p_eq = n * 1.380649e-23 * T
    assert p_w == pytest.approx(p_eq, rel=0.02)
    assert abs(q_w) / max(n * 1.380649e-23 * T * math.sqrt(260.0 * T), 1.0) < 0.05


def test_empty_incident_is_finite_zeros():
    out = evaluate_sample_disk(
        x_m=1.5,
        r_m=0.02,
        Tw_K=300.0,
        n_inf=0.0,
        T_inf=0.0,
        U_inf=0.0,
        R_specific=260.0,
        gamma=1.4,
        plume_mode="sudden_freeze",
        mole_fractions={"O2": 1.0},
    )
    assert out["regime"] == "kinetic"  # Kn_obj is huge
    assert math.isfinite(out["kn_obj"])
    assert out["p_w_Pa"] == 0.0
    assert out["q_w_W_m2"] == 0.0
    assert out["p_stag_Pa"] is None
    assert out["model"] == "khasawneh_diffuse"


def test_continuum_fixture_newtonian_billig_and_q_stag():
    r_m = 0.08
    out = evaluate_sample_disk(
        x_m=0.15,
        r_m=r_m,
        Tw_K=300.0,
        n_inf=_CONT_N,
        T_inf=_CONT_T,
        U_inf=_CONT_U,
        R_specific=_CONT_R,
        gamma=_CONT_G,
        plume_mode="sudden_freeze",
        mole_fractions={"O2": 1.0},
        mass_fractions={"O2": 1.0},
        h_tot_J_kg=0.5 * _CONT_U**2 + 1000.0 * _CONT_T,
    )
    assert out["kn_obj"] < KN_CRIT
    assert out["regime"] == "continuum"
    assert out["model"] == "newtonian_billig"
    assert out["q_stag_W_m2"] is not None and out["q_stag_W_m2"] > 0.0
    assert out["p_stag_Pa"] is not None and out["p_stag_Pa"] > 0.0
    assert math.isfinite(out["p_w_Pa"]) and out["p_w_Pa"] > 0.0
    assert math.isfinite(out["q_w_W_m2"]) and out["q_w_W_m2"] > 0.0
    assert out["q_w_W_m2"] == pytest.approx((2.0 / 3.0) * out["q_stag_W_m2"], rel=1e-9)
    assert isinstance(out["bow_xy"], list)
    if out["bow_xy"]:
        assert out["bow_xy"][0][1] == pytest.approx(0.0, abs=1e-12)
        ys = [p[1] for p in out["bow_xy"]]
        assert ys == sorted(ys)


def test_collisionless_mode_stays_kinetic_even_if_kn_low():
    out = evaluate_sample_disk(
        x_m=0.15,
        r_m=0.08,
        Tw_K=300.0,
        n_inf=_CONT_N,
        T_inf=_CONT_T,
        U_inf=_CONT_U,
        R_specific=_CONT_R,
        gamma=_CONT_G,
        plume_mode="collisionless",
        mole_fractions={"O2": 1.0},
    )
    assert out["kn_obj"] < KN_CRIT
    assert out["regime"] == "kinetic"
    assert out["model"] == "khasawneh_diffuse"
    assert out["p_stag_Pa"] is None
    assert out["q_stag_W_m2"] is None
    assert math.isfinite(out["p_w_Pa"]) and out["p_w_Pa"] > 0.0
    assert math.isfinite(out["q_w_W_m2"])


def test_rayleigh_pitot_subsonic_and_supersonic():
    p = 10.0
    g = 1.4
    p0_sub = rayleigh_pitot(p, 0.5, g)
    p0_super = rayleigh_pitot(p, 3.0, g)
    assert p0_sub > p
    assert p0_super > p0_sub


def test_sutton_graves_cold_wall_positive():
    q = sutton_graves_q(500.0, 0.02, 8.0e6, 3.0e5, {"O2": 1.0})
    assert q > 0.0
    q_hot = sutton_graves_q(500.0, 0.02, 3.0e5, 8.0e6, {"O2": 1.0})
    assert q_hot == 0.0


def test_attach_empty_payload_does_not_raise():
    payload = {
        "cea": {"exit": {"n0": 0.0, "T0": 0.0, "R": 260.0, "gamma": 1.4, "mole_fractions": {}}},
        "plume": {
            "mode": "sudden_freeze",
            "nx": 1,
            "ny": 1,
            "x": [0.0],
            "y": [0.0],
            "n_ratio": [0.0],
            "u": [0.0],
            "v": [0.0],
            "t_ratio": [0.0],
        },
    }
    out = attach_sample_disk(payload, 0.4, 20.0, 300.0)
    probe = out["plume"]["probe"]
    assert probe is not None
    assert math.isfinite(probe["kn_obj"])
    assert math.isfinite(probe["p_w_Pa"])
    assert math.isfinite(probe["q_w_W_m2"])


@pytest.fixture(scope="module")
def client():
    from app import app

    return TestClient(app)


def _solve_body(**extra):
    body = {
        "pinj_Pa": 100.0,
        "hinj_MJ_kg": 23.0,
        "mixture": {"O2": 1.0},
        "basis": "mole",
        "d_c_mm": 37.0,
        "d_t_mm": 20.0,
        "d_e_mm": 40.0,
        "xmax_m": 1.0,
        "ymax_m": 0.5,
        "nx": 17,
        "ny": 17,
        "plume_mode": "collisionless",
    }
    body.update(extra)
    return body


def test_api_no_probe_is_null(client):
    r = client.post("/api/solve", json=_solve_body())
    assert r.status_code == 200, r.text
    assert r.json()["plume"]["probe"] is None


def test_api_collisionless_probe_finite_kn_pw_qw(client):
    r = client.post(
        "/api/solve",
        json=_solve_body(probe_x_m=0.35, probe_r_mm=20.0, probe_Tw_K=300.0),
    )
    assert r.status_code == 200, r.text
    probe = r.json()["plume"]["probe"]
    assert probe is not None
    assert probe["model"] == "khasawneh_diffuse"
    assert probe["regime"] == "kinetic"
    assert math.isfinite(probe["kn_obj"]) and probe["kn_obj"] > 0.0
    assert math.isfinite(probe["p_w_Pa"]) and probe["p_w_Pa"] > 0.0
    assert math.isfinite(probe["q_w_W_m2"])
    assert probe["p_stag_Pa"] is None
    assert probe["q_stag_W_m2"] is None
    # bow not required in kinetic; if present it must be a list
    assert isinstance(probe.get("bow_xy"), list)


def test_api_continuum_low_kn_newtonian_billig(client):
    r = client.post(
        "/api/solve",
        json={
            "pinj_Pa": 8000.0,
            "hinj_MJ_kg": 8.0,
            "mixture": {"O2": 1.0},
            "d_c_mm": 37.0,
            "d_t_mm": 20.0,
            "d_e_mm": 40.0,
            "xmax_m": 0.6,
            "ymax_m": 0.4,
            "nx": 17,
            "ny": 17,
            "plume_mode": "sudden_freeze",
            "probe_x_m": 0.05,
            "probe_r_mm": 80.0,
            "probe_Tw_K": 300.0,
        },
    )
    assert r.status_code == 200, r.text
    probe = r.json()["plume"]["probe"]
    assert probe is not None
    assert probe["kn_obj"] < KN_CRIT
    assert probe["model"] == "newtonian_billig"
    assert probe["regime"] == "continuum"
    assert probe["q_stag_W_m2"] is not None and probe["q_stag_W_m2"] > 0.0
    assert probe["p_stag_Pa"] is not None and probe["p_stag_Pa"] > 0.0


def test_api_probe_far_downstream_does_not_500(client):
    r = client.post(
        "/api/solve",
        json=_solve_body(
            xmax_m=0.5,
            probe_x_m=0.49,
            plume_mode="collisionless",
        ),
    )
    assert r.status_code == 200, r.text
    probe = r.json()["plume"]["probe"]
    assert probe is not None
    assert math.isfinite(probe["kn_obj"])
    assert math.isfinite(probe["p_w_Pa"])
    assert math.isfinite(probe["q_w_W_m2"])


def test_api_probe_x_at_or_beyond_xmax_is_400(client):
    r = client.post("/api/solve", json=_solve_body(xmax_m=1.0, probe_x_m=1.0))
    assert r.status_code == 422
