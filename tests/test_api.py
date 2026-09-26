"""HTTP 接口测试：扫描、临界点、扫参、参考工况与 400 报错。"""

import pytest
from fastapi.testclient import TestClient

from streeter_phelps.api import app

client = TestClient(app)

VALID_PARAMS = {"k1": 0.2, "k2": 0.4, "u": 30.0, "l0": 10.0, "d0": 2.0, "csat": 10.0}


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_reference_scenario_endpoint_shows_sag():
    resp = client.get("/scenario/reference")
    assert resp.status_code == 200
    body = resp.json()
    assert body["exists"] is True
    assert body["critical_point"]["distance"] > 0.0
    assert body["critical_point"]["deficit"] > body["params"]["d0"]


def test_scan_endpoint_returns_curve_and_critical_point():
    resp = client.post("/scan", json={"params": VALID_PARAMS, "x_max": 150.0, "n_points": 61})
    assert resp.status_code == 200
    body = resp.json()
    assert body["has_critical_point"] is True
    assert len(body["points"]) == 61
    # 时间必须由流速折算：t = x / U
    for pt in body["points"]:
        assert pt["t"] == pytest.approx(pt["x"] / 30.0, rel=1e-9)
    analytic = body["analytic_critical"]
    assert analytic["distance"] == pytest.approx(70.5, abs=0.5)
    assert body["refined_critical"]["t_critical"] == pytest.approx(
        analytic["t_critical"], abs=1e-6
    )


def test_scan_endpoint_default_x_max():
    resp = client.post("/scan", json={"params": VALID_PARAMS})
    assert resp.status_code == 200
    body = resp.json()
    assert body["x_max"] > body["analytic_critical"]["distance"]


def test_critical_point_endpoint_monotonic_case():
    params = {**VALID_PARAMS, "l0": 5.0, "d0": 4.0, "k2": 0.5}
    resp = client.post("/critical-point", json=params)
    assert resp.status_code == 200
    body = resp.json()
    assert body["exists"] is False
    assert body["critical_point"] is None
    assert body["reason"] == "monotonic_reaeration"


def test_critical_point_endpoint_special_case():
    params = {**VALID_PARAMS, "k1": 0.3, "k2": 0.3, "d0": 1.0}
    resp = client.post("/critical-point", json=params)
    assert resp.status_code == 200
    body = resp.json()
    assert body["exists"] is True
    assert body["special_case"] is True
    assert body["critical_point"]["t_critical"] == pytest.approx(3.0, rel=1e-9)


def test_sweep_endpoint_reports_trend():
    resp = client.post(
        "/sweep",
        json={"base": VALID_PARAMS, "parameter": "k2", "start": 0.4, "stop": 0.8, "n": 5},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["trend"] == "decreasing"
    assert len(body["records"]) == 5
    deficits = [r["critical_deficit"] for r in body["records"]]
    assert all(b < a for a, b in zip(deficits, deficits[1:]))


@pytest.mark.parametrize(
    "field,bad",
    [
        ("k1", 0.0),
        ("k1", -0.2),
        ("k2", -0.4),
        ("u", 0.0),
        ("u", -30.0),
        ("csat", 0.0),
        ("l0", -5.0),
        ("d0", -1.0),
        ("d0", 99.0),  # 超过 Csat
    ],
)
def test_invalid_params_return_400(field, bad):
    resp = client.post("/scan", json={"params": {**VALID_PARAMS, field: bad}})
    assert resp.status_code == 400
    assert resp.json()["detail"]


def test_invalid_params_return_400_on_all_endpoints():
    bad_params = {**VALID_PARAMS, "k2": -1.0}
    assert client.post("/critical-point", json=bad_params).status_code == 400
    assert (
        client.post(
            "/sweep",
            json={"base": bad_params, "parameter": "k2", "start": 0.4, "stop": 0.8, "n": 5},
        ).status_code
        == 400
    )


def test_sweep_rejects_bad_range_with_400():
    resp = client.post(
        "/sweep",
        json={"base": VALID_PARAMS, "parameter": "k2", "start": -0.4, "stop": 0.8, "n": 5},
    )
    assert resp.status_code == 400


def test_missing_field_returns_422():
    incomplete = {k: v for k, v in VALID_PARAMS.items() if k != "k1"}
    resp = client.post("/scan", json={"params": incomplete})
    assert resp.status_code == 422
