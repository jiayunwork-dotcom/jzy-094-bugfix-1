"""k2 与 k1 近到数值上分不开时的回归测试。

现场系数常由单位换算而来、区间扫参的网格点也会落在距 k1 仅几个 ulp 的
k2 上（例如 0.2 到 0.4 取 3 点，中间一格就是 0.30000000000000004）。
此时一般式发生灾难性抵消，曾出现：临界距离 120 km（应为 90 km）、
临界亏氧 4.952（应为 6.0985）、数值复核 173.7 km 自相矛盾、沿程 DO 出现
负值、扫参曲线凭空凹口。这类输入必须按 k2 -> k1 的连续极限处理，与
k2 严格等于 k1 时一致；而 k2=0.3000003 这类有实际差别的输入仍按一般式
原样计算，不被抹平。
"""

import math

import pytest
from fastapi.testclient import TestClient

from streeter_phelps.api import app
from streeter_phelps.critical_point import critical_point, critical_point_numeric
from streeter_phelps.model import deficit
from streeter_phelps.parameters import SagParams
from streeter_phelps.scan import scan_along_river
from streeter_phelps.sweep import TREND_DECREASING, sweep_parameter

client = TestClient(app)

#: 演示工况：k1=0.3/day、u=30 km/day、L0=15、D0=1.5、Csat=10
BASE = dict(k1=0.3, u=30.0, l0=15.0, d0=1.5, csat=10.0)

#: k2 严格等于 k1 的参照工况（临界距离 90 km、临界亏氧约 6.0985 mg/L）
STRICT_EQUAL = SagParams(k2=0.3, **BASE)

#: 距 k1 仅几个 ulp、数值上分不开的 k2（曾全部算错）
NEAR_K2 = (
    0.30000000000000004,  # 0.3 上方 1 个 ulp，扫参网格中间点
    0.29999999999999993,  # 0.3 下方 1 个 ulp
    0.300000000000003,  # 相对差 1e-14
)

#: 相对差 1e-6、有实际差别的 k2，必须按一般式原样计算
DISTINCT_K2 = (0.3000003, 0.2999997)


def _params(k2: float) -> SagParams:
    return SagParams(k2=k2, **BASE)


def test_strictly_equal_anchor_values_unchanged():
    # k2 恰好等于 k1：t_c = 1/k - D0/(k·L0) = 3 day，临界距离 90 km
    result = critical_point(STRICT_EQUAL)
    assert result.exists
    assert result.special_case
    assert result.point.t_critical == pytest.approx(3.0, rel=1e-12)
    assert result.point.distance == pytest.approx(90.0, rel=1e-12)
    assert result.point.deficit == pytest.approx(6.0985, abs=1e-4)
    assert result.point.do == pytest.approx(3.9015, abs=1e-4)


@pytest.mark.parametrize("k2", NEAR_K2)
def test_near_equal_k2_critical_point_matches_strict_equality(k2):
    # 临界点必须与 k2 严格等于 k1 时一致（远严于要求的 1e-6）
    near = critical_point(_params(k2)).point
    strict = critical_point(STRICT_EQUAL).point
    assert near.t_critical == pytest.approx(strict.t_critical, rel=1e-9)
    assert near.distance == pytest.approx(strict.distance, rel=1e-9)
    assert near.deficit == pytest.approx(strict.deficit, rel=1e-9)
    assert near.do == pytest.approx(strict.do, rel=1e-9)


@pytest.mark.parametrize("k2", NEAR_K2)
def test_near_equal_k2_numeric_cross_check_agrees_with_analytic(k2):
    # 同一响应里解析结果与二分求根复核必须自洽（曾一个 120 km、一个 173.7 km）
    analytic = critical_point(_params(k2)).point
    numeric = critical_point_numeric(_params(k2))
    assert numeric is not None
    assert numeric.point.distance == pytest.approx(analytic.distance, abs=1e-6)
    assert numeric.point.deficit == pytest.approx(analytic.deficit, abs=1e-6)
    assert numeric.point.distance == pytest.approx(90.0, abs=1e-6)


@pytest.mark.parametrize("k2", NEAR_K2)
def test_near_equal_k2_deficit_curve_matches_strict_equality(k2):
    # 沿程亏氧曲线逐点与严格相等时一致
    for t in (0.0, 1.0 / 3.0, 1.0, 3.0, 10.0 / 3.0, 10.0):
        assert deficit(t, _params(k2)) == pytest.approx(
            deficit(t, STRICT_EQUAL), rel=1e-9
        )


def test_near_equal_k2_scan_has_no_negative_do_and_min_at_90km():
    # 曾复现的沿程扫描：x_max=300、301 个点，最低 DO 落在 10 km 且为 -0.357
    result = scan_along_river(_params(0.30000000000000004), x_max=300.0, n_points=301)
    assert all(pt.do >= 0.0 for pt in result.points)  # DO 不能出现负值
    # 网格最低点落在解析临界点 90 km 处（网格步长 1 km），DO 约 3.9015
    assert abs(result.grid_min.x - 90.0) <= 1.0 + 1e-9
    assert result.grid_min.do == pytest.approx(3.9015, abs=1e-3)
    # 90 km 处亏氧约 6.0985（曾错算成 9.61）
    at_90 = min(result.points, key=lambda p: abs(p.x - 90.0))
    assert at_90.deficit == pytest.approx(6.0985, abs=1e-3)
    # 解析与二分复核的临界点都应在 90 km
    assert result.analytic_critical.distance == pytest.approx(90.0, abs=1e-6)
    assert result.refined_critical.distance == pytest.approx(90.0, abs=1e-6)


def test_sweep_grid_point_landing_on_near_equal_k2():
    # k2 在 [0.2, 0.4] 取 3 点扫参：中间一格恰为 0.30000000000000004
    summary = sweep_parameter(_params(0.3), "k2", 0.2, 0.4, 3)
    assert summary.trend == TREND_DECREASING
    assert [r.value for r in summary.records] == [0.2, 0.30000000000000004, 0.4]
    # 两端结果保持原样
    assert summary.records[0].critical_distance == pytest.approx(111.8, abs=0.1)
    assert summary.records[0].critical_deficit == pytest.approx(7.356, abs=1e-3)
    assert summary.records[2].critical_distance == pytest.approx(76.1, abs=0.1)
    assert summary.records[2].critical_deficit == pytest.approx(5.254, abs=1e-3)
    # 中间点与严格相等时一致（曾算出 120 km / 4.952，曲线凭空凹口）
    strict = critical_point(STRICT_EQUAL).point
    assert summary.records[1].critical_distance == pytest.approx(
        strict.distance, rel=1e-9
    )
    assert summary.records[1].critical_deficit == pytest.approx(
        strict.deficit, rel=1e-9
    )
    # 复氧越强、谷底越浅：距离与亏氧都严格单调，不允许凹口
    distances = [r.critical_distance for r in summary.records]
    deficits = [r.critical_deficit for r in summary.records]
    assert all(b < a for a, b in zip(distances, distances[1:]))
    assert all(b < a for a, b in zip(deficits, deficits[1:]))


def test_results_continuous_as_k2_approaches_k1_from_both_sides():
    # k2 从两侧逼近 k1：与严格相等结果的偏差单调收窄、不跳，
    # 进入「数值上分不开」区间后严格落在连续极限上
    strict = critical_point(STRICT_EQUAL).point
    for sign in (+1.0, -1.0):
        deviations = []
        for delta in (3e-4, 3e-6, 3e-8, 3e-10, 3e-12, 3e-14):
            point = critical_point(_params(0.3 + sign * delta)).point
            deviations.append(abs(point.distance - strict.distance))
        assert all(b <= a for a, b in zip(deviations, deviations[1:]))
        assert deviations[-1] == 0.0
        assert deviations[-2] == 0.0


@pytest.mark.parametrize("k2", DISTINCT_K2)
def test_genuinely_different_k2_computed_as_is(k2):
    # 相对差 1e-6 的输入有实际意义：走一般式，数值不被抹平
    result = critical_point(_params(k2))
    assert not result.special_case
    # 与独立按一般式计算的结果逐位吻合（rel 1e-12）
    expected_t = math.log((k2 / 0.3) * (1.0 - 1.5 * (k2 - 0.3) / (0.3 * 15.0))) / (
        k2 - 0.3
    )
    expected_d = 0.3 * 15.0 / (k2 - 0.3) * (
        math.exp(-0.3 * expected_t) - math.exp(-k2 * expected_t)
    ) + 1.5 * math.exp(-k2 * expected_t)
    assert result.point.t_critical == pytest.approx(expected_t, rel=1e-12)
    assert result.point.deficit == pytest.approx(expected_d, rel=1e-12)
    # 结果贴近但绝不等于 k2=k1 的极限值（约 90 km / 6.0985）
    assert abs(result.point.distance - 90.0) < 0.01
    assert result.point.distance != pytest.approx(90.0, rel=1e-9)
    # 亏氧求值同样走一般式
    t = 2.5
    expected_deficit = 0.3 * 15.0 / (k2 - 0.3) * (
        math.exp(-0.3 * t) - math.exp(-k2 * t)
    ) + 1.5 * math.exp(-k2 * t)
    assert deficit(t, _params(k2)) == pytest.approx(expected_deficit, rel=1e-12)


def test_api_critical_point_near_equal_k2_self_consistent():
    resp = client.post("/critical-point", json={"k2": 0.30000000000000004, **BASE})
    assert resp.status_code == 200
    body = resp.json()
    assert body["exists"] is True
    assert body["special_case"] is True
    assert body["critical_point"]["distance"] == pytest.approx(90.0, abs=1e-6)
    assert body["critical_point"]["deficit"] == pytest.approx(6.0985, abs=1e-3)
    # 解析结果与附带的数值复核自洽
    assert body["numeric_cross_check"]["distance"] == pytest.approx(90.0, abs=1e-6)


def test_api_scan_near_equal_k2_no_negative_do():
    resp = client.post(
        "/scan",
        json={"params": {"k2": 0.30000000000000004, **BASE}, "x_max": 300.0, "n_points": 301},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert all(pt["do"] >= 0.0 for pt in body["points"])
    assert abs(body["grid_min"]["x"] - 90.0) <= 1.0 + 1e-9
    assert body["analytic_critical"]["distance"] == pytest.approx(90.0, abs=1e-6)
    assert body["refined_critical"]["distance"] == pytest.approx(90.0, abs=1e-6)


def test_api_sweep_near_equal_grid_point_monotone():
    resp = client.post(
        "/sweep",
        json={"base": {"k2": 0.3, **BASE}, "parameter": "k2", "start": 0.2, "stop": 0.4, "n": 3},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["trend"] == "decreasing"
    assert body["records"][1]["value"] == 0.30000000000000004
    assert body["records"][1]["critical_distance"] == pytest.approx(90.0, abs=1e-6)
    assert body["records"][1]["critical_deficit"] == pytest.approx(6.0985, abs=1e-3)
    deficits = [r["critical_deficit"] for r in body["records"]]
    assert all(b < a for a, b in zip(deficits, deficits[1:]))
