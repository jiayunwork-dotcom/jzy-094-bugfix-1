"""沿程扫描测试：时间-河程换算、网格最低点与临界点一致性、参考工况。"""

import pytest

from streeter_phelps.critical_point import critical_point
from streeter_phelps.model import deficit
from streeter_phelps.parameters import SagParams
from streeter_phelps.scan import scan_along_river
from streeter_phelps.scenarios import REFERENCE_SCENARIO
from streeter_phelps.validation import InvalidParameterError


def test_scan_converts_distance_to_time_via_velocity():
    # 每个采样点必须满足 t = x / U，绝不能把河长直接当时间
    p = SagParams(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)
    result = scan_along_river(p, x_max=150.0, n_points=51)
    assert len(result.points) == 51
    for pt in result.points:
        assert pt.t == pytest.approx(pt.x / 30.0, rel=1e-12)
        assert pt.deficit == pytest.approx(deficit(pt.x / 30.0, p), rel=1e-12)
        assert pt.do == pytest.approx(10.0 - pt.deficit, rel=1e-12)


def test_scan_grid_min_tracks_analytic_critical_point():
    p = SagParams(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)
    result = scan_along_river(p, x_max=200.0, n_points=401)
    analytic = critical_point(p).point
    # 网格最低点应贴近解析临界点（网格步长 0.5 km）
    assert abs(result.grid_min.x - analytic.distance) <= 0.5 + 1e-9
    assert result.grid_min.do == pytest.approx(analytic.do, abs=1e-2)
    # 二分复核的临界点与解析结果一致
    assert result.refined_critical is not None
    assert result.refined_critical.t_critical == pytest.approx(
        analytic.t_critical, abs=1e-6
    )


def test_scan_default_x_max_covers_the_sag():
    p = SagParams(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)
    result = scan_along_river(p)  # 不给 x_max，自动取
    analytic = critical_point(p).point
    assert result.x_max > analytic.distance  # 扫描范围必须覆盖谷底
    assert result.analytic_critical is not None


def test_scan_monotonic_case_reports_no_critical_point():
    p = SagParams(k1=0.2, k2=0.5, u=30.0, l0=5.0, d0=4.0, csat=10.0)
    result = scan_along_river(p, x_max=100.0, n_points=51)
    assert result.analytic_critical is None
    assert result.refined_critical is None
    # 单调复氧：网格最低点就是排污口
    assert result.grid_min.x == pytest.approx(0.0)


def test_scan_rejects_invalid_grid():
    p = SagParams(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)
    with pytest.raises(InvalidParameterError):
        scan_along_river(p, x_max=0.0)
    with pytest.raises(InvalidParameterError):
        scan_along_river(p, x_max=100.0, n_points=1)


def test_reference_scenario_has_positive_critical_distance():
    # 预置参考工况：必须出现明显氧垂、临界距离为正
    result = scan_along_river(REFERENCE_SCENARIO)
    analytic = result.analytic_critical
    assert analytic is not None
    assert analytic.distance > 0.0
    assert analytic.t_critical == pytest.approx(2.939, abs=0.01)
    assert analytic.distance == pytest.approx(88.2, abs=0.5)
    assert analytic.deficit == pytest.approx(4.167, abs=0.01)
    # 氧垂要"明显"：临界亏氧显著大于初始亏氧
    assert analytic.deficit > 2.0 * REFERENCE_SCENARIO.d0
    # 最低溶解氧低于饱和值但保持为正
    assert 0.0 < analytic.do < REFERENCE_SCENARIO.csat
