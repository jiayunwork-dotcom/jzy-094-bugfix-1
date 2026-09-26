"""k1、k2 仅差末位时的回归测试。

背景：现场系数常由换算或扫参网格产生（如 0.2 + 0.1 = 0.30000000000000004）。
这类输入此前因严格 `k1 == k2` 判定失配而落入一般式，极小分母把指数差的
舍入误差放大：临界点跑到 120 km / 4.952 mg/L、二分复核却给 173.7 km，
沿程扫描甚至出现负溶解氧。修复要求：

1. 两者近到数值上无法区分时，临界点、数值复核、沿程曲线、扫参记录都与
   k2 == k1 严格相等时一致（相对偏差 1e-6 以内，实际为逐位一致）；
2. k2 从两侧逼近 k1 时各项结果连续、不跳；
3. k2 = 0.3000003 这种有实际差别的输入仍走一般式，数值不被抹平；
4. 全程溶解氧不得为负。
"""

import pytest

from streeter_phelps.critical_point import critical_point, critical_point_numeric
from streeter_phelps.model import (
    COEFFICIENT_EQUAL_RELTOL,
    coefficients_effectively_equal,
    deficit,
)
from streeter_phelps.parameters import SagParams
from streeter_phelps.scan import scan_along_river
from streeter_phelps.sweep import TREND_DECREASING, sweep_parameter

#: 现场例题工况：k1=0.3/day、U=30 km/day、L0=15、D0=1.5、Csat=10
BASE = dict(k1=0.3, u=30.0, l0=15.0, d0=1.5, csat=10.0)

#: 仅差末位、浮点上无法与 0.3 区分的网格值（[0.2, 0.4] 取 3 点的中间格）
GRID_MIDDLE = 0.2 + 0.1  # 0.30000000000000004
NEAR_BELOW = 0.3 - 5.551115123125783e-17  # 0.29999999999999993
NEAR_ABOVE_TINY = 0.3 + 5.551115123125783e-15  # 0.30000000000000554

#: 修复前出错的精确数值锚点（k2 == k1 特解的正确结果）
EXPECTED_T_CRITICAL = 3.0
EXPECTED_DISTANCE_KM = 90.0
EXPECTED_DEFICIT = 6.098544896108987


def _params(k2: float) -> SagParams:
    return SagParams(k2=k2, **BASE)


# ---------------------------------------------------------------------------
# 临界点：末位误差输入必须与严格相等一致
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k2", [GRID_MIDDLE, NEAR_BELOW, NEAR_ABOVE_TINY])
def test_critical_point_near_equal_matches_exact_special_case(k2):
    exact = critical_point(_params(0.3))
    result = critical_point(_params(k2))
    assert result.exists
    assert result.special_case is True
    # 要求相对偏差 1e-6 以内；实际实现对末位误差输入逐位等同
    assert result.point.distance == pytest.approx(EXPECTED_DISTANCE_KM, rel=1e-6)
    assert result.point.deficit == pytest.approx(EXPECTED_DEFICIT, rel=1e-6)
    assert result.point.t_critical == pytest.approx(EXPECTED_T_CRITICAL, rel=1e-6)
    assert result.point.distance == pytest.approx(exact.point.distance, rel=1e-6)
    assert result.point.deficit == pytest.approx(exact.point.deficit, rel=1e-6)


def test_grid_middle_value_reproduces_demo_incident_values():
    # 现场演砸的那一格：修复前是 120 km / 4.952 mg/L
    result = critical_point(_params(GRID_MIDDLE))
    assert result.point.distance == pytest.approx(90.0, rel=1e-6)
    assert result.point.deficit == pytest.approx(6.098545, rel=1e-6)


def test_numeric_cross_check_agrees_for_near_equal_coefficients():
    # 修复前同一响应里数值复核自相矛盾地给出约 173.7 km
    analytic = critical_point(_params(GRID_MIDDLE))
    numeric = critical_point_numeric(_params(GRID_MIDDLE))
    assert numeric is not None
    assert numeric.point.distance == pytest.approx(EXPECTED_DISTANCE_KM, rel=1e-6)
    assert numeric.point.deficit == pytest.approx(EXPECTED_DEFICIT, rel=1e-6)
    assert numeric.point.distance == pytest.approx(analytic.point.distance, abs=1e-3)
    assert numeric.special_case is True


@pytest.mark.parametrize("k2", [GRID_MIDDLE, NEAR_BELOW, NEAR_ABOVE_TINY])
def test_deficit_curve_near_equal_is_bit_identical_to_special_formula(k2):
    p_exact, p_near = _params(0.3), _params(k2)
    for t in (0.0, 1.0 / 3.0, 1.0, 3.0, 10.0):
        assert deficit(t, p_near) == deficit(t, p_exact)


# ---------------------------------------------------------------------------
# 沿程扫描：不得出现负 DO，最低点、90 km 处取值正确
# ---------------------------------------------------------------------------


def test_scan_near_equal_has_no_negative_dissolved_oxygen():
    # 修复前：x_max=300、301 点，最低 DO 落在 10 km 处且为 -0.357 mg/L
    result = scan_along_river(_params(GRID_MIDDLE), x_max=300.0, n_points=301)
    assert all(pt.do >= 0.0 for pt in result.points)
    assert all(pt.deficit <= BASE["csat"] for pt in result.points)
    # 网格最低点必须在临界河程附近（步长 1 km），而不是 10 km 的假谷底
    assert result.grid_min.x == pytest.approx(EXPECTED_DISTANCE_KM, abs=1.0)
    assert result.grid_min.do == pytest.approx(10.0 - EXPECTED_DEFICIT, abs=1e-3)
    at_90 = next(pt for pt in result.points if pt.x == pytest.approx(90.0))
    assert at_90.deficit == pytest.approx(6.098545, rel=1e-6)
    # 解析临界点与二分复核一致，且都标记为特解
    assert result.analytic_critical.special_case is True
    assert result.refined_critical is not None
    assert result.refined_critical.distance == pytest.approx(EXPECTED_DISTANCE_KM, abs=1e-3)


def test_scan_below_side_near_equal_matches_exact_case():
    # 0.29999999999999993：修复前临界亏氧错成 9.452 mg/L
    result = scan_along_river(_params(NEAR_BELOW), x_max=300.0, n_points=301)
    assert result.analytic_critical.deficit == pytest.approx(EXPECTED_DEFICIT, rel=1e-6)
    assert result.analytic_critical.distance == pytest.approx(EXPECTED_DISTANCE_KM, rel=1e-6)
    assert all(pt.do >= 0.0 for pt in result.points)


# ---------------------------------------------------------------------------
# 区间扫参：中间网格值不再凹口，趋势单调
# ---------------------------------------------------------------------------


def test_sweep_demo_grid_has_smooth_middle_point_and_decreasing_trend():
    summary = sweep_parameter(_params(0.3), "k2", 0.2, 0.4, 3)
    values = [r.value for r in summary.records]
    assert values[1] == GRID_MIDDLE  # 确认扫到的正是出问题的网格值
    assert summary.trend == TREND_DECREASING
    # 两头保持原有结果
    assert summary.records[0].critical_distance == pytest.approx(111.8026, abs=1e-3)
    assert summary.records[0].critical_deficit == pytest.approx(7.3558, abs=1e-3)
    assert summary.records[2].critical_distance == pytest.approx(76.1342, abs=1e-3)
    assert summary.records[2].critical_deficit == pytest.approx(5.2542, abs=1e-3)
    # 中间格与 k2 == k1 一致，不再凭空凹口
    mid = summary.records[1]
    assert mid.critical_distance == pytest.approx(EXPECTED_DISTANCE_KM, rel=1e-6)
    assert mid.critical_deficit == pytest.approx(EXPECTED_DEFICIT, rel=1e-6)
    # 临界亏氧严格随 k2 增大而降低
    deficits = [r.critical_deficit for r in summary.records]
    assert deficits[0] > deficits[1] > deficits[2]


# ---------------------------------------------------------------------------
# 两侧逼近的连续性
# ---------------------------------------------------------------------------


def test_results_continuous_as_k2_approaches_k1_from_both_sides():
    exact = critical_point(_params(0.3)).point
    factors = (-1e-6, -1e-9, -1e-12, -1e-14, 0.0, 1e-14, 1e-12, 1e-9, 1e-6)
    distances, deficits = [], []
    for f in factors:
        k2 = 0.3 * (1.0 + f) if f != 0.0 else 0.3
        point = critical_point(_params(k2)).point
        distances.append(point.distance)
        deficits.append(point.deficit)
        # 越接近相等越贴近特解
        assert point.distance == pytest.approx(exact.distance, rel=max(2.0 * abs(f), 1e-12))
        assert point.deficit == pytest.approx(exact.deficit, rel=max(2.0 * abs(f), 1e-12))
    # 相邻取值之间不允许跳变（相对差不超过物理变化量的合理上界）
    for a, b in zip(distances, distances[1:]):
        assert abs(b - a) / exact.distance < 1e-6
    for a, b in zip(deficits, deficits[1:]):
        assert abs(b - a) / exact.deficit < 1e-6


# ---------------------------------------------------------------------------
# 有实际差别的输入不得被抹平
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k2", [0.3000003, 0.2999997])
def test_genuinely_different_k2_still_uses_general_formula(k2):
    # 相对差 1e-6，超过「数值上不可分」容差，必须照常按一般式计算
    assert not coefficients_effectively_equal(0.3, k2)
    result = critical_point(_params(k2))
    assert result.special_case is False
    assert result.point.distance == pytest.approx(EXPECTED_DISTANCE_KM, abs=1e-3)
    assert result.point.deficit == pytest.approx(EXPECTED_DEFICIT, abs=1e-4)
    # 与严格相等确有可分辨的差别，未被抹平
    exact = critical_point(_params(0.3)).point
    assert result.point.distance != exact.distance


def test_effective_equality_predicate_boundaries():
    assert coefficients_effectively_equal(0.3, GRID_MIDDLE)
    assert coefficients_effectively_equal(0.3, NEAR_BELOW)
    assert not coefficients_effectively_equal(0.3, 0.3 * (1.0 + 10.0 * COEFFICIENT_EQUAL_RELTOL))
    # 容差阈值本身按相对尺度判定
    assert coefficients_effectively_equal(300.0, 300.0 * (1.0 + 0.5 * COEFFICIENT_EQUAL_RELTOL))
