"""临界点定位测试：解析式、k1==k2 特解、无临界点判定、二分求根复核。"""

import math

import pytest

from streeter_phelps.critical_point import (
    REASON_MONOTONIC_REAERATION,
    critical_point,
    critical_point_numeric,
)
from streeter_phelps.model import deficit, deficit_rate
from streeter_phelps.parameters import SagParams

BASE = dict(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)


def test_critical_point_matches_textbook_formula():
    p = SagParams(**BASE)
    result = critical_point(p)
    assert result.exists
    assert not result.special_case
    # t_c = ln[(k2/k1)·(1 − D0·(k2−k1)/(k1·L0))] / (k2 − k1)
    expected_t = math.log((0.4 / 0.2) * (1.0 - 2.0 * (0.4 - 0.2) / (0.2 * 10.0))) / (0.4 - 0.2)
    assert result.point.t_critical == pytest.approx(expected_t, rel=1e-12)
    # 临界距离必须由流速把时间折算成河程
    assert result.point.distance == pytest.approx(30.0 * expected_t)
    # 临界点上 dD/dt = 0
    assert deficit_rate(result.point.t_critical, p) == pytest.approx(0.0, abs=1e-9)
    # 临界亏氧是全程最大亏氧
    assert result.point.deficit == pytest.approx(deficit(result.point.t_critical, p))


def test_critical_point_is_global_deficit_maximum():
    p = SagParams(**BASE)
    result = critical_point(p)
    t_c = result.point.t_critical
    for t in (0.0, 0.5 * t_c, 2.0 * t_c, 5.0 * t_c):
        assert deficit(t, p) <= result.point.deficit + 1e-12


def test_special_case_k1_equals_k2_uses_closed_form():
    # k1 == k2 = k：t_c = 1/k − D0/(k·L0)，不能套一般分母公式
    p = SagParams(k1=0.3, k2=0.3, u=25.0, l0=10.0, d0=1.0, csat=10.0)
    result = critical_point(p)
    assert result.exists
    assert result.special_case
    expected_t = 1.0 / 0.3 - 1.0 / (0.3 * 10.0)
    assert result.point.t_critical == pytest.approx(expected_t, rel=1e-12)
    assert result.point.distance == pytest.approx(25.0 * expected_t)
    assert deficit_rate(result.point.t_critical, p) == pytest.approx(0.0, abs=1e-9)


def test_no_critical_point_when_monotonic_reaeration():
    # k1·L0 <= k2·D0：全程单调复氧，必须报告无临界点，绝不能给负河程
    p = SagParams(k1=0.2, k2=0.5, u=30.0, l0=5.0, d0=4.0, csat=10.0)
    result = critical_point(p)
    assert not result.exists
    assert result.point is None
    assert result.reason == REASON_MONOTONIC_REAERATION


def test_no_critical_point_when_initial_deficit_dominates_equal_coefficients():
    # k1 == k2 且 D0 >= L0：特解公式会给出 t_c <= 0，同样必须判为无临界点
    p = SagParams(k1=0.3, k2=0.3, u=30.0, l0=2.0, d0=3.0, csat=10.0)
    result = critical_point(p)
    assert not result.exists
    assert result.point is None


def test_zero_initial_bod_means_no_sag():
    p = SagParams(k1=0.2, k2=0.4, u=30.0, l0=0.0, d0=1.0, csat=10.0)
    result = critical_point(p)
    assert not result.exists


def test_numeric_bisection_agrees_with_analytic():
    # 手写二分求根与解析式交叉复核（一般情形）
    p = SagParams(**BASE)
    analytic = critical_point(p)
    numeric = critical_point_numeric(p)
    assert numeric is not None
    assert numeric.point.t_critical == pytest.approx(analytic.point.t_critical, abs=1e-6)
    assert numeric.point.deficit == pytest.approx(analytic.point.deficit, abs=1e-6)


def test_numeric_bisection_agrees_with_analytic_special_case():
    # 手写二分求根与解析式交叉复核（k1 == k2 特解）
    p = SagParams(k1=0.3, k2=0.3, u=25.0, l0=10.0, d0=1.0, csat=10.0)
    analytic = critical_point(p)
    numeric = critical_point_numeric(p)
    assert numeric is not None
    assert numeric.point.t_critical == pytest.approx(analytic.point.t_critical, abs=1e-6)


def test_numeric_returns_none_when_no_critical_point():
    p = SagParams(k1=0.2, k2=0.5, u=30.0, l0=5.0, d0=4.0, csat=10.0)
    assert critical_point_numeric(p) is None
