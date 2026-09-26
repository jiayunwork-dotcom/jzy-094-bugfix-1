"""区间批量扫参测试：趋势汇报、无临界点记录、非法区间。"""

import pytest

from streeter_phelps.parameters import SagParams
from streeter_phelps.sweep import (
    TREND_DECREASING,
    TREND_INCREASING,
    TREND_MIXED,
    sweep_parameter,
)
from streeter_phelps.validation import InvalidParameterError

BASE = SagParams(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)


def test_k2_sweep_reports_decreasing_critical_deficit():
    # 复氧系数区间扫参：临界亏氧随 k2 增大而降低
    summary = sweep_parameter(BASE, "k2", 0.4, 0.8, 9)
    assert summary.parameter == "k2"
    assert summary.trend == TREND_DECREASING
    assert len(summary.records) == 9
    assert all(r.exists for r in summary.records)
    deficits = [r.critical_deficit for r in summary.records]
    assert all(b < a for a, b in zip(deficits, deficits[1:]))


def test_l0_sweep_reports_increasing_critical_deficit():
    # 初始 BOD 区间扫参：临界亏氧随 L0 增大而升高
    summary = sweep_parameter(BASE, "l0", 5.0, 15.0, 6)
    assert summary.trend == TREND_INCREASING
    deficits = [r.critical_deficit for r in summary.records]
    assert all(b > a for a, b in zip(deficits, deficits[1:]))


def test_u_sweep_scales_critical_distance_proportionally():
    # 流速扫参：临界距离与流速成正比，临界亏氧不变
    summary = sweep_parameter(BASE, "u", 20.0, 60.0, 5)
    assert all(r.exists for r in summary.records)
    distances = [r.critical_distance for r in summary.records]
    velocities = [r.value for r in summary.records]
    for v, x in zip(velocities, distances):
        assert x == pytest.approx(distances[0] * v / velocities[0], rel=1e-9)
    deficits = [r.critical_deficit for r in summary.records]
    assert all(d == pytest.approx(deficits[0], rel=1e-12) for d in deficits)


def test_sweep_records_no_critical_point_honestly():
    # 区间内部分取值 k1·L0 <= k2·D0：这些点必须如实记录无临界点
    summary = sweep_parameter(BASE, "l0", 0.0, 10.0, 6)  # L0 <= 4 时无氧垂
    assert not summary.records[0].exists
    assert summary.records[0].critical_deficit is None
    assert summary.records[0].critical_distance is None
    assert summary.records[-1].exists
    assert summary.trend == TREND_MIXED


def test_sweep_rejects_invalid_parameter_name():
    with pytest.raises(InvalidParameterError):
        sweep_parameter(BASE, "temperature", 0.1, 1.0, 5)


def test_sweep_rejects_invalid_range():
    with pytest.raises(InvalidParameterError):
        sweep_parameter(BASE, "k2", -0.1, 0.8, 5)  # k2 必须为正
    with pytest.raises(InvalidParameterError):
        sweep_parameter(BASE, "k2", 0.8, 0.4, 5)  # 终点小于起点
    with pytest.raises(InvalidParameterError):
        sweep_parameter(BASE, "k2", 0.4, 0.8, 1)  # 点数不足
