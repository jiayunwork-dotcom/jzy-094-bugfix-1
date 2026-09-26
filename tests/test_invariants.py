"""物理规律不变量测试：复氧增强氧垂变浅、BOD 加倍最大亏氧升高、
流速加倍临界距离按比例拉长。"""

import pytest

from streeter_phelps.critical_point import critical_point
from streeter_phelps.parameters import SagParams

BASE = dict(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)


def test_larger_reaeration_coefficient_shallows_the_sag():
    # 只把复氧系数明显调大：临界亏氧变小、氧垂变浅
    p_base = SagParams(**BASE)
    p_high_k2 = SagParams(**{**BASE, "k2": 0.8})
    d_base = critical_point(p_base).point.deficit
    d_high = critical_point(p_high_k2).point.deficit
    assert d_high < d_base


def test_reaeration_increase_monotonically_reduces_critical_deficit():
    deficits = []
    # k2 增大到 0.9 时 k1·L0 = 2.0 仍大于 k2·D0 = 1.8，氧垂始终存在
    for k2 in (0.4, 0.5, 0.6, 0.8, 0.9):
        p = SagParams(**{**BASE, "k2": k2})
        result = critical_point(p)
        assert result.exists
        deficits.append(result.point.deficit)
    assert deficits == sorted(deficits, reverse=True)
    assert all(b < a for a, b in zip(deficits, deficits[1:]))


def test_doubling_initial_bod_raises_max_deficit():
    # 只把初始 BOD 加倍：最大亏氧升高（csat 取大些避免 DO 穿零干扰判断）
    base = {**BASE, "csat": 12.0}
    d_normal = critical_point(SagParams(**base)).point.deficit
    d_doubled = critical_point(SagParams(**{**base, "l0": 20.0})).point.deficit
    assert d_doubled > d_normal


def test_doubling_velocity_doubles_critical_distance():
    # 流速加倍：同一临界时刻对应的河程加倍，临界距离按比例拉长
    p1 = SagParams(**BASE)
    p2 = SagParams(**{**BASE, "u": 60.0})
    r1 = critical_point(p1).point
    r2 = critical_point(p2).point
    # 临界时刻与流速无关
    assert r2.t_critical == pytest.approx(r1.t_critical, rel=1e-12)
    # 临界距离严格按流速比例放大
    assert r2.distance == pytest.approx(2.0 * r1.distance, rel=1e-12)
    # 临界亏氧不受流速影响
    assert r2.deficit == pytest.approx(r1.deficit, rel=1e-12)


def test_velocity_scales_distance_for_arbitrary_factor():
    p1 = SagParams(**BASE)
    p2 = SagParams(**{**BASE, "u": BASE["u"] * 3.5})
    r1 = critical_point(p1).point
    r2 = critical_point(p2).point
    assert r2.distance == pytest.approx(3.5 * r1.distance, rel=1e-12)
