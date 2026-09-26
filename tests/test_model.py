"""闭式解求值测试：一般式、特解、边界与非法输入。"""

import math

import pytest

from streeter_phelps.model import bod_remaining, deficit, deficit_rate, dissolved_oxygen
from streeter_phelps.parameters import SagParams
from streeter_phelps.validation import InvalidParameterError

BASE = dict(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)


def test_deficit_at_t_zero_equals_d0():
    p = SagParams(**BASE)
    assert deficit(0.0, p) == pytest.approx(2.0)
    assert dissolved_oxygen(0.0, p) == pytest.approx(8.0)


def test_bod_decays_exponentially():
    p = SagParams(**BASE)
    t = 3.7
    assert bod_remaining(t, p) == pytest.approx(10.0 * math.exp(-0.2 * t))


def test_deficit_matches_closed_form_general_case():
    p = SagParams(**BASE)
    t = 2.5
    expected = 0.2 * 10.0 / (0.4 - 0.2) * (math.exp(-0.2 * t) - math.exp(-0.4 * t)) + 2.0 * math.exp(
        -0.4 * t
    )
    assert deficit(t, p) == pytest.approx(expected, rel=1e-12)


def test_deficit_special_case_k1_equals_k2():
    # k1 == k2 必须走特解 D = (k·L0·t + D0)·exp(-k·t)，而非一般分母公式
    p = SagParams(k1=0.3, k2=0.3, u=30.0, l0=10.0, d0=1.0, csat=10.0)
    t = 1.8
    expected = (0.3 * 10.0 * t + 1.0) * math.exp(-0.3 * t)
    assert deficit(t, p) == pytest.approx(expected, rel=1e-12)


def test_special_case_is_continuous_limit_of_general_case():
    # k2 从两侧逼近 k1 时，一般式的值应收敛到特解
    t = 2.3
    common = dict(u=30.0, l0=8.0, d0=1.5, csat=10.0)
    d_special = deficit(t, SagParams(k1=0.25, k2=0.25, **common))
    for sign in (+1, -1):
        eps = 1e-7
        d_near = deficit(t, SagParams(k1=0.25, k2=0.25 + sign * eps, **common))
        assert d_near == pytest.approx(d_special, rel=1e-5)


def test_deficit_rate_satisfies_ode():
    # dD/dt 必须等于 k1·L - k2·D（与闭式解自洽）
    p = SagParams(**BASE)
    t = 1.3
    expected = p.k1 * bod_remaining(t, p) - p.k2 * deficit(t, p)
    assert deficit_rate(t, p) == pytest.approx(expected, rel=1e-12)


def test_negative_time_rejected():
    p = SagParams(**BASE)
    with pytest.raises(InvalidParameterError):
        deficit(-1.0, p)
