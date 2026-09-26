"""输入校验测试：耗氧系数、复氧系数、流速、饱和溶解氧、初始 BOD 任一非法即报错。"""

import pytest

from streeter_phelps.parameters import SagParams
from streeter_phelps.validation import InvalidParameterError

BASE = dict(k1=0.2, k2=0.4, u=30.0, l0=10.0, d0=2.0, csat=10.0)


def test_valid_params_accepted():
    p = SagParams(**BASE)
    assert p.k1 == 0.2 and p.k2 == 0.4 and p.u == 30.0


@pytest.mark.parametrize("field", ["k1", "k2", "u", "csat"])
@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf"), -float("inf")])
def test_positive_fields_reject_non_positive_and_non_finite(field, bad):
    with pytest.raises(InvalidParameterError):
        SagParams(**{**BASE, field: bad})


@pytest.mark.parametrize("bad", [-0.5, float("nan"), float("inf")])
def test_initial_bod_rejects_negative_and_non_finite(bad):
    with pytest.raises(InvalidParameterError):
        SagParams(**{**BASE, "l0": bad})


def test_initial_bod_zero_allowed():
    p = SagParams(**{**BASE, "l0": 0.0})
    assert p.l0 == 0.0


def test_initial_deficit_rejects_negative():
    with pytest.raises(InvalidParameterError):
        SagParams(**{**BASE, "d0": -0.1})


def test_initial_deficit_must_not_exceed_saturation():
    # D0 > Csat 意味着初始溶解氧为负，非物理
    with pytest.raises(InvalidParameterError):
        SagParams(**{**BASE, "d0": 10.5})


@pytest.mark.parametrize("field", ["k1", "k2", "u", "l0", "d0", "csat"])
@pytest.mark.parametrize("bad", ["0.2", None, True, [0.2]])
def test_non_numeric_values_rejected(field, bad):
    with pytest.raises(InvalidParameterError):
        SagParams(**{**BASE, field: bad})
