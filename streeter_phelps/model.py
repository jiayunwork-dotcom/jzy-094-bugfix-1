"""Streeter–Phelps 闭式解求值。

模型（碳质 BOD 单级耗氧 + 大气复氧）：
    L(t) = L0 · exp(-k1 · t)
    dD/dt = k1 · L - k2 · D

一般情形 k1 != k2：
    D(t) = k1·L0/(k2-k1) · (exp(-k1·t) - exp(-k2·t)) + D0 · exp(-k2·t)

特解 k1 == k2 = k（对一般式取 k2 -> k1 的极限）：
    D(t) = (k·L0·t + D0) · exp(-k·t)

实际溶解氧 DO(t) = Csat - D(t)。
"""

from __future__ import annotations

import math

from .parameters import SagParams
from .validation import InvalidParameterError


def bod_remaining(t: float, params: SagParams) -> float:
    """时刻 t 的剩余碳质 BOD：L = L0 · exp(-k1 · t)。"""
    _check_time(t)
    return params.l0 * math.exp(-params.k1 * t)


def deficit(t: float, params: SagParams) -> float:
    """时刻 t 的亏氧 D(t) [mg/L]，闭式解，k1 == k2 时自动走特解。"""
    _check_time(t)
    k1, k2, l0, d0 = params.k1, params.k2, params.l0, params.d0
    if k1 == k2:
        # 特解：D = (k·L0·t + D0)·exp(-k·t)
        return (k1 * l0 * t + d0) * math.exp(-k1 * t)
    return (
        k1 * l0 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))
        + d0 * math.exp(-k2 * t)
    )


def dissolved_oxygen(t: float, params: SagParams) -> float:
    """时刻 t 的实际溶解氧 DO = Csat - D(t)。"""
    return params.csat - deficit(t, params)


def deficit_rate(t: float, params: SagParams) -> float:
    """亏氧变化率 dD/dt = k1·L0·exp(-k1·t) - k2·D(t)，供求根与临界点判定用。"""
    _check_time(t)
    return params.k1 * params.l0 * math.exp(-params.k1 * t) - params.k2 * deficit(
        t, params
    )


def state_at(t: float, params: SagParams) -> dict[str, float]:
    """时刻 t 的全量状态，便于接口层一次性取用。"""
    d = deficit(t, params)
    return {
        "t": t,
        "bod": bod_remaining(t, params),
        "deficit": d,
        "do": params.csat - d,
        "deficit_rate": params.k1 * params.l0 * math.exp(-params.k1 * t) - params.k2 * d,
    }


def _check_time(t: float) -> None:
    if isinstance(t, bool) or not isinstance(t, (int, float)) or not math.isfinite(t):
        raise InvalidParameterError(f"时间 t 必须是有限数值，收到 {t!r}")
    if t < 0.0:
        raise InvalidParameterError(f"时间 t 不能为负，收到 {t}")
