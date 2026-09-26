"""Streeter–Phelps 闭式解求值。

模型（碳质 BOD 单级耗氧 + 大气复氧）：
    L(t) = L0 · exp(-k1 · t)
    dD/dt = k1 · L - k2 · D

一般情形 k1 != k2：
    D(t) = k1·L0/(k2-k1) · (exp(-k1·t) - exp(-k2·t)) + D0 · exp(-k2·t)

特解 k1 == k2 = k（对一般式取 k2 -> k1 的极限）：
    D(t) = (k·L0·t + D0) · exp(-k·t)

实际溶解氧 DO(t) = Csat - D(t)。

数值注意：现场系数常由换算或扫参网格产生，k2 与 k1 可能仅差最后几位
（如 0.2 + 0.1 = 0.30000000000000004）。这种差值在浮点上远小于
机器精度量级的「有效差别」，却会让一般式分母 k2-k1 极小、把
exp(-k1·t) - exp(-k2·t) 的舍入误差放大成完全错误的曲线（甚至负 DO）。
因此一律经 coefficients_effectively_equal 判定：相对差在
COEFFICIENT_EQUAL_RELTOL 以内即视为 k1 == k2 走特解；一般式本身也用
expm1 稳定求值，保证越过容差阈值时结果连续、不跳变。
"""

from __future__ import annotations

import math

from .parameters import SagParams
from .validation import InvalidParameterError

#: 判定 k1、k2 在数值上不可分的相对容差。
#: 取 1e-9：末位舍入误差（~1e-16）远在其内应走特解；
#: 而 0.3000003（相对差 1e-6）这类真实差别必须照常按一般式计算，不被抹平。
COEFFICIENT_EQUAL_RELTOL = 1e-9


def coefficients_effectively_equal(k1: float, k2: float) -> bool:
    """k1、k2 是否近到浮点上无法区分（相对差 <= COEFFICIENT_EQUAL_RELTOL）。

    以两者中较大值为尺度；k1、k2 经校验恒为正，无需处理零与非有限值。
    """
    return math.fabs(k2 - k1) <= COEFFICIENT_EQUAL_RELTOL * max(
        math.fabs(k1), math.fabs(k2)
    )


def bod_remaining(t: float, params: SagParams) -> float:
    """时刻 t 的剩余碳质 BOD：L = L0 · exp(-k1 · t)。"""
    _check_time(t)
    return params.l0 * math.exp(-params.k1 * t)


def deficit(t: float, params: SagParams) -> float:
    """时刻 t 的亏氧 D(t) [mg/L]，闭式解，k1 与 k2 数值上相等时自动走特解。"""
    _check_time(t)
    k1, k2, l0, d0 = params.k1, params.k2, params.l0, params.d0
    if coefficients_effectively_equal(k1, k2):
        # 特解（含 k2 仅差末位的输入；直接用 k1，结果与 k1 == k2 严格一致）：
        # D = (k·L0·t + D0)·exp(-k·t)
        return (k1 * l0 * t + d0) * math.exp(-k1 * t)
    # 一般式。exp(-k1·t) - exp(-k2·t) 在 k1、k2 接近时直接相减会发生
    # 灾难性抵消，改写成 expm1 的稳定形式：
    #   e^{-a} - e^{-b} = e^{-b}·expm1(b-a)  （b >= a）
    #                   = e^{-a}·expm1(a-b)  （a > b）
    delta = k2 - k1
    if delta > 0.0:
        coupled = k1 * l0 / delta * math.exp(-k2 * t) * math.expm1(delta * t)
    else:
        # delta < 0：分母为负，指数差同样为负，整体取正，
        # 用 -delta = k1-k2 > 0 写成正因子相除，避免符号出错。
        gap = -delta
        coupled = k1 * l0 / gap * math.exp(-k1 * t) * math.expm1(gap * t)
    return coupled + d0 * math.exp(-k2 * t)


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
