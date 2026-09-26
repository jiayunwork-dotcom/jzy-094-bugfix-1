"""Streeter–Phelps 闭式解求值。

模型（碳质 BOD 单级耗氧 + 大气复氧）：
    L(t) = L0 · exp(-k1 · t)
    dD/dt = k1 · L - k2 · D

一般情形 k1 != k2：
    D(t) = k1·L0/(k2-k1) · (exp(-k1·t) - exp(-k2·t)) + D0 · exp(-k2·t)

特解 k1 == k2 = k（对一般式取 k2 -> k1 的极限）：
    D(t) = (k·L0·t + D0) · exp(-k·t)

k1 与 k2 相近到相对差小于 NEARLY_EQUAL_REL_TOL 时，一般式中
exp(-k1·t) - exp(-k2·t) 与 1/(k2-k1) 发生灾难性抵消，双精度下结果完全
失真（亏氧可算出超过饱和值的假值）；此时一律按特解求值 —— 特解正是
一般式在 k2 -> k1 时的连续极限。该阈值远小于任何有实际意义的系数差
（1e-6 量级），不会抹平真实差别。

实际溶解氧 DO(t) = Csat - D(t)。
"""

from __future__ import annotations

import math

from .parameters import SagParams
from .validation import InvalidParameterError

#: k1 与 k2 的相对差小于该阈值时视为「数值上分不开」，统一走 k1 == k2 特解。
#: 取 1e-9：远大于双精度舍入噪声开始污染一般式的量级（~1e-13），又远小于
#: 现场换算系数的有意义差别（>= 1e-6 量级），两侧都留足余量。
NEARLY_EQUAL_REL_TOL = 1e-9


def coefficients_effectively_equal(k1: float, k2: float) -> bool:
    """k1 与 k2 是否在数值上分不开（含严格相等）。

    为 True 时一般式会遭遇灾难性抵消，所有计算必须改走 k1 == k2 特解
    （一般式的连续极限），并以 k1 为代表取值，保证与 k2 严格等于 k1
    时的结果逐位一致。
    """
    return math.isclose(k1, k2, rel_tol=NEARLY_EQUAL_REL_TOL, abs_tol=0.0)


def bod_remaining(t: float, params: SagParams) -> float:
    """时刻 t 的剩余碳质 BOD：L = L0 · exp(-k1 · t)。"""
    _check_time(t)
    return params.l0 * math.exp(-params.k1 * t)


def deficit(t: float, params: SagParams) -> float:
    """时刻 t 的亏氧 D(t) [mg/L]，闭式解；k1 与 k2 数值上分不开时走特解。"""
    _check_time(t)
    k1, k2, l0, d0 = params.k1, params.k2, params.l0, params.d0
    if coefficients_effectively_equal(k1, k2):
        # 特解：D = (k·L0·t + D0)·exp(-k·t)；以 k1 为代表取值，
        # 保证与 k2 严格等于 k1 时逐位一致
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
