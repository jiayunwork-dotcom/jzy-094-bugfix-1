"""输入参数校验。

所有工况参数（耗氧系数、复氧系数、流速、饱和溶解氧、初始 BOD、初始亏氧）
在进入任何计算之前都必须先过这里的校验，非法一律抛 InvalidParameterError。
"""

from __future__ import annotations

import math


class InvalidParameterError(ValueError):
    """工况参数非法（非有限数、非正、量纲矛盾等）。"""


def _as_finite_float(name: str, value: object) -> float:
    # bool 是 int 的子类，单独挡掉，避免 True 被当成 1.0 蒙混过关
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidParameterError(f"{name} 必须是数值，收到 {value!r}")
    v = float(value)
    if not math.isfinite(v):
        raise InvalidParameterError(f"{name} 必须是有限数值，收到 {value!r}")
    return v


def require_positive(name: str, value: object) -> float:
    """要求严格为正（> 0）。"""
    v = _as_finite_float(name, value)
    if v <= 0.0:
        raise InvalidParameterError(f"{name} 必须为正数，收到 {v}")
    return v


def require_non_negative(name: str, value: object) -> float:
    """要求非负（>= 0）。"""
    v = _as_finite_float(name, value)
    if v < 0.0:
        raise InvalidParameterError(f"{name} 必须为非负数，收到 {v}")
    return v


def validate_deoxygenation_coefficient(k1: object) -> float:
    """耗氧系数 k1 [1/day]，必须为正。"""
    return require_positive("耗氧系数 k1", k1)


def validate_reaeration_coefficient(k2: object) -> float:
    """复氧系数 k2 [1/day]，必须为正。"""
    return require_positive("复氧系数 k2", k2)


def validate_velocity(u: object) -> float:
    """河流流速 U [km/day]，必须为正；时间-河程换算 t = x / U 全靠它。"""
    return require_positive("流速 U", u)


def validate_saturation_do(csat: object) -> float:
    """饱和溶解氧 Csat [mg/L]，必须为正。"""
    return require_positive("饱和溶解氧 Csat", csat)


def validate_initial_bod(l0: object) -> float:
    """初始 BOD L0 [mg/L]，允许为 0（此时不会有氧垂），不允许为负。"""
    return require_non_negative("初始 BOD L0", l0)


def validate_initial_deficit(d0: object, csat: float) -> float:
    """初始亏氧 D0 [mg/L]，非负且不能超过饱和溶解氧（否则初始 DO 为负，非物理）。"""
    v = require_non_negative("初始亏氧 D0", d0)
    if v > csat:
        raise InvalidParameterError(
            f"初始亏氧 D0={v} 超过饱和溶解氧 Csat={csat}，初始溶解氧将为负，非物理"
        )
    return v
