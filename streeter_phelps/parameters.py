"""工况参数对象：一组 Streeter–Phelps 计算所需的全部输入。

单位约定（同一套时间基准即可，默认取「天」）：
- k1, k2 : 1/day
- U      : km/day（时间与河程经 t = x / U 换算，绝不能把河长直接当时间代入）
- L0, D0, Csat : mg/L
"""

from __future__ import annotations

from dataclasses import dataclass

from .validation import (
    validate_deoxygenation_coefficient,
    validate_initial_bod,
    validate_initial_deficit,
    validate_reaeration_coefficient,
    validate_saturation_do,
    validate_velocity,
)


@dataclass(frozen=True)
class SagParams:
    """一组完整工况。构造时即完成校验，非法参数直接抛 InvalidParameterError。"""

    k1: float  # 耗氧系数 [1/day]
    k2: float  # 复氧系数 [1/day]
    u: float  # 流速 [km/day]
    l0: float  # 初始 BOD [mg/L]
    d0: float  # 初始亏氧 [mg/L]
    csat: float  # 饱和溶解氧 [mg/L]

    def __post_init__(self) -> None:
        object.__setattr__(self, "k1", validate_deoxygenation_coefficient(self.k1))
        object.__setattr__(self, "k2", validate_reaeration_coefficient(self.k2))
        object.__setattr__(self, "u", validate_velocity(self.u))
        object.__setattr__(self, "l0", validate_initial_bod(self.l0))
        object.__setattr__(self, "csat", validate_saturation_do(self.csat))
        object.__setattr__(self, "d0", validate_initial_deficit(self.d0, self.csat))
