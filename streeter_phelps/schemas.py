"""接口层请求模型（Pydantic）。

这里只做结构与类型约束；物理合法性（正数、量纲矛盾等）统一由
validation 模块在构造 SagParams 时把关，保证核心计算不依赖接口层。
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class ParamsIn(BaseModel):
    """一组完整工况。单位：k1/k2 [1/day]，u [km/day]，其余 [mg/L]。"""

    k1: float = Field(description="耗氧系数 [1/day]")
    k2: float = Field(description="复氧系数 [1/day]")
    u: float = Field(description="流速 [km/day]")
    l0: float = Field(description="初始 BOD [mg/L]")
    d0: float = Field(description="初始亏氧 [mg/L]")
    csat: float = Field(description="饱和溶解氧 [mg/L]")


class ScanRequest(BaseModel):
    params: ParamsIn
    x_max: float | None = Field(
        default=None, description="扫描终点河程 [km]；缺省时按工况自动取"
    )
    n_points: int = Field(default=201, description="采样点数，>= 2")


class SweepRequest(BaseModel):
    base: ParamsIn
    parameter: str = Field(description="扫参系数名：k1/k2/u/l0/d0/csat")
    start: float
    stop: float
    n: int = Field(default=11, description="区间采样点数，>= 2")
