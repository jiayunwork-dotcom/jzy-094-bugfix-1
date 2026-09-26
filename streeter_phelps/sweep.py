"""单系数区间批量扫参：汇报临界亏氧随该系数的变化趋势。

在 [start, stop] 上均匀取 n 个系数值，逐个构造工况、解析定位临界点，
汇总临界亏氧 / 临界距离序列及整体变化趋势。某个取值下不存在氧垂时，
该点如实记录 exists=False，绝不编造临界值。
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from .critical_point import critical_point
from .parameters import SagParams
from .validation import (
    InvalidParameterError,
    require_non_negative,
    require_positive,
)

#: 允许扫参的系数
SWEEPABLE_PARAMETERS = ("k1", "k2", "u", "l0", "d0", "csat")

#: 趋势判定码
TREND_INCREASING = "increasing"  # 临界亏氧随系数增大而升高
TREND_DECREASING = "decreasing"  # 临界亏氧随系数增大而降低
TREND_MIXED = "mixed"  # 区间内部分取值无临界点，无法判定单一趋势


@dataclass(frozen=True)
class SweepRecord:
    """区间内一个系数取值下的临界结果。"""

    value: float  # 该点系数取值
    exists: bool  # 该工况下是否存在临界点
    t_critical: float | None
    critical_distance: float | None  # [km]
    critical_deficit: float | None  # [mg/L]
    critical_do: float | None  # [mg/L]
    reason: str


@dataclass(frozen=True)
class SweepSummary:
    """整个区间的扫参汇总。"""

    parameter: str
    trend: str  # TREND_* 之一
    records: list[SweepRecord]


def sweep_parameter(
    base: SagParams,
    parameter: str,
    start: float,
    stop: float,
    n: int,
) -> SweepSummary:
    """在 [start, stop] 区间上对单个系数扫参，其余工况参数保持不变。"""
    if parameter not in SWEEPABLE_PARAMETERS:
        raise InvalidParameterError(
            f"不支持对 {parameter!r} 扫参，可选：{list(SWEEPABLE_PARAMETERS)}"
        )
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise InvalidParameterError(f"扫参点数 n 必须是 >= 2 的整数，收到 {n!r}")

    if parameter in ("k1", "k2", "u", "csat"):
        start = require_positive(f"{parameter} 区间起点", start)
        stop = require_positive(f"{parameter} 区间终点", stop)
    else:  # l0, d0 允许从 0 起
        start = require_non_negative(f"{parameter} 区间起点", start)
        stop = require_non_negative(f"{parameter} 区间终点", stop)
    if stop < start:
        raise InvalidParameterError(f"区间终点 {stop} 小于起点 {start}")

    step = (stop - start) / (n - 1)
    records: list[SweepRecord] = []
    for i in range(n):
        value = start + i * step
        params = replace(base, **{parameter: value})  # replace 会重新触发校验
        result = critical_point(params)
        if result.exists and result.point is not None:
            p = result.point
            records.append(
                SweepRecord(
                    value=value,
                    exists=True,
                    t_critical=p.t_critical,
                    critical_distance=p.distance,
                    critical_deficit=p.deficit,
                    critical_do=p.do,
                    reason=result.reason,
                )
            )
        else:
            records.append(
                SweepRecord(
                    value=value,
                    exists=False,
                    t_critical=None,
                    critical_distance=None,
                    critical_deficit=None,
                    critical_do=None,
                    reason=result.reason,
                )
            )

    deficits = [r.critical_deficit for r in records if r.exists]
    if len(deficits) != len(records) or len(deficits) < 2:
        # 区间内有取值不存在氧垂，无法判定单一趋势
        trend = TREND_MIXED
    elif deficits[-1] > deficits[0]:
        trend = TREND_INCREASING
    elif deficits[-1] < deficits[0]:
        trend = TREND_DECREASING
    else:
        trend = TREND_MIXED

    return SweepSummary(parameter=parameter, trend=trend, records=records)
