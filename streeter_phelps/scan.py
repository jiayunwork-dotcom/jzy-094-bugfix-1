"""沿程扫描：把亏氧曲线沿河程扫出来，并定位网格上的最低溶解氧点。

河程 x 与时间 t 通过 t = x / U 换算 —— 扫描在河程坐标上进行，
代入闭式解之前必须先折成时间，绝不能把公里数直接当时间用。
"""

from __future__ import annotations

from dataclasses import dataclass

from . import model
from .critical_point import CriticalPoint, _bisect_deficit_rate_root, critical_point
from .parameters import SagParams
from .validation import InvalidParameterError, require_positive

#: 默认扫描点数
DEFAULT_N_POINTS = 201


@dataclass(frozen=True)
class ScanPoint:
    """沿程一个采样点。"""

    x: float  # 河程 [km]
    t: float  # 对应时间 t = x / U [day]
    bod: float  # 剩余 BOD [mg/L]
    deficit: float  # 亏氧 [mg/L]
    do: float  # 溶解氧 [mg/L]


@dataclass(frozen=True)
class ScanResult:
    """一次沿程扫描的完整结果。"""

    points: list[ScanPoint]
    grid_min: ScanPoint  # 网格上溶解氧最低（亏氧最大）的点
    analytic_critical: CriticalPoint | None  # 解析临界点（无氧垂时为 None）
    refined_critical: CriticalPoint | None  # 手写二分求根复核的临界点
    x_max: float  # 实际使用的扫描终点 [km]


def default_x_max(params: SagParams) -> float:
    """默认扫描终点：保证覆盖氧垂谷底及其后的复氧段。"""
    result = critical_point(params)
    t_c = result.point.t_critical if result.exists and result.point else 0.0
    t_max = max(3.0 * t_c, 5.0 / params.k2)
    return params.u * t_max


def scan_along_river(
    params: SagParams,
    x_max: float | None = None,
    n_points: int = DEFAULT_N_POINTS,
) -> ScanResult:
    """在 [0, x_max] 河程区间上均匀采样，返回亏氧/溶解氧曲线及临界点。"""
    if isinstance(n_points, bool) or not isinstance(n_points, int) or n_points < 2:
        raise InvalidParameterError(f"扫描点数 n_points 必须是 >= 2 的整数，收到 {n_points!r}")

    if x_max is None:
        x_max = default_x_max(params)
    else:
        x_max = require_positive("扫描终点 x_max", x_max)

    step = x_max / (n_points - 1)
    points: list[ScanPoint] = []
    for i in range(n_points):
        x = i * step
        t = x / params.u  # 关键换算：河程 -> 时间
        d = model.deficit(t, params)
        points.append(
            ScanPoint(
                x=x,
                t=t,
                bod=model.bod_remaining(t, params),
                deficit=d,
                do=params.csat - d,
            )
        )

    grid_min = min(points, key=lambda p: p.do)

    analytic = critical_point(params)
    analytic_point = analytic.point if analytic.exists else None

    refined_point: CriticalPoint | None = None
    t_refined = _bisect_deficit_rate_root(params)
    if t_refined is not None:
        d_refined = model.deficit(t_refined, params)
        refined_point = CriticalPoint(
            t_critical=t_refined,
            distance=params.u * t_refined,
            deficit=d_refined,
            do=params.csat - d_refined,
            special_case=params.k1 == params.k2,
        )

    return ScanResult(
        points=points,
        grid_min=grid_min,
        analytic_critical=analytic_point,
        refined_critical=refined_point,
        x_max=x_max,
    )
