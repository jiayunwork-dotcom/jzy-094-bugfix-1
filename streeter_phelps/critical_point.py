"""临界点（最大亏氧 / 最低溶解氧）定位。

临界时刻满足 dD/dt = 0，解析式：

    一般情形 k1 != k2：
        t_c = ln[ (k2/k1) · (1 - D0·(k2-k1)/(k1·L0)) ] / (k2 - k1)
    特解 k1 == k2 = k：
        t_c = 1/k - D0/(k·L0)

存在性：临界点存在的充要条件是初始亏氧增速为正，即 k1·L0 > k2·D0
（等价地，一般式中对数真数为正且求得的 t_c > 0）。
若 k1·L0 <= k2·D0，曲线从起点单调复氧，不存在氧垂 —— 此时必须明确报告
「无临界点」，绝不能把公式算出的负 t_c 当成结果返回。

除解析式外，本模块还用手写的二分求根在 dD/dt 上数值定位 t_c，
用于交叉复核解析结果（不依赖任何第三方数值库）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import model
from .parameters import SagParams
from .validation import InvalidParameterError

#: 临界点存在性判定码
REASON_OK = "ok"
REASON_MONOTONIC_REAERATION = "monotonic_reaeration"  # k1·L0 <= k2·D0，全程单调复氧
REASON_NONPOSITIVE_LOG_ARGUMENT = "nonpositive_log_argument"  # 对数真数非正，无实数解

#: 求根数值容差 [day]
ROOT_TOL = 1e-10


@dataclass(frozen=True)
class CriticalPoint:
    """临界点的完整信息。distance = u · t_critical，时间经流速折算成河程。"""

    t_critical: float  # 临界时刻 [day]
    distance: float  # 临界河程 x_c = U · t_c [km]
    deficit: float  # 临界亏氧 D_c [mg/L]
    do: float  # 临界溶解氧（全程最低 DO）[mg/L]
    special_case: bool  # 是否走了 k1 == k2 特解分支


@dataclass(frozen=True)
class CriticalPointResult:
    """临界点求解结果：要么给出临界点，要么明确报告无临界点。"""

    exists: bool
    point: CriticalPoint | None
    reason: str  # REASON_* 之一
    special_case: bool  # 是否命中 k1 == k2 特解（与是否存在临界点无关）


def critical_point(params: SagParams) -> CriticalPointResult:
    """解析求解临界点；不存在时 exists=False 且 point=None。"""
    k1, k2, l0, d0 = params.k1, params.k2, params.l0, params.d0
    special = k1 == k2

    # 存在性充要条件：初始时刻亏氧在增长（dD/dt|0 = k1·L0 - k2·D0 > 0）
    if k1 * l0 <= k2 * d0:
        return CriticalPointResult(
            exists=False,
            point=None,
            reason=REASON_MONOTONIC_REAERATION,
            special_case=special,
        )

    if special:
        # 特解分支：t_c = 1/k - D0/(k·L0)，此时必有 t_c > 0
        t_c = 1.0 / k1 - d0 / (k1 * l0)
    else:
        arg = (k2 / k1) * (1.0 - d0 * (k2 - k1) / (k1 * l0))
        if arg <= 0.0:
            return CriticalPointResult(
                exists=False,
                point=None,
                reason=REASON_NONPOSITIVE_LOG_ARGUMENT,
                special_case=False,
            )
        t_c = math.log(arg) / (k2 - k1)
        if t_c <= 0.0:  # 双保险：正常不会走到，绝不返回负河程
            return CriticalPointResult(
                exists=False,
                point=None,
                reason=REASON_MONOTONIC_REAERATION,
                special_case=False,
            )

    d_c = model.deficit(t_c, params)
    point = CriticalPoint(
        t_critical=t_c,
        distance=params.u * t_c,  # 时间 × 流速 = 河程
        deficit=d_c,
        do=params.csat - d_c,
        special_case=special,
    )
    return CriticalPointResult(exists=True, point=point, reason=REASON_OK, special_case=special)


def critical_point_numeric(
    params: SagParams,
    tol: float = ROOT_TOL,
    max_iter: int = 200,
) -> CriticalPointResult | None:
    """用二分求根在 dD/dt = 0 上数值定位临界点，与解析式交叉复核。

    返回 None 表示无法构造变号区间（即不存在临界点）。
    """
    t_c = _bisect_deficit_rate_root(params, tol=tol, max_iter=max_iter)
    if t_c is None:
        return None
    d_c = model.deficit(t_c, params)
    point = CriticalPoint(
        t_critical=t_c,
        distance=params.u * t_c,
        deficit=d_c,
        do=params.csat - d_c,
        special_case=params.k1 == params.k2,
    )
    return CriticalPointResult(
        exists=True, point=point, reason=REASON_OK, special_case=params.k1 == params.k2
    )


def _bisect_deficit_rate_root(
    params: SagParams,
    tol: float = ROOT_TOL,
    max_iter: int = 200,
) -> float | None:
    """手写二分法求 dD/dt = 0 的正根。

    仅当 f(0) > 0（初始耗氧强于复氧）时才存在正根；从 t=0 起成倍扩张
    上界直到 f(hi) < 0 形成变号区间，再二分收敛。f(0) <= 0 时曲线单调
    复氧，返回 None。
    """
    f = lambda t: model.deficit_rate(t, params)  # noqa: E731

    f_lo = f(0.0)
    if f_lo <= 0.0:
        return None

    lo = 0.0
    hi = 1.0 / params.k2  # 复氧特征时间起步
    f_hi = f(hi)
    expansions = 0
    while f_hi > 0.0:
        hi *= 2.0
        f_hi = f(hi)
        expansions += 1
        if expansions > 200:
            raise InvalidParameterError("二分求根无法为 dD/dt 构造变号区间")

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_mid == 0.0 or (hi - lo) <= 2.0 * tol:
            return mid
        if f_mid > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
