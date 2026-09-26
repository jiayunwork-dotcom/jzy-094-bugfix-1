"""FastAPI 接口层：仅提供 HTTP 计算接口，无页面。

接口一览：
    GET  /health            健康检查
    GET  /scenario/reference 预置参考工况及其临界点（启动后可直接核对氧垂）
    POST /scan              单一工况沿程扫描（亏氧曲线 + 临界点）
    POST /critical-point    单一工况临界点解析定位
    POST /sweep             单系数区间批量扫参

参数非法一律返回 400；未捕获异常返回 500。
"""

from __future__ import annotations

from dataclasses import asdict

from fastapi import FastAPI, HTTPException

from .critical_point import critical_point, critical_point_numeric
from .parameters import SagParams
from .scan import scan_along_river
from .scenarios import REFERENCE_SCENARIO
from .schemas import ParamsIn, ScanRequest, SweepRequest
from .sweep import sweep_parameter
from .validation import InvalidParameterError

app = FastAPI(
    title="Streeter–Phelps 氧垂分析计算服务",
    description="河流溶解氧沿程计算内核：亏氧闭式解、沿程扫描、临界点定位、区间扫参。",
    version="1.0.0",
)


def _build_params(p: ParamsIn) -> SagParams:
    return SagParams(k1=p.k1, k2=p.k2, u=p.u, l0=p.l0, d0=p.d0, csat=p.csat)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/scenario/reference")
def reference_scenario() -> dict:
    """预置参考工况 + 其临界点，便于启动后立即核对氧垂是否正常出现。"""
    result = critical_point(REFERENCE_SCENARIO)
    return {
        "params": asdict(REFERENCE_SCENARIO),
        "critical_point": asdict(result.point) if result.point else None,
        "exists": result.exists,
        "reason": result.reason,
    }


@app.post("/scan")
def scan(req: ScanRequest) -> dict:
    """单一工况沿程扫描：返回亏氧/溶解氧曲线与临界点（解析 + 二分复核）。"""
    try:
        params = _build_params(req.params)
        result = scan_along_river(params, x_max=req.x_max, n_points=req.n_points)
    except InvalidParameterError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {
        "params": asdict(params),
        "x_max": result.x_max,
        "n_points": len(result.points),
        "points": [asdict(p) for p in result.points],
        "grid_min": asdict(result.grid_min),
        "analytic_critical": (
            asdict(result.analytic_critical) if result.analytic_critical else None
        ),
        "refined_critical": (
            asdict(result.refined_critical) if result.refined_critical else None
        ),
        "has_critical_point": result.analytic_critical is not None,
    }


@app.post("/critical-point")
def critical_point_endpoint(req: ParamsIn) -> dict:
    """单一工况临界点定位；无氧垂时明确报告 exists=False，绝不返回负河程。"""
    try:
        params = _build_params(req)
        analytic = critical_point(params)
        numeric = critical_point_numeric(params)
    except InvalidParameterError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {
        "params": asdict(params),
        "exists": analytic.exists,
        "reason": analytic.reason,
        "special_case": analytic.special_case,
        "critical_point": asdict(analytic.point) if analytic.point else None,
        "numeric_cross_check": asdict(numeric.point) if numeric and numeric.point else None,
    }


@app.post("/sweep")
def sweep(req: SweepRequest) -> dict:
    """单系数区间批量扫参：汇报临界亏氧随该系数的变化趋势。"""
    try:
        params = _build_params(req.base)
        summary = sweep_parameter(
            params, req.parameter, req.start, req.stop, req.n
        )
    except InvalidParameterError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {
        "parameter": summary.parameter,
        "trend": summary.trend,
        "records": [asdict(r) for r in summary.records],
    }
