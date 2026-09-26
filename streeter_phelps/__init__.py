"""Streeter–Phelps 河流溶解氧氧垂分析计算内核。

模块划分：
- validation    输入参数校验
- parameters    工况参数对象
- model         亏氧 / BOD / 溶解氧闭式解求值
- critical_point 临界点（最大亏氧）解析定位 + 手写二分求根复核
- scan          沿程扫描（t = x / U 时间-河程换算）
- sweep         单系数区间批量扫参
- scenarios     预置参考工况
- api           FastAPI 接口层
"""
