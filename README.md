# Streeter–Phelps 氧垂分析计算服务

河流溶解氧沿程计算内核：把排污口下游的亏氧曲线扫出来，定位溶解氧最低的
临界位置。仅经 HTTP 提供，无页面。

## 模型

碳质 BOD 单级耗氧 + 大气复氧的教科书闭式解：

- `L(t) = L0·exp(−k1·t)`
- `dD/dt = k1·L − k2·D`
- 一般情形 `k1 ≠ k2`：`D(t) = k1·L0/(k2−k1)·(exp(−k1·t) − exp(−k2·t)) + D0·exp(−k2·t)`
- 特解 `k1 = k2 = k`：`D(t) = (k·L0·t + D0)·exp(−k·t)`
- `DO = Csat − D`

河程与时间经 **t = x / U** 换算（U 为流速），临界时刻：

- 一般：`t_c = ln[(k2/k1)·(1 − D0·(k2−k1)/(k1·L0))] / (k2 − k1)`
- 特解 `k1 = k2`：`t_c = 1/k − D0/(k·L0)`
- 存在条件：`k1·L0 > k2·D0`（否则单调复氧，服务明确报告「无临界点」）

## 单位约定

| 参数 | 含义 | 单位 |
| --- | --- | --- |
| `k1` | 耗氧系数 | 1/day |
| `k2` | 复氧系数 | 1/day |
| `u` | 流速 | km/day（如 0.35 m/s ≈ 30.24 km/day） |
| `l0` / `d0` / `csat` | 初始 BOD / 初始亏氧 / 饱和溶解氧 | mg/L |

`k1, k2, u, csat` 必须为正；`l0, d0` 非负且 `d0 ≤ csat`；非法一律返回 400。

## 模块划分

```
streeter_phelps/
  validation.py      输入校验（InvalidParameterError）
  parameters.py      工况参数对象 SagParams（构造即校验）
  model.py           亏氧 / BOD / DO 闭式解求值（含 k1=k2 特解）
  critical_point.py  临界点解析定位 + 手写二分求根复核
  scan.py            沿程扫描（t = x/U 换算，网格最低点 + 临界点）
  sweep.py           单系数区间批量扫参与趋势汇报
  scenarios.py       预置参考工况
  schemas.py         接口请求模型
  api.py             FastAPI 接口层
tests/               pytest 自动化测试
```

求根（二分法）与区间扫描均为标准库手写，不依赖 numpy/scipy。

## HTTP 接口

- `GET  /health` 健康检查
- `GET  /scenario/reference` 预置参考工况及其临界点（启动后可直接核对氧垂）
- `POST /scan` 沿程扫描：`{"params": {...}, "x_max": 150.0, "n_points": 201}`
  （`x_max` 缺省时按工况自动取，保证覆盖谷底）
- `POST /critical-point` 临界点定位：`{"k1":..., "k2":..., "u":..., "l0":..., "d0":..., "csat":...}`
- `POST /sweep` 区间扫参：`{"base": {...}, "parameter": "k2", "start": 0.4, "stop": 0.8, "n": 9}`

无临界点时响应 `exists=false` 与原因码（`monotonic_reaeration` 等），
绝不返回负河程。

## 运行

```bash
docker build -t oxygen-sag .
docker run --rm -p 8000:8000 oxygen-sag
# 服务监听固定端口 8000
curl http://localhost:8000/scenario/reference
```

## 测试

```bash
# 容器内执行
docker run --rm oxygen-sag pytest -q
# 或本地
pip install -r requirements.txt
pytest -q
```

测试覆盖：复氧系数增大氧垂变浅、初始 BOD 加倍最大亏氧升高、k1=k2 特解
（含与一般式的连续性极限）、流速加倍临界距离按比例拉长、无临界点判定、
非法参数报错，以及全部 HTTP 接口。
