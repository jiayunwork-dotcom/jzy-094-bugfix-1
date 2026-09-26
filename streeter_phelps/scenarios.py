"""预置参考工况：一个排污口下游应当出现明显氧垂的典型河流。

单位：k1/k2 [1/day]，U [km/day]，L0/D0/Csat [mg/L]。
该工况下 k1·L0 = 3.0 > k2·D0 = 0.6，初始耗氧明显强于复氧，
临界点解析存在且临界距离为正（约 88 km 量级，临界亏氧约 4.17 mg/L，
约为初始亏氧的 2.8 倍），服务启动后可直接扫出核对。
"""

from __future__ import annotations

from .parameters import SagParams

#: 参考工况：典型缓流河道，耗氧 0.2/day、复氧 0.4/day、流速 30 km/day
REFERENCE_SCENARIO = SagParams(
    k1=0.2,
    k2=0.4,
    u=30.0,
    l0=15.0,
    d0=1.5,
    csat=10.0,
)
