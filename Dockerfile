FROM python:3.12-slim

WORKDIR /app

# 仅安装轻量 Web 依赖；求根与区间扫描全部为标准库手写，不引入数值计算库
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY streeter_phelps ./streeter_phelps
COPY tests ./tests
COPY pytest.ini .

EXPOSE 8000

# 固定端口 8000 对外提供 HTTP 计算服务
CMD ["uvicorn", "streeter_phelps.api:app", "--host", "0.0.0.0", "--port", "8000"]
