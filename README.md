# Molyte Web

电解质分子动力学模拟和量子化学计算的 Web 平台。

## 项目简介

Molyte Web 是一个用于电解质体系分子动力学（MD）模拟和量子化学（QC）计算的 Web 应用平台。支持：

- 🧪 电解质体系的 MD 模拟（基于 LAMMPS）
- ⚛️ 分子的量子化学计算（基于 Gaussian）
- 📊 任务管理和结果可视化
- 🔄 异步任务处理和调度

## 技术栈

### 后端
- **FastAPI**: Python Web 框架
- **Celery**: 分布式任务队列
- **Redis**: 消息代理
- **SQLite**: 数据库
- **LAMMPS**: 分子动力学模拟引擎
- **Gaussian**: 量子化学计算软件

### 前端
- **React 18**: UI 框架
- **TypeScript**: 类型安全
- **Vite**: 构建工具
- **Ant Design**: UI 组件库
- **ECharts**: 数据可视化
- **Zustand**: 状态管理

## 项目结构

```
molyte_web/
├── backend/              # 后端代码
│   ├── app/
│   │   ├── api/         # API 路由
│   │   ├── core/        # 核心配置
│   │   ├── models/      # 数据模型
│   │   ├── schemas/     # Pydantic schemas
│   │   ├── tasks/       # Celery 任务
│   │   └── utils/       # 工具函数
│   └── requirements.txt
├── frontend/            # 前端代码
│   ├── src/
│   │   ├── api/        # API 调用
│   │   ├── components/ # React 组件
│   │   ├── pages/      # 页面
│   │   ├── stores/     # 状态管理
│   │   └── types/      # TypeScript 类型
│   └── package.json
├── data/               # 数据目录
│   ├── md_work/       # MD 任务工作目录
│   ├── qc_work/       # QC 任务工作目录
│   ├── initial_salts/ # 盐类配置文件
│   └── charges/       # 电荷计算结果
├── start_complete.sh  # 启动所有服务
└── stop_complete.sh   # 停止所有服务
```

## 快速开始

### 环境要求

- Python 3.7+
- Node.js 16+
- Redis
- Conda 环境（推荐）

### 安装依赖

**后端依赖**:
```bash
conda activate molyte
cd backend
pip install -r requirements.txt
```

**前端依赖**:
```bash
cd frontend
npm install
```

### 启动服务

**启动所有服务**（推荐）:
```bash
bash start_complete.sh
```

这将启动：
- Redis（消息队列）
- Celery Worker（任务处理）
- Celery Beat（定时任务）
- FastAPI 后端（端口 8000）
- React 前端（端口 3000）

**停止所有服务**:
```bash
bash stop_complete.sh
```

### 访问应用

- **前端应用**: http://localhost:3000
- **后端 API**: http://localhost:8000
- **API 文档**: http://localhost:8000/docs

## 功能特性

### MD 模拟
- 支持多种盐类（Li, K, Mg, Na, FSI, DFOB, BF4, ClO4, NO3, PF6, TFSI, FBS）
- 自动生成 LAMMPS 输入文件
- 支持 Slurm 集群任务提交
- 实时任务状态监控

### QC 计算
- Gaussian 量子化学计算
- 电荷计算和优化
- 结果自动解析和存储

### 任务管理
- 异步任务处理
- 任务队列管理
- 定时任务调度
- 任务状态追踪

## 配置说明

主要配置文件：`backend/app/core/config.py`

关键配置项：
- `MOLYTE_WORK_BASE_PATH`: MD 任务工作目录
- `QC_WORK_BASE_PATH`: QC 任务工作目录
- `MOLYTE_INITIAL_SALTS_PATH`: 盐类配置文件路径
- `MOLYTE_CHARGE_SAVE_PATH`: 电荷计算结果路径

## 开发指南

### 后端开发

启动后端开发服务器：
```bash
cd backend
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

### 前端开发

启动前端开发服务器：
```bash
cd frontend
npm run dev
```

### 添加新的 API 端点

1. 在 `backend/app/api/` 中创建新的路由文件
2. 在 `backend/app/schemas/` 中定义请求/响应模型
3. 在 `backend/app/api/api.py` 中注册路由

### 添加新的 Celery 任务

1. 在 `backend/app/tasks/` 中创建任务文件
2. 使用 `@celery_app.task` 装饰器定义任务
3. 在 API 中调用 `task.delay()` 异步执行

## 日志文件

- Redis: `/public/home/xiaoji/software/redis/redis.log`
- Celery Worker: `/tmp/celery_worker.log`
- Celery Beat: `/tmp/celery_beat.log`
- Backend: `/tmp/backend.log`
- Frontend: `/tmp/frontend.log`

## 故障排查

### 前端无法启动
- 检查 Node.js 是否安装：`node --version`
- 检查依赖是否安装：`cd frontend && npm install`
- 查看日志：`tail -f /tmp/frontend.log`

### 后端无法启动
- 检查 Python 环境：`conda activate molyte`
- 检查依赖是否安装：`pip install -r backend/requirements.txt`
- 查看日志：`tail -f /tmp/backend.log`

### Redis 连接失败
- 检查 Redis 是否运行：`redis-cli ping`
- 启动 Redis：`redis-server redis.conf --daemonize yes`

### Celery 任务不执行
- 检查 Celery Worker 是否运行：`ps aux | grep celery`
- 查看 Celery 日志：`tail -f /tmp/celery_worker.log`

## 许可证

MIT License

## 联系方式

如有问题，请提交 Issue。

