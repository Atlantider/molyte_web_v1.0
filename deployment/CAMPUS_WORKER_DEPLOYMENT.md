# 校园网集群 Worker 部署指南

本文档说明如何在校园网集群上部署轮询 Worker，实现混合云架构。

## 📋 架构概述

```
┌─────────────────────────────────────────────────────────────┐
│  腾讯云 (www.molyte.xyz)                                     │
│  ┌────────────────┐  ┌────────────────┐  ┌──────────────┐  │
│  │  前端 (React)  │  │  后端 (FastAPI)│  │  数据库 (PG) │  │
│  └────────────────┘  └────────────────┘  └──────────────┘  │
│           ▲                  ▲                              │
└───────────┼──────────────────┼──────────────────────────────┘
            │                  │
            │ HTTPS            │ HTTPS (轮询)
            │                  │
┌───────────┼──────────────────┼──────────────────────────────┐
│  校园网集群                   │                              │
│  ┌────────────────────────────┼──────────────────────────┐  │
│  │  Polling Worker            │                          │  │
│  │  - 每 30 秒轮询一次        │                          │  │
│  │  - 获取待处理任务          │                          │  │
│  │  - 提交到 Slurm            │                          │  │
│  │  - 上传结果到 COS          │                          │  │
│  └────────────────────────────┼──────────────────────────┘  │
│                               ▼                              │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  Slurm 集群                                            │ │
│  │  - LAMMPS (MD 计算)                                    │ │
│  │  - Gaussian (QC 计算)                                  │ │
│  └────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

## 🚀 部署步骤

### 步骤 1：在腾讯云上更新后端代码

```bash
# SSH 连接到腾讯云服务器
ssh root@42.193.107.130

# 进入项目目录
cd /path/to/molyte_web

# 拉取最新代码
git pull origin main

# 重启后端服务
# 如果使用 systemd
sudo systemctl restart molyte-backend

# 或者如果使用 PM2
pm2 restart molyte-backend

# 或者手动重启
pkill -f "uvicorn"
source venv/bin/activate
nohup uvicorn backend.app.main:app --host 0.0.0.0 --port 8000 &
```

### 步骤 2：验证 Worker API 可用

```bash
# 测试 Worker API（在校园网集群上执行）
curl -k -X GET "https://www.molyte.xyz/api/v1/workers/jobs/pending?job_type=MD&limit=10" \
  -H "Authorization: Bearer YOUR_WORKER_TOKEN" \
  -H "Content-Type: application/json"

# 应该返回空数组 [] 或任务列表，而不是 401/403 错误
```

### 步骤 3：在校园网集群上启动 Worker

```bash
# 进入项目目录
cd /public/home/xiaoji/molyte_web

# 启动 Worker
bash deployment/start_polling_worker.sh tencent

# 查看日志
tail -f /tmp/polling_worker_stdout.log
```

### 步骤 4：验证 Worker 正常运行

```bash
# 检查 Worker 进程
ps aux | grep polling_worker

# 查看最新日志
tail -30 /tmp/polling_worker_stdout.log

# 应该看到类似的日志：
# 2025-11-30 23:29:31 | INFO | Worker 'campus-worker-01' 已启动
# 2025-11-30 23:29:31 | INFO | 开始轮询...
# 2025-11-30 23:29:31 | INFO | 获取到 0 个待处理的 MD 任务
```

## 🔧 配置说明

### 配置文件位置

`deployment/polling_worker_config_tencent.yaml`

### 关键配置项

```yaml
api:
  base_url: "https://www.molyte.xyz/api/v1"
  worker_token: "YOUR_WORKER_TOKEN"  # 从腾讯云后端获取
  poll_interval: 30                   # 轮询间隔（秒）
  verify_ssl: false                   # SSL 验证（生产环境应为 true）

cos:
  secret_id: "YOUR_SECRET_ID"         # 腾讯云 COS 密钥
  secret_key: "YOUR_SECRET_KEY"
  region: "ap-beijing"
  bucket: "molyte-results-1308567295"

local:
  md_work_dir: "/public/home/xiaoji/molyte_web/data/md_work"
  qc_work_dir: "/public/home/xiaoji/molyte_web/data/qc_work"
  lammps_bin: "/public/software/lammps/bin/lmp"
  gaussian_bin: "/public/software/gaussian/g16"
  conda_env: "molyte"

worker:
  name: "campus-worker-01"
  max_concurrent_jobs: 3
  heartbeat_interval: 300
```

## 🐛 故障排查

### 问题 1：SSL 证书验证失败

**错误信息**：
```
SSLError: certificate verify failed: certificate has expired
```

**解决方案**：
在配置文件中设置 `verify_ssl: false`（仅用于开发/测试）

### 问题 2：401 认证失败

**错误信息**：
```
Could not validate credentials
```

**解决方案**：
1. 在腾讯云服务器上运行 `python deployment/create_worker_user.py`
2. 复制生成的 Token 到配置文件的 `api.worker_token`
3. 重启 Worker

### 问题 3：403 权限不足

**错误信息**：
```
Only admin/worker users can fetch pending jobs
```

**解决方案**：
确保 Worker 用户的 `role` 是 `ADMIN`：

```python
# 在腾讯云服务器上
from backend.app.core.database import SessionLocal
from backend.app.models.user import User, UserRole

db = SessionLocal()
worker = db.query(User).filter(User.username == "worker").first()
worker.role = UserRole.ADMIN
db.commit()
db.close()
```

### 问题 4：Worker 无法连接到腾讯云

**检查步骤**：
```bash
# 1. 测试网络连通性
ping www.molyte.xyz

# 2. 测试 HTTPS 连接
curl -k https://www.molyte.xyz/

# 3. 测试 API 端点
curl -k https://www.molyte.xyz/api/v1/workers/jobs/pending?job_type=MD&limit=10 \
  -H "Authorization: Bearer YOUR_TOKEN"
```

## 📊 监控和维护

### 查看 Worker 状态

```bash
# 查看进程
ps aux | grep polling_worker

# 查看日志
tail -f /tmp/polling_worker_stdout.log

# 查看最近的错误
grep ERROR /tmp/polling_worker_stdout.log | tail -20
```

### 重启 Worker

```bash
# 停止
bash deployment/stop_polling_worker.sh

# 启动
bash deployment/start_polling_worker.sh tencent
```

### 设置开机自启动

创建 systemd 服务文件 `/etc/systemd/system/molyte-worker.service`：

```ini
[Unit]
Description=Molyte Polling Worker
After=network.target

[Service]
Type=forking
User=xiaoji
WorkingDirectory=/public/home/xiaoji/molyte_web
ExecStart=/bin/bash /public/home/xiaoji/molyte_web/deployment/start_polling_worker.sh tencent
ExecStop=/bin/bash /public/home/xiaoji/molyte_web/deployment/stop_polling_worker.sh
Restart=on-failure
RestartSec=10

[Install]
WantedBy=multi-user.target
```

启用服务：

```bash
sudo systemctl daemon-reload
sudo systemctl enable molyte-worker
sudo systemctl start molyte-worker
sudo systemctl status molyte-worker
```

## ✅ 验证端到端流程

### 1. 在 Web 界面提交任务

访问 https://www.molyte.xyz，登录后提交一个 MD 任务。

### 2. 观察 Worker 日志

```bash
tail -f /tmp/polling_worker_stdout.log
```

应该看到：
```
INFO | 获取到 1 个待处理的 MD 任务
INFO | 开始处理 MD 任务 #123
INFO | 创建工作目录: /public/home/xiaoji/molyte_web/data/md_work/job_123
INFO | 提交 Slurm 任务: sbatch job_123.sh
INFO | Slurm 任务 ID: 456789
```

### 3. 检查 Slurm 任务

```bash
squeue -u xiaoji
```

### 4. 等待任务完成

Worker 会自动：
1. 监控 Slurm 任务状态
2. 任务完成后提取最后一帧
3. 上传结果到腾讯云 COS
4. 更新任务状态为 COMPLETED

### 5. 在 Web 界面查看结果

刷新任务列表，应该看到任务状态变为"已完成"，可以下载结果文件。

## 📝 注意事项

1. **SSL 证书**：生产环境应该配置有效的 SSL 证书，并设置 `verify_ssl: true`
2. **Token 安全**：Worker Token 应该妥善保管，不要提交到公开仓库
3. **资源限制**：根据集群负载调整 `max_concurrent_jobs`
4. **日志轮转**：定期清理日志文件，避免占用过多磁盘空间
5. **网络稳定性**：确保校园网到腾讯云的网络连接稳定

## 🎯 下一步

- [ ] 配置 SSL 证书（Let's Encrypt）
- [ ] 设置日志轮转
- [ ] 配置监控告警
- [ ] 优化轮询间隔
- [ ] 添加任务优先级支持

