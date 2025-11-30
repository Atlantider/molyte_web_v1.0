# 混合云部署检查清单

## 📋 腾讯云部署检查清单

### 1. 服务器基础配置

- [ ] 服务器已购买并运行（42.193.107.130）
- [ ] 域名已解析到服务器 IP（www.molyte.xyz → 42.193.107.130）
- [ ] 安全组已配置（开放 80, 443, 22 端口）
- [ ] SSH 可以正常连接

### 2. 代码部署

- [ ] Git 已安装
- [ ] 代码已克隆：`git clone https://github.com/Atlantider/molyte_web_v1.0.git`
- [ ] 代码已更新到最新版本：`git pull origin main`

### 3. 后端部署

- [ ] Python 3.7+ 已安装
- [ ] 虚拟环境已创建：`python3 -m venv venv`
- [ ] 依赖已安装：`pip install -r backend/requirements.txt`
- [ ] PostgreSQL 数据库已创建
- [ ] 数据库连接配置正确（`backend/.env`）
- [ ] 数据库表已创建：`alembic upgrade head`
- [ ] Worker 用户已创建：`python deployment/create_worker_user.py`
- [ ] Worker Token 已记录
- [ ] 后端服务已启动（端口 8000）
- [ ] 后端 API 可访问：`curl http://localhost:8000/api/v1/docs`

### 4. 前端部署

- [ ] Node.js 16+ 已安装
- [ ] 依赖已安装：`cd frontend && npm install`
- [ ] 环境变量已配置（`frontend/.env.production`）
- [ ] 前端已构建：`npm run build`
- [ ] Nginx 已安装并配置
- [ ] 前端文件已部署到 Nginx 目录
- [ ] Nginx 已重启：`sudo systemctl restart nginx`

### 5. SSL 证书配置

- [ ] Certbot 已安装
- [ ] SSL 证书已申请：`certbot --nginx -d www.molyte.xyz`
- [ ] 证书自动续期已配置
- [ ] HTTPS 可正常访问：`curl https://www.molyte.xyz`

### 6. 腾讯云 COS 配置

- [ ] COS Bucket 已创建（molyte-results-1308567295）
- [ ] 访问密钥已创建（SecretId, SecretKey）
- [ ] Bucket 权限已配置（私有读写）
- [ ] CORS 已配置（允许前端访问）

### 7. 服务管理

- [ ] 后端服务已配置为 systemd 服务或 PM2 管理
- [ ] 服务已设置开机自启动
- [ ] 日志轮转已配置

---

## 📋 校园网集群部署检查清单

### 1. 环境准备

- [ ] Conda 环境已创建（molyte）
- [ ] Python 依赖已安装：
  - [ ] requests
  - [ ] pyyaml
  - [ ] cos-python-sdk-v5
- [ ] 代码已克隆或更新

### 2. 软件配置

- [ ] LAMMPS 已安装并可用
- [ ] Gaussian 已安装并可用（如需 QC 计算）
- [ ] Slurm 集群可正常使用
- [ ] 工作目录已创建：
  - [ ] `/public/home/xiaoji/molyte_web/data/md_work`
  - [ ] `/public/home/xiaoji/molyte_web/data/qc_work`

### 3. 网络连接测试

- [ ] 可以 ping 通腾讯云：`ping www.molyte.xyz`
- [ ] 可以访问 HTTPS：`curl -k https://www.molyte.xyz`
- [ ] 可以访问 API：`curl -k https://www.molyte.xyz/api/v1/workers/jobs/pending?job_type=MD&limit=10 -H "Authorization: Bearer TOKEN"`

### 4. Worker 配置

- [ ] 配置文件已更新：`deployment/polling_worker_config_tencent.yaml`
  - [ ] API URL 正确
  - [ ] Worker Token 正确
  - [ ] COS 密钥正确
  - [ ] 本地路径正确
  - [ ] Conda 环境名称正确
- [ ] 配置文件权限已设置（避免密钥泄露）：`chmod 600 deployment/polling_worker_config_tencent.yaml`

### 5. Worker 启动

- [ ] Worker 已启动：`bash deployment/start_polling_worker.sh tencent`
- [ ] Worker 进程正在运行：`ps aux | grep polling_worker`
- [ ] 日志无错误：`tail -f /tmp/polling_worker_stdout.log`
- [ ] Worker 可以成功轮询任务（日志中看到 "获取到 X 个待处理的任务"）

### 6. 开机自启动（可选）

- [ ] systemd 服务文件已创建：`/etc/systemd/system/molyte-worker.service`
- [ ] 服务已启用：`sudo systemctl enable molyte-worker`
- [ ] 服务已启动：`sudo systemctl start molyte-worker`
- [ ] 服务状态正常：`sudo systemctl status molyte-worker`

---

## 📋 端到端测试检查清单

### 1. 用户注册和登录

- [ ] 可以访问 https://www.molyte.xyz
- [ ] 可以注册新用户
- [ ] 可以登录
- [ ] 可以查看个人信息

### 2. MD 任务提交

- [ ] 可以选择阳离子
- [ ] 可以选择阴离子
- [ ] 可以选择溶剂
- [ ] 可以设置浓度
- [ ] 可以提交任务
- [ ] 任务状态显示为 "PENDING"

### 3. Worker 处理任务

- [ ] Worker 日志显示获取到任务
- [ ] Worker 创建工作目录
- [ ] Worker 提交 Slurm 任务
- [ ] Slurm 任务正在运行：`squeue -u xiaoji`
- [ ] 任务状态更新为 "RUNNING"

### 4. 任务完成

- [ ] Slurm 任务完成
- [ ] Worker 检测到任务完成
- [ ] Worker 提取最后一帧
- [ ] Worker 上传结果到 COS
- [ ] 任务状态更新为 "COMPLETED"

### 5. 结果查看

- [ ] 在 Web 界面可以看到任务完成
- [ ] 可以下载结果文件
- [ ] 可以查看分析图表
- [ ] 可以查看最后一帧结构

---

## 🐛 常见问题快速检查

### 问题：Worker 无法连接到腾讯云

**检查步骤**：
```bash
# 1. 网络连通性
ping www.molyte.xyz

# 2. HTTPS 连接
curl -k https://www.molyte.xyz

# 3. API 可用性
curl -k https://www.molyte.xyz/api/v1/workers/jobs/pending?job_type=MD&limit=10 \
  -H "Authorization: Bearer YOUR_TOKEN"
```

### 问题：Worker 认证失败（401）

**检查步骤**：
1. 确认 Worker 用户已创建
2. 确认 Token 正确
3. 确认 Token 未过期

**解决方案**：
```bash
# 在腾讯云服务器上重新生成 Token
cd /path/to/molyte_web
source venv/bin/activate
python deployment/create_worker_user.py
```

### 问题：Worker 权限不足（403）

**检查步骤**：
确认 Worker 用户的 role 是 ADMIN

**解决方案**：
```python
# 在腾讯云服务器上
from backend.app.core.database import SessionLocal
from backend.app.models.user import User, UserRole

db = SessionLocal()
worker = db.query(User).filter(User.username == "worker").first()
print(f"Current role: {worker.role}")
worker.role = UserRole.ADMIN
db.commit()
print("✅ Worker role updated to ADMIN")
db.close()
```

### 问题：SSL 证书验证失败

**临时解决方案**：
在配置文件中设置 `verify_ssl: false`

**永久解决方案**：
在腾讯云服务器上配置有效的 SSL 证书（Let's Encrypt）

### 问题：任务一直处于 PENDING 状态

**检查步骤**：
1. Worker 是否正在运行
2. Worker 日志是否有错误
3. Worker 是否能成功轮询到任务

**解决方案**：
```bash
# 查看 Worker 状态
ps aux | grep polling_worker

# 查看日志
tail -50 /tmp/polling_worker_stdout.log

# 重启 Worker
bash deployment/stop_polling_worker.sh
bash deployment/start_polling_worker.sh tencent
```

---

## 📊 部署完成验证

全部检查完成后，执行以下验证：

```bash
# 1. 在浏览器访问
https://www.molyte.xyz

# 2. 提交一个测试任务

# 3. 在校园网集群查看 Worker 日志
tail -f /tmp/polling_worker_stdout.log

# 4. 查看 Slurm 任务
squeue -u xiaoji

# 5. 等待任务完成（约 5-10 分钟）

# 6. 在 Web 界面查看结果
```

如果以上步骤全部成功，恭喜！🎉 混合云架构部署完成！

