# 🚀 Molyte Web 快速参考

## 📦 当前版本

**版本**: `v1.0.0-hybrid-polling`  
**分支**: `deploy/hybrid`  
**部署方式**: 混合云部署（校园网计算 + 腾讯云前端）

---

## 🌿 分支管理

### 查看所有分支
```bash
git branch -a
```

### 切换分支
```bash
# 切换到混合部署分支（当前推荐）
git checkout deploy/hybrid

# 切换到主开发分支
git checkout main

# 切换到校园网部署分支（未来）
git checkout deploy/campus

# 切换到腾讯云部署分支（未来）
git checkout deploy/tencent
```

### 查看当前版本
```bash
# 查看当前分支
git branch

# 查看当前标签
git describe --tags

# 查看所有标签
git tag -l

# 查看标签详细信息
git show v1.0.0-hybrid-polling
```

---

## 🔄 代码同步

### 使用同步脚本（推荐）

```bash
# 查看帮助
bash deployment/sync_code.sh

# 检查状态
bash deployment/sync_code.sh status

# 推送到 GitHub
bash deployment/sync_code.sh push

# 从 GitHub 拉取
bash deployment/sync_code.sh pull

# 部署到腾讯云
bash deployment/sync_code.sh deploy

# 完整同步（本地 → GitHub → 腾讯云）
bash deployment/sync_code.sh sync
```

### 手动同步

```bash
# 1. 在校园网集群修改代码
vim backend/app/api/v1/worker.py

# 2. 提交到 GitHub
git add .
git commit -m "Fix: xxx"
git push origin deploy/hybrid

# 3. 在腾讯云服务器同步
ssh root@42.193.107.130
cd /root/molyte_web
git checkout deploy/hybrid
git pull origin deploy/hybrid
sudo systemctl restart molyte-backend
```

---

## 🏷️ 版本发布

### 创建新版本

```bash
# 1. 确保在正确的分支
git checkout deploy/hybrid

# 2. 创建标签
git tag -a v1.1.0-hybrid-polling -m "版本说明"

# 3. 推送标签
git push origin v1.1.0-hybrid-polling

# 4. 推送分支
git push origin deploy/hybrid
```

### 版本号规范

```
v{major}.{minor}.{patch}-{deployment}-{description}

示例:
v1.0.0-hybrid-polling       # 混合部署轮询版本
v1.1.0-hybrid-websocket     # 混合部署 WebSocket 版本
v1.0.1-hybrid-hotfix        # 热修复版本
v2.0.0-production           # 生产版本
```

---

## 🚀 部署流程

### 校园网集群部署 Worker

```bash
# 1. 克隆仓库（首次）
git clone https://github.com/Atlantider/molyte_web_v1.0.git
cd molyte_web_v1.0

# 2. 切换到混合部署分支
git checkout deploy/hybrid

# 3. 配置环境
cp deployment/polling_worker_config_tencent.yaml.template \
   deployment/polling_worker_config_tencent.yaml
vim deployment/polling_worker_config_tencent.yaml

# 4. 启动 Worker
bash deployment/start_polling_worker.sh tencent

# 5. 检查状态
bash deployment/check_worker_status.sh

# 6. 查看日志
tail -f /tmp/polling_worker_stdout.log
```

### 腾讯云服务器部署

```bash
# 1. 克隆仓库（首次）
git clone https://github.com/Atlantider/molyte_web_v1.0.git
cd molyte_web_v1.0

# 2. 切换到混合部署分支
git checkout deploy/hybrid

# 3. 配置环境
cp backend/.env.template backend/.env
vim backend/.env

# 4. 安装依赖
cd backend && pip install -r requirements.txt
cd ../frontend && npm install

# 5. 启动服务
bash scripts/start_all.sh

# 6. 检查服务
curl http://localhost:8000/api/v1/health
curl http://localhost:3000
```

---

## 🔧 常用命令

### Worker 管理

```bash
# 启动 Worker
bash deployment/start_polling_worker.sh tencent

# 停止 Worker
bash deployment/stop_polling_worker.sh

# 重启 Worker
bash deployment/stop_polling_worker.sh && \
sleep 2 && \
bash deployment/start_polling_worker.sh tencent

# 查看 Worker 状态
ps aux | grep polling_worker

# 查看 Worker 日志
tail -f /tmp/polling_worker_stdout.log
tail -f /tmp/polling_worker_stderr.log
```

### 服务管理（腾讯云）

```bash
# 启动所有服务
bash scripts/start_all.sh

# 停止所有服务
bash scripts/stop_all.sh

# 重启后端
sudo systemctl restart molyte-backend

# 查看后端日志
sudo journalctl -u molyte-backend -f

# 查看前端日志
pm2 logs molyte-frontend
```

---

## 🔍 故障排查

### Worker 无法连接到腾讯云

```bash
# 1. 检查网络连接
ping www.molyte.xyz
curl -k https://www.molyte.xyz/api/v1/health

# 2. 检查 Worker Token
# 编辑配置文件，确认 worker_token 正确
vim deployment/polling_worker_config_tencent.yaml

# 3. 检查 SSL 证书
# 如果证书过期，设置 verify_ssl: false
vim deployment/polling_worker_config_tencent.yaml

# 4. 查看详细日志
tail -100 /tmp/polling_worker_stdout.log
```

### 腾讯云服务无响应

```bash
# 1. 检查服务状态
sudo systemctl status molyte-backend
pm2 status

# 2. 检查端口占用
netstat -tunlp | grep 8000
netstat -tunlp | grep 3000

# 3. 重启服务
sudo systemctl restart molyte-backend
pm2 restart molyte-frontend

# 4. 查看日志
sudo journalctl -u molyte-backend -n 100
pm2 logs molyte-frontend --lines 100
```

### Git 冲突解决

```bash
# 1. 查看冲突文件
git status

# 2. 手动解决冲突
vim <冲突文件>

# 3. 标记为已解决
git add <冲突文件>

# 4. 完成合并
git commit -m "Resolve merge conflict"

# 5. 推送
git push origin deploy/hybrid
```

---

## 📊 监控和日志

### Worker 监控

```bash
# 实时查看 Worker 日志
tail -f /tmp/polling_worker_stdout.log

# 查看最近 100 行日志
tail -100 /tmp/polling_worker_stdout.log

# 搜索错误日志
grep -i error /tmp/polling_worker_stdout.log

# 查看 Worker 进程
ps aux | grep polling_worker
```

### 服务监控（腾讯云）

```bash
# 后端健康检查
curl https://www.molyte.xyz/api/v1/health

# 查看系统资源
htop
df -h
free -h

# 查看数据库连接
psql -U molyte -d molyte -c "SELECT count(*) FROM pg_stat_activity;"
```

---

## 🔒 安全检查

### 检查敏感文件

```bash
# 确保配置文件不在 Git 中
git status | grep -E "(config|secret|\.env)"

# 检查 .gitignore 是否生效
git check-ignore -v deployment/polling_worker_config_tencent.yaml

# 扫描仓库中的密钥（需要安装 gitleaks）
gitleaks detect --source . --verbose
```

### 更新密钥

```bash
# 1. 在腾讯云控制台撤销旧密钥
# 2. 生成新密钥
# 3. 更新配置文件
vim deployment/polling_worker_config_tencent.yaml

# 4. 重启 Worker
bash deployment/stop_polling_worker.sh
bash deployment/start_polling_worker.sh tencent
```

---

## 📚 相关文档

- [版本管理策略](./VERSION_MANAGEMENT.md) - 完整的版本管理指南
- [校园网部署指南](./CAMPUS_WORKER_DEPLOYMENT.md) - Worker 部署详细步骤
- [安全配置指南](./SECURITY_GUIDE.md) - 密钥管理和安全最佳实践
- [文件上传策略](./FILE_UPLOAD_STRATEGY.md) - 智能文件上传说明
- [代码同步工具](./sync_code.sh) - 自动化同步脚本

---

## 🆘 紧急联系

如果遇到紧急问题：

1. **查看日志**: `tail -f /tmp/polling_worker_stdout.log`
2. **检查状态**: `bash deployment/sync_code.sh status`
3. **重启服务**: `bash deployment/stop_polling_worker.sh && bash deployment/start_polling_worker.sh tencent`
4. **回滚版本**: `git checkout v1.0.0-hybrid-polling`

---

**最后更新**: 2024-12-01  
**维护者**: xiaoji  
**当前版本**: v1.0.0-hybrid-polling

