# 腾讯云部署指南

本指南将帮助您在腾讯云上部署 Molyte Web 平台（www.molyte.xyz）。

## 前置准备

### 1. 腾讯云账号
- ✅ 已有腾讯云账号
- 完成实名认证
- 充值（建议至少 500 元）

### 2. 域名
- ✅ 已有域名：www.molyte.xyz
- 需要在腾讯云 DNS 解析中添加 A 记录指向服务器 IP
- 如果服务器在中国大陆，需要完成域名备案

## 资源购买

### 1. CVM 云服务器

**推荐配置**：
- 实例规格：标准型 S5.MEDIUM4（2核4GB）或 S5.LARGE8（4核8GB）
- 操作系统：Ubuntu 22.04 LTS
- 系统盘：50GB SSD 云硬盘
- 带宽：5 Mbps（按量付费）
- 地域：根据用户分布选择（推荐：北京、上海、广州）

**购买步骤**：
1. 登录腾讯云控制台：https://console.cloud.tencent.com/
2. 进入云服务器 CVM 产品页面
3. 点击"新建"创建实例
4. 选择配置并完成购买
5. 记录公网 IP 地址

### 2. TencentDB for PostgreSQL 数据库

**推荐配置**：
- 实例规格：1核2GB 或 2核4GB
- 数据库版本：PostgreSQL 14
- 存储空间：20GB SSD
- 地域：与 CVM 相同

**购买步骤**：
1. 进入云数据库 TencentDB 产品页面
2. 选择 PostgreSQL 引擎
3. 选择配置并完成购买
4. 创建数据库 `molyte_db`
5. 创建用户并授权
6. 记录内网连接地址（如：`xxx.postgres.tencentcdb.com:5432`）

### 3. COS 对象存储

**配置**：
- 存储类型：标准存储
- 访问权限：私有读写
- 地域：与 CVM 相同

**购买步骤**：
1. 进入对象存储 COS 产品页面
2. 创建存储桶（如 `molyte-results-1234567890`）
3. 设置访问权限为"私有读写"
4. 创建 API 密钥（SecretId 和 SecretKey）
   - 访问：https://console.cloud.tencent.com/cam/capi
   - 点击"新建密钥"
5. 记录存储桶名称、地域、SecretId 和 SecretKey

### 4. 安全组配置

在 CVM 安全组中开放以下端口：
- 22（SSH）
- 80（HTTP）
- 443（HTTPS）
- 8000（后端 API，可选，建议通过 Nginx 代理）

**配置步骤**：
1. 进入 CVM 控制台
2. 点击实例 ID 进入详情页
3. 点击"安全组"标签
4. 编辑入站规则，添加上述端口

## 部署步骤

### 阶段 1：CVM 服务器初始化

#### 1.1 连接到 CVM

```bash
ssh ubuntu@<CVM_PUBLIC_IP>
# 或者使用腾讯云控制台的"登录"按钮
```

#### 1.2 更新系统

```bash
sudo apt update
sudo apt upgrade -y
```

#### 1.3 安装基础软件

```bash
# 安装 Python 3.10+
sudo apt install -y python3 python3-pip python3-venv

# 安装 Node.js 18
curl -fsSL https://deb.nodesource.com/setup_18.x | sudo bash -
sudo apt install -y nodejs

# 安装 Nginx
sudo apt install -y nginx

# 安装 PostgreSQL 客户端
sudo apt install -y postgresql-client

# 安装 Git
sudo apt install -y git

# 安装其他工具
sudo apt install -y vim curl wget htop
```

### 阶段 2：部署后端

#### 2.1 克隆代码

```bash
cd /opt
sudo git clone https://github.com/Atlantider/molyte_web_v1.0.git
sudo chown -R ubuntu:ubuntu molyte_web_v1.0
cd molyte_web_v1.0
```

#### 2.2 创建 Python 虚拟环境

```bash
cd backend
python3 -m venv venv
source venv/bin/activate
```

#### 2.3 安装依赖

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

#### 2.4 配置环境变量

创建 `.env` 文件：

```bash
vim .env
```

编辑内容：

```bash
# 数据库配置
DATABASE_URL=postgresql://用户名:密码@TencentDB内网地址:5432/molyte_db

# JWT 密钥（生成随机字符串）
SECRET_KEY=your-secret-key-here-change-in-production-$(openssl rand -hex 32)

# COS 配置
COS_SECRET_ID=your-secret-id
COS_SECRET_KEY=your-secret-key
COS_REGION=ap-guangzhou  # 或 ap-beijing, ap-shanghai
COS_BUCKET=molyte-results-1234567890

# Redis 配置（如果使用）
REDIS_HOST=localhost
REDIS_PORT=6379

# 应用配置
ENVIRONMENT=production
DEBUG=False
```

#### 2.5 初始化数据库

```bash
# 连接到 TencentDB
psql -h <TencentDB_HOST> -U <USERNAME> -d molyte_db

# 如果有初始化脚本，执行它
# \i ../database/init_db.sql

# 或者使用 Alembic 迁移
# alembic upgrade head

# 退出
\q
```

#### 2.6 测试后端

```bash
# 启动后端（测试）
uvicorn app.main:app --host 0.0.0.0 --port 8000

# 在另一个终端测试
curl http://localhost:8000/
```

按 `Ctrl+C` 停止测试。

#### 2.7 使用 Systemd 管理后端服务

创建服务文件：

```bash
sudo vim /etc/systemd/system/molyte-backend.service
```

内容：

```ini
[Unit]
Description=Molyte Web Backend API
After=network.target

[Service]
Type=simple
User=ubuntu
WorkingDirectory=/opt/molyte_web_v1.0/backend
Environment="PATH=/opt/molyte_web_v1.0/backend/venv/bin"
ExecStart=/opt/molyte_web_v1.0/backend/venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8000
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

启动服务：

```bash
sudo systemctl daemon-reload
sudo systemctl enable molyte-backend
sudo systemctl start molyte-backend
sudo systemctl status molyte-backend
```

### 阶段 3：部署前端

#### 3.1 构建前端

```bash
cd /opt/molyte_web_v1.0/frontend

# 安装依赖
npm install

# 修改 API 地址
vim src/api/client.ts
# 将 baseURL 改为: https://www.molyte.xyz/api/v1

# 构建
npm run build
```

#### 3.2 配置 Nginx

```bash
sudo vim /etc/nginx/sites-available/molyte
```

内容：

```nginx
server {
    listen 80;
    server_name www.molyte.xyz molyte.xyz;

    # 前端静态文件
    location / {
        root /opt/molyte_web_v1.0/frontend/dist;
        try_files $uri $uri/ /index.html;
        
        # 缓存静态资源
        location ~* \.(js|css|png|jpg|jpeg|gif|ico|svg|woff|woff2|ttf|eot)$ {
            expires 1y;
            add_header Cache-Control "public, immutable";
        }
    }

    # 后端 API 代理
    location /api {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # 超时设置
        proxy_connect_timeout 60s;
        proxy_send_timeout 60s;
        proxy_read_timeout 60s;
    }

    # API 文档
    location /docs {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    # WebSocket 支持（如果需要）
    location /ws {
        proxy_pass http://localhost:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
    }
}
```

启用站点：

```bash
sudo ln -s /etc/nginx/sites-available/molyte /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

#### 3.3 配置 DNS 解析

1. 登录腾讯云控制台
2. 进入 DNS 解析 DNSPod
3. 找到域名 `molyte.xyz`
4. 添加 A 记录：
   - 主机记录：`www`
   - 记录类型：`A`
   - 记录值：`<CVM_PUBLIC_IP>`
   - TTL：`600`
5. 添加 A 记录（可选，用于裸域名）：
   - 主机记录：`@`
   - 记录类型：`A`
   - 记录值：`<CVM_PUBLIC_IP>`
   - TTL：`600`

等待 DNS 生效（通常 5-10 分钟）。

### 阶段 4：配置 HTTPS（推荐）

#### 4.1 安装 Certbot

```bash
sudo apt install -y certbot python3-certbot-nginx
```

#### 4.2 获取 SSL 证书

```bash
sudo certbot --nginx -d www.molyte.xyz -d molyte.xyz
```

按提示输入邮箱并同意服务条款。

#### 4.3 自动续期

Certbot 会自动配置续期，可以测试：

```bash
sudo certbot renew --dry-run
```

### 阶段 5：创建 Worker 用户和 API Token

#### 5.1 创建管理员用户

访问 `https://www.molyte.xyz`，注册第一个用户（会自动成为管理员）。

#### 5.2 生成 Worker API Token

在后端数据库中创建一个专用的 Worker 用户，或者使用管理员账号生成 API Token。

临时方案（使用管理员 Token）：
1. 登录系统
2. 打开浏览器开发者工具（F12）
3. 在 Application/Storage -> Local Storage 中找到 `token`
4. 复制 Token 值

## 校园网 Worker 配置

### 1. 编辑配置文件

在校园网服务器上：

```bash
cd /public/home/xiaoji/molyte_web

# 编辑配置
vim deployment/polling_worker_config.yaml
```

修改以下内容：

```yaml
# API 配置
api:
  base_url: "https://www.molyte.xyz/api/v1"  # 腾讯云后端地址
  worker_token: "your-worker-token-here"      # 从腾讯云后端获取的 Token
  poll_interval: 10                            # 轮询间隔（秒）
  timeout: 30                                  # 请求超时（秒）

# COS 配置（腾讯云对象存储）
cos:
  secret_id: "your-secret-id"                  # 腾讯云 API 密钥 ID
  secret_key: "your-secret-key"                # 腾讯云 API 密钥 Key
  region: "ap-guangzhou"                       # COS 地域
  bucket: "molyte-results-1234567890"          # COS 存储桶名称
  result_prefix: "results/"                    # 结果文件前缀

# 本地环境配置（校园网集群）
local:
  work_base_path: "/public/home/xiaoji/molyte_web/data/md_work"
  initial_salts_path: "/public/home/xiaoji/molyte_web/data/initial_salts"
  ligpargen_path: "/public/home/xiaoji/molyte_web/ligpargen"
  packmol_path: "/public/software/packmol-20.16.0/packmol"
  ltemplify_path: "/public/home/xiaoji/molyte_web/ltemplify.py"
  moltemplate_path: "/public/software/moltemplate-2.20.22/moltemplate/scripts/moltemplate.sh"
  charge_save_path: "/public/home/xiaoji/molyte_web/data/charge_save"

# Worker 配置
worker:
  name: "campus-worker-01"
  max_concurrent_jobs: 5
  heartbeat_interval: 60
  log_file: "/tmp/polling_worker.log"
  log_level: "INFO"

# 上传配置
upload:
  result_files:
    - "*.data"
    - "*.in.*"
    - "*.log"
    - "*.out"
    - "*.xyz"
    - "slurm-*.out"
  max_file_size: 100  # MB

# Slurm 配置
slurm:
  partition: "cpu"
  account: "your-account"
```

### 2. 安装 COS Python SDK

```bash
source /public/software/anaconda3/bin/activate molyte
pip install cos-python-sdk-v5
```

### 3. 修改 Worker 代码以支持腾讯云 COS

需要修改 `deployment/polling_worker.py` 中的 OSS 部分，改为使用腾讯云 COS SDK。

### 4. 启动 Worker

```bash
bash deployment/start_polling_worker.sh
```

### 5. 查看日志

```bash
tail -f /tmp/polling_worker.log
```

## 测试

### 1. 访问前端

浏览器访问：`https://www.molyte.xyz`

### 2. 注册用户并登录

### 3. 提交测试任务

选择阴阳离子、溶剂等参数，提交一个 MD 任务。

### 4. 检查 Worker 日志

在校园网服务器上：

```bash
tail -f /tmp/polling_worker.log
```

确认 Worker 接收到任务并开始处理。

### 5. 检查任务状态

在前端界面查看任务状态变化：
- PENDING → PROCESSING → RUNNING → COMPLETED

## 监控和维护

### 1. 查看服务状态

```bash
# 后端服务
sudo systemctl status molyte-backend

# Nginx
sudo systemctl status nginx

# Worker（在校园网）
ps aux | grep polling_worker
```

### 2. 查看日志

```bash
# 后端日志
sudo journalctl -u molyte-backend -f

# Nginx 访问日志
sudo tail -f /var/log/nginx/access.log

# Nginx 错误日志
sudo tail -f /var/log/nginx/error.log

# Worker 日志（校园网）
tail -f /tmp/polling_worker.log
```

### 3. 数据库备份

```bash
# 手动备份
pg_dump -h <TencentDB_HOST> -U <USERNAME> molyte_db > backup_$(date +%Y%m%d).sql

# 上传到 COS（可选）
# 使用 COSCMD 工具
```

### 4. 重启服务

```bash
# 重启后端
sudo systemctl restart molyte-backend

# 重启 Nginx
sudo systemctl restart nginx

# 重启 Worker（校园网）
bash deployment/stop_polling_worker.sh
bash deployment/start_polling_worker.sh
```

## 成本估算

### 月度成本（约 300-500 元/月）

1. **CVM 云服务器**：
   - 2核4GB：约 150-200 元/月
   - 4核8GB：约 300-400 元/月

2. **TencentDB PostgreSQL**：
   - 1核2GB：约 100-150 元/月

3. **COS 对象存储**：
   - 存储费用：0.118 元/GB/月
   - 流量费用：0.5 元/GB（外网下行）
   - 预计：20-50 元/月

4. **带宽**：
   - 5 Mbps：约 50-80 元/月

### 优化建议

1. **使用包年优惠**：CVM 和数据库包年可享受 6-7 折优惠
2. **COS 生命周期**：设置自动删除 30 天前的旧文件
3. **CDN 加速**：用户量大时启用 CDN，降低带宽成本
4. **按量付费**：初期使用按量付费，稳定后改为包年

## 故障排查

### 1. 后端无法启动

```bash
# 查看日志
sudo journalctl -u molyte-backend -n 50

# 检查端口占用
sudo netstat -tlnp | grep 8000

# 检查数据库连接
psql -h <TencentDB_HOST> -U <USERNAME> -d molyte_db
```

### 2. Worker 无法连接腾讯云

```bash
# 测试网络连接
curl https://www.molyte.xyz/api/v1/

# 检查 Token
# 确认配置文件中的 Token 正确

# 测试 DNS 解析
nslookup www.molyte.xyz
```

### 3. COS 上传失败

```bash
# 检查 SecretId 和 SecretKey
# 确认存储桶名称和地域正确
# 检查网络连接

# 测试 COS 连接
python3 -c "
from qcloud_cos import CosConfig, CosS3Client
config = CosConfig(Region='ap-guangzhou', SecretId='xxx', SecretKey='xxx')
client = CosS3Client(config)
print(client.list_buckets())
"
```

### 4. HTTPS 证书问题

```bash
# 检查证书状态
sudo certbot certificates

# 手动续期
sudo certbot renew

# 重新申请
sudo certbot --nginx -d www.molyte.xyz -d molyte.xyz --force-renewal
```

## 下一步

- [ ] 配置监控告警（腾讯云监控）
- [ ] 设置自动备份
- [ ] 优化性能（Redis 缓存、数据库索引）
- [ ] 添加日志分析（腾讯云 CLS）
- [ ] 配置 CDN 加速
- [ ] 设置定时任务清理旧数据

