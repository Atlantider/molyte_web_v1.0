# 阿里云部署指南

本指南将帮助您在阿里云上部署 Molyte Web 平台。

## 前置准备

### 1. 阿里云账号
- ✅ 已有阿里云账号
- 确保账号已完成实名认证
- 建议充值至少 500 元

### 2. 域名
- ✅ 已有域名：www.molyte.xyz
- 确保域名已解析到阿里云 ECS 公网 IP
- 如果服务器在中国大陆，需要完成域名备案

## 资源购买

### 1. ECS 云服务器

**推荐配置**：
- 实例规格：ecs.c6.xlarge（4核8GB）
- 操作系统：Ubuntu 22.04 LTS
- 系统盘：40GB SSD
- 带宽：5 Mbps（按量付费）
- 地域：根据用户分布选择（如华东2-上海）

**购买步骤**：
1. 登录阿里云控制台
2. 进入 ECS 产品页面
3. 点击"创建实例"
4. 选择配置并完成购买
5. 记录公网 IP 地址

### 2. RDS PostgreSQL 数据库

**推荐配置**：
- 实例规格：rds.pg.s2.large（2核4GB）
- 数据库版本：PostgreSQL 14
- 存储空间：20GB SSD
- 地域：与 ECS 相同

**购买步骤**：
1. 进入 RDS 产品页面
2. 选择 PostgreSQL 引擎
3. 选择配置并完成购买
4. 创建数据库 `molyte_db`
5. 创建用户并授权
6. 记录内网连接地址

### 3. OSS 对象存储

**配置**：
- 存储类型：标准存储
- 读写权限：私有
- 地域：与 ECS 相同

**购买步骤**：
1. 进入 OSS 产品页面
2. 创建 Bucket（如 `molyte-results`）
3. 设置访问权限为"私有"
4. 创建 AccessKey（用于 API 访问）
5. 记录 Endpoint 和 AccessKey

### 4. 安全组配置

在 ECS 安全组中开放以下端口：
- 22（SSH）
- 80（HTTP）
- 443（HTTPS）
- 8000（后端 API，可选，建议通过 Nginx 代理）

## 部署步骤

### 阶段 1：ECS 服务器初始化

#### 1.1 连接到 ECS

```bash
ssh root@<ECS_PUBLIC_IP>
```

#### 1.2 更新系统

```bash
apt update
apt upgrade -y
```

#### 1.3 安装基础软件

```bash
# 安装 Python 3.9+
apt install -y python3.9 python3-pip python3-venv

# 安装 Node.js 18
curl -fsSL https://deb.nodesource.com/setup_18.x | bash -
apt install -y nodejs

# 安装 Nginx
apt install -y nginx

# 安装 PostgreSQL 客户端
apt install -y postgresql-client

# 安装 Git
apt install -y git

# 安装其他工具
apt install -y vim curl wget htop
```

### 阶段 2：部署后端

#### 2.1 克隆代码

```bash
cd /opt
git clone https://github.com/Atlantider/molyte_web_v1.0.git
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

```bash
cp .env.example .env
vim .env
```

编辑 `.env` 文件：

```bash
# 数据库配置
DATABASE_URL=postgresql://用户名:密码@RDS内网地址:5432/molyte_db

# JWT 密钥（生成随机字符串）
SECRET_KEY=your-secret-key-here-change-in-production

# OSS 配置（可选，如果使用 OSS）
OSS_ENDPOINT=oss-cn-shanghai.aliyuncs.com
OSS_ACCESS_KEY_ID=your-access-key-id
OSS_ACCESS_KEY_SECRET=your-access-key-secret
OSS_BUCKET_NAME=molyte-results

# Redis 配置（如果使用）
REDIS_HOST=localhost
REDIS_PORT=6379
```

#### 2.5 初始化数据库

```bash
# 连接到 RDS
psql -h <RDS_HOST> -U <USERNAME> -d molyte_db

# 执行初始化脚本
\i ../database/init_db.sql

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

#### 2.7 使用 Systemd 管理后端服务

创建服务文件：

```bash
vim /etc/systemd/system/molyte-backend.service
```

内容：

```ini
[Unit]
Description=Molyte Web Backend API
After=network.target

[Service]
Type=simple
User=root
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
systemctl daemon-reload
systemctl enable molyte-backend
systemctl start molyte-backend
systemctl status molyte-backend
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
vim /etc/nginx/sites-available/molyte
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
    }

    # 后端 API 代理
    location /api {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # 文档
    location /docs {
        proxy_pass http://localhost:8000;
    }

    # WebSocket 支持（如果需要）
    location /ws {
        proxy_pass http://localhost:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

启用站点：

```bash
ln -s /etc/nginx/sites-available/molyte /etc/nginx/sites-enabled/
nginx -t
systemctl restart nginx
```

### 阶段 4：配置 HTTPS（推荐）

#### 4.1 安装 Certbot

```bash
apt install -y certbot python3-certbot-nginx
```

#### 4.2 获取 SSL 证书

```bash
certbot --nginx -d www.molyte.xyz -d molyte.xyz
```

按提示完成配置，选择自动重定向 HTTP 到 HTTPS。

### 阶段 5：创建 Worker 用户

在阿里云后端创建一个专用的 Worker 用户：

```bash
# 连接到数据库
psql -h <RDS_HOST> -U <USERNAME> -d molyte_db

# 创建 Worker 用户（需要在应用中实现）
# 或者使用管理员账号，生成 API Token
```

## 校园网 Worker 配置

### 1. 配置 Worker

在校园网服务器上：

```bash
cd /public/home/xiaoji/molyte_web

# 复制配置文件
cp deployment/polling_worker_config.yaml.example deployment/polling_worker_config.yaml

# 编辑配置
vim deployment/polling_worker_config.yaml
```

填入：
- 阿里云 API 地址：`https://your-domain.com/api/v1`
- Worker Token：从阿里云后端获取
- OSS 配置：AccessKey、Endpoint、Bucket

### 2. 启动 Worker

```bash
bash deployment/start_polling_worker.sh
```

### 3. 查看日志

```bash
tail -f /tmp/polling_worker.log
```

## 测试

### 1. 访问前端

浏览器访问：`https://your-domain.com`

### 2. 注册用户并登录

### 3. 提交测试任务

### 4. 检查 Worker 日志

确认 Worker 接收到任务并开始处理。

## 监控和维护

### 1. 查看服务状态

```bash
# 后端服务
systemctl status molyte-backend

# Nginx
systemctl status nginx

# Worker（在校园网）
ps aux | grep polling_worker
```

### 2. 查看日志

```bash
# 后端日志
journalctl -u molyte-backend -f

# Nginx 日志
tail -f /var/log/nginx/access.log
tail -f /var/log/nginx/error.log

# Worker 日志
tail -f /tmp/polling_worker.log
```

### 3. 数据库备份

RDS 会自动备份，也可以手动备份：

```bash
pg_dump -h <RDS_HOST> -U <USERNAME> molyte_db > backup_$(date +%Y%m%d).sql
```

## 成本优化

1. **使用按量付费**：初期用户少时使用按量付费
2. **购买资源包**：用户增长后购买资源包更划算
3. **OSS 生命周期**：设置 OSS 生命周期规则，自动删除旧文件
4. **CDN 加速**：用户量大时启用 CDN

## 故障排查

### 后端无法启动

```bash
# 查看日志
journalctl -u molyte-backend -n 50

# 检查端口占用
netstat -tlnp | grep 8000

# 检查数据库连接
psql -h <RDS_HOST> -U <USERNAME> -d molyte_db
```

### Worker 无法连接

```bash
# 测试网络连接
curl https://your-domain.com/api/v1/

# 检查 Token
# 在配置文件中确认 Token 正确
```

### OSS 上传失败

```bash
# 检查 AccessKey 权限
# 确认 Bucket 存在
# 检查网络连接
```

## 下一步

- 配置监控告警
- 设置自动备份
- 优化性能
- 添加日志分析

