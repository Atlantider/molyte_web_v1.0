# 腾讯云快速部署指南

本指南帮助您快速在腾讯云上部署 Molyte Web 平台（www.molyte.xyz）。

## 📋 前置条件

- ✅ 腾讯云账号（已有）
- ✅ 域名：www.molyte.xyz
- 💰 账户余额：建议至少 500 元

## 🚀 快速部署步骤

### 第一步：购买云资源

#### 1. 购买 CVM 云服务器

访问：https://console.cloud.tencent.com/cvm/instance

**推荐配置**：
- 地域：广州/北京/上海（选一个）
- 实例：标准型 S5.MEDIUM4（2核4GB）
- 镜像：Ubuntu 22.04 LTS
- 系统盘：50GB SSD
- 带宽：5 Mbps
- 安全组：开放 22, 80, 443 端口

**记录**：公网 IP 地址 `_______________`

#### 2. 购买 TencentDB for PostgreSQL

访问：https://console.cloud.tencent.com/postgres

**推荐配置**：
- 地域：与 CVM 相同
- 规格：1核2GB
- 版本：PostgreSQL 14
- 存储：20GB

**记录**：
- 内网地址：`_______________`
- 端口：`5432`
- 用户名：`_______________`
- 密码：`_______________`

#### 3. 创建 COS 存储桶

访问：https://console.cloud.tencent.com/cos

**操作**：
1. 点击"创建存储桶"
2. 名称：`molyte-results`（会自动加上 APPID 后缀）
3. 地域：与 CVM 相同
4. 访问权限：私有读写

**记录**：
- 存储桶名称：`molyte-results-_______________`
- 地域：`ap-guangzhou`（或其他）

#### 4. 创建 API 密钥

访问：https://console.cloud.tencent.com/cam/capi

**操作**：
1. 点击"新建密钥"
2. 复制 SecretId 和 SecretKey

**记录**：
- SecretId：`_______________`
- SecretKey：`_______________`

### 第二步：配置 DNS 解析

访问：https://console.cloud.tencent.com/cns

**操作**：
1. 找到域名 `molyte.xyz`
2. 添加记录：
   - 主机记录：`www`
   - 记录类型：`A`
   - 记录值：`<CVM公网IP>`
   - TTL：`600`

等待 5-10 分钟生效。

### 第三步：部署后端和前端

#### 1. SSH 连接到 CVM

```bash
ssh ubuntu@<CVM公网IP>
```

#### 2. 一键部署脚本

复制以下脚本并执行：

```bash
# 更新系统
sudo apt update && sudo apt upgrade -y

# 安装基础软件
sudo apt install -y python3 python3-pip python3-venv nodejs npm nginx postgresql-client git vim curl wget

# 克隆代码
cd /opt
sudo git clone https://github.com/Atlantider/molyte_web_v1.0.git
sudo chown -R ubuntu:ubuntu molyte_web_v1.0
cd molyte_web_v1.0

# 部署后端
cd backend
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

# 创建 .env 文件
cat > .env << 'EOF'
DATABASE_URL=postgresql://用户名:密码@TencentDB地址:5432/molyte_db
SECRET_KEY=$(openssl rand -hex 32)
COS_SECRET_ID=your-secret-id
COS_SECRET_KEY=your-secret-key
COS_REGION=ap-guangzhou
COS_BUCKET=molyte-results-1234567890
ENVIRONMENT=production
DEBUG=False
EOF

# 编辑 .env 文件，填入真实的数据库和 COS 配置
vim .env

# 初始化数据库（如果需要）
# psql -h <TencentDB地址> -U <用户名> -d molyte_db < ../database/init_db.sql

# 测试后端
uvicorn app.main:app --host 0.0.0.0 --port 8000
# 按 Ctrl+C 停止

# 创建 systemd 服务
sudo tee /etc/systemd/system/molyte-backend.service > /dev/null << 'EOF'
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
EOF

# 启动后端服务
sudo systemctl daemon-reload
sudo systemctl enable molyte-backend
sudo systemctl start molyte-backend
sudo systemctl status molyte-backend

# 部署前端
cd /opt/molyte_web_v1.0/frontend
npm install

# 修改 API 地址
sed -i 's|http://localhost:8000/api/v1|https://www.molyte.xyz/api/v1|g' src/api/client.ts

# 构建
npm run build

# 配置 Nginx
sudo tee /etc/nginx/sites-available/molyte > /dev/null << 'EOF'
server {
    listen 80;
    server_name www.molyte.xyz molyte.xyz;

    location / {
        root /opt/molyte_web_v1.0/frontend/dist;
        try_files $uri $uri/ /index.html;
    }

    location /api {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location /docs {
        proxy_pass http://localhost:8000;
    }
}
EOF

# 启用站点
sudo ln -s /etc/nginx/sites-available/molyte /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx

# 配置 HTTPS
sudo apt install -y certbot python3-certbot-nginx
sudo certbot --nginx -d www.molyte.xyz -d molyte.xyz

echo "✅ 部署完成！访问 https://www.molyte.xyz"
```

### 第四步：配置校园网 Worker

在校园网服务器上：

#### 1. 编辑配置文件

```bash
cd /public/home/xiaoji/molyte_web

# 编辑腾讯云配置
vim deployment/polling_worker_config_tencent.yaml
```

填入以下信息：

```yaml
api:
  base_url: "https://www.molyte.xyz/api/v1"
  worker_token: "从腾讯云后端获取的 Token"

cos:
  secret_id: "你的 SecretId"
  secret_key: "你的 SecretKey"
  region: "ap-guangzhou"  # 或其他地域
  bucket: "molyte-results-1234567890"  # 你的存储桶名称
```

#### 2. 安装依赖

```bash
source /public/software/anaconda3/bin/activate molyte
pip install cos-python-sdk-v5 pyyaml
```

#### 3. 启动 Worker

```bash
bash deployment/start_polling_worker.sh tencent
```

#### 4. 查看日志

```bash
tail -f /tmp/polling_worker.log
```

### 第五步：测试

#### 1. 访问网站

浏览器打开：https://www.molyte.xyz

#### 2. 注册用户

第一个注册的用户会自动成为管理员。

#### 3. 获取 Worker Token

临时方案：
1. 登录系统
2. 按 F12 打开开发者工具
3. Application -> Local Storage -> token
4. 复制 Token 值
5. 填入校园网 Worker 配置文件

#### 4. 提交测试任务

- 选择阴阳离子（如 Li + TFSI）
- 选择溶剂（如 DME）
- 提交任务

#### 5. 检查任务状态

- 前端：查看任务列表，状态应该从 PENDING → PROCESSING → RUNNING → COMPLETED
- 校园网：`tail -f /tmp/polling_worker.log` 查看 Worker 日志

## 📊 验证清单

- [ ] CVM 服务器可以访问
- [ ] 数据库连接正常
- [ ] COS 存储桶创建成功
- [ ] DNS 解析生效（`ping www.molyte.xyz`）
- [ ] 后端服务运行正常（`sudo systemctl status molyte-backend`）
- [ ] Nginx 运行正常（`sudo systemctl status nginx`）
- [ ] HTTPS 证书配置成功
- [ ] 前端页面可以访问
- [ ] 用户可以注册登录
- [ ] 校园网 Worker 启动成功
- [ ] 可以提交任务
- [ ] 任务可以正常执行

## 🔧 常见问题

### 1. 无法访问网站

```bash
# 检查 Nginx
sudo systemctl status nginx
sudo nginx -t

# 检查防火墙
sudo ufw status
```

### 2. 后端 500 错误

```bash
# 查看后端日志
sudo journalctl -u molyte-backend -n 50

# 检查数据库连接
psql -h <TencentDB地址> -U <用户名> -d molyte_db
```

### 3. Worker 无法连接

```bash
# 测试网络
curl https://www.molyte.xyz/api/v1/

# 检查配置
cat deployment/polling_worker_config_tencent.yaml
```

### 4. COS 上传失败

```bash
# 测试 COS 连接
python3 << 'EOF'
from qcloud_cos import CosConfig, CosS3Client

config = CosConfig(
    Region='ap-guangzhou',
    SecretId='your-secret-id',
    SecretKey='your-secret-key'
)
client = CosS3Client(config)
print(client.list_buckets())
EOF
```

## 📞 获取帮助

如果遇到问题：
1. 查看日志文件
2. 检查配置文件
3. 参考完整部署指南：`deployment/TENCENT_CLOUD_DEPLOYMENT_GUIDE.md`

## 🎉 完成

恭喜！您已经成功部署了 Molyte Web 平台！

现在您可以：
- 通过 https://www.molyte.xyz 访问平台
- 提交 MD 和 QC 计算任务
- 在校园网集群上运行计算
- 在腾讯云上存储和管理结果

