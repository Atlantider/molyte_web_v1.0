# Molyte Web 部署文档

本目录包含 Molyte Web 平台的混合云部署方案和相关配置文件。

现在是腾讯云+校园网的混合架构，只有校园网能够访问腾讯云，腾讯云无法访问校园网。所以方法deployment里面。现在你是在腾讯云架构上面，我们要拓展一些新的功能，然后我通过git来与校园网同步，每次修改都要记录到git，然后推送，网络不好推送一次不成功就先进入下一步，而不是持续推送。我的Nginx要运行到molyte.xyz看效果。
腾讯云服务器：运行后端 API + 本地 PostgreSQL 数据库
腾讯云 COS：存储文件（轨迹、结果等）
校园网集群：运行 polling worker + Slurm 计算任务

## 📚 文档索引

### 快速开始

- **[腾讯云快速部署指南](QUICK_START_TENCENT.md)** ⭐ 推荐
  - 适用于域名 www.molyte.xyz
  - 包含完整的一键部署脚本
  - 适合快速上线

### 完整部署指南

- **[腾讯云部署指南](TENCENT_CLOUD_DEPLOYMENT_GUIDE.md)**
  - 详细的腾讯云部署步骤
  - 包含资源购买、配置、监控等完整流程
  - 成本估算和优化建议

- **[阿里云部署指南](ALIYUN_DEPLOYMENT_GUIDE.md)**
  - 阿里云部署方案
  - 适用于已有阿里云资源的用户

- **[混合云架构设计](../docs/DEPLOYMENT_HYBRID_CLOUD.md)**
  - 混合云架构原理
  - 轮询方案设计
  - 数据流和安全考虑

- **[文件上传策略说明](FILE_UPLOAD_STRATEGY.md)** 🔥 重要
  - 智能选择性上传策略
  - 轨迹文件处理方案
  - 存储成本优化（节省 95%）

## 🏗️ 架构概述

```
┌─────────────────────────────────────────────────────────────┐
│                    云端（腾讯云/阿里云）                        │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐      │
│  │   前端 Web   │  │  后端 API    │  │  PostgreSQL  │      │
│  │   (Nginx)    │  │  (FastAPI)   │  │   数据库     │      │
│  └──────────────┘  └──────────────┘  └──────────────┘      │
│                           │                                  │
│                    ┌──────────────┐                         │
│                    │  COS/OSS     │                         │
│                    │  对象存储    │                         │
│                    └──────────────┘                         │
└─────────────────────────────────────────────────────────────┘
                            ↕ HTTPS API
┌─────────────────────────────────────────────────────────────┐
│                    校园网计算集群                             │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐      │
│  │ 轮询 Worker  │→ │    Slurm     │→ │   计算节点   │      │
│  │              │  │   调度器     │  │ LAMMPS/G16   │      │
│  └──────────────┘  └──────────────┘  └──────────────┘      │
└─────────────────────────────────────────────────────────────┘
```

**工作流程**：
1. 用户通过云端 Web 界面提交任务
2. 任务存储在云端数据库（状态：PENDING）
3. 校园网 Worker 定期轮询云端 API，获取待处理任务
4. Worker 下载任务数据，生成输入文件，提交到 Slurm
5. 计算完成后，Worker 上传结果到 COS/OSS
6. Worker 更新云端数据库中的任务状态（COMPLETED）
7. 用户在 Web 界面查看结果

## 📁 文件说明

### 配置文件

- **`polling_worker_config.yaml`** - 阿里云 OSS 配置模板
- **`polling_worker_config_tencent.yaml`** - 腾讯云 COS 配置模板 ⭐

### 脚本文件

- **`polling_worker.py`** - 轮询 Worker 主程序
  - 支持腾讯云 COS 和阿里云 OSS
  - 自动检测配置文件中的云服务商
  - 处理 MD 和 QC 任务

- **`start_polling_worker.sh`** - Worker 启动脚本
  - 用法：`bash start_polling_worker.sh tencent` （腾讯云）
  - 用法：`bash start_polling_worker.sh` （阿里云，默认）
  - 自动检测和安装依赖

- **`stop_polling_worker.sh`** - Worker 停止脚本

## 🚀 快速部署

### 云端部署（腾讯云）

1. **购买资源**
   - CVM 云服务器（2核4GB）
   - TencentDB PostgreSQL（1核2GB）
   - COS 对象存储

2. **部署应用**
   ```bash
   # SSH 连接到 CVM
   ssh ubuntu@<CVM_IP>
   
   # 克隆代码
   cd /opt
   sudo git clone https://github.com/Atlantider/molyte_web_v1.0.git
   
   # 按照 QUICK_START_TENCENT.md 中的脚本执行
   ```

3. **配置 DNS**
   - 将 www.molyte.xyz 解析到 CVM 公网 IP

4. **配置 HTTPS**
   ```bash
   sudo certbot --nginx -d www.molyte.xyz
   ```

### 校园网 Worker 部署

1. **编辑配置文件**
   ```bash
   cd /public/home/xiaoji/molyte_web
   vim deployment/polling_worker_config_tencent.yaml
   ```
   
   填入：
   - API 地址：`https://www.molyte.xyz/api/v1`
   - Worker Token（从云端获取）
   - COS 配置（SecretId、SecretKey、Bucket）

2. **安装依赖**
   ```bash
   source /public/software/anaconda3/bin/activate molyte
   pip install cos-python-sdk-v5 pyyaml
   ```

3. **启动 Worker**
   ```bash
   bash deployment/start_polling_worker.sh tencent
   ```

4. **查看日志**
   ```bash
   tail -f /tmp/polling_worker.log
   ```

## 🔧 配置说明

### API 配置

```yaml
api:
  base_url: "https://www.molyte.xyz/api/v1"  # 云端 API 地址
  worker_token: "your-token"                  # Worker 认证 Token
  poll_interval: 30                           # 轮询间隔（秒）
  timeout: 60                                 # 请求超时（秒）
```

### 腾讯云 COS 配置

```yaml
cos:
  secret_id: "AKIDxxxxx"                      # API 密钥 ID
  secret_key: "xxxxx"                         # API 密钥 Key
  region: "ap-guangzhou"                      # 地域
  bucket: "molyte-results-1234567890"         # 存储桶名称
  result_prefix: "results/"                   # 文件前缀
```

### 阿里云 OSS 配置

```yaml
oss:
  endpoint: "oss-cn-shanghai.aliyuncs.com"    # OSS Endpoint
  access_key_id: "xxxxx"                      # AccessKey ID
  access_key_secret: "xxxxx"                  # AccessKey Secret
  bucket_name: "molyte-results"               # Bucket 名称
  result_prefix: "results/"                   # 文件前缀
```

## 📊 监控和维护

### 查看服务状态

**云端（CVM）**：
```bash
# 后端服务
sudo systemctl status molyte-backend

# Nginx
sudo systemctl status nginx

# 查看日志
sudo journalctl -u molyte-backend -f
```

**校园网（Worker）**：
```bash
# 检查进程
ps aux | grep polling_worker

# 查看日志
tail -f /tmp/polling_worker.log

# 重启 Worker
bash deployment/stop_polling_worker.sh
bash deployment/start_polling_worker.sh tencent
```

### 常用命令

```bash
# 重启后端
sudo systemctl restart molyte-backend

# 重启 Nginx
sudo systemctl restart nginx

# 查看 Worker 状态
ps aux | grep polling_worker

# 停止 Worker
bash deployment/stop_polling_worker.sh

# 启动 Worker（腾讯云）
bash deployment/start_polling_worker.sh tencent

# 启动 Worker（阿里云）
bash deployment/start_polling_worker.sh
```

## 💰 成本估算

### 腾讯云（约 300-500 元/月）

- CVM 2核4GB：150-200 元/月
- TencentDB 1核2GB：100-150 元/月
- COS 存储+流量：20-50 元/月
- 带宽 5Mbps：50-80 元/月

### 优化建议

1. 使用包年优惠（6-7 折）
2. 设置 COS 生命周期规则，自动删除旧文件
3. 初期使用按量付费，稳定后改为包年
4. 用户量大时启用 CDN

## 🔒 安全建议

1. **定期更新系统**
   ```bash
   sudo apt update && sudo apt upgrade
   ```

2. **配置防火墙**
   ```bash
   sudo ufw allow 22,80,443/tcp
   sudo ufw enable
   ```

3. **定期备份数据库**
   ```bash
   pg_dump -h <DB_HOST> -U <USER> molyte_db > backup.sql
   ```

4. **使用强密码**
   - 数据库密码
   - API Token
   - COS/OSS 密钥

5. **启用 HTTPS**
   - 使用 Let's Encrypt 免费证书
   - 自动续期

## 📞 故障排查

### 问题 1：Worker 无法连接云端

**检查**：
```bash
# 测试网络
curl https://www.molyte.xyz/api/v1/

# 检查 Token
grep worker_token deployment/polling_worker_config_tencent.yaml
```

### 问题 2：COS/OSS 上传失败

**检查**：
```bash
# 测试 COS 连接
python3 -c "
from qcloud_cos import CosConfig, CosS3Client
config = CosConfig(Region='ap-guangzhou', SecretId='xxx', SecretKey='xxx')
client = CosS3Client(config)
print(client.list_buckets())
"
```

### 问题 3：任务一直 PENDING

**检查**：
```bash
# Worker 是否运行
ps aux | grep polling_worker

# 查看 Worker 日志
tail -50 /tmp/polling_worker.log

# 检查配置
cat deployment/polling_worker_config_tencent.yaml
```

## 📖 更多资源

- [GitHub 仓库](https://github.com/Atlantider/molyte_web_v1.0)
- [腾讯云文档](https://cloud.tencent.com/document)
- [阿里云文档](https://help.aliyun.com/)

## 🎯 下一步

部署完成后，您可以：

1. ✅ 访问 https://www.molyte.xyz
2. ✅ 注册用户并登录
3. ✅ 提交 MD/QC 计算任务
4. ✅ 在校园网集群上运行计算
5. ✅ 查看和下载结果

祝您使用愉快！🎉

