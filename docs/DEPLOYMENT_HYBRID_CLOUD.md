# 混合云部署架构文档

## 架构概述

本项目采用**混合云轮询架构**，将 Web 服务部署在阿里云（公网可访问），计算任务在校园网集群执行。

```
┌─────────────────────────────────────────────────────────────────┐
│                    阿里云（公网可访问）                            │
│                                                                   │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐              │
│  │   前端       │  │  后端 API    │  │  PostgreSQL │              │
│  │  (Nginx)    │  │  (FastAPI)  │  │   (RDS)     │              │
│  └─────────────┘  └─────────────┘  └─────────────┘              │
│                                                                   │
│  ┌─────────────┐  ┌─────────────┐                               │
│  │   Redis     │  │  OSS 存储    │                               │
│  │  (可选)     │  │  (结果文件)  │                               │
│  └─────────────┘  └─────────────┘                               │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
                            ↑ HTTPS API 轮询
                            │
┌─────────────────────────────────────────────────────────────────┐
│                      校园网集群（内网）                            │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │  轮询 Worker (polling_worker.py)                         │    │
│  │  - 每 30 秒轮询阿里云 API                                 │    │
│  │  - 获取 PENDING 状态的任务                                │    │
│  │  - 下载输入文件                                           │    │
│  │  - 提交到 Slurm                                          │    │
│  │  - 监控任务状态                                           │    │
│  │  - 上传结果到 OSS                                        │    │
│  │  - 更新任务状态                                           │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │  Slurm 集群                                              │    │
│  │  - LAMMPS (MD 模拟)                                      │    │
│  │  - Gaussian (QC 计算)                                    │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

## 数据流

### 1. 用户提交任务
```
用户 → 阿里云前端 → 阿里云后端 API → PostgreSQL (状态: PENDING)
```

### 2. 任务执行
```
校园网 Worker 轮询 → 发现 PENDING 任务 → 下载输入数据 → 
生成 LAMMPS/Gaussian 输入文件 → 提交到 Slurm → 
监控任务状态 → 任务完成
```

### 3. 结果上传
```
校园网 Worker → 上传结果文件到 OSS → 更新数据库状态 (COMPLETED) →
用户从前端查看结果 → 从 OSS 下载文件
```

## 阿里云资源清单

### 必需资源

1. **ECS 云服务器**
   - 规格：2核4GB（起步）/ 4核8GB（推荐）
   - 系统：Ubuntu 22.04 / CentOS 7
   - 带宽：5 Mbps（按量付费）
   - 用途：运行前端 + 后端 API

2. **RDS PostgreSQL**
   - 规格：1核2GB（起步）/ 2核4GB（推荐）
   - 版本：PostgreSQL 14
   - 存储：20GB SSD
   - 用途：存储用户数据、任务元数据

3. **OSS 对象存储**
   - 存储类型：标准存储
   - 容量：按需（预估 100GB/月起）
   - 用途：存储计算结果、可视化图表

4. **域名 + SSL 证书**
   - 域名：example.com
   - SSL：免费证书（Let's Encrypt）或阿里云 SSL

### 可选资源

5. **Redis（可选）**
   - 规格：256MB（起步）
   - 用途：缓存、会话管理
   - 备注：也可以在 ECS 上自建

6. **CDN（可选）**
   - 用途：加速前端静态资源
   - 备注：用户量大时启用

## 成本估算（按月）

| 资源 | 规格 | 价格（元/月） |
|------|------|--------------|
| ECS | 2核4GB | ~100 |
| RDS PostgreSQL | 1核2GB | ~200 |
| OSS | 100GB 存储 + 流量 | ~50 |
| 带宽 | 5 Mbps | ~150 |
| 域名 | .com | ~60/年 |
| SSL 证书 | 免费 | 0 |
| **总计** | | **~500/月** |

*注：价格仅供参考，实际以阿里云官网为准*

## 部署步骤

### 阶段 1：阿里云部署

1. **购买资源**
   - ECS、RDS、OSS
   - 配置安全组（开放 80、443、22 端口）

2. **部署后端**
   ```bash
   # 安装依赖
   sudo apt update
   sudo apt install python3.9 python3-pip postgresql-client nginx
   
   # 克隆代码
   git clone https://github.com/Atlantider/molyte_web_v1.0.git
   cd molyte_web_v1.0/backend
   
   # 安装 Python 依赖
   pip3 install -r requirements.txt
   
   # 配置环境变量
   cp .env.example .env
   # 编辑 .env，配置数据库连接、OSS 密钥等
   
   # 初始化数据库
   psql -h <RDS_HOST> -U <USERNAME> -d molyte_db < ../database/init_db.sql
   
   # 启动后端
   uvicorn app.main:app --host 0.0.0.0 --port 8000
   ```

3. **部署前端**
   ```bash
   cd ../frontend
   npm install
   npm run build
   
   # 配置 Nginx
   sudo cp dist/* /var/www/html/
   ```

4. **配置 Nginx**
   ```nginx
   server {
       listen 80;
       server_name your-domain.com;
       
       # 前端
       location / {
           root /var/www/html;
           try_files $uri $uri/ /index.html;
       }
       
       # 后端 API
       location /api {
           proxy_pass http://localhost:8000;
           proxy_set_header Host $host;
           proxy_set_header X-Real-IP $remote_addr;
       }
   }
   ```

### 阶段 2：校园网 Worker 部署

1. **安装依赖**
   ```bash
   cd /public/home/xiaoji/molyte_web
   source /public/software/anaconda3/bin/activate molyte
   pip install requests oss2
   ```

2. **配置 Worker**
   ```bash
   # 编辑 deployment/polling_worker_config.yaml
   # 填入阿里云 API 地址、OSS 密钥等
   ```

3. **启动 Worker**
   ```bash
   python deployment/polling_worker.py
   ```

4. **设置开机自启**
   ```bash
   # 使用 systemd 或 supervisor
   ```

## 安全考虑

1. **API 认证**
   - Worker 使用专用 API Token
   - Token 存储在环境变量中

2. **数据传输**
   - 全程 HTTPS 加密
   - OSS 使用签名 URL

3. **访问控制**
   - RDS 仅允许 ECS 访问
   - OSS 设置访问策略

4. **数据备份**
   - RDS 自动备份（7天）
   - OSS 版本控制

## 监控和运维

1. **日志**
   - 阿里云：使用云监控
   - 校园网：本地日志文件

2. **告警**
   - Worker 离线告警
   - 任务失败告警
   - 资源使用告警

3. **性能优化**
   - 数据库索引优化
   - OSS CDN 加速
   - Worker 并发控制

## 扩展性

### 水平扩展
- 增加多个 Worker 节点
- 使用任务锁避免重复执行

### 垂直扩展
- 升级 ECS/RDS 规格
- 增加 OSS 存储容量

## 故障恢复

1. **Worker 故障**
   - 任务状态回滚到 PENDING
   - 其他 Worker 自动接管

2. **网络中断**
   - Worker 自动重连
   - 任务状态持久化

3. **计算节点故障**
   - Slurm 自动重调度
   - Worker 监控任务状态

## 下一步

1. 实现轮询 Worker 脚本
2. 配置 OSS SDK
3. 修改后端 API 支持轮询模式
4. 测试完整流程

