# 🔒 安全配置指南

## ⚠️ 重要提醒

**永远不要将包含敏感信息的配置文件提交到 Git！**

敏感信息包括：
- API Token / JWT Token
- 数据库密码
- 云服务密钥（SecretId, SecretKey, Access Key 等）
- 私钥文件
- 任何认证凭据

## 📋 配置文件安全检查清单

### 1. 使用配置模板

所有包含敏感信息的配置文件都应该有对应的 `.template` 文件：

```bash
# ✅ 正确：提交模板文件
deployment/polling_worker_config_tencent.yaml.template

# ❌ 错误：提交真实配置
deployment/polling_worker_config_tencent.yaml
```

### 2. 配置 .gitignore

确保 `.gitignore` 包含以下规则：

```gitignore
# Configuration files with secrets
deployment/polling_worker_config_tencent.yaml
deployment/polling_worker_config_aliyun.yaml
deployment/polling_worker_config.yaml
backend/.env
backend/.env.production
frontend/.env.production
*_config_secret.yaml
*_secret.yaml
*.secret.yaml
```

### 3. 从本地模板创建配置文件

```bash
# 复制模板
cp deployment/polling_worker_config_tencent.yaml.template \
   deployment/polling_worker_config_tencent.yaml

# 编辑配置文件，填入真实的密钥
vim deployment/polling_worker_config_tencent.yaml

# 设置文件权限（仅所有者可读写）
chmod 600 deployment/polling_worker_config_tencent.yaml
```

### 4. 验证配置文件不会被提交

```bash
# 检查 Git 状态
git status

# 如果看到配置文件出现在 "Changes to be committed" 或 "Untracked files"
# 说明 .gitignore 配置有问题！

# 正确的输出应该是：
# nothing to commit, working tree clean
```

## 🚨 如果密钥已经泄露怎么办？

### 立即行动清单

#### 1. **立即撤销泄露的密钥**

**腾讯云 COS 密钥：**
```
1. 登录腾讯云控制台
2. 访问管理 -> API密钥管理
3. 禁用或删除泄露的密钥
4. 创建新的密钥对
```

**Worker Token：**
```bash
# 在腾讯云服务器上
cd /path/to/molyte_web
source venv/bin/activate

# 删除旧的 Worker 用户
python -c "
from backend.app.core.database import SessionLocal
from backend.app.models.user import User

db = SessionLocal()
worker = db.query(User).filter(User.username == 'worker').first()
if worker:
    db.delete(worker)
    db.commit()
    print('✅ Old worker user deleted')
db.close()
"

# 创建新的 Worker 用户
python deployment/create_worker_user.py
```

#### 2. **从 Git 历史中删除敏感文件**

**方法 1：使用 git filter-branch（彻底删除）**

```bash
# ⚠️ 警告：这会重写 Git 历史！

# 删除文件
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch deployment/polling_worker_config_tencent.yaml" \
  --prune-empty --tag-name-filter cat -- --all

# 清理引用
git for-each-ref --format="delete %(refname)" refs/original | git update-ref --stdin
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# 强制推送（覆盖远程仓库）
git push origin --force --all
git push origin --force --tags
```

**方法 2：使用 BFG Repo-Cleaner（推荐，更快）**

```bash
# 安装 BFG
# macOS: brew install bfg
# Linux: 下载 https://rtyley.github.io/bfg-repo-cleaner/

# 删除文件
bfg --delete-files polling_worker_config_tencent.yaml

# 清理
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# 强制推送
git push origin --force --all
```

#### 3. **通知 GitHub 删除缓存**

即使从 Git 历史中删除，GitHub 可能仍有缓存：

1. 访问 https://github.com/Atlantider/molyte_web_v1.0/settings
2. 滚动到底部，点击 "Delete this repository"（如果仓库是公开的且密钥已泄露）
3. 或者联系 GitHub Support 请求清除缓存

#### 4. **检查是否有未授权访问**

**腾讯云：**
```
1. 云审计 -> 操作记录
2. 检查是否有异常的 API 调用
3. 检查 COS Bucket 是否有异常访问
```

**服务器：**
```bash
# 检查登录日志
sudo grep "Accepted" /var/log/auth.log | tail -50

# 检查运行中的进程
ps aux | grep -v "$(whoami)"

# 检查网络连接
netstat -tunlp
```

## 🛡️ 最佳实践

### 1. 使用环境变量

**后端配置（backend/.env）：**
```bash
# 数据库
DATABASE_URL=postgresql://user:password@localhost/molyte

# JWT
SECRET_KEY=your-secret-key-here
ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30

# 腾讯云
TENCENT_SECRET_ID=your-secret-id
TENCENT_SECRET_KEY=your-secret-key
```

**在代码中读取：**
```python
import os
from dotenv import load_dotenv

load_dotenv()

SECRET_ID = os.getenv("TENCENT_SECRET_ID")
SECRET_KEY = os.getenv("TENCENT_SECRET_KEY")
```

### 2. 使用密钥管理服务

**腾讯云 KMS（密钥管理系统）：**
- 存储敏感配置
- 自动轮换密钥
- 访问审计

**HashiCorp Vault：**
- 开源密钥管理
- 动态密钥生成
- 细粒度访问控制

### 3. 定期轮换密钥

```bash
# 每 90 天轮换一次密钥
# 设置日历提醒或使用自动化脚本
```

### 4. 最小权限原则

**腾讯云 COS 权限：**
```json
{
  "version": "2.0",
  "statement": [
    {
      "effect": "allow",
      "action": [
        "name/cos:PutObject",
        "name/cos:GetObject"
      ],
      "resource": [
        "qcs::cos:ap-beijing:uid/1308567295:molyte-results-1308567295/*"
      ]
    }
  ]
}
```

只授予必要的权限，不要使用管理员密钥。

### 5. 监控和告警

**设置告警：**
- API 调用异常
- 流量异常
- 费用异常
- 登录异常

## 📝 配置文件管理流程

### 开发环境

```bash
# 1. 克隆仓库
git clone https://github.com/Atlantider/molyte_web_v1.0.git
cd molyte_web_v1.0

# 2. 从模板创建配置
cp deployment/polling_worker_config_tencent.yaml.template \
   deployment/polling_worker_config_tencent.yaml

# 3. 填写配置（使用开发环境密钥）
vim deployment/polling_worker_config_tencent.yaml

# 4. 验证不会被提交
git status  # 不应该看到 polling_worker_config_tencent.yaml
```

### 生产环境

```bash
# 1. 使用环境变量或密钥管理服务
# 2. 配置文件存储在安全位置（不在代码仓库中）
# 3. 使用配置管理工具（Ansible, Chef, Puppet）
```

## 🔍 安全审计

### 定期检查

```bash
# 1. 检查 Git 历史中是否有敏感信息
git log --all --full-history --source -- "*secret*" "*password*" "*token*"

# 2. 使用工具扫描
# 安装 gitleaks
brew install gitleaks  # macOS
# 或从 https://github.com/gitleaks/gitleaks/releases 下载

# 扫描仓库
gitleaks detect --source . --verbose

# 3. 检查 .gitignore 是否生效
git check-ignore -v deployment/polling_worker_config_tencent.yaml
# 应该输出：.gitignore:XX:pattern    deployment/polling_worker_config_tencent.yaml
```

## 📞 紧急联系方式

如果发现安全问题：

1. **立即撤销密钥**
2. **联系管理员**
3. **检查访问日志**
4. **评估影响范围**
5. **制定修复计划**

## ✅ 安全检查清单

部署前检查：

- [ ] 所有敏感配置文件都有 `.template` 版本
- [ ] `.gitignore` 包含所有敏感文件
- [ ] 真实配置文件不在 Git 仓库中
- [ ] 文件权限设置正确（600 或 400）
- [ ] 使用最小权限原则
- [ ] 密钥定期轮换
- [ ] 设置监控和告警
- [ ] 团队成员了解安全规范

---

**记住：安全是一个持续的过程，不是一次性的任务！**

