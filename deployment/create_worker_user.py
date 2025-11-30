#!/usr/bin/env python3
"""
创建 Worker 用户并生成 Token

用法:
    python deployment/create_worker_user.py
"""

import sys
import os

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from backend.app.core.database import SessionLocal
from backend.app.models.user import User
from backend.app.core.security import get_password_hash, create_access_token
from sqlalchemy.orm import Session


def create_worker_user(db: Session):
    """创建 Worker 用户"""
    
    # 检查是否已存在 worker 用户
    existing_user = db.query(User).filter(User.username == "worker").first()
    
    if existing_user:
        print("✅ Worker 用户已存在")
        print(f"   用户名: {existing_user.username}")
        print(f"   邮箱: {existing_user.email}")
        print(f"   是否管理员: {existing_user.is_admin}")
        
        # 生成新 Token
        token = create_access_token(data={"sub": existing_user.username})
        print(f"\n🔑 Worker Token:")
        print(f"   {token}")
        print(f"\n请将此 Token 复制到 deployment/polling_worker_config_tencent.yaml 中的 api.worker_token")
        return token
    
    # 创建新的 worker 用户
    print("创建新的 Worker 用户...")
    
    worker_user = User(
        username="worker",
        email="worker@molyte.local",
        hashed_password=get_password_hash("worker_password_change_me"),
        is_admin=True,  # Worker 需要管理员权限来访问所有任务
        is_active=True
    )
    
    db.add(worker_user)
    db.commit()
    db.refresh(worker_user)
    
    print("✅ Worker 用户创建成功")
    print(f"   用户名: {worker_user.username}")
    print(f"   邮箱: {worker_user.email}")
    print(f"   密码: worker_password_change_me")
    
    # 生成 Token
    token = create_access_token(data={"sub": worker_user.username})
    print(f"\n🔑 Worker Token:")
    print(f"   {token}")
    print(f"\n请将此 Token 复制到 deployment/polling_worker_config_tencent.yaml 中的 api.worker_token")
    
    return token


def main():
    print("=" * 60)
    print("  创建 Worker 用户并生成 Token")
    print("=" * 60)
    print()
    
    # 创建数据库会话
    db = SessionLocal()
    
    try:
        token = create_worker_user(db)
        
        print("\n" + "=" * 60)
        print("  完成！")
        print("=" * 60)
        print()
        print("下一步:")
        print("1. 将上面的 Token 复制到配置文件:")
        print("   deployment/polling_worker_config_tencent.yaml")
        print()
        print("2. 重启 Worker:")
        print("   bash deployment/stop_polling_worker.sh")
        print("   bash deployment/start_polling_worker.sh tencent")
        print()
        
    except Exception as e:
        print(f"\n❌ 错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        db.close()


if __name__ == "__main__":
    main()

