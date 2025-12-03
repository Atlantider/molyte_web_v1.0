#!/usr/bin/env python3
"""
测试数据库枚举类型修复
"""
import sys
sys.path.insert(0, 'backend')

from app.database import SessionLocal
from app.models.job import PostprocessJob, PostprocessType, JobStatus
from sqlalchemy import text

def test_enum_values():
    """测试枚举值是否正确"""
    db = SessionLocal()
    try:
        # 1. 测试从数据库读取枚举值
        print("1. 测试数据库枚举值...")
        result = db.execute(text("SELECT enum_range(NULL::postprocesstype)"))
        enum_values = result.fetchone()[0]
        print(f"   数据库枚举值: {enum_values}")
        
        if 'DESOLVATION_ENERGY' in enum_values:
            print("   ✅ DESOLVATION_ENERGY 存在于数据库枚举中")
        else:
            print("   ❌ DESOLVATION_ENERGY 不存在于数据库枚举中")
            return False
        
        # 2. 测试 Python 枚举
        print("\n2. 测试 Python 枚举...")
        python_enum_values = [e.value for e in PostprocessType]
        print(f"   Python 枚举值: {python_enum_values}")
        
        if 'DESOLVATION_ENERGY' in python_enum_values:
            print("   ✅ DESOLVATION_ENERGY 存在于 Python 枚举中")
        else:
            print("   ❌ DESOLVATION_ENERGY 不存在于 Python 枚举中")
            return False
        
        # 3. 测试创建 PostprocessJob（不提交）
        print("\n3. 测试创建 PostprocessJob...")
        test_job = PostprocessJob(
            md_job_id=1,
            job_type=PostprocessType.DESOLVATION_ENERGY,
            status=JobStatus.CREATED,
            config={"test": True}
        )
        db.add(test_job)
        db.flush()  # 不提交，只验证
        print(f"   ✅ 成功创建测试任务 (ID: {test_job.id})")
        
        # 4. 回滚测试数据
        db.rollback()
        print("   ✅ 测试数据已回滚")
        
        print("\n" + "="*50)
        print("✅ 所有测试通过！枚举类型修复成功！")
        print("="*50)
        return True
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        db.close()

if __name__ == "__main__":
    success = test_enum_values()
    sys.exit(0 if success else 1)

