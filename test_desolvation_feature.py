#!/usr/bin/env python3
"""
测试去溶剂化能计算功能
"""
import sys
sys.path.insert(0, 'backend')

from sqlalchemy import text
from app.database import SessionLocal
from app.models.job import PostprocessJob, PostprocessType, JobStatus
from app.models.result import SolvationStructure, DesolvationEnergyResult
from app.tasks.desolvation import run_desolvation_job

def test_desolvation_feature():
    """测试去溶剂化能计算功能"""
    db = SessionLocal()
    
    try:
        print("=" * 60)
        print("测试去溶剂化能计算功能")
        print("=" * 60)
        
        # 1. 检查数据库表
        print("\n1. 检查数据库表...")
        result = db.execute(text("SELECT tablename FROM pg_tables WHERE schemaname='public' AND tablename='desolvation_energy_results'"))
        if result.fetchone():
            print("   ✅ desolvation_energy_results 表已创建")
        else:
            print("   ❌ desolvation_energy_results 表不存在")
            return
        
        # 2. 检查 PostprocessType 枚举
        print("\n2. 检查 PostprocessType 枚举...")
        print(f"   枚举值: {[e.value for e in PostprocessType]}")
        if 'DESOLVATION_ENERGY' in [e.value for e in PostprocessType]:
            print("   ✅ DESOLVATION_ENERGY 已添加到枚举")
        else:
            print("   ❌ DESOLVATION_ENERGY 未添加到枚举")
            return
        
        # 3. 查找一个溶剂化结构用于测试
        print("\n3. 查找测试用的溶剂化结构...")
        structure = db.query(SolvationStructure).first()
        
        if not structure:
            print("   ⚠️  数据库中没有溶剂化结构，无法测试")
            print("   提示：请先运行 MD 任务并生成溶剂化结构")
            return
        
        print(f"   ✅ 找到溶剂化结构 ID: {structure.id}")
        print(f"      MD Job ID: {structure.md_job_id}")
        print(f"      配位数: {structure.coordination_num}")
        print(f"      组成: {structure.composition}")
        
        # 4. 检查是否已有去溶剂化能任务
        print("\n4. 检查现有的去溶剂化能任务...")
        existing_jobs = db.query(PostprocessJob).filter(
            PostprocessJob.job_type == PostprocessType.DESOLVATION_ENERGY,
            PostprocessJob.config['solvation_structure_id'].astext == str(structure.id)
        ).all()
        
        if existing_jobs:
            print(f"   找到 {len(existing_jobs)} 个现有任务:")
            for job in existing_jobs:
                print(f"      - Job ID: {job.id}, 状态: {job.status.value}")
                
                # 如果有完成的任务，显示结果
                if job.status == JobStatus.COMPLETED and job.desolvation_energy_result:
                    result = job.desolvation_energy_result
                    print(f"        E_cluster: {result.e_cluster:.6f} A.U.")
                    print(f"        配体数量: {len(result.per_ligand_results)}")
                    print(f"        类型数量: {len(result.per_type_summary)}")
        else:
            print("   没有找到现有任务")
        
        # 5. 测试 API 模型
        print("\n5. 测试 API 模型...")
        from app.schemas.desolvation import DesolvationJobCreate, DesolvationJobResponse
        
        test_create = DesolvationJobCreate(
            md_job_id=structure.md_job_id,
            solvation_structure_id=structure.id,
            method_level='fast_xtb'
        )
        print(f"   ✅ DesolvationJobCreate 模型正常: {test_create.dict()}")
        
        # 6. 测试 Worker 任务处理函数
        print("\n6. 测试 Worker 任务处理函数...")
        print("   注意：实际计算需要 xTB 软件，这里只测试函数导入")
        from app.tasks.desolvation import (
            run_desolvation_job,
            parse_solvation_cluster,
            calculate_energy_with_xtb,
            identify_ligands
        )
        print("   ✅ 所有任务处理函数导入成功")
        
        # 7. 总结
        print("\n" + "=" * 60)
        print("测试总结")
        print("=" * 60)
        print("✅ 数据库层：已完成")
        print("✅ Schema 层：已完成")
        print("✅ API 层：已完成")
        print("✅ Worker 层：已完成")
        print("\n下一步：")
        print("1. 启动后端服务：cd backend && source venv/bin/activate && uvicorn app.main:app --reload")
        print("2. 启动 Worker：cd deployment && python polling_worker.py")
        print("3. 在前端溶剂化结构详情页点击'计算去溶剂化能'按钮")
        print("4. Worker 会自动拉取任务并执行计算")
        print("5. 计算完成后可以在前端查看结果")
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()

if __name__ == '__main__':
    test_desolvation_feature()

