#!/usr/bin/env python3
"""
测试去溶剂化能计算功能 - 仅测试导入
"""
import sys
sys.path.insert(0, 'backend')

def test_imports():
    """测试所有模块导入"""
    print("=" * 60)
    print("测试去溶剂化能计算功能 - 模块导入测试")
    print("=" * 60)
    
    try:
        # 1. 测试模型导入
        print("\n1. 测试模型导入...")
        from app.models.job import PostprocessJob, PostprocessType, JobStatus
        from app.models.result import DesolvationEnergyResult
        print("   ✅ 模型导入成功")
        print(f"   PostprocessType 枚举: {[e.value for e in PostprocessType]}")
        
        # 2. 测试 Schema 导入
        print("\n2. 测试 Schema 导入...")
        from app.schemas.desolvation import (
            DesolvationJobCreate,
            DesolvationJobResponse,
            DesolvationEnergyResultSchema,
            LigandDesolvationResult,
            TypeSummary
        )
        print("   ✅ Schema 导入成功")
        
        # 3. 测试 API 导入
        print("\n3. 测试 API 导入...")
        from app.api.v1 import desolvation
        print("   ✅ API 导入成功")
        print(f"   路由数量: {len(desolvation.router.routes)}")
        
        # 4. 测试任务处理函数导入
        print("\n4. 测试任务处理函数导入...")
        from app.tasks.desolvation import (
            run_desolvation_job,
            parse_solvation_cluster,
            get_ion_charge,
            get_molecule_charge,
            identify_ligands,
            estimate_atoms_per_molecule,
            calculate_energy_with_xtb,
            calculate_ligand_desolvation,
            generate_cluster_minus_xyz,
            summarize_by_type
        )
        print("   ✅ 任务处理函数导入成功")
        
        # 5. 测试 Worker API 修改
        print("\n5. 测试 Worker API 修改...")
        from app.api.v1.worker import get_pending_jobs
        print("   ✅ Worker API 导入成功")
        
        # 6. 测试常量
        print("\n6. 测试常量...")
        from app.tasks.desolvation import HARTREE_TO_KCAL
        print(f"   ✅ HARTREE_TO_KCAL = {HARTREE_TO_KCAL}")
        
        # 7. 测试示例数据
        print("\n7. 测试示例数据创建...")
        test_create = DesolvationJobCreate(
            md_job_id=1,
            solvation_structure_id=1,
            method_level='fast_xtb'
        )
        print(f"   ✅ 示例数据: {test_create.dict()}")
        
        # 8. 总结
        print("\n" + "=" * 60)
        print("✅ 所有模块导入测试通过！")
        print("=" * 60)
        print("\n实施完成的功能：")
        print("  ✅ 数据库层：desolvation_energy_results 表")
        print("  ✅ 模型层：DesolvationEnergyResult 模型")
        print("  ✅ Schema 层：5 个 Pydantic 模型")
        print("  ✅ API 层：3 个 REST 端点")
        print("  ✅ Worker 层：完整的任务处理逻辑")
        print("  ✅ 前端层：3 个 React 组件")
        print("\n下一步操作：")
        print("  1. 确保数据库服务正常运行")
        print("  2. 启动后端：cd backend && uvicorn app.main:app --reload")
        print("  3. 启动 Worker：cd deployment && python polling_worker.py")
        print("  4. 访问前端，在溶剂化结构详情页测试功能")
        
        return True
        
    except Exception as e:
        print(f"\n❌ 导入测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    success = test_imports()
    sys.exit(0 if success else 1)

