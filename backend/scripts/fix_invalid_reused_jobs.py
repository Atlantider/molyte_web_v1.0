#!/usr/bin/env python3
"""
修复无效的复用 QC 任务

问题描述：
1. 很多复用任务（is_reused=True）没有对应的 qc_results 记录
2. 复用链条过长，形成 A->B->C->D... 的链条
3. 链条末端的任务可能没有实际结果

修复策略：
1. 对于能追溯到有效根任务的复用任务：复制结果
2. 对于无法追溯的复用任务：重置为 SUBMITTED 状态，等待重新计算
3. 更新 reused_from_job_id 直接指向根任务（扁平化链条）

使用方法：
    # 仅检测（不修改数据库）
    python fix_invalid_reused_jobs.py --dry-run
    
    # 执行修复
    python fix_invalid_reused_jobs.py --fix
"""
import sys
import os
import argparse
import logging

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.utils.qc_reuse import batch_fix_invalid_reused_jobs

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='修复无效的复用 QC 任务')
    parser.add_argument('--dry-run', action='store_true', 
                       help='仅检测，不修改数据库')
    parser.add_argument('--fix', action='store_true',
                       help='执行修复')
    args = parser.parse_args()
    
    if not args.dry_run and not args.fix:
        print("请指定 --dry-run 或 --fix 参数")
        print("  --dry-run: 仅检测，不修改数据库")
        print("  --fix: 执行修复")
        sys.exit(1)
    
    db = SessionLocal()
    try:
        dry_run = not args.fix
        
        if dry_run:
            logger.info("=== 检测模式（不修改数据库）===")
        else:
            logger.info("=== 修复模式 ===")
        
        stats = batch_fix_invalid_reused_jobs(db, dry_run=dry_run)
        
        print("\n" + "=" * 50)
        print("修复统计:")
        print("=" * 50)
        print(f"  总复用任务数: {stats['total_reused_jobs']}")
        print(f"  有效任务数: {stats['valid_jobs']}")
        print(f"  可通过复制结果修复: {stats['fixed_by_copy']}")
        print(f"  需要重置为待计算: {stats['reset_to_submitted']}")
        
        if stats['errors']:
            print(f"\n  错误数: {len(stats['errors'])}")
            for err in stats['errors'][:10]:
                print(f"    - Job {err['job_id']}: {err['error']}")
        
        if dry_run:
            print("\n[提示] 这是检测模式，数据库未被修改。")
            print("       使用 --fix 参数执行实际修复。")
        else:
            print("\n[完成] 数据库已修复。")
        
    finally:
        db.close()


if __name__ == "__main__":
    main()

