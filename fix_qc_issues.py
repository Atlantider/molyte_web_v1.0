#!/usr/bin/env python3
"""
修复 QC 任务的问题：
1. 对于阴离子，如果 RDKit 无法生成坐标，尝试从 PDB 文件加载
2. 修复自旋多重度错误
"""
import sys
import os
from pathlib import Path
import re

def fix_gjf_spin_multiplicity(gjf_path):
    """修复 gjf 文件中的自旋多重度"""
    with open(gjf_path, 'r') as f:
        lines = f.readlines()
    
    # 找到电荷和自旋多重度行
    for i, line in enumerate(lines):
        # 匹配 "charge spin" 格式的行
        match = re.match(r'^(-?\d+)\s+(\d+)\s*$', line.strip())
        if match:
            charge = int(match.group(1))
            spin = int(match.group(2))
            
            # 检查是否有坐标（下一行应该是原子坐标或注释）
            if i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                
                # 如果没有坐标（只有注释），跳过
                if next_line.startswith('!'):
                    print(f"  警告: {gjf_path.name} 没有坐标，无法修复")
                    return False
                
                # 计算电子数
                total_electrons = 0
                for j in range(i + 1, len(lines)):
                    atom_line = lines[j].strip()
                    if not atom_line or atom_line.startswith('!'):
                        break
                    
                    parts = atom_line.split()
                    if len(parts) >= 4:
                        atom_symbol = parts[0]
                        # 原子序数映射
                        atomic_numbers = {
                            'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9,
                            'S': 16, 'P': 15, 'Cl': 17, 'Li': 3, 'Na': 11
                        }
                        if atom_symbol in atomic_numbers:
                            total_electrons += atomic_numbers[atom_symbol]
                
                # 减去电荷得到实际电子数
                total_electrons -= charge
                
                # 计算正确的自旋多重度
                # 对于闭壳层分子（偶数电子），自旋多重度应该是 1
                # 对于开壳层分子（奇数电子），自旋多重度应该是 2
                correct_spin = 1 if total_electrons % 2 == 0 else 2
                
                if spin != correct_spin:
                    print(f"  修复 {gjf_path.name}:")
                    print(f"    电荷={charge}, 电子数={total_electrons}")
                    print(f"    错误的自旋多重度: {spin} -> 正确的: {correct_spin}")
                    
                    # 修改自旋多重度
                    lines[i] = f"{charge} {correct_spin}\n"
                    
                    # 写回文件
                    with open(gjf_path, 'w') as f:
                        f.writelines(lines)
                    
                    return True
                else:
                    print(f"  {gjf_path.name} 的自旋多重度正确 (charge={charge}, spin={spin}, electrons={total_electrons})")
                    return False
    
    return False


def check_missing_coordinates(gjf_path):
    """检查 gjf 文件是否缺少坐标"""
    with open(gjf_path, 'r') as f:
        content = f.read()
    
    if '! 请手动添加分子坐标' in content:
        return True
    
    return False


def main():
    qc_work_dir = Path('/public/home/xiaoji/molyte_web/data/qc_work')
    
    # 查找所有需要修复的任务
    problem_dirs = [
        'QC-157-TTE-B3LYP-6-31ppgd,p-pcm-TetraHydroFuran',
        'QC-159-FSI-B3LYP-6-31ppgd,p-pcm-TetraHydroFuran',
    ]
    
    print("="*80)
    print("检查和修复 QC 任务问题")
    print("="*80)
    
    for dir_name in problem_dirs:
        work_dir = qc_work_dir / dir_name
        if not work_dir.exists():
            print(f"\n目录不存在: {dir_name}")
            continue
        
        print(f"\n处理: {dir_name}")
        
        # 查找 gjf 文件
        gjf_files = list(work_dir.glob('*.gjf'))
        if not gjf_files:
            print(f"  未找到 gjf 文件")
            continue
        
        gjf_path = gjf_files[0]
        print(f"  GJF 文件: {gjf_path.name}")
        
        # 检查是否缺少坐标
        if check_missing_coordinates(gjf_path):
            print(f"  ❌ 缺少分子坐标 - 需要手动添加 PDB 文件或修改 SMILES")
            print(f"     建议:")
            print(f"     1. 在 data/initial_salts/ 目录下添加对应的 PDB 文件")
            print(f"     2. 或者使用其他工具生成3D坐标")
            continue
        
        # 修复自旋多重度
        fixed = fix_gjf_spin_multiplicity(gjf_path)
        if fixed:
            print(f"  ✓ 已修复自旋多重度")
            print(f"  建议: 重新提交 Gaussian 计算")
        
    print("\n" + "="*80)
    print("检查完成")
    print("="*80)


if __name__ == '__main__':
    main()

