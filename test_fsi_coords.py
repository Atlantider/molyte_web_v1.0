#!/usr/bin/env python3
"""测试 FSI 坐标生成"""
from pathlib import Path

# 模拟 worker 的查找逻辑
initial_salts_path = Path("/public/home/xiaoji/molyte_web/data/initial_salts")

# FSI 的可能名称
molecule_names = [
    "FSI-",
    "FSI",
    "FSI-B3LYP-6-31ppgd,p-pcm-TetraHydroFuran",
]

for molecule_name in molecule_names:
    print(f"\n测试分子名称: {molecule_name}")
    
    # 清理名称
    clean_name = molecule_name.replace("+", "").replace("-", "").strip()
    print(f"  清理后的名称: {clean_name}")
    
    # 可能的路径
    possible_paths = [
        initial_salts_path / f"{clean_name}.pdb",
        initial_salts_path / f"{molecule_name}.pdb",
    ]
    
    for pdb_path in possible_paths:
        print(f"  检查: {pdb_path}")
        if pdb_path.exists():
            print(f"    ✓ 找到!")
        else:
            print(f"    ✗ 不存在")

# 列出实际存在的 PDB 文件
print(f"\n实际存在的 PDB 文件:")
for pdb_file in initial_salts_path.glob("*.pdb"):
    print(f"  {pdb_file.name}")

