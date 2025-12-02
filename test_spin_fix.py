#!/usr/bin/env python3
"""测试自旋多重度修复"""
import sys
sys.path.insert(0, '/public/home/xiaoji/molyte_web/backend')

from rdkit import Chem

def calculate_spin_OLD(smiles: str, charge: int) -> int:
    """旧版本（有BUG）"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 1
    
    # BUG: 没有添加氢原子
    total_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms())
    total_electrons -= charge
    
    if total_electrons % 2 == 0:
        return 1
    else:
        return 2


def calculate_spin_NEW(smiles: str, charge: int) -> int:
    """新版本（已修复）"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 1
    
    # 添加氢原子以获得完整的电子数
    mol_with_h = Chem.AddHs(mol)
    
    # 计算总电子数
    total_electrons = sum(atom.GetAtomicNum() for atom in mol_with_h.GetAtoms())
    total_electrons -= charge
    
    # 检查是否有显式自由基
    num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    
    if num_radical_electrons > 0:
        spin_multiplicity = num_radical_electrons + 1
    else:
        if total_electrons % 2 == 0:
            spin_multiplicity = 1
        else:
            spin_multiplicity = 2
    
    return spin_multiplicity, total_electrons


# 测试用例
test_cases = [
    ("TTE", "FC(F)(F)COC(C(F)(F)F)C(F)(F)F", 0, 1),  # 应该是 1
    ("FSI-", "O=S(=O)(NS(=O)(F)=O)F", -1, 1),  # 应该是 1
    ("EC", "C1COC(=O)O1", 0, 1),  # 应该是 1
    ("EMC", "CCOC(=O)OC", 0, 1),  # 应该是 1
    ("Li+", "[Li+]", 1, 1),  # 应该是 1
    ("PF6-", "F[P-](F)(F)(F)(F)F", -1, 1),  # 应该是 1
]

print("="*100)
print("测试自旋多重度计算修复")
print("="*100)

errors_found = []

for name, smiles, charge, expected_spin in test_cases:
    print(f"\n{name}: {smiles} (charge={charge})")
    
    spin_old = calculate_spin_OLD(smiles, charge)
    spin_new, electrons = calculate_spin_NEW(smiles, charge)
    
    print(f"  旧版本: spin={spin_old}")
    print(f"  新版本: spin={spin_new}, electrons={electrons}")
    print(f"  预期值: spin={expected_spin}")
    
    if spin_old != expected_spin:
        print(f"  ❌ 旧版本错误！")
        errors_found.append((name, "旧版本", spin_old, expected_spin))
    else:
        print(f"  ✓ 旧版本正确")
        
    if spin_new != expected_spin:
        print(f"  ❌ 新版本错误！")
        errors_found.append((name, "新版本", spin_new, expected_spin))
    else:
        print(f"  ✓ 新版本正确")

print("\n" + "="*100)
print("总结")
print("="*100)

if errors_found:
    print(f"\n发现 {len(errors_found)} 个错误:")
    for name, version, got, expected in errors_found:
        print(f"  - {name} ({version}): 得到 {got}, 预期 {expected}")
else:
    print("\n✓ 所有测试通过！")

print("\n修复说明:")
print("  1. 旧版本没有添加氢原子，导致电子数计算不准确")
print("  2. 新版本使用 Chem.AddHs() 添加氢原子后再计算电子数")
print("  3. 这样可以正确判断分子的自旋多重度")

