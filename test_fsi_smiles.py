#!/usr/bin/env python3
"""测试 FSI SMILES"""
from rdkit import Chem

# 不同的 FSI SMILES 表示
fsi_smiles_variants = [
    ("无电荷", "O=S(=O)(NS(=O)(F)=O)F"),
    ("带电荷1", "O=S(=O)([N-]S(=O)(F)=O)F"),
    ("带电荷2", "[O-]S(=O)(NS(=O)(F)=O)F"),
    ("标准", "FS(=O)(=O)[N-]S(F)(=O)=O"),
]

print("="*100)
print("测试 FSI 的不同 SMILES 表示")
print("="*100)

for name, smiles in fsi_smiles_variants:
    print(f"\n{name}: {smiles}")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("  ❌ 无法解析")
        continue
    
    # 计算形式电荷
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # 添加氢原子
    mol_with_h = Chem.AddHs(mol)
    
    # 计算总电子数
    total_electrons = sum(atom.GetAtomicNum() for atom in mol_with_h.GetAtoms())
    
    # 考虑电荷
    electrons_with_charge = total_electrons - formal_charge
    
    print(f"  形式电荷: {formal_charge}")
    print(f"  总电子数: {total_electrons}")
    print(f"  考虑电荷后: {electrons_with_charge}")
    print(f"  奇偶性: {'偶数' if electrons_with_charge % 2 == 0 else '奇数'}")
    print(f"  自旋多重度: {1 if electrons_with_charge % 2 == 0 else 2}")
    
    # 打印原子信息
    print(f"  原子: ", end="")
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        charge = atom.GetFormalCharge()
        if charge != 0:
            print(f"{symbol}({charge:+d}) ", end="")
        else:
            print(f"{symbol} ", end="")
    print()

print("\n" + "="*100)
print("结论:")
print("="*100)
print("FSI- 是一个阴离子，应该有 -1 的电荷")
print("如果 SMILES 中没有包含电荷信息，需要在计算时手动指定 charge=-1")
print("正确的做法是:")
print("  1. 使用带电荷的 SMILES: FS(=O)(=O)[N-]S(F)(=O)=O")
print("  2. 或者在计算时指定 charge=-1")

