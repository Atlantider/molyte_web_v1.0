#!/usr/bin/env python3
"""测试自旋多重度计算逻辑"""
from rdkit import Chem

def get_charge_and_spin_from_smiles_OLD(smiles: str):
    """旧版本：从SMILES获取电荷和自旋多重度（有BUG）"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0, 1

        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        unpaired_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        spin_multiplicity = unpaired_electrons + 1

        return total_charge, spin_multiplicity
    except Exception as e:
        print(f"Failed to parse SMILES {smiles}: {e}")
        return 0, 1


def get_charge_and_spin_from_smiles_NEW(smiles: str):
    """新版本：正确计算自旋多重度"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0, 1

        # 添加氢原子以获得完整的电子数
        mol_with_h = Chem.AddHs(mol)
        
        # 计算总电荷
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        
        # 计算总电子数
        total_electrons = sum(atom.GetAtomicNum() for atom in mol_with_h.GetAtoms())
        total_electrons -= total_charge  # 减去电荷
        
        # 检查是否有显式自由基
        unpaired_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        
        if unpaired_electrons > 0:
            # 有显式自由基
            spin_multiplicity = unpaired_electrons + 1
        else:
            # 根据电子数判断
            # 偶数电子 -> 闭壳层 -> 自旋多重度 = 1
            # 奇数电子 -> 开壳层 -> 自旋多重度 = 2
            if total_electrons % 2 == 0:
                spin_multiplicity = 1
            else:
                spin_multiplicity = 2

        return total_charge, spin_multiplicity, total_electrons, unpaired_electrons
    except Exception as e:
        print(f"Failed to parse SMILES {smiles}: {e}")
        return 0, 1, 0, 0


# 测试用例
test_cases = [
    ("TTE", "FC(F)(F)COC(C(F)(F)F)C(F)(F)F", 0),  # 中性分子
    ("FSI-", "O=S(=O)(NS(=O)(F)=O)F", -1),  # 阴离子
    ("Li+", "[Li+]", 1),  # 阳离子
    ("EC", "C1COC(=O)O1", 0),  # 碳酸乙烯酯
    ("EMC", "CCOC(=O)OC", 0),  # 碳酸甲乙酯
    ("PF6-", "F[P-](F)(F)(F)(F)F", -1),  # 六氟磷酸根
]

print("="*100)
print("测试自旋多重度计算")
print("="*100)

for name, smiles, expected_charge in test_cases:
    print(f"\n{name}: {smiles}")
    print(f"  预期电荷: {expected_charge}")
    
    # 旧版本
    charge_old, spin_old = get_charge_and_spin_from_smiles_OLD(smiles)
    print(f"  旧版本: charge={charge_old}, spin={spin_old}")
    
    # 新版本
    charge_new, spin_new, electrons, radicals = get_charge_and_spin_from_smiles_NEW(smiles)
    print(f"  新版本: charge={charge_new}, spin={spin_new}, electrons={electrons}, radicals={radicals}")
    
    # 验证
    if charge_new != expected_charge:
        print(f"  ⚠️  电荷不匹配！预期 {expected_charge}，得到 {charge_new}")
    
    if electrons % 2 == 0 and spin_new != 1:
        print(f"  ❌ 错误！偶数电子 ({electrons}) 应该是自旋多重度 1，但得到 {spin_new}")
    elif electrons % 2 == 1 and spin_new != 2 and radicals == 0:
        print(f"  ❌ 错误！奇数电子 ({electrons}) 应该是自旋多重度 2，但得到 {spin_new}")
    else:
        print(f"  ✓ 正确")

print("\n" + "="*100)
print("结论:")
print("="*100)
print("旧版本的问题:")
print("  1. 没有添加氢原子，导致电子数计算不准确")
print("  2. 只依赖 GetNumRadicalElectrons()，对于大多数分子返回 0")
print("  3. 导致 spin_multiplicity = 0 + 1 = 1，即使对于奇数电子的分子")
print("\n新版本的改进:")
print("  1. 添加氢原子后计算总电子数")
print("  2. 根据电子数奇偶性判断自旋多重度")
print("  3. 偶数电子 -> 自旋多重度 1 (闭壳层)")
print("  4. 奇数电子 -> 自旋多重度 2 (开壳层)")

