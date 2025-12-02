#!/usr/bin/env python3
"""为 FSI 重新生成 gjf 文件（使用 PDB 坐标）"""
from pathlib import Path

def parse_pdb_coordinates(pdb_path):
    """解析 PDB 文件坐标"""
    coords = []
    try:
        with open(pdb_path, 'r', encoding='latin1') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # PDB 格式: ATOM/HETATM, 序号, 原子名, 残基名, 链, 残基号, x, y, z
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    # 提取元素符号（通常是原子名的第一个字符）
                    element = atom_name[0]
                    if len(atom_name) > 1 and atom_name[1].isalpha() and atom_name[1].islower():
                        element = atom_name[:2]
                    
                    coords.append((element, x, y, z))
        
        return coords
    except Exception as e:
        print(f"解析 PDB 文件失败: {e}")
        return None


def generate_gjf_with_pdb_coords():
    """为 FSI 生成包含 PDB 坐标的 gjf 文件"""
    # 路径
    pdb_path = Path("/public/home/xiaoji/molyte_web/data/initial_salts/FSI.pdb")
    work_dir = Path("/public/home/xiaoji/molyte_web/data/qc_work/QC-159-FSI-B3LYP-6-31ppgd,p-pcm-TetraHydroFuran")
    gjf_path = work_dir / "FSI-B3LYP-6-31ppgd_p-pcm-TetraHydroFuran.gjf"
    
    # 解析 PDB 坐标
    print(f"从 PDB 文件读取坐标: {pdb_path}")
    coords = parse_pdb_coordinates(pdb_path)
    
    if not coords:
        print("❌ 无法读取 PDB 坐标")
        return False
    
    print(f"✓ 读取到 {len(coords)} 个原子")
    
    # 生成 gjf 内容
    gjf_content = """%nprocshared=16
%mem=8GB
%chk=FSI-B3LYP-6-31ppgd_p-pcm-TetraHydroFuran.chk
# opt freq B3LYP/6-31++g(d,p) em=gd3bj scrf=(pcm,solvent=TetraHydroFuran)

FSI-B3LYP-6-31ppgd,p-pcm-TetraHydroFuran

-1 1
"""
    
    # 添加坐标
    for atom, x, y, z in coords:
        gjf_content += f" {atom:<2}  {x:>12.8f}  {y:>12.8f}  {z:>12.8f}\n"
    
    gjf_content += "\n\n"
    
    # 备份原文件
    if gjf_path.exists():
        backup_path = gjf_path.with_suffix('.gjf.bak')
        gjf_path.rename(backup_path)
        print(f"✓ 备份原文件到: {backup_path}")
    
    # 写入新文件
    with open(gjf_path, 'w') as f:
        f.write(gjf_content)
    
    print(f"✓ 生成新的 gjf 文件: {gjf_path}")
    print(f"\n预览前10行:")
    print("="*60)
    for i, line in enumerate(gjf_content.split('\n')[:15], 1):
        print(f"{i:3d}: {line}")
    print("="*60)
    
    return True


if __name__ == '__main__':
    print("="*80)
    print("为 FSI 重新生成 gjf 文件")
    print("="*80)
    
    if generate_gjf_with_pdb_coords():
        print("\n✓ 成功！现在可以重新提交 Gaussian 计算")
    else:
        print("\n❌ 失败")

