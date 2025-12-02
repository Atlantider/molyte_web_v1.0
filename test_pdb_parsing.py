#!/usr/bin/env python3
"""测试 PDB 文件解析"""
from pathlib import Path

def parse_pdb_coordinates(pdb_path: Path):
    """解析 PDB 文件坐标"""
    coords = []
    
    # 尝试不同的编码
    encodings = ['utf-8', 'latin1', 'gbk', 'gb2312']
    
    for encoding in encodings:
        try:
            with open(pdb_path, 'r', encoding=encoding) as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        try:
                            # PDB 格式: ATOM/HETATM, atom_num, atom_name, ..., x, y, z
                            # 标准 PDB 格式的列位置
                            atom_name = line[12:16].strip()
                            
                            # 提取元素符号
                            # 优先使用 76-78 列的元素符号（如果存在）
                            element = ""
                            if len(line) > 76:
                                element = line[76:78].strip()
                                if element:
                                    # 去掉电荷符号（如 Mg2+）
                                    element = ''.join(c for c in element if c.isalpha())
                            
                            # 如果没有元素符号列，从原子名称提取
                            if not element or not element.strip():
                                element = ''.join(c for c in atom_name if c.isalpha())
                            
                            if not element:
                                element = atom_name[0] if atom_name else 'C'
                            
                            # 提取坐标（标准 PDB 格式）
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            
                            coords.append((element, x, y, z))
                        except (ValueError, IndexError) as e:
                            # 跳过无法解析的行
                            print(f"  跳过无法解析的行: {line.strip()}, 错误: {e}")
                            continue
            
            # 如果成功读取到坐标，返回结果
            if coords:
                print(f"  ✓ 使用 {encoding} 编码成功解析")
                return coords
                
        except UnicodeDecodeError:
            # 尝试下一个编码
            continue
        except Exception as e:
            print(f"  ✗ 使用 {encoding} 编码失败: {e}")
            continue
    
    # 所有编码都失败
    if not coords:
        print(f"  ✗ 无法解析（尝试了所有编码）")
    
    return coords if coords else None


# 测试文件
test_files = [
    "data/initial_salts/FSI.pdb",
    "data/initial_salts/Mg.pdb",
    "data/initial_salts/ClO4.pdb",
    "data/initial_salts/Li.pdb",
]

print("="*100)
print("测试 PDB 文件解析")
print("="*100)

for pdb_file in test_files:
    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        print(f"\n{pdb_file}: 文件不存在")
        continue
    
    print(f"\n{pdb_file}:")
    coords = parse_pdb_coordinates(pdb_path)
    
    if coords:
        print(f"  成功解析 {len(coords)} 个原子:")
        for element, x, y, z in coords:
            print(f"    {element:>2s}  {x:8.3f}  {y:8.3f}  {z:8.3f}")
    else:
        print(f"  解析失败")

print("\n" + "="*100)

