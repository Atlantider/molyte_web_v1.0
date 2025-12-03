#!/usr/bin/env python3
"""
测试密度提取修复

验证新的密度提取方法是否正确
"""

def test_parse_thermo_line():
    """测试密度、盒子尺寸、温度和能量提取"""

    # LAMMPS thermo 输出格式: Step CPU CPULeft Temp Density Lx Ly Lz TotEng KinEng PotEng ...
    # 列索引:                  0    1   2       3    4       5  6  7  8      9       10
    def parse_thermo_line(parts):
        density_val = None
        box_vals = []
        temp_val = None
        energy_vals = {}
        try:
            # 根据列位置提取密度（第5列，索引4）
            if len(parts) > 4:
                density_val = float(parts[4])

            # 根据列位置提取盒子尺寸（第6-8列，索引5-7）
            if len(parts) > 7:
                box_vals = [float(parts[5]), float(parts[6]), float(parts[7])]

            # 根据列位置提取温度（第4列，索引3）
            if len(parts) > 3:
                temp_val = float(parts[3])

            # 根据列位置提取能量值（第9-11列，索引8-10）
            if len(parts) > 10:
                energy_vals = {
                    'total_energy': float(parts[8]),
                    'kinetic_energy': float(parts[9]),
                    'potential_energy': float(parts[10])
                }
        except (ValueError, IndexError):
            # 如果按列位置提取失败，返回 None
            density_val = None
            box_vals = []
            temp_val = None
            energy_vals = {}

        return density_val, box_vals, temp_val, energy_vals
    
    # 测试数据（来自实际的 LAMMPS log 文件）
    test_cases = [
        # 初始状态：密度 0.50212718，盒子 40×40×40，温度 0，能量 19131.958
        {
            'line': '0 0 0 0 0.50212718 40 40 40 19131.958 0 19131.958 14198.18 7122.4663 -7119.4818 3144.3554 1118.0087 668.42896 0',
            'expected_density': 0.50212718,
            'expected_box': [40.0, 40.0, 40.0],
            'expected_temp': 0.0,
            'expected_energy': {'total_energy': 19131.958, 'kinetic_energy': 0.0, 'potential_energy': 19131.958},
            'description': '初始状态'
        },
        # 中间状态：密度 0.53554349，盒子 39.15011×39.15011×39.15011，温度 24.083998
        {
            'line': '500 0.77026475 7701.9396 24.083998 0.53554349 39.15011 39.15011 39.15011 -5932.0803 198.71423 -6130.7945 -809.81366 2391.6909 -8831.1829 14.253354 520.79861 583.45927 0',
            'expected_density': 0.53554349,
            'expected_box': [39.15011, 39.15011, 39.15011],
            'expected_temp': 24.083998,
            'expected_energy': {'total_energy': -5932.0803, 'kinetic_energy': 198.71423, 'potential_energy': -6130.7945},
            'description': '中间状态'
        },
        # 最终状态：密度 1.665438，盒子 32.505703×32.505703×32.505703，温度 298.45175
        {
            'line': '10000000 14698.757 0 298.45175 1.665438 32.505703 32.505703 32.505703 -1981.8958 2219.6218 -4201.5176 -1070.1361 8232.0376 -12318.317 750.85511 1229.3299 -1025.2876 0',
            'expected_density': 1.665438,
            'expected_box': [32.505703, 32.505703, 32.505703],
            'expected_temp': 298.45175,
            'expected_energy': {'total_energy': -1981.8958, 'kinetic_energy': 2219.6218, 'potential_energy': -4201.5176},
            'description': '最终状态'
        },
    ]
    
    print("=" * 80)
    print("测试密度、盒子尺寸、温度和能量提取")
    print("=" * 80)

    all_passed = True
    for i, test_case in enumerate(test_cases, 1):
        parts = test_case['line'].split()
        density, box, temp, energy = parse_thermo_line(parts)

        print(f"\n测试 {i}: {test_case['description']}")
        print(f"  输入行: {test_case['line'][:80]}...")

        # 检查密度
        density_ok = abs(density - test_case['expected_density']) < 1e-6
        print(f"  密度: {density} (期望: {test_case['expected_density']}) {'✓' if density_ok else '✗'}")

        # 检查盒子尺寸
        box_ok = all(abs(box[i] - test_case['expected_box'][i]) < 1e-5 for i in range(3))
        print(f"  盒子: {box} (期望: {test_case['expected_box']}) {'✓' if box_ok else '✗'}")

        # 检查温度
        temp_ok = abs(temp - test_case['expected_temp']) < 1e-5
        print(f"  温度: {temp} K (期望: {test_case['expected_temp']} K) {'✓' if temp_ok else '✗'}")

        # 检查能量值
        energy_ok = all(
            abs(energy.get(key, 0) - test_case['expected_energy'][key]) < 1e-3
            for key in test_case['expected_energy']
        )
        print(f"  能量:")
        print(f"    总能量: {energy.get('total_energy')} (期望: {test_case['expected_energy']['total_energy']}) {'✓' if abs(energy.get('total_energy', 0) - test_case['expected_energy']['total_energy']) < 1e-3 else '✗'}")
        print(f"    动能: {energy.get('kinetic_energy')} (期望: {test_case['expected_energy']['kinetic_energy']}) {'✓' if abs(energy.get('kinetic_energy', 0) - test_case['expected_energy']['kinetic_energy']) < 1e-3 else '✗'}")
        print(f"    势能: {energy.get('potential_energy')} (期望: {test_case['expected_energy']['potential_energy']}) {'✓' if abs(energy.get('potential_energy', 0) - test_case['expected_energy']['potential_energy']) < 1e-3 else '✗'}")

        if not (density_ok and box_ok and temp_ok and energy_ok):
            all_passed = False

    print("\n" + "=" * 80)
    if all_passed:
        print("✓ 所有测试通过！")
    else:
        print("✗ 有测试失败")
    print("=" * 80)

    return all_passed


def test_concentration_calculation():
    """测试浓度计算"""
    print("\n" + "=" * 70)
    print("测试浓度计算")
    print("=" * 70)
    
    AVOGADRO = 6.022e23
    
    # 测试场景：初始密度 1.0434 g/cm³，最终密度 0.8155 g/cm³
    # 初始盒子：55.88 × 55.88 × 55.88 Å
    # 最终盒子：22.68 × 59.71 × 59.71 Å
    
    # 初始浓度计算
    init_box_volume = 55.88 * 55.88 * 55.88  # Å³
    init_volume_L = init_box_volume * 1e-27  # L
    cation_count = 10  # 假设有 10 个阳离子
    init_conc = (cation_count / AVOGADRO) / init_volume_L
    
    # 最终浓度计算
    final_box_volume = 22.68 * 59.71 * 59.71  # Å³
    final_volume_L = final_box_volume * 1e-27  # L
    final_conc = (cation_count / AVOGADRO) / final_volume_L
    
    print(f"\n初始状态:")
    print(f"  盒子体积: {init_box_volume:.2f} Å³ = {init_volume_L:.2e} L")
    print(f"  阳离子数: {cation_count}")
    print(f"  初始浓度: {init_conc:.4f} mol/L")
    
    print(f"\n最终状态:")
    print(f"  盒子体积: {final_box_volume:.2f} Å³ = {final_volume_L:.2e} L")
    print(f"  阳离子数: {cation_count}")
    print(f"  最终浓度: {final_conc:.4f} mol/L")
    
    print(f"\n分析:")
    print(f"  体积变化: {init_box_volume:.2f} → {final_box_volume:.2f} Å³")
    print(f"  体积比: {final_box_volume / init_box_volume:.4f}")
    print(f"  浓度比: {final_conc / init_conc:.4f}")
    
    if final_conc > init_conc:
        print(f"  ✓ 体积减小，浓度增大（符合物理规律）")
    else:
        print(f"  ✗ 体积减小，但浓度没有增大（不符合物理规律）")
    
    print("=" * 70)


if __name__ == "__main__":
    test_parse_thermo_line()
    test_concentration_calculation()

