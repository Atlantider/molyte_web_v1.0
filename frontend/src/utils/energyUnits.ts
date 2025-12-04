/**
 * 能量单位转换工具
 * Energy unit conversion utilities
 */

export type EnergyUnit = 'kcal/mol' | 'kJ/mol' | 'eV' | 'Hartree';

// 转换因子 (从 kcal/mol 到其他单位)
const CONVERSION_FACTORS: Record<EnergyUnit, number> = {
  'kcal/mol': 1,
  'kJ/mol': 4.184,        // 1 kcal = 4.184 kJ
  'eV': 0.0433641,        // 1 kcal/mol = 0.0433641 eV
  'Hartree': 0.00159362,  // 1 kcal/mol = 0.00159362 Hartree
};

// 单位显示名称
export const UNIT_LABELS: Record<EnergyUnit, string> = {
  'kcal/mol': 'kcal/mol',
  'kJ/mol': 'kJ/mol',
  'eV': 'eV',
  'Hartree': 'A.U.',
};

// 单位精度
export const UNIT_PRECISION: Record<EnergyUnit, number> = {
  'kcal/mol': 2,
  'kJ/mol': 2,
  'eV': 4,
  'Hartree': 6,
};

/**
 * 将能量值从 kcal/mol 转换到目标单位
 */
export function convertEnergy(value: number, toUnit: EnergyUnit): number {
  return value * CONVERSION_FACTORS[toUnit];
}

/**
 * 格式化能量值显示
 */
export function formatEnergy(value: number, unit: EnergyUnit): string {
  const converted = convertEnergy(value, unit);
  const precision = UNIT_PRECISION[unit];
  return converted.toFixed(precision);
}

/**
 * 获取单位选项列表
 */
export function getUnitOptions(): { value: EnergyUnit; label: string }[] {
  return [
    { value: 'kcal/mol', label: 'kcal/mol' },
    { value: 'kJ/mol', label: 'kJ/mol' },
    { value: 'eV', label: 'eV' },
    { value: 'Hartree', label: 'Hartree (A.U.)' },
  ];
}

