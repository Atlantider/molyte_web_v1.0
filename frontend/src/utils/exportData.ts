/**
 * 数据导出工具
 * Data export utilities for CSV and Excel
 */
import type { DesolvationEnergyResult, LigandDesolvationResult, TypeSummary } from '../types/desolvation';
import type { EnergyUnit } from './energyUnits';
import { convertEnergy, UNIT_LABELS } from './energyUnits';

/**
 * 将去溶剂化能结果导出为 CSV
 */
export function exportDesolvationToCSV(
  result: DesolvationEnergyResult,
  unit: EnergyUnit,
  filename: string = 'desolvation_energy'
): void {
  const lines: string[] = [];
  
  // 头部信息
  lines.push('# Desolvation Energy Results');
  lines.push(`# Method Level: ${result.method_level}`);
  lines.push(`# Functional: ${result.functional || 'N/A'}`);
  lines.push(`# Basis Set: ${result.basis_set || 'N/A'}`);
  lines.push(`# E_cluster: ${result.e_cluster} A.U.`);
  lines.push(`# Energy Unit: ${UNIT_LABELS[unit]}`);
  lines.push('');
  
  // Per-ligand results
  lines.push('## Per Ligand Results');
  lines.push(`Ligand,Type,Delta_E (${UNIT_LABELS[unit]}),E_ligand (A.U.),E_cluster_minus (A.U.)`);
  
  for (const lig of result.per_ligand_results) {
    const deltaE = convertEnergy(lig.delta_e, unit);
    lines.push(`${lig.ligand_label},${lig.ligand_type},${deltaE.toFixed(4)},${lig.e_ligand.toFixed(6)},${lig.e_cluster_minus.toFixed(6)}`);
  }
  
  lines.push('');
  
  // Per-type summary
  lines.push('## Per Type Summary');
  lines.push(`Type,Count,Avg_Delta_E (${UNIT_LABELS[unit]}),Std_Delta_E,Min_Delta_E,Max_Delta_E`);
  
  for (const type of result.per_type_summary) {
    const avg = convertEnergy(type.avg_delta_e, unit);
    const std = convertEnergy(type.std_delta_e, unit);
    const min = convertEnergy(type.min_delta_e, unit);
    const max = convertEnergy(type.max_delta_e, unit);
    lines.push(`${type.ligand_type},${type.count},${avg.toFixed(4)},${std.toFixed(4)},${min.toFixed(4)},${max.toFixed(4)}`);
  }
  
  // 下载文件
  downloadFile(lines.join('\n'), `${filename}.csv`, 'text/csv;charset=utf-8');
}

/**
 * 批量导出多个结果到 CSV
 */
export function exportBatchDesolvationToCSV(
  results: Array<{ compositionKey: string; result: DesolvationEnergyResult }>,
  unit: EnergyUnit,
  filename: string = 'desolvation_batch'
): void {
  const lines: string[] = [];
  
  lines.push(`# Batch Desolvation Energy Results`);
  lines.push(`# Energy Unit: ${UNIT_LABELS[unit]}`);
  lines.push('');
  
  // 汇总表头
  lines.push('Composition,Ligand,Type,Delta_E,E_cluster,E_ligand,E_cluster_minus');
  
  for (const { compositionKey, result } of results) {
    for (const lig of result.per_ligand_results) {
      const deltaE = convertEnergy(lig.delta_e, unit);
      lines.push(`${compositionKey},${lig.ligand_label},${lig.ligand_type},${deltaE.toFixed(4)},${result.e_cluster.toFixed(6)},${lig.e_ligand.toFixed(6)},${lig.e_cluster_minus.toFixed(6)}`);
    }
  }
  
  downloadFile(lines.join('\n'), `${filename}.csv`, 'text/csv;charset=utf-8');
}

/**
 * 导出类型汇总到 CSV
 */
export function exportTypeSummaryToCSV(
  summaries: Array<{ compositionKey: string; summary: TypeSummary[] }>,
  unit: EnergyUnit,
  filename: string = 'desolvation_type_summary'
): void {
  const lines: string[] = [];
  
  lines.push(`Composition,Ligand_Type,Count,Avg_Delta_E (${UNIT_LABELS[unit]}),Std_Delta_E,Min,Max`);
  
  for (const { compositionKey, summary } of summaries) {
    for (const type of summary) {
      const avg = convertEnergy(type.avg_delta_e, unit);
      const std = convertEnergy(type.std_delta_e, unit);
      const min = convertEnergy(type.min_delta_e, unit);
      const max = convertEnergy(type.max_delta_e, unit);
      lines.push(`${compositionKey},${type.ligand_type},${type.count},${avg.toFixed(4)},${std.toFixed(4)},${min.toFixed(4)},${max.toFixed(4)}`);
    }
  }
  
  downloadFile(lines.join('\n'), `${filename}.csv`, 'text/csv;charset=utf-8');
}

/**
 * 下载文件
 */
function downloadFile(content: string, filename: string, mimeType: string): void {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

