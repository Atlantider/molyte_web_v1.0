/**
 * Desolvation energy calculation types
 * 去溶剂化能计算相关类型定义
 */

/**
 * 去溶剂化模式
 * - stepwise: 逐级去溶剂（每次去掉一个配体）
 * - full: 全部去溶剂（直接计算中心离子单独能量）
 */
export type DesolvationMode = 'stepwise' | 'full';

export interface DesolvationJobCreate {
  md_job_id: number;
  solvation_structure_id: number;
  method_level: string;
  desolvation_mode?: DesolvationMode;  // 默认为 'stepwise'
}

export interface LigandDesolvationResult {
  ligand_id: number;
  ligand_type: string;
  ligand_label: string;
  e_ligand: number;
  e_cluster_minus: number;
  delta_e: number;
}

export interface TypeSummary {
  ligand_type: string;
  avg_delta_e: number;
  std_delta_e: number;
  count: number;
  min_delta_e: number;
  max_delta_e: number;
}

export interface DesolvationEnergyResult {
  id: number;
  postprocess_job_id: number;
  solvation_structure_id: number;
  method_level: string;
  basis_set?: string;
  functional?: string;
  e_cluster: number;
  per_ligand_results: LigandDesolvationResult[];
  per_type_summary: TypeSummary[];
  created_at: string;
}

export interface DesolvationJobResponse {
  job_id: number;
  status: string;
  method_level: string;
  desolvation_mode: DesolvationMode;
  created_at: string;
  started_at?: string;
  finished_at?: string;
  elapsed_seconds?: number;
  error_message?: string;
  result?: DesolvationEnergyResult;
}

