/**
 * Desolvation energy calculation types
 * 去溶剂化能计算相关类型定义
 */

export interface DesolvationJobCreate {
  md_job_id: number;
  solvation_structure_id: number;
  method_level: string;
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
  created_at: string;
  started_at?: string;
  finished_at?: string;
  elapsed_seconds?: number;
  error_message?: string;
  result?: DesolvationEnergyResult;
}

