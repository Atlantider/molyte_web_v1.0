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

/**
 * 溶剂模型类型
 * - gas: 气相（无溶剂效应）
 * - pcm: PCM 极化连续介质模型
 * - smd: SMD 溶剂密度模型
 * - custom: 自定义溶剂参数
 */
export type SolventModel = 'gas' | 'pcm' | 'smd' | 'custom';

/**
 * 溶剂配置
 */
export interface SolventConfig {
  model: SolventModel;
  solvent_name?: string;  // 预定义溶剂名称（如 water, acetonitrile）
  // 自定义溶剂参数（用于 custom 模型）
  eps?: number;           // 介电常数 ε
  eps_inf?: number;       // 光学介电常数 n²
  hbond_acidity?: number; // Abraham氢键酸度 α
  hbond_basicity?: number; // Abraham氢键碱度 β
  surface_tension?: number; // 表面张力 γ (cal/mol·Å²)
  carbon_aromaticity?: number; // 芳香碳原子比例 φ
  halogenicity?: number;  // 卤素原子比例 ψ
}

export interface DesolvationJobCreate {
  md_job_id: number;
  solvation_structure_id: number;
  method_level: string;
  desolvation_mode?: DesolvationMode;  // 默认为 'stepwise'
  solvent_config?: SolventConfig;      // 溶剂配置（可选，默认气相）
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

export interface QCProgress {
  total: number;
  completed: number;
  running: number;
  failed: number;
  progress_percent: number;
}

export interface DesolvationJobResponse {
  job_id: number;
  status: string;
  method_level: string;
  desolvation_mode: DesolvationMode;
  solvent_config?: SolventConfig;  // 溶剂模型配置
  created_at: string;
  started_at?: string;
  finished_at?: string;
  elapsed_seconds?: number;
  error_message?: string;
  result?: DesolvationEnergyResult;
  // 溯源信息
  solvation_structure_id?: number;
  composition_key?: string;  // 如 "Li-EC2-DMC1-PF6_1"
  md_job_id?: number;
  electrolyte_name?: string;
  qc_progress?: QCProgress;
}

export interface BatchDesolvationJobCreate {
  md_job_id: number;
  structure_ids: number[];
  method_level: string;
  desolvation_mode?: DesolvationMode;
  solvent_config?: SolventConfig;
  // Slurm 资源配置
  slurm_partition?: string;
  slurm_cpus?: number;
  slurm_time?: number;
}

export interface BatchDesolvationJobResponse {
  created_count: number;
  skipped_count: number;
  jobs: DesolvationJobResponse[];
}

export interface DesolvationOverviewResponse {
  md_job_id: number;
  electrolyte_name?: string;
  total_jobs: number;
  status_summary: Record<string, number>;
  jobs: DesolvationJobResponse[];
}

