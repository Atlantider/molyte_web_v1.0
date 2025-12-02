/**
 * QC量子化学计算相关类型定义
 */

// QC任务状态
export enum QCJobStatus {
  CREATED = 'CREATED',
  QUEUED = 'QUEUED',
  RUNNING = 'RUNNING',
  POSTPROCESSING = 'POSTPROCESSING',
  COMPLETED = 'COMPLETED',
  FAILED = 'FAILED',
  CANCELLED = 'CANCELLED',
}

// 分子类型
export enum MoleculeType {
  SOLVENT = 'solvent',
  CATION = 'cation',
  ANION = 'anion',
  CUSTOM = 'custom',
}

// QC任务
export interface QCJob {
  id: number;
  user_id: number;
  md_job_id?: number;
  molecule_name: string;
  smiles: string;
  molecule_type: string;
  basis_set: string;
  functional: string;
  charge: number;
  spin_multiplicity: number;
  status: string;
  slurm_job_id?: string;
  progress: number;
  work_dir?: string;
  log_file?: string;
  error_message?: string;
  config: Record<string, any>;
  created_at: string;
  updated_at?: string;
  started_at?: string;
  finished_at?: string;
  results?: QCResult[];
  // 额外字段
  accuracy_level?: string;
  solvent_config?: SolventConfig;
  slurm_partition?: string;
  slurm_cpus?: number;
  slurm_time?: number;
  // 复用已有计算结果
  is_reused?: boolean;
  reused_from_job_id?: number;
}

// 精度等级
export enum QCAccuracyLevel {
  FAST = 'fast',
  STANDARD = 'standard',
  ACCURATE = 'accurate',
  CUSTOM = 'custom',
}

// 溶剂模型
export enum SolventModel {
  GAS = 'gas',
  PCM = 'pcm',
  SMD = 'smd',
  CUSTOM = 'custom',
}

// 溶剂配置
export interface SolventConfig {
  model: string;
  solvent_name?: string;
  // 自定义溶剂参数
  eps?: number;
  eps_inf?: number;
  hbond_acidity?: number;
  hbond_basicity?: number;
  surface_tension?: number;
  carbon_aromaticity?: number;
  halogenicity?: number;
}

// QC任务创建参数
export interface QCJobCreate {
  molecule_name: string;
  smiles: string;
  molecule_type?: string;
  basis_set?: string;
  functional?: string;
  charge?: number;
  spin_multiplicity?: number;
  md_job_id?: number;
  accuracy_level?: string;
  solvent_config?: SolventConfig;
  auto_spin?: boolean;
  config?: Record<string, any>;
  // Slurm 资源配置
  slurm_partition?: string;
  slurm_cpus?: number;
  slurm_time?: number;
}

// 批量创建QC任务
export interface QCJobBatchCreate {
  molecules: QCJobCreate[];
  md_job_id?: number;
  basis_set?: string;
  functional?: string;
}

// QC计算结果
export interface QCResult {
  id: number;
  qc_job_id: number;
  smiles: string;
  energy_au?: number;
  homo?: number;
  lumo?: number;
  homo_lumo_gap?: number;
  esp_min_kcal?: number;
  esp_max_kcal?: number;
  dipole_moment?: number;
  polarizability?: number;
  esp_image_path?: string;
  homo_image_path?: string;
  lumo_image_path?: string;
  fchk_file_path?: string;
  log_file_path?: string;
  cube_density_path?: string;
  cube_esp_path?: string;
  cube_homo_path?: string;
  cube_lumo_path?: string;
  additional_properties: Record<string, any>;
  created_at: string;
  // 计算属性
  homo_ev?: number;
  lumo_ev?: number;
}

// 分子QC缓存
export interface MoleculeQCCache {
  id: number;
  smiles: string;
  molecule_name?: string;
  basis_set?: string;
  functional?: string;
  energy_au?: number;
  homo_ev?: number;
  lumo_ev?: number;
  homo_lumo_gap_ev?: number;
  esp_min_kcal?: number;
  esp_max_kcal?: number;
  esp_image_path?: string;
  homo_image_path?: string;
  lumo_image_path?: string;
  preferred_qc_result_id?: number;
  calculation_count: number;
  created_at: string;
  updated_at?: string;
}

// QC任务列表响应
export interface QCJobListResponse {
  total: number;
  jobs: QCJob[];
}

// QC搜索参数
export interface QCSearchParams {
  smiles?: string;
  molecule_name?: string;
  basis_set?: string;
  lumo_min?: number;
  lumo_max?: number;
  homo_min?: number;
  homo_max?: number;
  has_esp?: boolean;
  skip?: number;
  limit?: number;
}

// 基组选项
export interface BasisSetOption {
  value: string;
  label: string;
  description?: string;
}

// 泛函选项
export interface FunctionalOption {
  value: string;
  label: string;
  description?: string;
}
