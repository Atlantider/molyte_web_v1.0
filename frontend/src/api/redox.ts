/**
 * 热力学循环计算氧化还原电位 API
 *
 * ⚠️ 高风险警告：
 * - 结果对方法/基组/溶剂模型/构型高度敏感
 * - 计算量大，经常不收敛
 * - 数值可能存在数百 mV 的系统性偏差
 * - 仅供研究参考，不应作为定量预测
 */
import apiClient from './client';

// ============================================================================
// 物理常数
// ============================================================================

export const PhysicalConstants = {
  FARADAY_C_MOL: 96485.33,
  FARADAY_KCAL_MOL_V: 23.061,
  HARTREE_TO_KCAL: 627.509,
  HARTREE_TO_EV: 27.2114,
  LI_ABSOLUTE_POTENTIAL_VS_SHE: -3.04,  // V
  SHE_ABSOLUTE_POTENTIAL: 4.44,  // V vs vacuum
};

// ============================================================================
// 枚举类型
// ============================================================================

export type RedoxCalculationMode = 'cheap' | 'standard' | 'heavy';
export type RedoxJobStatus = 'CREATED' | 'SUBMITTED' | 'RUNNING' | 'COMPLETED' | 'FAILED';
export type RedoxType = 'oxidation' | 'reduction';
export type ReorgEnergyJobStatus = 'CREATED' | 'SUBMITTED' | 'RUNNING' | 'COMPLETED' | 'FAILED';

// ============================================================================
// 热力学循环类型定义
// ============================================================================

export interface SpeciesConfig {
  name: string;
  qc_job_id?: number;  // 基于已有 QC 任务（推荐）
  smiles?: string;
  xyz_content?: string;
  charge: number;
  multiplicity: number;
  redox_type: RedoxType;
}

// ============================================================================
// 可用 Cluster 类型定义
// ============================================================================

export interface ClusterTypeInfo {
  type_name: string;
  count: number;
  charge: number;
  multiplicity: number;
  example_smiles: string;
  molecule_type: string;
  energy_mean_au?: number;
  energy_std_au?: number;
}

export interface ClusterInfo {
  qc_job_id: number;
  molecule_name: string;
  smiles: string;
  charge: number;
  multiplicity: number;
  functional: string;
  basis_set: string;
  molecule_type: string;
  energy_au?: number;
  homo_ev?: number;
  lumo_ev?: number;
  xyz_content?: string;
}

export interface AvailableClustersForRedoxResponse {
  md_job_id: number;
  total_clusters: number;
  cluster_types: ClusterTypeInfo[];
  clusters_by_type?: Record<string, ClusterInfo[]>;
  note: string;
}

export interface RedoxJobConfig {
  species_list: SpeciesConfig[];
  mode: RedoxCalculationMode;
  functional: string;
  basis_set: string;
  solvent_model: string;
  solvent: string;
  use_dispersion: boolean;
  scf_max_cycles: number;
  opt_max_cycles: number;
  li_reference_potential: number;
  reuse_existing_qc: boolean;  // 是否复用已有 QC 结果
}

export interface SpeciesRedoxResult {
  name: string;
  redox_type: RedoxType;
  e_neutral_gas?: number;
  e_charged_gas?: number;
  e_neutral_sol?: number;
  e_charged_sol?: number;
  dg_gas_kcal?: number;
  dg_solv_neutral_kcal?: number;
  dg_solv_charged_kcal?: number;
  dg_sol_kcal?: number;
  e_abs_v?: number;
  e_vs_li_v?: number;
  converged: boolean;
  warnings: string[];
  qc_job_ids: number[];
}

export interface RedoxJobResult {
  species_results: SpeciesRedoxResult[];
  oxidation_potentials_v: number[];
  reduction_potentials_v: number[];
  oxidation_limit_v?: number;
  reduction_limit_v?: number;
  electrochemical_window_v?: number;
  global_warnings: string[];
  calculation_mode: string;
  reference_note: string;
}

export interface RedoxJobResponse {
  id: number;
  md_job_id?: number;
  status: RedoxJobStatus;
  progress: number;
  error_message?: string;
  config?: RedoxJobConfig;
  result?: RedoxJobResult;
  created_at: string;
  updated_at?: string;
  started_at?: string;
  finished_at?: string;
}

export interface RedoxJobListResponse {
  total: number;
  jobs: RedoxJobResponse[];
}

// ============================================================================
// 重组能类型定义
// ============================================================================

export interface ReorgSpeciesConfig {
  name: string;
  qc_job_id?: number;  // 基于已有 QC 任务（推荐）
  smiles?: string;
  xyz_content?: string;
  charge_neutral: number;
  charge_oxidized: number;
  multiplicity_neutral: number;
  multiplicity_oxidized: number;
}

export interface ReorgEnergyJobConfig {
  species_list: ReorgSpeciesConfig[];
  functional: string;
  basis_set: string;
  use_dispersion: boolean;
  scf_max_cycles: number;
  opt_max_cycles: number;
  reuse_existing_qc: boolean;  // 是否复用已有 QC 结果
}

export interface SpeciesReorgResult {
  name: string;
  e_neutral_at_neutral_geom?: number;
  e_oxidized_at_neutral_geom?: number;
  e_neutral_at_oxidized_geom?: number;
  e_oxidized_at_oxidized_geom?: number;
  lambda_ox_ev?: number;
  lambda_red_ev?: number;
  lambda_total_ev?: number;
  converged: boolean;
  warnings: string[];
  qc_job_ids: number[];
}

export interface ReorgEnergyJobResult {
  species_results: SpeciesReorgResult[];
  lambda_ox_mean_ev?: number;
  lambda_red_mean_ev?: number;
  lambda_total_mean_ev?: number;
  global_warnings: string[];
  reference_note: string;
}

export interface ReorgEnergyJobResponse {
  id: number;
  md_job_id?: number;
  status: ReorgEnergyJobStatus;
  progress: number;
  error_message?: string;
  config?: ReorgEnergyJobConfig;
  result?: ReorgEnergyJobResult;
  created_at: string;
  updated_at?: string;
  started_at?: string;
  finished_at?: string;
}

// ============================================================================
// 热力学循环 API 函数
// ============================================================================

/**
 * 创建热力学循环任务
 */
export async function createRedoxJob(data: {
  md_job_id?: number;
  config: Partial<RedoxJobConfig>;
}): Promise<RedoxJobResponse> {
  const config: RedoxJobConfig = {
    species_list: data.config.species_list || [],
    mode: data.config.mode || 'cheap',
    functional: data.config.functional || 'B3LYP',
    basis_set: data.config.basis_set || '6-31G*',
    solvent_model: data.config.solvent_model || 'SMD',
    solvent: data.config.solvent || 'water',
    use_dispersion: data.config.use_dispersion ?? true,
    scf_max_cycles: data.config.scf_max_cycles || 200,
    opt_max_cycles: data.config.opt_max_cycles || 100,
    li_reference_potential: data.config.li_reference_potential || PhysicalConstants.LI_ABSOLUTE_POTENTIAL_VS_SHE,
    reuse_existing_qc: data.config.reuse_existing_qc ?? true,
  };

  const response = await apiClient.post('/redox/jobs', {
    md_job_id: data.md_job_id,
    config,
  });
  return response.data;
}

/**
 * 获取热力学循环任务列表
 */
export async function listRedoxJobs(params?: {
  md_job_id?: number;
  status?: RedoxJobStatus;
  skip?: number;
  limit?: number;
}): Promise<RedoxJobListResponse> {
  const response = await apiClient.get('/redox/jobs', { params });
  return response.data;
}

/**
 * 获取热力学循环任务详情
 */
export async function getRedoxJob(jobId: number): Promise<RedoxJobResponse> {
  const response = await apiClient.get(`/redox/jobs/${jobId}`);
  return response.data;
}

/**
 * 提交热力学循环任务
 */
export async function submitRedoxJob(jobId: number): Promise<RedoxJobResponse> {
  const response = await apiClient.post(`/redox/jobs/${jobId}/submit`);
  return response.data;
}

/**
 * 删除热力学循环任务
 */
export async function deleteRedoxJob(jobId: number): Promise<{ message: string; job_id: number }> {
  const response = await apiClient.delete(`/redox/jobs/${jobId}`);
  return response.data;
}

// ============================================================================
// 重组能 API 函数
// ============================================================================

/**
 * 获取重组能任务列表
 */
export async function listReorgEnergyJobs(params?: {
  md_job_id?: number;
  status?: ReorgEnergyJobStatus;
  skip?: number;
  limit?: number;
}): Promise<ReorgEnergyJobResponse[]> {
  const response = await apiClient.get('/redox/reorganization-energy/jobs', { params });
  return response.data;
}

/**
 * 创建重组能任务
 */
export async function createReorgEnergyJob(data: {
  md_job_id?: number;
  config: Partial<ReorgEnergyJobConfig>;
}): Promise<ReorgEnergyJobResponse> {
  const config: ReorgEnergyJobConfig = {
    species_list: data.config.species_list || [],
    functional: data.config.functional || 'B3LYP',
    basis_set: data.config.basis_set || '6-31G*',
    use_dispersion: data.config.use_dispersion ?? true,
    scf_max_cycles: data.config.scf_max_cycles || 200,
    opt_max_cycles: data.config.opt_max_cycles || 150,
    reuse_existing_qc: data.config.reuse_existing_qc ?? true,
  };

  const response = await apiClient.post('/redox/reorganization-energy/jobs', {
    md_job_id: data.md_job_id,
    config,
  });
  return response.data;
}

/**
 * 获取重组能任务详情
 */
export async function getReorgEnergyJob(jobId: number): Promise<ReorgEnergyJobResponse> {
  const response = await apiClient.get(`/redox/reorganization-energy/jobs/${jobId}`);
  return response.data;
}

/**
 * 提交重组能任务
 */
export async function submitReorgEnergyJob(jobId: number): Promise<ReorgEnergyJobResponse> {
  const response = await apiClient.post(`/redox/reorganization-energy/jobs/${jobId}/submit`);
  return response.data;
}

/**
 * 删除重组能任务
 */
export async function deleteReorgEnergyJob(jobId: number): Promise<{ message: string; job_id: number }> {
  const response = await apiClient.delete(`/redox/reorganization-energy/jobs/${jobId}`);
  return response.data;
}

/**
 * 获取物理常数
 */
export async function getPhysicalConstants(): Promise<{
  faraday_c_mol: number;
  faraday_kcal_mol_v: number;
  hartree_to_kcal: number;
  hartree_to_ev: number;
  li_absolute_potential_vs_she: number;
  she_absolute_potential: number;
  notes: {
    li_reference: string;
    she_reference: string;
    warning: string;
  };
}> {
  const response = await apiClient.get('/redox/constants');
  return response.data;
}

/**
 * 获取可用于 Redox/重组能计算的 Cluster 列表
 * @param mdJobId MD 任务 ID
 * @param includeXyz 是否包含 XYZ 结构内容
 */
export async function getAvailableClustersForRedox(
  mdJobId: number,
  includeXyz: boolean = false
): Promise<AvailableClustersForRedoxResponse> {
  const response = await apiClient.get(`/qc/available-clusters-for-redox/${mdJobId}`, {
    params: { include_xyz: includeXyz },
  });
  return response.data;
}
