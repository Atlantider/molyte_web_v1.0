/**
 * Cluster Analysis API - 统一的 Cluster 高级计算规划
 */
import apiClient from './client';

// ============================================================================
// 类型定义
// ============================================================================

export type ClusterCalcType =
  | 'BINDING_TOTAL'
  | 'BINDING_PAIRWISE'
  | 'DESOLVATION_STEPWISE'
  | 'DESOLVATION_FULL'
  | 'REDOX'
  | 'REORGANIZATION';

export type AdvancedClusterJobStatus =
  | 'CREATED'
  | 'SUBMITTED'
  | 'RUNNING'
  | 'WAITING_QC'
  | 'CALCULATING'
  | 'COMPLETED'
  | 'FAILED'
  | 'CANCELLED';

export interface QCConfig {
  functional: string;
  basis_set: string;
  solvent_model?: string;
  solvent?: string;
  use_dispersion: boolean;
  charge_cluster: number;
  charge_ion: number;
}

export interface PlannedQCTask {
  task_type: string;
  description: string;
  smiles?: string;
  structure_id?: number;
  charge: number;
  multiplicity: number;
  calc_mode: 'opt' | 'sp';  // 计算模式: opt (几何优化) / sp (单点能量)
  status: 'new' | 'reused' | 'local_reused' | 'pending';
  existing_qc_job_id?: number;
  existing_energy?: number;
}

export interface CalcTypeRequirements {
  calc_type: ClusterCalcType;
  description: string;
  required_qc_tasks: PlannedQCTask[];
  new_tasks_count: number;
  reused_tasks_count: number;
}

export interface RedoxOptions {
  include_molecule: boolean;
  include_dimer: boolean;
  include_cluster: boolean;
}

export interface ReorganizationOptions {
  include_molecule: boolean;
  include_cluster: boolean;
}

export interface ClusterAnalysisPlanRequest {
  md_job_id: number;
  solvation_structure_ids?: number[];
  composition_keys?: string[];
  calc_types: ClusterCalcType[];
  redox_options?: RedoxOptions;
  reorganization_options?: ReorganizationOptions;
  qc_config?: Partial<QCConfig>;
}

export interface ClusterAnalysisPlanResponse {
  md_job_id: number;
  selected_structures_count: number;
  selected_structure_ids: number[];
  calc_requirements: CalcTypeRequirements[];
  total_new_qc_tasks: number;
  total_reused_qc_tasks: number;
  estimated_compute_hours: number;
  warnings: string[];
}

export interface AdvancedClusterJob {
  id: number;
  md_job_id: number;
  user_id: number;
  status: AdvancedClusterJobStatus;
  progress: number;
  calc_types: string[];
  selected_structures: {
    solvation_structure_ids: number[];
    count: number;
  };
  qc_config: Record<string, unknown>;
  qc_task_plan: {
    planned_qc_tasks: PlannedQCTask[];
    reused_qc_jobs: number[];
    new_qc_jobs: number[];
    total_qc_tasks: number;
    completed_qc_tasks: number;
  };
  results: Record<string, unknown>;
  error_message?: string;
  created_at: string;
  updated_at: string;
  started_at?: string;
  finished_at?: string;
}

export interface AddCalcTypePlanResponse {
  job_id: number;
  existing_calc_types: string[];
  additional_calc_types: string[];
  new_qc_tasks_required: number;
  reused_from_existing: number;
  details: CalcTypeRequirements[];
}

// ============================================================================
// 计算类型描述
// ============================================================================

export const CALC_TYPE_INFO: Record<ClusterCalcType, {
  label: string;
  description: string;
  formula: string;
  riskLevel: 'low' | 'medium' | 'high';
  icon: string;
}> = {
  BINDING_TOTAL: {
    label: '总 Binding Energy',
    description: '计算整个溶剂化簇的总脱溶剂化能',
    formula: 'E_bind = E_cluster - (E_ion + Σ n_j × E_ligand_j)',
    riskLevel: 'low',
    icon: '🔗',
  },
  BINDING_PAIRWISE: {
    label: '分子-Li Binding',
    description: '计算单个分子与 Li+ 的 binding energy',
    formula: 'E_bind = E(Li-X) - E(Li) - E(X)',
    riskLevel: 'low',
    icon: '⚛️',
  },
  DESOLVATION_STEPWISE: {
    label: '逐级去溶剂化',
    description: '逐个移除配体，计算每步的去溶剂化能',
    formula: 'ΔE_i = E_cluster - (E_minus_i + E_ligand_i)',
    riskLevel: 'medium',
    icon: '📉',
  },
  DESOLVATION_FULL: {
    label: '完全去溶剂化',
    description: '计算从完整簇到裸离子的总去溶剂化能',
    formula: 'ΔE = E_cluster - (E_ion + Σ E_ligand_i)',
    riskLevel: 'low',
    icon: '🎯',
  },
  REDOX: {
    label: '氧化还原电位',
    description: '热力学循环计算氧化还原电位',
    formula: 'E° = -ΔG(sol) / nF',
    riskLevel: 'high',
    icon: '⚡',
  },
  REORGANIZATION: {
    label: 'Marcus 重组能',
    description: 'Marcus 理论 4 点方案计算重组能',
    formula: 'λ = (λ_ox + λ_red) / 2',
    riskLevel: 'high',
    icon: '🔄',
  },
};

// ============================================================================
// API 函数
// ============================================================================

export interface RecommendedSolventResponse {
  recommended_solvent: string;
  recommended_dielectric: number;
  average_dielectric: number;
  composition_analyzed: Record<string, { count: number; dielectric: number }>;
  reason: string;
}

export async function recommendPCMSolvent(
  mdJobId: number
): Promise<RecommendedSolventResponse> {
  const response = await apiClient.get(`/cluster-analysis/recommend-solvent/${mdJobId}`);
  return response.data;
}

export async function planClusterAnalysis(
  request: ClusterAnalysisPlanRequest
): Promise<ClusterAnalysisPlanResponse> {
  const response = await apiClient.post('/cluster-analysis/plan', request);
  return response.data;
}

export async function submitClusterAnalysis(
  request: ClusterAnalysisPlanRequest
): Promise<AdvancedClusterJob> {
  const response = await apiClient.post('/cluster-analysis/submit', request);
  return response.data;
}

export async function listClusterAnalysisJobs(
  mdJobId?: number,
  skip = 0,
  limit = 50
): Promise<AdvancedClusterJob[]> {
  const params: Record<string, unknown> = { skip, limit };
  if (mdJobId) params.md_job_id = mdJobId;
  const response = await apiClient.get('/cluster-analysis/jobs', { params });
  return response.data;
}

export async function getClusterAnalysisJob(jobId: number): Promise<AdvancedClusterJob> {
  const response = await apiClient.get(`/cluster-analysis/jobs/${jobId}`);
  return response.data;
}

export async function cancelClusterAnalysisJob(jobId: number): Promise<{
  message: string;
  job_id: number;
  cancelled_qc_count: number;
}> {
  const response = await apiClient.post(`/cluster-analysis/jobs/${jobId}/cancel`);
  return response.data;
}

export async function planAddCalcTypes(
  jobId: number,
  additionalCalcTypes: ClusterCalcType[]
): Promise<AddCalcTypePlanResponse> {
  const response = await apiClient.post(`/cluster-analysis/jobs/${jobId}/add-calc-types`, {
    job_id: jobId,
    additional_calc_types: additionalCalcTypes,
  });
  return response.data;
}

// ============================================================================
// 结果查询 API
// ============================================================================

export interface QCTaskInfo {
  task_type: string;
  name: string;
  description: string;
  smiles: string;
  charge: number;
  multiplicity: number;
  status: 'new' | 'reused' | 'local_reused';
  qc_job_id: number | null;
  qc_status: string | null;
  slurm_job_id?: string | null;
}

export interface QCStatus {
  job_id: number;
  total_qc_jobs: number;
  completed: number;
  running: number;
  pending: number;
  failed: number;
  all_completed: boolean;
  calc_types?: string[];
  tasks_by_calc_type?: Record<string, QCTaskInfo[]>;
  qc_jobs?: Array<{
    id: number;
    status: string;
    molecule_name: string;
    task_type?: string;
  }>;
}

export interface ClusterAnalysisResults {
  job_id: number;
  status: string;
  progress: number;
  calc_types: string[];
  results: Record<string, unknown>;
  qc_task_plan: Record<string, unknown>;
}

// Binding Total 结果
export interface BindingTotalResult {
  e_cluster_au: number;
  e_ion_au: number;
  ligand_energies_au: Record<string, number>;
  e_bind_au: number;
  e_bind_ev: number;
  e_bind_kcal_mol: number;
  error?: string;
}

// Binding Pairwise 结果
export interface BindingPairwiseResult {
  pairwise_bindings: Array<{
    ligand: string;
    e_dimer_au: number;
    e_ligand_au: number;
    e_ion_au: number;
    e_bind_au: number;
    e_bind_ev: number;
    e_bind_kcal_mol: number;
  }>;
}

// Desolvation Stepwise 结果
export interface DesolvationStepwiseResult {
  stepwise_desolvation: Array<{
    ligand: string;
    e_cluster_au: number;
    e_minus_au: number;
    e_ligand_au: number;
    delta_e_au: number;
    delta_e_ev: number;
    delta_e_kcal_mol: number;
  }>;
}

// Redox 结果
export interface RedoxResult {
  redox_potentials: Array<{
    smiles: string;
    e_neutral_gas_au: number;
    e_charged_gas_au: number;
    e_neutral_sol_au: number;
    e_charged_sol_au: number;
    delta_g_sol_ev: number;
    oxidation_potential_v: number;
  }>;
}

export async function getClusterAnalysisResults(jobId: number): Promise<ClusterAnalysisResults> {
  const response = await apiClient.get(`/cluster-analysis/jobs/${jobId}/results`);
  return response.data;
}

export async function getClusterAnalysisQCStatus(jobId: number): Promise<QCStatus> {
  const response = await apiClient.get(`/cluster-analysis/jobs/${jobId}/qc-status`);
  return response.data;
}

export async function calculateClusterAnalysisResults(jobId: number): Promise<{ status: string; results: Record<string, unknown> }> {
  const response = await apiClient.post(`/cluster-analysis/jobs/${jobId}/calculate`);
  return response.data;
}
