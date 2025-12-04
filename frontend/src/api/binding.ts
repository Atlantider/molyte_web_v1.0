/**
 * Binding Analysis API
 * Li-配体 Binding Energy 分析相关的 API
 */
import apiClient from './client';

// ============================================================================
// Types
// ============================================================================

export type BindingAnalysisStatus = 'CREATED' | 'SUBMITTED' | 'RUNNING' | 'COMPLETED' | 'FAILED';

export interface BindingAnalysisConfig {
  composition_keys?: string[];
  functional?: string;
  basis_set?: string;
  solvent_model?: string;
  solvent_name?: string;
  reuse_existing_qc?: boolean;
}

export interface ClusterBindingResult {
  composition_key: string;
  cluster_energy_au?: number;
  center_ion_energy_au?: number;
  ligand_energies_au?: Record<string, number>;
  ligand_counts?: Record<string, number>;
  binding_energy_au?: number;
  binding_energy_kcal?: number;
  per_ligand_binding_kcal?: Record<string, number>;
  cluster_qc_job_id?: number;
  center_ion_qc_job_id?: number;
  ligand_qc_job_ids?: Record<string, number>;
  converged: boolean;
  warnings?: string[];
}

export interface BindingAnalysisSummary {
  total_clusters: number;
  completed_clusters: number;
  failed_clusters: number;
  mean_total_binding_kcal?: number;
  std_total_binding_kcal?: number;
  per_type_statistics?: Array<{
    ligand_type: string;
    count: number;
    mean_binding_kcal: number;
    std_binding_kcal: number;
    min_binding_kcal: number;
    max_binding_kcal: number;
  }>;
  warnings?: string[];
}

export interface BindingAnalysisResult {
  per_cluster_results: ClusterBindingResult[];
  summary: BindingAnalysisSummary;
}

export interface BindingAnalysisJob {
  id: number;
  md_job_id: number;
  user_id: number;
  status: BindingAnalysisStatus;
  progress: number;
  error_message?: string;
  config?: BindingAnalysisConfig;
  result?: BindingAnalysisResult;
  qc_job_ids?: number[];
  created_at: string;
  updated_at: string;
  started_at?: string;
  finished_at?: string;
}

export interface AvailableClustersInfo {
  md_job_id: number;
  total_solvation_structures: number;
  composition_keys: string[];
  clusters_by_composition: Record<string, Array<{
    id: number;
    frame: number;
    center_ion: string;
    ligand_types: Record<string, number>;
    has_xyz: boolean;
  }>>;
  existing_qc_by_type: Record<string, {
    count: number;
    examples: string[];
  }>;
  note: string;
}

// ============================================================================
// API Functions
// ============================================================================

/**
 * 创建 Binding 分析任务
 */
export async function createBindingAnalysisJob(data: {
  md_job_id: number;
  config?: BindingAnalysisConfig;
}): Promise<BindingAnalysisJob> {
  const response = await apiClient.post('/binding/jobs', data);
  return response.data;
}

/**
 * 获取 Binding 分析任务列表
 */
export async function getBindingAnalysisJobs(params?: {
  md_job_id?: number;
  status?: BindingAnalysisStatus;
  skip?: number;
  limit?: number;
}): Promise<{ items: BindingAnalysisJob[]; total: number }> {
  const response = await apiClient.get('/binding/jobs', { params });
  return response.data;
}

/**
 * 获取 Binding 分析任务详情
 */
export async function getBindingAnalysisJob(jobId: number): Promise<BindingAnalysisJob> {
  const response = await apiClient.get(`/binding/jobs/${jobId}`);
  return response.data;
}

/**
 * 提交 Binding 分析任务
 */
export async function submitBindingAnalysisJob(jobId: number): Promise<{ message: string; job_id: number; status: string }> {
  const response = await apiClient.post(`/binding/jobs/${jobId}/submit`);
  return response.data;
}

/**
 * 删除 Binding 分析任务
 */
export async function deleteBindingAnalysisJob(jobId: number): Promise<{ message: string; job_id: number }> {
  const response = await apiClient.delete(`/binding/jobs/${jobId}`);
  return response.data;
}

/**
 * 获取可用于 Binding 分析的 Cluster 列表
 */
export async function getAvailableClustersForBinding(mdJobId: number): Promise<AvailableClustersInfo> {
  const response = await apiClient.get(`/binding/available-clusters/${mdJobId}`);
  return response.data;
}

