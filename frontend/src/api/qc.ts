/**
 * QC量子化学计算API服务
 */
import apiClient from './client';
import type {
  QCJob,
  QCJobCreate,
  QCJobBatchCreate,
  QCResult,
  MoleculeQCCache,
  QCJobListResponse,
  QCSearchParams,
  BasisSetOption,
  FunctionalOption,
} from '../types/qc';

// ============================================================================
// QC Job APIs
// ============================================================================

/**
 * 获取QC任务列表
 */
export async function getQCJobs(params?: {
  status?: string;
  md_job_id?: number;
  molecule_name?: string;
  smiles?: string;
  functional?: string;
  basis_set?: string;
  visibility?: string;
  skip?: number;
  limit?: number;
}): Promise<QCJobListResponse> {
  const response = await apiClient.get('/qc/jobs', { params });
  return response.data;
}

/**
 * 获取QC任务详情
 */
export async function getQCJob(jobId: number): Promise<QCJob> {
  const response = await apiClient.get(`/qc/jobs/${jobId}`);
  return response.data;
}

/**
 * 获取QC任务状态（轻量级轮询）
 */
export async function getQCJobStatus(jobId: number): Promise<{
  id: number;
  status: string;
  progress: number;
  error_message?: string;
  slurm_job_id?: string;
  updated_at?: string;
}> {
  const response = await apiClient.get(`/qc/jobs/${jobId}/status`);
  return response.data;
}

/**
 * 创建QC任务
 */
export async function createQCJob(data: QCJobCreate): Promise<QCJob> {
  const response = await apiClient.post('/qc/jobs', data);
  return response.data;
}

/**
 * 批量创建QC任务
 */
export async function createQCJobsBatch(data: QCJobBatchCreate): Promise<QCJob[]> {
  const response = await apiClient.post('/qc/jobs/batch', data);
  return response.data;
}

/**
 * 提交QC任务到集群
 */
export async function submitQCJob(jobId: number): Promise<{
  message: string;
  job_id: number;
  celery_task_id: string;
}> {
  const response = await apiClient.post(`/qc/jobs/${jobId}/submit`);
  return response.data;
}

/**
 * 编辑QC任务（仅CREATED状态可编辑）
 */
export async function updateQCJob(jobId: number, data: Partial<QCJobCreate>): Promise<QCJob> {
  const response = await apiClient.put(`/qc/jobs/${jobId}`, data);
  return response.data;
}

/**
 * 删除/取消QC任务
 */
export async function deleteQCJob(jobId: number): Promise<{
  message: string;
  id: number;
}> {
  const response = await apiClient.delete(`/qc/jobs/${jobId}`);
  return response.data;
}

/**
 * 批量删除QC任务
 */
export async function batchDeleteQCJobs(ids: number[]): Promise<{
  deleted_count: number;
  cancelled_count: number;
  failed_ids: number[];
  message: string;
}> {
  const response = await apiClient.delete('/qc/jobs/batch/delete', { data: ids });
  return response.data;
}

/**
 * 批量提交QC任务
 */
export async function batchSubmitQCJobs(ids: number[]): Promise<{
  success_count: number;
  failed_count: number;
  errors: Array<{ job_id: number; error: string }>;
  message: string;
}> {
  const response = await apiClient.post('/qc/jobs/batch/submit', ids);
  return response.data;
}

/**
 * 批量取消QC任务
 */
export async function batchCancelQCJobs(ids: number[]): Promise<{
  success_count: number;
  failed_count: number;
  errors: Array<{ job_id: number; error: string }>;
  message: string;
}> {
  const response = await apiClient.post('/qc/jobs/batch/cancel', ids);
  return response.data;
}

/**
 * 重新计算QC任务（基于已有任务创建新任务）
 */
export async function recalculateQCJob(
  jobId: number,
  params: {
    functional?: string;
    basis_set?: string;
    solvent_config?: {
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
    };
    accuracy_level?: string;
    slurm_partition?: string;
    slurm_cpus?: number;
    slurm_time?: number;
  }
): Promise<QCJob> {
  const response = await apiClient.post(`/qc/jobs/${jobId}/recalculate`, params);
  return response.data;
}

// ============================================================================
// QC Results APIs
// ============================================================================

/**
 * 获取QC任务的计算结果
 */
export async function getQCResults(jobId: number): Promise<QCResult[]> {
  const response = await apiClient.get(`/qc/results/${jobId}`);
  return response.data;
}

/**
 * 根据SMILES查询QC结果
 */
export async function getQCResultsBySmiles(
  smiles: string,
  basisSet?: string
): Promise<QCResult[]> {
  const response = await apiClient.get('/qc/results/by-smiles', {
    params: { smiles, basis_set: basisSet },
  });
  return response.data;
}

// ============================================================================
// Molecule QC Cache APIs
// ============================================================================

/**
 * 获取分子的QC缓存数据
 */
export async function getMoleculeQCCache(smiles: string): Promise<MoleculeQCCache> {
  const response = await apiClient.get(`/qc/cache/${encodeURIComponent(smiles)}`);
  return response.data;
}

/**
 * 搜索分子QC缓存
 */
export async function searchMoleculeQCCache(params: QCSearchParams): Promise<{
  total: number;
  data: MoleculeQCCache[];
}> {
  const response = await apiClient.get('/qc/cache', { params });
  return response.data;
}

// ============================================================================
// Utility APIs
// ============================================================================

/**
 * 获取可用的基组列表
 */
export async function getBasisSets(): Promise<{ basis_sets: BasisSetOption[] }> {
  const response = await apiClient.get('/qc/config/basis-sets');
  return response.data;
}

/**
 * 获取可用的泛函列表
 */
export async function getFunctionals(): Promise<{ functionals: FunctionalOption[] }> {
  const response = await apiClient.get('/qc/config/functionals');
  return response.data;
}

/**
 * 获取精度等级列表
 */
export async function getAccuracyLevels(): Promise<{
  levels: Array<{
    value: string;
    label: string;
    functional: string | null;
    basis_set: string | null;
    description: string;
    estimated_time: string;
  }>;
}> {
  const response = await apiClient.get('/qc/config/accuracy-levels');
  return response.data;
}

/**
 * 获取溶剂模型列表
 */
export async function getSolventModels(): Promise<{
  models: Array<{
    value: string;
    label: string;
    description: string;
  }>;
}> {
  const response = await apiClient.get('/qc/config/solvent-models');
  return response.data;
}

/**
 * 获取可用溶剂列表
 */
export async function getSolvents(): Promise<{
  solvents: Array<{
    value: string;
    label: string;
    eps: number;
    description: string;
  }>;
}> {
  const response = await apiClient.get('/qc/config/solvents');
  return response.data;
}

/**
 * 获取常用分子列表
 */
export async function getCommonMolecules(): Promise<{
  categories: Array<{
    name: string;
    molecules: Array<{
      name: string;
      smiles: string;
      charge: number;
    }>;
  }>;
}> {
  const response = await apiClient.get('/qc/config/common-molecules');
  return response.data;
}

/**
 * 获取自定义溶剂参数说明
 */
export async function getCustomSolventParamsInfo(): Promise<{
  description: string;
  parameters: Array<{
    name: string;
    label: string;
    description: string;
    example: number;
    unit: string;
    range: string;
  }>;
  example_solvents: Record<string, Record<string, number>>;
}> {
  const response = await apiClient.get('/qc/config/custom-solvent-params');
  return response.data;
}

/**
 * 计算自旋多重度
 */
export async function calculateSpinMultiplicity(
  smiles: string,
  charge: number = 0
): Promise<{
  smiles: string;
  charge: number;
  total_electrons: number;
  num_radical_electrons: number;
  spin_multiplicity: number;
  description: string;
}> {
  const response = await apiClient.post('/qc/calculate-spin', null, {
    params: { smiles, charge },
  });
  return response.data;
}

// ============================================================================
// 重复计算检查 APIs
// ============================================================================

export interface MoleculeCheckRequest {
  smiles: string;
  molecule_name?: string;
  functional: string;
  basis_set: string;
  solvent_model: string;
  solvent_name?: string;
  charge: number;
  spin_multiplicity: number;
}

export interface MoleculeCheckResult {
  smiles: string;
  molecule_name?: string;
  has_existing_result: boolean;
  existing_qc_job_id?: number;
  existing_result_id?: number;
  functional?: string;
  basis_set?: string;
  solvent_model?: string;
  solvent_name?: string;
  energy_au?: number;
  homo_ev?: number;
  lumo_ev?: number;
  homo_lumo_gap_ev?: number;
  completed_at?: string;
}

export interface DuplicateCheckResponse {
  total_molecules: number;
  existing_count: number;
  new_count: number;
  results: MoleculeCheckResult[];
}

/**
 * 检查分子是否已有相同参数的QC计算结果（全局共享）
 */
export async function checkDuplicateCalculations(
  molecules: MoleculeCheckRequest[]
): Promise<DuplicateCheckResponse> {
  const response = await apiClient.post('/qc/check-duplicates', { molecules });
  return response.data;
}


// ============================================================================
// ESP Image APIs
// ============================================================================

/**
 * 获取ESP图片（返回Blob）
 */
export async function getESPImage(resultId: number): Promise<Blob> {
  const response = await apiClient.get(`/qc/esp-image/${resultId}`, {
    responseType: 'blob',
  });
  return response.data;
}

/**
 * 获取ESP图片的URL（带认证token）
 */
export function getESPImageUrl(resultId: number): string {
  const token = localStorage.getItem('token');
  return `/api/v1/qc/esp-image/${resultId}?token=${token}`;
}

/**
 * 下载ESP图片
 */
export async function downloadESPImage(resultId: number, filename?: string): Promise<void> {
  try {
    const blob = await getESPImage(resultId);
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename || `esp_${resultId}.png`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
  } catch (error) {
    console.error('下载ESP图片失败:', error);
    throw error;
  }
}

/**
 * 获取HOMO轨道图片（返回Blob）
 */
export async function getHOMOImage(resultId: number): Promise<Blob> {
  const response = await apiClient.get(`/qc/homo-image/${resultId}`, {
    responseType: 'blob',
  });
  return response.data;
}

/**
 * 获取LUMO轨道图片（返回Blob）
 */
export async function getLUMOImage(resultId: number): Promise<Blob> {
  const response = await apiClient.get(`/qc/lumo-image/${resultId}`, {
    responseType: 'blob',
  });
  return response.data;
}

/**
 * 下载HOMO轨道图片
 */
export async function downloadHOMOImage(resultId: number, filename?: string): Promise<void> {
  try {
    const blob = await getHOMOImage(resultId);
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename || `homo_${resultId}.png`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
  } catch (error) {
    console.error('下载HOMO图片失败:', error);
    throw error;
  }
}

/**
 * 下载LUMO轨道图片
 */
export async function downloadLUMOImage(resultId: number, filename?: string): Promise<void> {
  try {
    const blob = await getLUMOImage(resultId);
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename || `lumo_${resultId}.png`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
  } catch (error) {
    console.error('下载LUMO图片失败:', error);
    throw error;
  }
}

// ============================================================================
// Cluster Statistics APIs
// ============================================================================

/**
 * 能级统计数据
 */
export interface OrbitalStatistics {
  count: number;
  mean: number;
  std: number;
  min: number;
  max: number;
  percentile_5: number;
  percentile_95: number;
  values: number[];
}

/**
 * 电化学窗口估计
 */
export interface ElectrochemicalWindowEstimate {
  oxidation_limit_ev: number;
  reduction_limit_ev: number;
  window_ev: number;
  note: string;
}

/**
 * Cluster 统计响应
 */
export interface ClusterStatisticsResponse {
  md_job_id: number | null;
  total_qc_jobs: number;
  total_with_orbital_data?: number;
  message?: string;

  // 整体统计
  homo_statistics: OrbitalStatistics | null;
  lumo_statistics: OrbitalStatistics | null;
  gap_statistics: OrbitalStatistics | null;

  // 按类型分组
  per_type_statistics?: Record<string, {
    count: number;
    homo: OrbitalStatistics | null;
    lumo: OrbitalStatistics | null;
    gap: OrbitalStatistics | null;
  }>;

  // 按是否含 Li 分组
  with_li_statistics?: {
    homo: OrbitalStatistics | null;
    lumo: OrbitalStatistics | null;
    gap: OrbitalStatistics | null;
  } | null;
  without_li_statistics?: {
    homo: OrbitalStatistics | null;
    lumo: OrbitalStatistics | null;
    gap: OrbitalStatistics | null;
  } | null;

  // 电化学窗口估计
  electrochemical_window_estimate: ElectrochemicalWindowEstimate | null;
}

/**
 * 获取 Cluster QC 统计
 */
export async function getClusterStatistics(params?: {
  md_job_id?: number;
  include_single_molecule?: boolean;
}): Promise<ClusterStatisticsResponse> {
  const response = await apiClient.get('/qc/cluster-statistics', { params });
  return response.data;
}
