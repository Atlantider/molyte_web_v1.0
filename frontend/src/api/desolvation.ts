/**
 * Desolvation energy calculation API client
 * 去溶剂化能计算 API 客户端
 */
import client from './client';
import type {
  DesolvationJobCreate,
  DesolvationJobResponse,
  BatchDesolvationJobCreate,
  BatchDesolvationJobResponse,
  DesolvationOverviewResponse
} from '../types/desolvation';

/**
 * 创建去溶剂化能任务
 */
export const createDesolvationJob = async (
  data: DesolvationJobCreate
): Promise<DesolvationJobResponse> => {
  const response = await client.post('/desolvation/jobs', data);
  return response.data;
};

/**
 * 批量创建去溶剂化能任务
 */
export const batchCreateDesolvationJobs = async (
  data: BatchDesolvationJobCreate
): Promise<BatchDesolvationJobResponse> => {
  const response = await client.post('/desolvation/batch', data);
  return response.data;
};

/**
 * 获取去溶剂化能任务详情
 */
export const getDesolvationJob = async (
  jobId: number
): Promise<DesolvationJobResponse> => {
  const response = await client.get(`/desolvation/jobs/${jobId}`);
  return response.data;
};

/**
 * 获取某个 cluster 的所有去溶剂化能任务
 */
export const listClusterDesolvationJobs = async (
  clusterId: number
): Promise<DesolvationJobResponse[]> => {
  const response = await client.get(`/desolvation/cluster/${clusterId}/jobs`);
  return response.data;
};

/**
 * 获取某个 MD 任务下所有去溶剂化计算的总览
 */
export const getDesolvationOverview = async (
  mdJobId: number
): Promise<DesolvationOverviewResponse> => {
  const response = await client.get(`/desolvation/md/${mdJobId}/overview`);
  return response.data;
};

/**
 * QC子任务信息
 */
export interface QCTaskInfo {
  id: number;
  molecule_name: string;
  task_type: 'cluster' | 'cluster_minus' | 'ligand';
  status: string;
  progress: number;
  charge: number;
  spin_multiplicity: number;
  basis_set: string;
  functional: string;
  is_reused: boolean;
  reused_from_job_id?: number;
  slurm_job_id?: string;
  error_message?: string;
  created_at?: string;
  started_at?: string;
  finished_at?: string;
}

/**
 * 去溶剂化任务QC子任务响应
 */
export interface DesolvationQCTasksResponse {
  job_id: number;
  composition_key?: string;
  total: number;
  completed: number;
  running: number;
  failed: number;
  queued: number;
  reused: number;
  qc_tasks: QCTaskInfo[];
}

/**
 * 获取某个去溶剂化任务的 QC 子任务列表
 */
export const getDesolvationQCTasks = async (
  jobId: number
): Promise<DesolvationQCTasksResponse> => {
  const response = await client.get(`/desolvation/jobs/${jobId}/qc-tasks`);
  return response.data;
};

/**
 * 结构预览数据类型
 */
export interface StructurePreview {
  name: string;
  xyz_content: string;
  atom_count: number;
  charge: number;
}

export interface LigandPreview {
  ligand_id: string;
  ligand_type: string;
  ligand_label: string;
  xyz_content: string;
  atom_count: number;
  charge: number;
}

export interface ClusterMinusPreview {
  name: string;
  removed_ligand: string;
  removed_ligand_type: string;
  xyz_content: string;
  atom_count: number;
  charge: number;
}

export interface DimerPreview {
  name: string;
  ligand_type: string;
  xyz_content: string;
  atom_count: number;
  charge: number;
  source_ligand_id: number;
}

export interface DesolvationPreviewResponse {
  structure_id: number;
  cluster_name: string;
  center_ion: string;
  total_charge: number;
  composition: Record<string, number>;
  cluster: StructurePreview;
  ligands: LigandPreview[];
  cluster_minus_structures: ClusterMinusPreview[];
  dimer_structures?: DimerPreview[];  // Li-配体 dimer（用于 pairwise binding）
  center_ion_structure: StructurePreview;
}

/**
 * 预览去溶剂化结构（不创建任务）
 */
export const previewDesolvationStructures = async (
  structureId: number
): Promise<DesolvationPreviewResponse> => {
  const response = await client.get(`/desolvation/preview/${structureId}`);
  return response.data;
};

/**
 * Binding Energy 汇总
 */
export interface LigandBindingData {
  ligand_id: string;
  ligand_type: string;
  ligand_label: string;
  binding_energy_kcal: number | null;
  e_ligand_au: number | null;
  e_cluster_minus_au: number | null;
}

export interface TypeBindingStats {
  count: number;
  mean: number;
  std: number;
  min: number;
  max: number;
  values: number[];
  percentile_25?: number;
  percentile_75?: number;
}

export interface BindingSummaryResponse {
  job_id: number;
  method_level: string;
  e_cluster_au: number;
  solvation_info: {
    center_ion: string;
    coordination_num: number;
    composition_key: string;
  } | null;
  per_ligand_binding: LigandBindingData[];
  per_type_stats: Record<string, TypeBindingStats>;
  last_layer_binding: {
    ligand_type: string;
    ligand_label: string;
    binding_energy_kcal: number;
  } | null;
  total_ligands: number;
}

/**
 * 获取单个去溶剂化任务的 Binding Energy 汇总
 */
export const getBindingSummary = async (
  jobId: number
): Promise<BindingSummaryResponse> => {
  const response = await client.get(`/desolvation/jobs/${jobId}/binding-summary`);
  return response.data;
};

/**
 * Binding 统计总览
 */
export interface BindingStatisticsOverview {
  md_job_id: number;
  total_completed: number;
  total_ligand_bindings: number;
  per_type_statistics: Record<string, TypeBindingStats>;
  last_layer_statistics: TypeBindingStats & {
    details: Array<{
      job_id: number;
      composition_key: string;
      ligand_type: string;
      binding_energy_kcal: number;
    }>;
  } | null;
  message?: string;
}

/**
 * 获取 MD 任务下所有去溶剂化结果的 Binding 统计
 */
export const getBindingStatisticsOverview = async (
  mdJobId: number
): Promise<BindingStatisticsOverview> => {
  const response = await client.get(`/desolvation/overview/${mdJobId}/binding-statistics`);
  return response.data;
};
