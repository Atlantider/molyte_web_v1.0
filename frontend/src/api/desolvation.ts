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
