/**
 * 计算任务管理 API
 */
import client from './client';
import type { MDJob, MDJobCreate, PostprocessJob, MDJobQCJobsResponse } from '../types';

/**
 * 获取所有 MD 任务
 */
export const getMDJobs = async (): Promise<MDJob[]> => {
  const response = await client.get('/jobs/');
  return response.data;
};

/**
 * 获取单个 MD 任务
 */
export const getMDJob = async (id: number): Promise<MDJob> => {
  const response = await client.get(`/jobs/${id}`);
  return response.data;
};

/**
 * 创建 MD 任务
 */
export const createMDJob = async (data: MDJobCreate): Promise<MDJob> => {
  const response = await client.post('/jobs/', data);
  return response.data;
};

/**
 * 批量创建 MD 任务
 */
export const batchCreateMDJobs = async (systemIds: number[], jobData: Omit<MDJobCreate, 'system_id'>): Promise<any> => {
  const response = await client.post('/jobs/batch', {
    system_ids: systemIds,
    ...jobData
  });
  return response.data;
};

/**
 * 更新 MD 任务配置
 */
export const updateMDJobConfig = async (id: number, config: any): Promise<MDJob> => {
  const response = await client.put(`/jobs/${id}/config`, config);
  return response.data;
};

/**
 * 提交 MD 任务到集群
 */
export const submitJobToCluster = async (id: number): Promise<MDJob> => {
  const response = await client.post(`/jobs/${id}/submit`);
  return response.data;
};

/**
 * 取消 MD 任务
 */
export const cancelMDJob = async (id: number): Promise<MDJob> => {
  const response = await client.post(`/jobs/${id}/cancel`);
  return response.data;
};

/**
 * 重新提交失败的 MD 任务
 */
export const resubmitMDJob = async (id: number): Promise<MDJob> => {
  const response = await client.post(`/jobs/${id}/resubmit`);
  return response.data;
};

/**
 * 获取所有后处理任务
 */
export const getPostprocessJobs = async (): Promise<PostprocessJob[]> => {
  const response = await client.get('/jobs/postprocess');
  return response.data;
};

/**
 * 获取单个后处理任务
 */
export const getPostprocessJob = async (id: number): Promise<PostprocessJob> => {
  const response = await client.get(`/jobs/postprocess/${id}`);
  return response.data;
};

/**
 * 检查任务创建配额
 */
export const checkJobQuota = async (): Promise<{
  can_create: boolean;
  current_count: number;
  limit: number;
  remaining: number;
}> => {
  const response = await client.get('/jobs/quota/check');
  return response.data;
};

/**
 * 获取任务的 Slurm 状态
 */
export interface SlurmJobStatus {
  job_id: number;
  slurm_job_id: string | null;
  status: string;
  raw_state?: string;
  exit_code?: string;
  start_time?: string;
  end_time?: string;
  elapsed?: string;
  cpu_time?: string;
  job_status?: string;
  message?: string;
}

export const getJobSlurmStatus = async (id: number): Promise<SlurmJobStatus> => {
  const response = await client.get(`/jobs/${id}/slurm_status`);
  return response.data;
};

/**
 * 同步任务状态（从 Slurm 更新到数据库）
 */
export const syncJobStatus = async (id: number): Promise<{
  job_id: number;
  slurm_job_id: string;
  slurm_status: string;
  job_status: string;
  progress: number;
  updated: boolean;
}> => {
  const response = await client.post(`/jobs/${id}/sync_status`);
  return response.data;
};

/**
 * 获取结构信息（密度、浓度等）
 */
export interface StructureInfo {
  available: boolean;
  message?: string;
  sample_name?: string;
  box_dimensions?: string;
  density?: number;
  concentration?: number;
  initial_density?: number;
  initial_concentration?: number;
  initial_box_dimensions?: string;
}

export const getStructureInfo = async (id: number): Promise<StructureInfo> => {
  const response = await client.get(`/jobs/${id}/structure_info`);
  return response.data;
};

// ============== 溶剂化结构 API ==============

/**
 * 溶剂化结构数据类型
 */
export interface SolvationStructure {
  id: number;
  center_ion: string;
  structure_type: string;
  coordination_num: number;
  composition: Record<string, number>;
  file_path: string | null;
  snapshot_frame: number;
  description: string;
  created_at: string;
}

/**
 * 溶剂化结构统计数据类型
 */
export interface SolvationStatistics {
  total_count: number;
  average_coordination_number: number;
  coordination_distribution: Record<string, number>;
  composition_distribution: Record<string, number>;
  molecule_counts: Record<string, number>;
  anion_coordination_distribution?: Record<string, number>;
}

/**
 * 体系结构数据类型
 */
export interface SystemStructure {
  frame_index: number;
  total_frames: number;
  atom_count: number;
  box: number[];
  xyz_content: string;
}

/**
 * 溶剂化结构XYZ内容类型
 */
export interface SolvationStructureContent {
  id: number;
  center_ion: string;
  coordination_num: number;
  composition: Record<string, number>;
  xyz_content: string;
  filename: string;
}

/**
 * 获取溶剂化结构列表
 */
export const getSolvationStructures = async (jobId: number): Promise<SolvationStructure[]> => {
  const response = await client.get(`/jobs/${jobId}/solvation`);
  return response.data;
};

/**
 * 刷新/重新计算溶剂化结构
 */
export const refreshSolvationStructures = async (jobId: number, cutoff: number = 3.0): Promise<{
  success: boolean;
  count: number;
  statistics?: SolvationStatistics;
  message: string;
}> => {
  const response = await client.post(`/jobs/${jobId}/solvation/refresh`, null, {
    params: { cutoff }
  });
  return response.data;
};

/**
 * 获取溶剂化结构统计信息
 */
export const getSolvationStatistics = async (jobId: number): Promise<SolvationStatistics> => {
  const response = await client.get(`/jobs/${jobId}/solvation/statistics`);
  return response.data;
};

/**
 * 获取溶剂化结构文件下载 URL
 */
export const getSolvationStructureFileUrl = (jobId: number, structureId: number): string => {
  return `/api/v1/jobs/${jobId}/solvation/structure/${structureId}`;
};

/**
 * 获取溶剂化结构 XYZ 内容（用于3D可视化）
 */
export const getSolvationStructureContent = async (jobId: number, structureId: number): Promise<SolvationStructureContent> => {
  const response = await client.get(`/jobs/${jobId}/solvation/structure/${structureId}`, {
    params: { format: 'content' }
  });
  return response.data;
};

/**
 * 获取整个体系的结构（用于3D可视化）
 */
export const getSystemStructure = async (jobId: number, frame: number = -1): Promise<SystemStructure> => {
  const response = await client.get(`/jobs/${jobId}/solvation/system-structure`, {
    params: { frame }
  });
  return response.data;
};

/**
 * 获取轨迹帧数
 */
export const getFrameCount = async (jobId: number): Promise<{ frame_count: number }> => {
  const response = await client.get(`/jobs/${jobId}/solvation/frame-count`);
  return response.data;
};

/**
 * 导出溶剂化数据
 */
export const exportSolvationData = async (jobId: number, format: 'json' | 'csv' = 'json'): Promise<Blob | any> => {
  const response = await client.get(`/jobs/${jobId}/solvation/export-data`, {
    params: { format },
    responseType: format === 'csv' ? 'blob' : 'json'
  });
  return response.data;
};

// ============== QC任务关联 API ==============

/**
 * 获取MD任务关联的QC任务列表和状态汇总
 */
export const getMDJobQCJobs = async (jobId: number): Promise<MDJobQCJobsResponse> => {
  const response = await client.get(`/jobs/${jobId}/qc-jobs`);
  return response.data;
};
