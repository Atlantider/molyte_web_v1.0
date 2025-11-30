/**
 * 数据可见性控制 API
 */
import client from './client';

// 可见性枚举
export enum DataVisibility {
  PRIVATE = 'PRIVATE',
  DELAYED = 'DELAYED',
  PUBLIC = 'PUBLIC',
  ADMIN_ONLY = 'ADMIN_ONLY',
}

// 可见性更新请求
export interface VisibilityUpdate {
  visibility: DataVisibility;
  delay_days?: number;
  anonymous_public?: boolean;
  allow_download?: boolean;
  reason?: string;
}

// 任务可见性信息
export interface JobVisibility {
  id: number;
  name: string;
  visibility: DataVisibility;
  visibility_delay_until: string | null;
  anonymous_public: boolean;
  allow_download: boolean;
  view_count: number;
  download_count: number;
  reward_claimed: boolean;
  is_free_quota: boolean;
  created_at: string;
  // 管理员视图额外字段
  user_id?: number;
  username?: string;
}

// 可见性统计
export interface VisibilityStats {
  total: number;
  public: number;
  private: number;
  delayed: number;
  private_quota_used?: number;
  private_quota_limit?: number;
  max_delay_years?: number;
  contribution_points?: number;
  public_data_count?: number;
  balance_cpu_hours?: number;
  // 管理员统计额外字段
  admin_only?: number;
  soon_public?: number;
}

// 积分兑换比例信息
export interface ExchangeRateInfo {
  ratio: number;
  description: string;
  current_points: number;
  max_exchangeable_cpu_hours: number;
  current_balance: number;
}

// 积分兑换结果
export interface ExchangeResult {
  message: string;
  points_used: number;
  cpu_hours_gained: number;
  remaining_points: number;
  new_balance: number;
}

// 分页响应
export interface PaginatedResponse<T> {
  total: number;
  page: number;
  page_size: number;
  items: T[];
}

/**
 * 获取当前用户的任务可见性列表
 */
export const getMyJobsVisibility = async (
  visibility?: DataVisibility,
  page: number = 1,
  pageSize: number = 20
): Promise<PaginatedResponse<JobVisibility>> => {
  const params = new URLSearchParams();
  if (visibility) params.append('visibility', visibility);
  params.append('page', page.toString());
  params.append('page_size', pageSize.toString());
  
  const response = await client.get(`/visibility/my-jobs?${params.toString()}`);
  return response.data;
};

/**
 * 更新任务可见性
 */
export const updateJobVisibility = async (
  jobId: number,
  update: VisibilityUpdate
): Promise<{ message: string; visibility: string }> => {
  const response = await client.put(`/visibility/job/${jobId}`, update);
  return response.data;
};

/**
 * 领取公开数据奖励
 */
export const claimPublicReward = async (
  jobId: number
): Promise<{ message: string; reward_hours: number; new_balance: number }> => {
  const response = await client.post(`/visibility/job/${jobId}/claim-reward`);
  return response.data;
};

/**
 * 获取当前用户的可见性统计
 */
export const getMyVisibilityStats = async (): Promise<VisibilityStats> => {
  const response = await client.get('/visibility/my-stats');
  return response.data;
};

// ============ 管理员 API ============

/**
 * 管理员获取所有任务可见性列表
 */
export const adminGetAllJobsVisibility = async (
  visibility?: DataVisibility,
  userId?: number,
  page: number = 1,
  pageSize: number = 20
): Promise<PaginatedResponse<JobVisibility>> => {
  const params = new URLSearchParams();
  if (visibility) params.append('visibility', visibility);
  if (userId) params.append('user_id', userId.toString());
  params.append('page', page.toString());
  params.append('page_size', pageSize.toString());
  
  const response = await client.get(`/visibility/admin/all-jobs?${params.toString()}`);
  return response.data;
};

/**
 * 管理员更新任务可见性
 */
export const adminUpdateJobVisibility = async (
  jobId: number,
  update: VisibilityUpdate
): Promise<{ message: string; visibility: string }> => {
  const response = await client.put(`/visibility/admin/job/${jobId}`, update);
  return response.data;
};

/**
 * 管理员批量更新可见性
 */
export const adminBatchUpdateVisibility = async (
  jobIds: number[],
  visibility: DataVisibility,
  delayDays?: number,
  reason?: string
): Promise<{ message: string; updated_count: number }> => {
  const response = await client.post('/visibility/admin/batch-update', {
    job_ids: jobIds,
    visibility,
    delay_days: delayDays,
    reason,
  });
  return response.data;
};

/**
 * 管理员获取全局可见性统计
 */
export const adminGetVisibilityStats = async (): Promise<VisibilityStats> => {
  const response = await client.get('/visibility/admin/stats');
  return response.data;
};

// ============ 积分兑换 API ============

/**
 * 获取积分兑换比例信息
 */
export const getExchangeRate = async (): Promise<ExchangeRateInfo> => {
  const response = await client.get('/visibility/exchange-rate');
  return response.data;
};

/**
 * 积分兑换核时
 */
export const exchangePointsForCpuHours = async (
  points: number
): Promise<ExchangeResult> => {
  const response = await client.post('/visibility/exchange-points', { points });
  return response.data;
};
