/**
 * 认证相关 API
 */
import client from './client';
import type { LoginRequest, LoginResponse, RegisterRequest, User } from '../types/api';

/**
 * 用户登录
 */
export const login = async (data: LoginRequest): Promise<LoginResponse> => {
  // FastAPI OAuth2PasswordRequestForm 需要 form-urlencoded 格式
  const formData = new URLSearchParams();
  formData.append('username', data.username);
  formData.append('password', data.password);
  
  const response = await client.post<LoginResponse>('/auth/login', formData, {
    headers: {
      'Content-Type': 'application/x-www-form-urlencoded',
    },
  });
  
  return response.data;
};

/**
 * 用户注册
 */
export const register = async (data: RegisterRequest): Promise<User> => {
  const response = await client.post<User>('/auth/register', data);
  return response.data;
};

/**
 * 获取当前用户信息
 */
export const getCurrentUser = async (): Promise<User> => {
  const response = await client.get<User>('/auth/me');
  return response.data;
};

/**
 * 退出登录
 */
export const logout = (): void => {
  localStorage.removeItem('access_token');
};

/**
 * 修改密码
 */
export const changePassword = async (oldPassword: string, newPassword: string): Promise<{ message: string }> => {
  const response = await client.post<{ message: string }>('/auth/change-password', {
    old_password: oldPassword,
    new_password: newPassword,
  });
  return response.data;
};

/**
 * 用户完整资料
 */
export interface UserProfile {
  id: number;
  username: string;
  email: string;
  role: string;
  user_type: string;
  organization: string | null;
  department: string | null;
  email_verified: boolean;
  is_active: boolean;
  created_at: string | null;
  last_login_at: string | null;
  // 配额配置
  daily_job_limit: number;
  concurrent_job_limit: number;
  storage_quota_gb: number;
  allowed_partitions: string[] | null;
  // 核时余额系统
  balance_cpu_hours: number;        // 账户余额
  frozen_cpu_hours: number;         // 冻结核时
  available_cpu_hours: number;      // 可用余额
  debt_cpu_hours: number;           // 欠费核时
  // 核时统计
  free_cpu_hours_granted: number;   // 初始赠送
  total_recharged: number;          // 总充值
  total_consumed: number;           // 总消耗
  total_cpu_hours: number;          // 总核时额度
  // 兼容旧字段
  used_cpu_hours: number;
  remaining_cpu_hours: number;
  // 贡献统计
  public_data_count: number;
  contribution_points: number;
  private_quota_used: number;
  private_quota_limit: number;
  // 使用情况
  today_jobs: number;
  running_jobs: number;
  // 任务统计
  total_jobs: number;
  completed_jobs: number;
  failed_jobs: number;
  // 使用率
  cpu_hours_usage_percent: number;
  daily_jobs_usage_percent: number;
}

/**
 * 获取当前用户完整资料
 */
export const getUserProfile = async (): Promise<UserProfile> => {
  const response = await client.get<UserProfile>('/auth/me/profile');
  return response.data;
};

