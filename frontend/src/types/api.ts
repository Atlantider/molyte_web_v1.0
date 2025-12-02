/**
 * API 类型定义
 */

// 用户角色
export enum UserRole {
  USER = 'USER',
  ADMIN = 'ADMIN',
  PREMIUM = 'PREMIUM',
  GUEST = 'GUEST',
}

// 用户类型
export enum UserType {
  STUDENT = 'STUDENT',
  RESEARCHER = 'RESEARCHER',
  COMPANY = 'COMPANY',
}

// 数据可见性
export enum DataVisibility {
  PRIVATE = 'PRIVATE',
  DELAYED = 'DELAYED',
  PUBLIC = 'PUBLIC',
  ADMIN_ONLY = 'ADMIN_ONLY',
}

// 用户信息
export interface User {
  id: number;
  email: string;
  username: string;
  role: UserRole;
  user_type?: UserType;
  organization?: string;
  department?: string;
  email_verified?: boolean;
  balance_cpu_hours?: number;
  frozen_cpu_hours?: number;
  public_data_count?: number;
  contribution_points?: number;
  created_at: string;
  updated_at: string;
}

// 登录请求
export interface LoginRequest {
  username: string;
  password: string;
}

// 登录响应
export interface LoginResponse {
  access_token: string;
  token_type: string;
}

// 注册请求
export interface RegisterRequest {
  email: string;
  username: string;
  password: string;
  role?: UserRole;
  user_type?: string;
  organization?: string;
  department?: string;
}

// API 响应
export interface ApiResponse<T = any> {
  data?: T;
  message?: string;
  error?: string;
}

// 项目信息
export interface Project {
  id: number;
  name: string;
  description?: string;
  user_id: number;
  created_at: string;
  updated_at: string;
}

// 创建项目请求
export interface ProjectCreate {
  name: string;
  description?: string;
}

// 更新项目请求
export interface ProjectUpdate {
  name?: string;
  description?: string;
}

// 分子规格（旧版本，保留兼容性）
export interface MoleculeSpec {
  name: string;
  smiles: string;
  number: number;
}

// 离子规格（新版本 - 使用浓度）
export interface IonSpec {
  name: string;
  concentration: number; // mol/L
}

// 溶剂规格（新版本 - 使用摩尔比）
export interface SolventSpec {
  name: string;
  smiles: string;
  molar_ratio: number; // 相对于第一种阳离子的摩尔比
}

// 离子信息（包含电荷）
export interface IonInfo {
  name: string;
  charge: number;
}

// 可用离子列表
export interface AvailableIons {
  cations: IonInfo[];
  anions: IonInfo[];
}

// 电解质系统
export interface ElectrolyteSystem {
  id: number;
  project_id: number;
  hash_key: string;
  name: string;
  cations: MoleculeSpec[];
  anions: MoleculeSpec[];
  solvents?: MoleculeSpec[];
  temperature: number;
  pressure?: number;
  density?: number;
  concentration?: number;
  box_size?: number;
  nsteps_npt?: number;
  nsteps_nvt?: number;
  timestep?: number;
  force_field?: string;
  created_at: string;
  // 兼容旧的 composition 字段
  composition?: {
    cations: MoleculeSpec[];
    anions: MoleculeSpec[];
    solvents: MoleculeSpec[];
  };
}

// 创建电解质系统请求
export interface ElectrolyteSystemCreate {
  project_id: number;
  name: string;
  cations: MoleculeSpec[];
  anions: MoleculeSpec[];
  solvents?: MoleculeSpec[];
  temperature?: number;
  pressure?: number;
  density?: number;
  concentration?: number;
  box_size?: number;
  nsteps_npt?: number;
  nsteps_nvt?: number;
  timestep?: number;
  force_field?: string;
}

// 更新电解质系统请求
export interface ElectrolyteSystemUpdate {
  name?: string;
  temperature?: number;
  pressure?: number;
  density?: number;
  concentration?: number;
  box_size?: number;
  nsteps_npt?: number;
  nsteps_nvt?: number;
  timestep?: number;
  force_field?: string;
}

// 任务状态
export enum JobStatus {
  CREATED = 'CREATED',
  QUEUED = 'QUEUED',
  RUNNING = 'RUNNING',
  POSTPROCESSING = 'POSTPROCESSING',
  COMPLETED = 'COMPLETED',
  FAILED = 'FAILED',
  CANCELLED = 'CANCELLED',
}

// 后处理类型
export enum PostprocessType {
  RDF = 'RDF',
  MSD = 'MSD',
  SOLVATION = 'SOLVATION',
}

// MD 任务
export interface MDJob {
  id: number;
  system_id: number;
  user_id: number;
  status: JobStatus;
  slurm_job_id?: string;
  progress: number;
  work_dir?: string;
  log_file?: string;
  error_message?: string;
  config?: {
    job_name?: string;
    nsteps_npt?: number;
    nsteps_nvt?: number;
    timestep?: number;
    temperature?: number;
    pressure?: number;
    freq_trj_npt?: number;
    freq_trj_nvt?: number;
    thermo_freq?: number;
    // Slurm 资源配置
    slurm_partition?: string;
    slurm_nodes?: number;
    slurm_ntasks?: number;
    slurm_cpus_per_task?: number;
    slurm_time?: number;
    slurm_job_id?: string;  // Slurm 任务 ID（兼容旧数据）
    // QC 计算配置
    qc_enabled?: boolean;
    qc_accuracy_level?: string;
    qc_basis_set?: string;
    qc_functional?: string;
    qc_solvent_model?: string;
    qc_solvent_name?: string;
    qc_use_recommended_params?: boolean;
    // 自定义溶剂参数
    qc_custom_eps?: number;
    qc_custom_eps_inf?: number;
    qc_custom_solvent_name?: string;
  };
  created_at: string;
  updated_at: string;
  started_at?: string;
  finished_at?: string;
  // 用户信息（仅管理员可见）
  username?: string;
  user_email?: string;
  // 计费相关
  cpu_cores?: number;
  estimated_cpu_hours?: number;
  actual_cpu_hours?: number;
  result_locked?: boolean;
  locked_reason?: string;
  billed?: boolean;
}

// 精度等级
export type AccuracyLevel = 'fast' | 'standard' | 'accurate' | 'custom';

// 创建 MD 任务请求
// QC计算选项（用于MD任务附带QC计算）
export interface MDJobQCOptions {
  enabled: boolean;
  molecules?: string[];  // SMILES列表，为空则自动从电解质配方提取
  accuracy_level?: AccuracyLevel;  // 精度等级
  basis_set?: string;
  functional?: string;
  solvent_model?: string;  // 溶剂模型: gas, pcm, smd
  solvent_name?: string;  // 溶剂名称
  use_recommended_params?: boolean;  // 是否使用推荐参数
}

export interface MDJobCreate {
  system_id: number;
  job_name?: string;
  accuracy_level?: AccuracyLevel;  // 精度等级
  nsteps_npt?: number;
  nsteps_nvt?: number;
  timestep?: number;
  temperature?: number;
  pressure?: number;
  freq_trj_npt?: number;
  freq_trj_nvt?: number;
  thermo_freq?: number;
  submit_to_cluster?: boolean;  // 是否提交到集群

  // Slurm 资源配置
  slurm_partition?: string;  // 队列/分区
  slurm_nodes?: number;  // 节点数
  slurm_ntasks?: number;  // 任务数
  slurm_cpus_per_task?: number;  // 每个任务的 CPU 核心数
  slurm_time?: number;  // 最大运行时间（分钟）

  // QC计算选项
  qc_options?: MDJobQCOptions;
}

// 后处理任务
export interface PostprocessJob {
  id: number;
  md_job_id: number;
  status: JobStatus;
  job_type: string;
  created_at: string;
  updated_at: string;
  completed_at?: string;
  error_message?: string;
}

// QC结果简要信息（用于MD任务关联展示）
export interface QCResultSummary {
  energy_au?: number;  // 能量 (A.U.)
  homo_ev?: number;  // HOMO能量 (eV)
  lumo_ev?: number;  // LUMO能量 (eV)
  homo_lumo_gap?: number;  // HOMO-LUMO能隙 (eV)
  esp_min_kcal?: number;  // ESP最小值 (kcal/mol)
  esp_max_kcal?: number;  // ESP最大值 (kcal/mol)
  dipole_moment?: number;  // 偶极矩
  has_esp_image?: boolean;  // 是否有ESP图片
  has_homo_image?: boolean;  // 是否有HOMO图片
  has_lumo_image?: boolean;  // 是否有LUMO图片
}

// QC任务简要信息（用于MD任务关联展示）
export interface QCJobSummary {
  id: number;
  molecule_name: string;
  smiles: string;
  molecule_type: string;
  status: string;
  progress: number;
  basis_set: string;
  functional: string;
  charge?: number;
  spin_multiplicity?: number;
  solvent_model?: string;  // gas, pcm, smd
  solvent_name?: string;  // 隐式溶剂名称
  accuracy_level?: string;  // fast, standard, accurate, custom
  is_reused?: boolean;  // 是否复用已有计算结果
  reused_from_job_id?: number;  // 如果是复用，原始任务ID
  slurm_job_id?: string;
  work_dir?: string;
  created_at: string;
  started_at?: string;
  finished_at?: string;
  error_message?: string;
  result?: QCResultSummary;  // 计算结果（仅已完成任务有）
}

// QC任务状态汇总
export interface QCJobsStatusSummary {
  total: number;
  created: number;
  queued: number;
  running: number;
  postprocessing: number;
  completed: number;
  failed: number;
  cancelled: number;
}

// MD任务关联的QC任务响应
export interface MDJobQCJobsResponse {
  md_job_id: number;
  qc_jobs: QCJobSummary[];
  status_summary: QCJobsStatusSummary;
  qc_enabled: boolean;
}

