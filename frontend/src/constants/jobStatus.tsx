/**
 * 统一的任务状态配置
 * 用于 MD 和 QC 任务的状态显示
 */
import React from 'react';
import {
  ClockCircleOutlined,
  SyncOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  StopOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import { JobStatus } from '../types/api';

// QC 任务状态与 MD 任务状态相同
export type QCJobStatus = JobStatus;

// MD 任务状态配置
export const MD_STATUS_CONFIG: Record<JobStatus, {
  color: string;
  text: string;
  icon: React.ReactNode;
}> = {
  [JobStatus.CREATED]: {
    color: 'default',
    text: '已创建',
    icon: <ClockCircleOutlined />,
  },
  [JobStatus.QUEUED]: {
    color: 'blue',
    text: '排队中',
    icon: <ClockCircleOutlined />,
  },
  [JobStatus.RUNNING]: {
    color: 'processing',
    text: '运行中',
    icon: <SyncOutlined spin />,
  },
  [JobStatus.POSTPROCESSING]: {
    color: 'cyan',
    text: '后处理中',
    icon: <SyncOutlined spin />,
  },
  [JobStatus.COMPLETED]: {
    color: 'success',
    text: '已完成',
    icon: <CheckCircleOutlined />,
  },
  [JobStatus.FAILED]: {
    color: 'error',
    text: '失败',
    icon: <CloseCircleOutlined />,
  },
  [JobStatus.CANCELLED]: {
    color: 'warning',
    text: '已取消',
    icon: <StopOutlined />,
  },
};

// QC 任务状态配置
export const QC_STATUS_CONFIG: Record<string, {
  color: string;
  text: string;
  icon?: React.ReactNode;
}> = {
  CREATED: {
    color: 'default',
    text: '已创建',
    icon: <ClockCircleOutlined />,
  },
  SUBMITTED: {
    color: 'blue',
    text: '已提交',
    icon: <ClockCircleOutlined />,
  },
  QUEUED: {
    color: 'processing',
    text: '排队中',
    icon: <ClockCircleOutlined />,
  },
  RUNNING: {
    color: 'processing',
    text: '运行中',
    icon: <SyncOutlined spin />,
  },
  POSTPROCESSING: {
    color: 'processing',
    text: '后处理中',
    icon: <SyncOutlined spin />,
  },
  COMPLETED: {
    color: 'success',
    text: '已完成',
    icon: <CheckCircleOutlined />,
  },
  FAILED: {
    color: 'error',
    text: '失败',
    icon: <CloseCircleOutlined />,
  },
  CANCELLED: {
    color: 'warning',
    text: '已取消',
    icon: <StopOutlined />,
  },
};

// 获取 MD 状态配置
export const getMDStatusConfig = (status: JobStatus) => {
  return MD_STATUS_CONFIG[status] || {
    color: 'default',
    text: status,
    icon: <ExperimentOutlined />,
  };
};

// 获取 QC 状态配置
export const getQCStatusConfig = (status: string) => {
  return QC_STATUS_CONFIG[status] || {
    color: 'default',
    text: status,
    icon: <ExperimentOutlined />,
  };
};

// 状态颜色映射（用于进度条等）
export const getStatusColor = (status: string): string => {
  switch (status) {
    case 'COMPLETED':
      return 'success';
    case 'RUNNING':
    case 'POSTPROCESSING':
      return 'processing';
    case 'QUEUED':
    case 'SUBMITTED':
      return 'warning';
    case 'FAILED':
    case 'CANCELLED':
      return 'error';
    default:
      return 'default';
  }
};

