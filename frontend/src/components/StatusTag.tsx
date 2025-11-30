/**
 * 任务状态标签组件
 */
import { Tag } from 'antd';
import {
  ClockCircleOutlined,
  SyncOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  StopOutlined,
} from '@ant-design/icons';
import { JobStatus } from '../types';

interface StatusTagProps {
  status: JobStatus;
}

export default function StatusTag({ status }: StatusTagProps) {
  const statusConfig = {
    [JobStatus.CREATED]: {
      color: 'default',
      icon: <ClockCircleOutlined />,
      text: '已创建',
    },
    [JobStatus.QUEUED]: {
      color: 'blue',
      icon: <ClockCircleOutlined />,
      text: '已提交',
    },
    [JobStatus.RUNNING]: {
      color: 'processing',
      icon: <SyncOutlined spin />,
      text: '运行中',
    },
    [JobStatus.POSTPROCESSING]: {
      color: 'cyan',
      icon: <SyncOutlined spin />,
      text: '后处理中',
    },
    [JobStatus.COMPLETED]: {
      color: 'success',
      icon: <CheckCircleOutlined />,
      text: '已完成',
    },
    [JobStatus.FAILED]: {
      color: 'error',
      icon: <CloseCircleOutlined />,
      text: '失败',
    },
    [JobStatus.CANCELLED]: {
      color: 'warning',
      icon: <StopOutlined />,
      text: '已取消',
    },
  };

  const config = statusConfig[status] || statusConfig[JobStatus.CREATED];

  return (
    <Tag color={config.color} icon={config.icon}>
      {config.text}
    </Tag>
  );
}

