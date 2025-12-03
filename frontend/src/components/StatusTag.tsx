/**
 * 任务状态标签组件
 */
import { Tag } from 'antd';
import { JobStatus } from '../types';
import { getMDStatusConfig } from '../constants/jobStatus';

interface StatusTagProps {
  status: JobStatus;
}

export default function StatusTag({ status }: StatusTagProps) {
  const config = getMDStatusConfig(status);

  return (
    <Tag color={config.color} icon={config.icon}>
      {config.text}
    </Tag>
  );
}

