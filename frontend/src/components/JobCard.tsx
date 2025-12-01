/**
 * 任务卡片组件
 */
import { Card, Space, Button, Popconfirm, Typography, Descriptions, Steps, message, Progress, Tag, Tooltip } from 'antd';
import {
  EyeOutlined,
  StopOutlined,
  RocketOutlined,
  SettingOutlined,
  CopyOutlined,
  RedoOutlined,
  CalendarOutlined,
  ExperimentOutlined,
  DeleteOutlined,
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import type { MDJob, ElectrolyteSystem, MDJobCreate } from '../types';
import { JobStatus } from '../types';
import StatusTag from './StatusTag';
import dayjs from 'dayjs';
import { createMDJob } from '../api/jobs';

const { Text } = Typography;

interface JobCardProps {
  job: MDJob;
  electrolyte?: ElectrolyteSystem;
  onCancel: (id: number) => void;
  onResubmit?: (job: MDJob) => void;
  onDelete?: (id: number) => void;
}

export default function JobCard({ job, electrolyte, onCancel, onResubmit, onDelete }: JobCardProps) {
  const navigate = useNavigate();

  // 判断是否可以取消（只有已提交到集群的任务才能取消）
  const canCancel = job.status === JobStatus.QUEUED || job.status === JobStatus.RUNNING;

  // 判断是否可以重新提交（失败、取消或已完成的任务都可以重新提交）
  const canResubmit = job.status === JobStatus.FAILED ||
                      job.status === JobStatus.CANCELLED ||
                      job.status === JobStatus.COMPLETED;

  // 判断是否可以配置（CREATED 和 CANCELLED 状态可以配置）
  const canConfigure = job.status === JobStatus.CREATED || job.status === JobStatus.CANCELLED;

  // 判断是否可以删除（非运行中和非排队中的任务可以删除）
  const canDelete = job.status !== JobStatus.QUEUED && job.status !== JobStatus.RUNNING;

  // 判断是否启用了QC计算
  const hasQCEnabled = job.config?.qc_enabled === true;

  // 生成任务名称
  const getJobTitle = () => {
    const taskType = hasQCEnabled ? 'MD+QC 任务' : 'MD 任务';
    if (electrolyte) {
      return `${electrolyte.name} - ${taskType}`;
    }
    return `${taskType} #${job.id}`;
  };

  // 处理按钮点击
  const handleConfigClick = () => {
    navigate(`/workspace/jobs/${job.id}/submit`);
  };

  const handleDetailClick = () => {
    navigate(`/workspace/jobs/${job.id}/detail`);
  };

  // 复制任务配置
  const handleCopyConfig = async () => {
    try {
      const newJobData: MDJobCreate = {
        system_id: job.system_id,
        nsteps_npt: job.config?.nsteps_npt || 100000,
        nsteps_nvt: job.config?.nsteps_nvt || 500000,
        timestep: job.config?.timestep || 1.0,
        // 复制资源配置
        slurm_partition: job.config?.slurm_partition || 'cpu',
        slurm_nodes: job.config?.slurm_nodes || 1,
        slurm_ntasks: job.config?.slurm_ntasks || 8,
        slurm_cpus_per_task: job.config?.slurm_cpus_per_task || 8,
        slurm_time: job.config?.slurm_time || 7200,
      };

      const newJob = await createMDJob(newJobData);
      message.success('任务配置已复制，请配置参数后提交');

      // 跳转到新任务的配置页面
      navigate(`/workspace/jobs/${newJob.id}/submit`);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '复制任务失败');
    }
  };

  // 处理卡片点击 - 整个卡片可点击跳转
  const handleCardClick = (e: React.MouseEvent) => {
    // 如果点击的是按钮或其子元素，不触发卡片跳转
    const target = e.target as HTMLElement;
    if (target.closest('button') || target.closest('.ant-popconfirm') || target.closest('.ant-card-actions')) {
      return;
    }
    // 根据状态决定跳转目标
    if (job.status === JobStatus.CREATED) {
      handleConfigClick();
    } else {
      handleDetailClick();
    }
  };

  // 获取状态对应的渐变色
  const getStatusGradient = () => {
    switch (job.status) {
      case JobStatus.RUNNING:
      case JobStatus.QUEUED:
      case JobStatus.POSTPROCESSING:
        return 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)';
      case JobStatus.COMPLETED:
        return 'linear-gradient(135deg, #11998e 0%, #38ef7d 100%)';
      case JobStatus.FAILED:
      case JobStatus.CANCELLED:
        return 'linear-gradient(135deg, #ff416c 0%, #ff4b2b 100%)';
      default:
        return 'linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)';
    }
  };

  return (
    <Card
      hoverable
      onClick={handleCardClick}
      style={{
        height: '100%',
        cursor: 'pointer',
        transition: 'all 0.3s ease',
        border: 'none',
        borderRadius: 12,
        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
      }}
      styles={{
        body: { padding: '20px' },
      }}
      actions={[
        canConfigure && (
          <Button
            key="config"
            type="link"
            icon={<SettingOutlined />}
            onClick={(e) => { e.stopPropagation(); handleConfigClick(); }}
          >
            {job.status === JobStatus.CANCELLED ? '重新配置' : '配置参数'}
          </Button>
        ),
        job.status !== JobStatus.CREATED && (
          <Button
            key="detail"
            type="link"
            icon={<EyeOutlined />}
            onClick={(e) => { e.stopPropagation(); handleDetailClick(); }}
          >
            查看详情
          </Button>
        ),
        <Button
          key="copy"
          type="link"
          icon={<CopyOutlined />}
          onClick={(e) => { e.stopPropagation(); handleCopyConfig(); }}
        >
          复制配置
        </Button>,
        canCancel ? (
          <Popconfirm
            key="cancel"
            title="确定要取消这个任务吗？"
            description="取消后任务将停止运行。"
            onConfirm={() => onCancel(job.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button type="link" danger icon={<StopOutlined />} onClick={(e) => e.stopPropagation()}>
              取消任务
            </Button>
          </Popconfirm>
        ) : null,
        canResubmit && onResubmit ? (
          <Button
            key="resubmit"
            type="link"
            icon={<RedoOutlined />}
            onClick={(e) => { e.stopPropagation(); onResubmit(job); }}
          >
            重新提交
          </Button>
        ) : null,
        canDelete && onDelete ? (
          <Popconfirm
            key="delete"
            title="确定要删除这个任务吗？"
            description="删除后不可恢复。"
            onConfirm={() => onDelete(job.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button type="link" danger icon={<DeleteOutlined />} onClick={(e) => e.stopPropagation()}>
              删除
            </Button>
          </Popconfirm>
        ) : null,
      ].filter(Boolean)}
    >
      <div style={{ display: 'flex', gap: 16 }}>
        {/* 左侧图标 */}
        <div style={{
          width: 48,
          height: 48,
          borderRadius: 12,
          background: getStatusGradient(),
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          flexShrink: 0,
        }}>
          <RocketOutlined style={{ fontSize: 24, color: '#fff' }} />
        </div>

        {/* 右侧内容 */}
        <div style={{ flex: 1, minWidth: 0 }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: 8 }}>
            <Text strong style={{ fontSize: 15 }}>
              {getJobTitle()}
            </Text>
            <StatusTag status={job.status} />
          </div>

          <div style={{ marginBottom: 8, display: 'flex', alignItems: 'center', gap: 8 }}>
            <Text type="secondary" style={{ fontSize: 12 }}>任务 ID: #{job.id}</Text>
            {job.slurm_job_id && (
              <Text type="secondary" style={{ fontSize: 12 }}>
                Slurm: {job.slurm_job_id}
              </Text>
            )}
            {job.config?.qc_enabled && (
              <Tooltip title={`QC计算: ${job.config.qc_functional || 'B3LYP'}/${job.config.qc_basis_set || '6-31++G(d,p)'}`}>
                <Tag
                  color="purple"
                  style={{ fontSize: 11, padding: '0 4px', lineHeight: '16px', marginLeft: 4 }}
                  icon={<ExperimentOutlined />}
                >
                  QC
                </Tag>
              </Tooltip>
            )}
          </div>

          {job.config && (
            <div style={{ marginBottom: 8 }}>
              <Text type="secondary" style={{ fontSize: 12 }}>
                NPT: {job.config.nsteps_npt?.toLocaleString()} | NVT: {job.config.nsteps_nvt?.toLocaleString()}
              </Text>
              {job.config.slurm_ntasks && job.config.slurm_cpus_per_task && (
                <Text type="secondary" style={{ fontSize: 12, display: 'block' }}>
                  资源: {job.config.slurm_ntasks * job.config.slurm_cpus_per_task} 核
                  {job.config.slurm_partition && ` (${job.config.slurm_partition})`}
                </Text>
              )}
            </div>
          )}

          {/* 进度条 - 仅运行中显示 */}
          {(job.status === JobStatus.RUNNING || job.status === JobStatus.POSTPROCESSING) && (
            <Progress
              percent={job.progress || 0}
              size="small"
              status="active"
              strokeColor={{
                '0%': '#108ee9',
                '100%': '#87d068',
              }}
              style={{ marginBottom: 8 }}
            />
          )}

          <Space size={4}>
            <CalendarOutlined style={{ color: '#bfbfbf', fontSize: 12 }} />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {dayjs(job.created_at).format('YYYY-MM-DD HH:mm')}
            </Text>
          </Space>
        </div>
      </div>

      {job.error_message && (
        <div style={{
          marginTop: 12,
          padding: 8,
          background: '#fff2f0',
          borderRadius: 8,
          border: '1px solid #ffccc7'
        }}>
          <Text type="danger" style={{ fontSize: 12 }}>
            错误: {job.error_message}
          </Text>
        </div>
      )}
    </Card>
  );
}

