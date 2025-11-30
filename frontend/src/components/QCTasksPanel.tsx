/**
 * QC任务状态面板组件
 * 用于在MD任务详情页展示关联的QC任务状态
 */
import React, { useState, useEffect, useCallback } from 'react';
import {
  Card,
  Row,
  Col,
  Tag,
  Progress,
  Typography,
  Space,
  Button,
  Tooltip,
  Spin,
  Empty,
  Badge,
  Statistic,
  List,
  message,
} from 'antd';
import {
  ExperimentOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  SyncOutlined,
  ClockCircleOutlined,
  ThunderboltOutlined,
  EyeOutlined,
  ReloadOutlined,
  RightOutlined,
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import { getMDJobQCJobs } from '../api/jobs';
import type { QCJobSummary, QCJobsStatusSummary } from '../types';

const { Text, Title } = Typography;

interface QCTasksPanelProps {
  mdJobId: number;
  refreshInterval?: number; // 自动刷新间隔（毫秒），默认10秒
}

// 状态颜色映射
const getStatusColor = (status: string): string => {
  switch (status) {
    case 'COMPLETED':
      return 'success';
    case 'RUNNING':
    case 'POSTPROCESSING':
      return 'processing';
    case 'QUEUED':
      return 'warning';
    case 'FAILED':
    case 'CANCELLED':
      return 'error';
    default:
      return 'default';
  }
};

// 状态文字映射
const getStatusText = (status: string): string => {
  switch (status) {
    case 'CREATED':
      return '已创建';
    case 'QUEUED':
      return '排队中';
    case 'RUNNING':
      return '计算中';
    case 'POSTPROCESSING':
      return '后处理';
    case 'COMPLETED':
      return '已完成';
    case 'FAILED':
      return '失败';
    case 'CANCELLED':
      return '已取消';
    default:
      return status;
  }
};

// 状态图标映射
const getStatusIcon = (status: string) => {
  switch (status) {
    case 'COMPLETED':
      return <CheckCircleOutlined style={{ color: '#52c41a' }} />;
    case 'RUNNING':
    case 'POSTPROCESSING':
      return <SyncOutlined spin style={{ color: '#1890ff' }} />;
    case 'QUEUED':
      return <ClockCircleOutlined style={{ color: '#faad14' }} />;
    case 'FAILED':
    case 'CANCELLED':
      return <CloseCircleOutlined style={{ color: '#ff4d4f' }} />;
    default:
      return <ExperimentOutlined style={{ color: '#d9d9d9' }} />;
  }
};

// 分子类型标签
const getMoleculeTypeTag = (type: string) => {
  switch (type) {
    case 'cation':
      return <Tag color="red">阳离子</Tag>;
    case 'anion':
      return <Tag color="blue">阴离子</Tag>;
    case 'solvent':
      return <Tag color="green">溶剂</Tag>;
    default:
      return <Tag>自定义</Tag>;
  }
};

// 溶剂模型显示
const getSolventModelText = (model?: string): string => {
  switch (model) {
    case 'gas':
      return '气相';
    case 'pcm':
    case 'PCM':
      return 'PCM';
    case 'smd':
    case 'SMD':
      return 'SMD';
    default:
      return model || '气相';
  }
};

// 精度等级显示
const getAccuracyLevelText = (level?: string): string => {
  switch (level) {
    case 'fast':
      return '快速';
    case 'standard':
      return '标准';
    case 'accurate':
      return '精确';
    case 'custom':
      return '自定义';
    default:
      return level || '标准';
  }
};

export default function QCTasksPanel({ mdJobId, refreshInterval = 10000 }: QCTasksPanelProps) {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(true);
  const [qcJobs, setQcJobs] = useState<QCJobSummary[]>([]);
  const [statusSummary, setStatusSummary] = useState<QCJobsStatusSummary | null>(null);
  const [qcEnabled, setQcEnabled] = useState(false);

  // 加载QC任务数据
  const loadQCJobs = useCallback(async () => {
    try {
      const data = await getMDJobQCJobs(mdJobId);
      setQcJobs(data.qc_jobs);
      setStatusSummary(data.status_summary);
      setQcEnabled(data.qc_enabled);
    } catch (error: any) {
      console.error('加载QC任务失败:', error);
      // 如果是404，说明没有关联的QC任务
      if (error.response?.status !== 404) {
        message.error('加载QC任务状态失败');
      }
    } finally {
      setLoading(false);
    }
  }, [mdJobId]);

  // 初始加载和自动刷新
  useEffect(() => {
    loadQCJobs();

    // 只有在有运行中的任务时才自动刷新
    const interval = setInterval(() => {
      if (statusSummary && (statusSummary.queued > 0 || statusSummary.running > 0 || statusSummary.postprocessing > 0)) {
        loadQCJobs();
      }
    }, refreshInterval);

    return () => clearInterval(interval);
  }, [loadQCJobs, refreshInterval, statusSummary]);

  // 手动刷新
  const handleRefresh = () => {
    setLoading(true);
    loadQCJobs();
  };

  // 如果没有启用QC计算，不显示面板
  if (!loading && !qcEnabled) {
    return null;
  }

  // 计算总体进度
  const calculateOverallProgress = () => {
    if (!statusSummary || statusSummary.total === 0) return 0;
    const completed = statusSummary.completed + statusSummary.failed + statusSummary.cancelled;
    return Math.round((completed / statusSummary.total) * 100);
  };

  // 判断是否全部完成
  const isAllCompleted = statusSummary && 
    (statusSummary.completed + statusSummary.failed + statusSummary.cancelled) === statusSummary.total;

  // 判断是否有运行中的任务
  const hasRunning = statusSummary && 
    (statusSummary.queued > 0 || statusSummary.running > 0 || statusSummary.postprocessing > 0);

  return (
    <Card
      title={
        <Space>
          <ExperimentOutlined style={{ color: '#722ed1' }} />
          <span>量子化学计算</span>
          {hasRunning && (
            <Badge status="processing" text="计算中" />
          )}
          {isAllCompleted && statusSummary && statusSummary.completed > 0 && (
            <Badge status="success" text="已完成" />
          )}
        </Space>
      }
      extra={
        <Button
          type="text"
          icon={<ReloadOutlined spin={loading} />}
          onClick={handleRefresh}
          disabled={loading}
        >
          刷新
        </Button>
      }
      style={{ marginBottom: 16 }}
    >
      {loading ? (
        <div style={{ textAlign: 'center', padding: 40 }}>
          <Spin tip="加载QC任务状态..." />
        </div>
      ) : qcJobs.length === 0 ? (
        <Empty description="暂无QC计算任务" />
      ) : (
        <>
          {/* 状态统计卡片 */}
          <Row gutter={16} style={{ marginBottom: 16 }}>
            <Col span={6}>
              <Statistic
                title="总任务数"
                value={statusSummary?.total || 0}
                prefix={<ExperimentOutlined />}
              />
            </Col>
            <Col span={6}>
              <Statistic
                title="已完成"
                value={statusSummary?.completed || 0}
                valueStyle={{ color: '#52c41a' }}
                prefix={<CheckCircleOutlined />}
              />
            </Col>
            <Col span={6}>
              <Statistic
                title="运行中"
                value={(statusSummary?.queued || 0) + (statusSummary?.running || 0) + (statusSummary?.postprocessing || 0)}
                valueStyle={{ color: '#1890ff' }}
                prefix={<SyncOutlined spin={!!hasRunning} />}
              />
            </Col>
            <Col span={6}>
              <Statistic
                title="失败"
                value={statusSummary?.failed || 0}
                valueStyle={{ color: statusSummary?.failed ? '#ff4d4f' : undefined }}
                prefix={<CloseCircleOutlined />}
              />
            </Col>
          </Row>

          {/* 总体进度条 */}
          {statusSummary && statusSummary.total > 0 && (
            <Progress
              percent={calculateOverallProgress()}
              status={isAllCompleted ? (statusSummary.failed > 0 ? 'exception' : 'success') : 'active'}
              strokeColor={isAllCompleted && statusSummary.failed === 0 ? '#52c41a' : undefined}
              style={{ marginBottom: 16 }}
            />
          )}

          {/* QC任务列表 */}
          <List
            size="small"
            dataSource={qcJobs}
            renderItem={(job) => (
              <List.Item
                key={job.id}
                style={{
                  padding: '12px 16px',
                  marginBottom: 8,
                  background: job.is_reused ? '#f6ffed' : '#fafafa',
                  borderRadius: 8,
                  border: job.is_reused ? '1px solid #b7eb8f' : '1px solid #f0f0f0'
                }}
              >
                <div style={{ width: '100%' }}>
                  {/* 标题行 */}
                  <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
                    <Space>
                      {getStatusIcon(job.status)}
                      <Text strong style={{ fontSize: 14 }}>{job.molecule_name}</Text>
                      {getMoleculeTypeTag(job.molecule_type)}
                      <Tag color={getStatusColor(job.status)}>
                        {getStatusText(job.status)}
                      </Tag>
                      {job.is_reused && (
                        <Tooltip title={`复用已有计算结果 (来自任务 #${job.reused_from_job_id})`}>
                          <Tag color="success" icon={<CheckCircleOutlined />}>
                            复用结果
                          </Tag>
                        </Tooltip>
                      )}
                    </Space>
                    <Tooltip title="查看详情">
                      <Button
                        type="link"
                        size="small"
                        icon={<EyeOutlined />}
                        onClick={() => navigate(`/workspace/qc-jobs/${job.id}`, { state: { fromMDJob: mdJobId } })}
                      >
                        详情
                      </Button>
                    </Tooltip>
                  </div>

                  {/* 计算参数行 - 第一行：基本参数 */}
                  <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>
                    <Space split={<span style={{ color: '#d9d9d9' }}>|</span>} wrap>
                      <span>泛函: <Text code style={{ fontSize: 11 }}>{job.functional}</Text></span>
                      <span>基组: <Text code style={{ fontSize: 11 }}>{job.basis_set}</Text></span>
                      <span>电荷: <Text code style={{ fontSize: 11 }}>{job.charge ?? 0}</Text></span>
                      <span>自旋多重度: <Text code style={{ fontSize: 11 }}>{job.spin_multiplicity ?? 1}</Text></span>
                    </Space>
                  </div>

                  {/* 计算参数行 - 第二行：溶剂模型和精度 */}
                  <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>
                    <Space split={<span style={{ color: '#d9d9d9' }}>|</span>} wrap>
                      <span>溶剂模型: <Text code style={{ fontSize: 11 }}>{getSolventModelText(job.solvent_model)}</Text></span>
                      {job.solvent_model && job.solvent_model !== 'gas' && job.solvent_name && (
                        <span>隐式溶剂: <Text code style={{ fontSize: 11 }}>{job.solvent_name}</Text></span>
                      )}
                      {job.accuracy_level && (
                        <span>精度: <Text code style={{ fontSize: 11 }}>{getAccuracyLevelText(job.accuracy_level)}</Text></span>
                      )}
                    </Space>
                  </div>

                  {/* SMILES 行 */}
                  {job.smiles && (
                    <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>
                      <Tooltip title={job.smiles}>
                        <span>SMILES: <Text type="secondary" style={{ fontSize: 11 }}>
                          {job.smiles.length > 50 ? job.smiles.substring(0, 50) + '...' : job.smiles}
                        </Text></span>
                      </Tooltip>
                    </div>
                  )}

                  {/* Slurm 信息行 */}
                  <div style={{ fontSize: 12, color: '#888' }}>
                    <Space split={<span style={{ color: '#d9d9d9' }}>|</span>} wrap>
                      {job.slurm_job_id && (
                        <span>Slurm ID: <Text code style={{ fontSize: 11, color: '#1890ff' }}>{job.slurm_job_id}</Text></span>
                      )}
                      {(job.status === 'RUNNING' || job.status === 'POSTPROCESSING') && (
                        <span>进度: <Text style={{ fontSize: 11, color: '#52c41a' }}>{job.progress}%</Text></span>
                      )}
                      {job.started_at && (
                        <span>开始: <Text type="secondary" style={{ fontSize: 11 }}>
                          {new Date(job.started_at).toLocaleString()}
                        </Text></span>
                      )}
                      {job.finished_at && (
                        <span>完成: <Text type="secondary" style={{ fontSize: 11 }}>
                          {new Date(job.finished_at).toLocaleString()}
                        </Text></span>
                      )}
                    </Space>
                  </div>

                  {/* 错误信息 */}
                  {job.error_message && (
                    <div style={{ marginTop: 6, padding: '4px 8px', background: '#fff2f0', borderRadius: 4, border: '1px solid #ffccc7' }}>
                      <Text type="danger" style={{ fontSize: 11 }}>
                        ⚠️ {job.error_message}
                      </Text>
                    </div>
                  )}
                </div>
              </List.Item>
            )}
          />

          {/* 查看全部链接 */}
          {qcJobs.length > 0 && (
            <div style={{ textAlign: 'center', marginTop: 16 }}>
              <Button
                type="link"
                onClick={() => navigate('/workspace/qc-jobs')}
                icon={<RightOutlined />}
              >
                查看全部QC任务
              </Button>
            </div>
          )}
        </>
      )}
    </Card>
  );
}

