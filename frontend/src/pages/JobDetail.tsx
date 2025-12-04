/**
 * 任务详情页面
 */
import { useState, useEffect, useRef, useCallback } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import {
  Card,
  Tabs,
  Button,
  Space,
  message,
  Spin,
  Progress,
  Typography,
  Statistic,
  Row,
  Col,
  Alert,
  Tooltip,
  Modal,
  theme,
} from 'antd';
import {
  ArrowLeftOutlined,
  ReloadOutlined,
  ThunderboltOutlined,
  LineChartOutlined,
  ExperimentOutlined,
  FileTextOutlined,
  RedoOutlined,
  SyncOutlined,
  LockOutlined,
  WalletOutlined,
} from '@ant-design/icons';
import type { MDJob, ElectrolyteSystem } from '../types';
import { JobStatus } from '../types';
import StatusTag from '../components/StatusTag';
import MoleculeViewer from '../components/MoleculeViewer';
import RDFCalculatorNature from '../components/RDFCalculatorNature';
import MSDCalculatorNature from '../components/MSDCalculatorNature';
import SolvationStructureNature from '../components/SolvationStructureNature';
import JobBasicInfo from '../components/JobBasicInfo';
import QCTasksPanel from '../components/QCTasksPanel';
import DesolvationBatchPanel from '../components/DesolvationBatchPanel';
import { getMDJob, resubmitMDJob, getJobSlurmStatus, syncJobStatus, type SlurmJobStatus } from '../api/jobs';
import { getElectrolyte } from '../api/electrolytes';
import { translateError } from '../utils/errorTranslator';
import { useThemeStore } from '../stores/themeStore';
import dayjs from 'dayjs';
import duration from 'dayjs/plugin/duration';

dayjs.extend(duration);

const { Title, Text } = Typography;

export default function JobDetail() {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const [loading, setLoading] = useState(false);
  const [job, setJob] = useState<MDJob | null>(null);
  const [electrolyte, setElectrolyte] = useState<ElectrolyteSystem | null>(null);
  const [slurmStatus, setSlurmStatus] = useState<SlurmJobStatus | null>(null);
  const [lastRefresh, setLastRefresh] = useState<Date>(new Date());
  const pollingRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const [activeTab, setActiveTab] = useState('info');

  // 检查任务是否处于活跃状态（需要轮询）
  const isJobActive = useCallback((jobData: MDJob | null) => {
    if (!jobData) return false;
    return [JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.POSTPROCESSING].includes(jobData.status);
  }, []);

  // 加载任务详情
  const loadJobDetail = useCallback(async () => {
    if (!id) return;

    setLoading(true);
    try {
      console.log(`[JobDetail] Loading job ${id}...`);

      // 加载任务数据
      let jobData = await getMDJob(Number(id));
      console.log(`[JobDetail] Job data loaded (before sync):`, {
        id: jobData.id,
        status: jobData.status,
        progress: jobData.progress,
        started_at: jobData.started_at,
        finished_at: jobData.finished_at,
      });

      // 如果任务处于活跃状态，先同步 Slurm 状态以更新进度
      if ([JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.POSTPROCESSING].includes(jobData.status)) {
        if (jobData.slurm_job_id || jobData.config?.slurm_job_id) {
          try {
            console.log(`[JobDetail] Syncing status for active job...`);
            const syncResult = await syncJobStatus(Number(id));
            console.log(`[JobDetail] Sync result:`, {
              job_id: syncResult.job_id,
              slurm_job_id: syncResult.slurm_job_id,
              slurm_status: syncResult.slurm_status,
              job_status: syncResult.job_status,
              progress: syncResult.progress,
              updated: syncResult.updated
            });

            // 重新加载任务数据以获取更新后的进度
            jobData = await getMDJob(Number(id));
            console.log(`[JobDetail] Job data loaded (after sync):`, {
              id: jobData.id,
              status: jobData.status,
              progress: jobData.progress,
              started_at: jobData.started_at,
              finished_at: jobData.finished_at,
            });
          } catch (e) {
            console.error('同步 Slurm 状态失败:', e);
          }
        }
      }

      setJob(jobData);
      setLastRefresh(new Date());

      // 加载配方数据
      const electrolyteData = await getElectrolyte(jobData.system_id);
      setElectrolyte(electrolyteData);

      // 加载 Slurm 状态详情
      if (jobData.slurm_job_id || jobData.config?.slurm_job_id) {
        try {
          const status = await getJobSlurmStatus(Number(id));
          console.log(`[JobDetail] Slurm status:`, status);
          setSlurmStatus(status);
        } catch (e) {
          console.error('获取 Slurm 状态失败:', e);
        }
      }

    } catch (error: any) {
      console.error(`[JobDetail] Error loading job ${id}:`, error);
      message.error('加载任务详情失败: ' + (error.response?.data?.detail || error.message));
      navigate('/workspace/jobs');
    } finally {
      setLoading(false);
    }
  }, [id, navigate]);

  // 同步 Slurm 状态
  const handleSyncStatus = async () => {
    if (!id) return;
    try {
      const result = await syncJobStatus(Number(id));
      if (result.updated) {
        message.success(`状态已更新: ${result.slurm_status}`);
        await loadJobDetail();
      } else {
        message.info(`状态未变化: ${result.slurm_status}`);
      }
    } catch (error: any) {
      message.error('同步状态失败: ' + (error.response?.data?.detail || error.message));
    }
  };

  useEffect(() => {
    loadJobDetail();
  }, [loadJobDetail]);

  // 智能轮询：只有在任务活跃时才轮询
  useEffect(() => {
    // 清除之前的轮询
    if (pollingRef.current) {
      clearInterval(pollingRef.current);
      pollingRef.current = null;
    }

    // 如果任务处于活跃状态，启动轮询（每 10 秒刷新一次）
    if (isJobActive(job)) {
      console.log(`[JobDetail] Starting polling for job ${job?.id}, status: ${job?.status}`);
      pollingRef.current = setInterval(() => {
        console.log(`[JobDetail] Polling job ${job?.id}...`);
        loadJobDetail();
      }, 10000);
    } else {
      console.log(`[JobDetail] Job ${job?.id} is not active (status: ${job?.status}), polling stopped`);
    }

    // 清理轮询
    return () => {
      if (pollingRef.current) {
        console.log(`[JobDetail] Cleaning up polling for job ${job?.id}`);
        clearInterval(pollingRef.current);
      }
    };
  }, [job?.id, job?.status, isJobActive, loadJobDetail]);

  // 计算运行时间
  const getRunningTime = (jobData: MDJob) => {
    if (jobData.started_at) {
      const end = jobData.finished_at ? dayjs(jobData.finished_at) : dayjs();
      const start = dayjs(jobData.started_at);
      const diff = end.diff(start);
      const dur = dayjs.duration(diff);
      return `${Math.floor(dur.asDays())}天 ${dur.hours()}小时 ${dur.minutes()}分钟`;
    }
    return '-';
  };

  // 重新提交任务
  const handleResubmit = async () => {
    if (!id) return;

    Modal.confirm({
      title: '确认重新提交',
      content: '确定要重新提交这个任务吗？这将重新生成输入文件并提交到集群。',
      okText: '确认',
      cancelText: '取消',
      onOk: async () => {
        try {
          setLoading(true);
          await resubmitMDJob(Number(id));
          message.success('任务已重新提交到集群');
          await loadJobDetail();
        } catch (error: any) {
          message.error('重新提交失败: ' + (error.response?.data?.detail || error.message));
        } finally {
          setLoading(false);
        }
      },
    });
  };

  if (loading) {
    return (
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        height: 'calc(100vh - 64px)',
        background: token.colorBgLayout,
      }}>
        <Spin size="large" />
      </div>
    );
  }

  if (!job && !loading) {
    return (
      <div style={{
        padding: '100px 24px',
        background: token.colorBgLayout,
        minHeight: 'calc(100vh - 64px)',
      }}>
        <Alert message="任务不存在" type="error" style={{ borderRadius: 8 }} />
      </div>
    );
  }

  if (!job || !electrolyte) {
    return null;
  }

  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)', transition: 'background 0.3s' }}>
      {/* 页面头部 - 优化布局 */}
      <Card
        style={{
          marginBottom: 16,
          borderRadius: 12,
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
          border: `1px solid ${token.colorBorder}`,
        }}
      >
        <Row justify="space-between" align="middle">
          <Col>
            <Space size="large">
              <Button
                icon={<ArrowLeftOutlined />}
                onClick={() => navigate('/workspace/jobs')}
                size="large"
              >
                返回
              </Button>
              <div>
                <Title level={3} style={{ margin: 0, marginBottom: 4 }}>
                  <ThunderboltOutlined style={{ color: '#1890ff' }} /> {electrolyte.name || `任务 #${job.id}`}
                </Title>
                <Text type="secondary" style={{ fontSize: 12 }}>
                  最后更新: {lastRefresh.toLocaleTimeString()}
                  {isJobActive(job) && (
                    <span style={{ color: '#52c41a', marginLeft: 8 }}>
                      <SyncOutlined spin /> 自动刷新中
                    </span>
                  )}
                </Text>
              </div>
            </Space>
          </Col>
          <Col>
            <Space>
              <Tooltip title="手动刷新任务详情">
                <Button icon={<ReloadOutlined />} onClick={loadJobDetail}>
                  刷新
                </Button>
              </Tooltip>
              {(job.slurm_job_id || job.config?.slurm_job_id) && (
                <Tooltip title="从 Slurm 同步最新状态">
                  <Button icon={<SyncOutlined />} onClick={handleSyncStatus}>
                    同步状态
                  </Button>
                </Tooltip>
              )}
              {(job.status === JobStatus.FAILED || job.status === JobStatus.CANCELLED) && (
                <Button
                  type="primary"
                  icon={<RedoOutlined />}
                  onClick={handleResubmit}
                  loading={loading}
                >
                  重新提交
                </Button>
              )}
            </Space>
          </Col>
        </Row>
      </Card>

      {/* 任务状态卡片 - 优化布局 */}
      <Card
        style={{
          marginBottom: 16,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
        }}
      >
        <Row gutter={[24, 16]}>
          {/* 任务状态 */}
          <Col xs={12} sm={8} md={4}>
            <div style={{ textAlign: 'center' }}>
              <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8 }}>
                任务状态
              </Text>
              <StatusTag status={job.status} />
            </div>
          </Col>

          {/* 任务进度 */}
          <Col xs={12} sm={8} md={6}>
            <div>
              <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8 }}>
                任务进度
              </Text>
              {(() => {
                // 根据状态计算显示的进度值
                const getProgressByStatus = () => {
                  switch (job.status) {
                    case JobStatus.CREATED:
                      return { percent: 0, status: 'normal' as const, text: '待配置' };
                    case JobStatus.QUEUED:
                      return { percent: 15, status: 'active' as const, text: '排队中' };
                    case JobStatus.RUNNING:
                      return { percent: job.progress || 50, status: 'active' as const, text: '运行中' };
                    case JobStatus.POSTPROCESSING:
                      return { percent: 90, status: 'active' as const, text: '后处理' };
                    case JobStatus.COMPLETED:
                      return { percent: 100, status: 'success' as const, text: '已完成' };
                    case JobStatus.FAILED:
                      return { percent: 100, status: 'exception' as const, text: '失败' };
                    case JobStatus.CANCELLED:
                      // 防止已完成的任务显示"已取消"
                      // 如果进度是 100%，说明任务实际上已完成，不应该显示为取消
                      if (job.progress === 100) {
                        return { percent: 100, status: 'success' as const, text: '已完成' };
                      }
                      return { percent: job.progress || 0, status: 'exception' as const, text: '已取消' };
                    default:
                      return { percent: 0, status: 'normal' as const, text: '' };
                  }
                };
                const progressInfo = getProgressByStatus();
                return (
                  <Progress
                    percent={progressInfo.percent}
                    status={progressInfo.status}
                    strokeColor={
                      progressInfo.status === 'active'
                        ? { '0%': '#108ee9', '100%': '#87d068' }
                        : undefined
                    }
                    size="small"
                    format={() => `${progressInfo.percent}%`}
                  />
                );
              })()}
            </div>
          </Col>

          {/* 运行时间 */}
          <Col xs={12} sm={8} md={4}>
            <div style={{ textAlign: 'center' }}>
              <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8 }}>
                运行时间
              </Text>
              <Text strong style={{ fontSize: 14 }}>
                {getRunningTime(job)}
              </Text>
            </div>
          </Col>

          {/* Slurm Job ID */}
          <Col xs={12} sm={8} md={4}>
            <div style={{ textAlign: 'center' }}>
              <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8 }}>
                Slurm Job ID
              </Text>
              <Text strong style={{ fontSize: 14 }}>
                {job.slurm_job_id || job.config?.slurm_job_id || '-'}
              </Text>
            </div>
          </Col>

          {/* Slurm 状态 */}
          <Col xs={24} sm={16} md={6}>
            <div style={{ textAlign: 'center' }}>
              <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8 }}>
                Slurm 状态
              </Text>
              {slurmStatus ? (
                <div>
                  <Text
                    strong
                    style={{
                      fontSize: 14,
                      color: slurmStatus.status === 'RUNNING' ? '#1890ff' :
                             slurmStatus.status === 'COMPLETED' ? '#52c41a' :
                             slurmStatus.status === 'FAILED' ? '#ff4d4f' : '#666'
                    }}
                  >
                    {slurmStatus.status}
                  </Text>
                  {slurmStatus.elapsed && (
                    <Text type="secondary" style={{ fontSize: 12, marginLeft: 8, whiteSpace: 'nowrap' }}>
                      (已运行: {slurmStatus.elapsed})
                    </Text>
                  )}
                </div>
              ) : (
                <Text strong style={{ fontSize: 14 }}>-</Text>
              )}
            </div>
          </Col>
        </Row>
      </Card>

      {/* 错误信息 - 紧凑显示 */}
      {job.error_message && job.status !== JobStatus.COMPLETED && (() => {
        const translatedError = translateError(job.error_message);
        return translatedError ? (
          <Alert
            message={
              <span style={{ fontSize: 13 }}>
                {translatedError.title}
                <Text type="secondary" style={{ marginLeft: 12, fontSize: 12 }}>
                  {translatedError.description}
                </Text>
              </span>
            }
            description={
              <span style={{ fontSize: 12 }}>
                💡 {translatedError.suggestion}
                {translatedError.originalError && (
                  <Tooltip title={<pre style={{ fontSize: 11, margin: 0, maxWidth: 400, overflow: 'auto' }}>{translatedError.originalError}</pre>}>
                    <Text type="secondary" style={{ marginLeft: 8, cursor: 'pointer', fontSize: 11 }}>
                      查看详情
                    </Text>
                  </Tooltip>
                )}
              </span>
            }
            type={translatedError.severity}
            showIcon
            style={{ marginBottom: 12, borderRadius: 8, padding: '8px 12px' }}
          />
        ) : null;
      })()}

      {/* 结果锁定警告 */}
      {job.result_locked && (
        <Alert
          message={<><LockOutlined /> 结果已锁定</>}
          description={
            <div>
              <p>{job.locked_reason || '由于账户欠费，任务结果已被锁定。'}</p>
              <Button
                type="primary"
                icon={<WalletOutlined />}
                onClick={() => navigate('/workspace/recharge')}
              >
                前往充值
              </Button>
            </div>
          }
          type="warning"
          showIcon
          style={{ marginBottom: 16, borderRadius: 8 }}
        />
      )}

      {/* 详细信息 Tabs - 优化样式 */}
      <Card
        style={{
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
          overflow: 'hidden'
        }}
        styles={{ body: { padding: 0 } }}
      >
        <Tabs
          activeKey={activeTab}
          onChange={setActiveTab}
          size="large"
          tabBarStyle={{
            margin: 0,
            padding: '0 24px',
            background: token.colorBgContainer,
            borderBottom: `2px solid ${token.colorBorder}`
          }}
          style={{ minHeight: '500px' }}
          items={[
            {
              key: 'info',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <FileTextOutlined style={{ marginRight: 6 }} />
                  基本信息
                </span>
              ),
              children: (
                <div style={{ padding: 0 }}>
                  <JobBasicInfo job={job} electrolyte={electrolyte} slurmStatus={slurmStatus} />
                </div>
              ),
            },
            {
              key: 'molecule_structure',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <ExperimentOutlined style={{ marginRight: 6 }} />
                  分子结构
                </span>
              ),
              children: (
                <div style={{ padding: 0 }}>
                  <MoleculeViewer jobId={job.id} />
                </div>
              ),
            },
            {
              key: 'rdf',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <LineChartOutlined style={{ marginRight: 6 }} />
                  径向分布函数 (RDF)
                </span>
              ),
              children: job.status === JobStatus.COMPLETED || job.status === JobStatus.POSTPROCESSING ? (
                <div style={{ padding: 0 }}>
                  <RDFCalculatorNature jobId={job.id} />
                </div>
              ) : (
                <div style={{ padding: 24 }}>
                  <Alert
                    message="任务未完成"
                    description="请等待 MD 任务完成后再进行 RDF 计算"
                    type="info"
                    showIcon
                  />
                </div>
              ),
            },
            {
              key: 'msd',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <LineChartOutlined style={{ marginRight: 6 }} />
                  均方位移 (MSD)
                </span>
              ),
              children: job.status === JobStatus.COMPLETED || job.status === JobStatus.POSTPROCESSING ? (
                <div style={{ padding: 0 }}>
                  <MSDCalculatorNature jobId={job.id} />
                </div>
              ) : (
                <div style={{ padding: 24 }}>
                  <Alert
                    message="任务未完成"
                    description="请等待 MD 任务完成后再进行 MSD 计算"
                    type="info"
                    showIcon
                  />
                </div>
              ),
            },
            {
              key: 'solvation',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <ExperimentOutlined style={{ marginRight: 6 }} />
                  溶剂化结构
                </span>
              ),
              children: job.status === JobStatus.COMPLETED || job.status === JobStatus.POSTPROCESSING ? (
                <div style={{ padding: 0 }}>
                  <SolvationStructureNature
                    jobId={job.id}
                    onGoToDesolvation={() => setActiveTab('desolvation')}
                  />
                </div>
              ) : (
                <div style={{ padding: 24 }}>
                  <Alert
                    message="任务未完成"
                    description="请等待 MD 任务完成后再进行溶剂化结构分析"
                    type="info"
                    showIcon
                  />
                </div>
              ),
            },
            {
              key: 'desolvation',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <ThunderboltOutlined style={{ marginRight: 6, color: '#1890ff' }} />
                  去溶剂化能
                </span>
              ),
              children: job.status === JobStatus.COMPLETED || job.status === JobStatus.POSTPROCESSING ? (
                <div style={{ padding: 24 }}>
                  <DesolvationBatchPanel jobId={job.id} />
                </div>
              ) : (
                <div style={{ padding: 24 }}>
                  <Alert
                    message="任务未完成"
                    description="请等待 MD 任务完成后再进行去溶剂化能计算"
                    type="info"
                    showIcon
                  />
                </div>
              ),
            },
            {
              key: 'qc',
              label: (
                <span style={{ fontSize: 14, fontWeight: 500 }}>
                  <ThunderboltOutlined style={{ marginRight: 6, color: '#722ed1' }} />
                  量子化学计算
                </span>
              ),
              children: (
                <div style={{ padding: 24 }}>
                  <QCTasksPanel mdJobId={job.id} />
                </div>
              ),
            },
          ]}
        />
      </Card>
    </div>
  );
}


