/**
 * QC任务详情页面
 */
import { useState, useEffect, useCallback } from 'react';
import { useParams, useNavigate, useLocation } from 'react-router-dom';
import {
  Card,
  Descriptions,
  Tag,
  Button,
  Space,
  Spin,
  message,
  Row,
  Col,
  Typography,
  Divider,
  Statistic,
  Progress,
  Alert,
  Popconfirm,
  Timeline,
  Tabs,
  theme,
} from 'antd';
import {
  ArrowLeftOutlined,
  ReloadOutlined,
  PlayCircleOutlined,
  DeleteOutlined,
  ExperimentOutlined,
  ThunderboltOutlined,
  ClockCircleOutlined,
  FieldTimeOutlined,
  CheckCircleOutlined,
  SyncOutlined,
  CloseCircleOutlined,
  HourglassOutlined,
  SettingOutlined,
  DatabaseOutlined,
} from '@ant-design/icons';
import { getQCJob, getQCJobStatus, submitQCJob, deleteQCJob, getQCResults } from '../api/qc';
import type { QCJob, QCResult } from '../types/qc';
import QCResultsPanel from '../components/QCResultsPanel';
import { useThemeStore } from '../stores/themeStore';
import dayjs from 'dayjs';
import duration from 'dayjs/plugin/duration';

dayjs.extend(duration);

const { Title, Text, Paragraph } = Typography;

// Dashboard 样式常量（与其他页面保持一致）
const DASHBOARD_STYLES = {
  cardBorderRadius: 12,
  cardPadding: 24,
  gutter: 24,
  titleFontSize: 16,
  titleFontWeight: 600,
};

// 响应式CSS样式
const RESPONSIVE_STYLES = `
  .qc-stats-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 16px;
    margin-bottom: 16px;
  }
  @media (max-width: 1200px) {
    .qc-stats-grid {
      grid-template-columns: repeat(2, 1fr);
    }
  }
  @media (max-width: 768px) {
    .qc-stats-grid {
      grid-template-columns: 1fr;
    }
  }
  .dashboard-card {
    transition: all 0.3s ease;
  }
  .dashboard-card:hover {
    box-shadow: 0 8px 24px rgba(15, 23, 42, 0.12);
    transform: translateY(-2px);
  }
`;

// 状态映射
const statusMap: Record<string, { color: string; text: string }> = {
  CREATED: { color: 'default', text: '已创建' },
  QUEUED: { color: 'processing', text: '排队中' },
  RUNNING: { color: 'processing', text: '运行中' },
  POSTPROCESSING: { color: 'processing', text: '后处理中' },
  COMPLETED: { color: 'success', text: '已完成' },
  FAILED: { color: 'error', text: '失败' },
  CANCELLED: { color: 'warning', text: '已取消' },
};

export default function QCJobDetail() {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const location = useLocation();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [job, setJob] = useState<QCJob | null>(null);
  const [results, setResults] = useState<QCResult[]>([]);
  const [loading, setLoading] = useState(true);
  const [submitting, setSubmitting] = useState(false);

  // 获取来源页面信息（从MD任务详情跳转过来时携带）
  const fromMDJob = location.state?.fromMDJob as number | undefined;

  // 返回按钮处理：如果从MD任务详情跳转过来，则返回MD任务详情
  const handleGoBack = () => {
    if (fromMDJob) {
      navigate(`/workspace/jobs/${fromMDJob}/detail`);
    } else {
      navigate('/workspace/qc-jobs');
    }
  };

  // 加载任务详情
  const loadJob = useCallback(async () => {
    if (!id) return;
    try {
      const data = await getQCJob(parseInt(id));
      setJob(data);
      if (data.results) {
        setResults(data.results);
      }
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载任务详情失败');
    } finally {
      setLoading(false);
    }
  }, [id]);

  // 轮询状态
  useEffect(() => {
    loadJob();
    
    const interval = setInterval(async () => {
      if (!id || !job) return;
      if (['QUEUED', 'RUNNING', 'POSTPROCESSING'].includes(job.status)) {
        try {
          const status = await getQCJobStatus(parseInt(id));
          setJob(prev => prev ? { ...prev, ...status } : null);
          if (status.status === 'COMPLETED') {
            loadJob(); // 完成后重新加载完整数据
          }
        } catch (error) {
          console.error('轮询状态失败:', error);
        }
      }
    }, 5000);

    return () => clearInterval(interval);
  }, [id, job?.status, loadJob]);

  // 提交任务
  const handleSubmit = async () => {
    if (!id) return;
    setSubmitting(true);
    try {
      await submitQCJob(parseInt(id));
      message.success('任务已提交到计算集群');
      loadJob();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '提交失败');
    } finally {
      setSubmitting(false);
    }
  };

  // 删除任务
  const handleDelete = async () => {
    if (!id) return;
    try {
      await deleteQCJob(parseInt(id));
      message.success('任务已删除');
      handleGoBack();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  if (loading) {
    return (
      <div style={{ textAlign: 'center', padding: 100 }}>
        <Spin size="large" />
      </div>
    );
  }

  if (!job) {
    return (
      <div style={{ padding: 24 }}>
        <Alert type="error" message="任务不存在" showIcon />
        <Button style={{ marginTop: 16 }} onClick={handleGoBack}>
          返回
        </Button>
      </div>
    );
  }

  const { color, text } = statusMap[job.status] || { color: 'default', text: job.status };

  // 计算运行时间
  const getRunningTime = () => {
    if (job.started_at) {
      const end = job.finished_at ? dayjs(job.finished_at) : dayjs();
      const start = dayjs(job.started_at);
      const diff = end.diff(start);
      const dur = dayjs.duration(diff);

      const days = Math.floor(dur.asDays());
      const hours = dur.hours();
      const minutes = dur.minutes();
      const seconds = dur.seconds();

      if (days > 0) {
        return `${days}天 ${hours}小时 ${minutes}分钟`;
      } else if (hours > 0) {
        return `${hours}小时 ${minutes}分钟 ${seconds}秒`;
      } else if (minutes > 0) {
        return `${minutes}分钟 ${seconds}秒`;
      } else {
        return `${seconds}秒`;
      }
    }
    return '-';
  };

  // 计算CPU核时（core-hours）
  const getCoreHours = () => {
    if (!job.started_at) return '-';

    const cpus = job.config?.slurm_cpus || 16;  // 默认16核

    const end = job.finished_at ? dayjs(job.finished_at) : dayjs();
    const start = dayjs(job.started_at);
    const diff = end.diff(start);
    const hours = diff / (1000 * 60 * 60);  // 转换为小时

    const coreHours = cpus * hours;
    return coreHours.toFixed(2);
  };

  // 获取状态时间线
  const getStatusTimeline = () => {
    const items = [];

    items.push({
      color: 'green',
      dot: <CheckCircleOutlined />,
      children: (
        <div>
          <Text strong>任务创建</Text>
          <br />
          <Text type="secondary" style={{ fontSize: 12 }}>
            {dayjs(job.created_at).format('YYYY-MM-DD HH:mm:ss')}
          </Text>
        </div>
      ),
    });

    if (job.started_at) {
      items.push({
        color: job.status === 'RUNNING' ? 'blue' : 'green',
        dot: job.status === 'RUNNING' ? <SyncOutlined spin /> : <CheckCircleOutlined />,
        children: (
          <div>
            <Text strong>开始计算</Text>
            <br />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {dayjs(job.started_at).format('YYYY-MM-DD HH:mm:ss')}
            </Text>
          </div>
        ),
      });
    }

    if (job.finished_at) {
      items.push({
        color: job.status === 'COMPLETED' ? 'green' : 'red',
        dot: job.status === 'COMPLETED' ? <CheckCircleOutlined /> : <CloseCircleOutlined />,
        children: (
          <div>
            <Text strong>{job.status === 'COMPLETED' ? '计算完成' : '计算失败'}</Text>
            <br />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {dayjs(job.finished_at).format('YYYY-MM-DD HH:mm:ss')}
            </Text>
          </div>
        ),
      });
    }

    return items;
  };

  // 卡片样式
  const dashboardCardStyle: React.CSSProperties = {
    background: token.colorBgContainer,
    borderRadius: DASHBOARD_STYLES.cardBorderRadius,
    boxShadow: mode === 'dark' ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(15, 23, 42, 0.08)',
    border: `1px solid ${token.colorBorder}`,
    transition: 'all 0.3s ease',
  };

  return (
    <div style={{ padding: '16px 24px', background: token.colorBgLayout, minHeight: '100vh', transition: 'background 0.3s' }}>
      <style>{RESPONSIVE_STYLES}</style>
      {/* 顶部导航栏 */}
      <div style={{
        marginBottom: 16,
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
        background: token.colorBgContainer,
        padding: '12px 20px',
        borderRadius: 8,
        boxShadow: mode === 'dark' ? '0 1px 2px rgba(0,0,0,0.2)' : '0 1px 2px rgba(0,0,0,0.03)'
      }}>
        <Space size="middle">
          <Button icon={<ArrowLeftOutlined />} onClick={handleGoBack}>
            {fromMDJob ? '返回MD任务详情' : '返回列表'}
          </Button>
          <Divider type="vertical" style={{ height: 24 }} />
          <ExperimentOutlined style={{ fontSize: 22, color: '#1890ff' }} />
          <Title level={4} style={{ margin: 0 }}>{job.molecule_name}</Title>
          <Tag color={color} style={{ fontSize: 13, padding: '2px 10px' }}>{text}</Tag>
          {job.slurm_job_id && (
            <Tag color="default">Slurm #{job.slurm_job_id}</Tag>
          )}
        </Space>
        <Space>
          <Button icon={<ReloadOutlined />} onClick={loadJob}>刷新</Button>
          {job.status === 'CREATED' && (
            <Button
              type="primary"
              icon={<PlayCircleOutlined />}
              loading={submitting}
              onClick={handleSubmit}
            >
              提交计算
            </Button>
          )}
          <Popconfirm
            title="确定要删除这个任务吗？"
            description="删除后将无法恢复"
            onConfirm={handleDelete}
            okText="确定"
            cancelText="取消"
          >
            <Button danger icon={<DeleteOutlined />}>删除</Button>
          </Popconfirm>
        </Space>
      </div>

      {/* 进度条 - 仅运行中显示 */}
      {['QUEUED', 'RUNNING', 'POSTPROCESSING'].includes(job.status) && (
        <Card size="small" style={{ marginBottom: 16 }}>
          <Row align="middle" gutter={16}>
            <Col flex="100px">
              <Space>
                <SyncOutlined spin style={{ color: '#1890ff' }} />
                <Text>
                  {job.status === 'QUEUED' ? '排队中' :
                   job.status === 'RUNNING' ? '计算中' : '后处理中'}
                </Text>
              </Space>
            </Col>
            <Col flex="auto">
              <Progress
                percent={job.progress || 0}
                status="active"
                strokeColor={{ from: '#108ee9', to: '#87d068' }}
              />
            </Col>
          </Row>
        </Card>
      )}

      {/* 错误信息 */}
      {job.error_message && (
        <Alert
          type="error"
          message="计算失败"
          description={job.error_message}
          style={{ marginBottom: 16 }}
          showIcon
        />
      )}

      {/* 统计卡片网格 */}
      <div className="qc-stats-grid">
        {/* 运行时长 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #1890ff 0%, #096dd9 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <ClockCircleOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>运行时长 (Runtime)</div>
                <div style={{ fontSize: 22, fontWeight: 700, color: '#1890ff' }}>
                  {getRunningTime()}
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* CPU核时 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #fa8c16 0%, #d46b08 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <FieldTimeOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>CPU核时 (Core-hours)</div>
                <div style={{ fontSize: 22, fontWeight: 700, color: '#fa8c16' }}>
                  {getCoreHours()} <span style={{ fontSize: 14, fontWeight: 400 }}>核时</span>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* CPU核数 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #52c41a 0%, #389e0d 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <SettingOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>CPU核数 (Cores)</div>
                <div style={{ fontSize: 22, fontWeight: 700, color: '#52c41a' }}>
                  {job.config?.slurm_cpus || 16} <span style={{ fontSize: 14, fontWeight: 400 }}>核</span>
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* 计算分区 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #722ed1 0%, #531dab 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <DatabaseOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>计算分区 (Partition)</div>
                <div style={{ fontSize: 22, fontWeight: 700, color: '#722ed1' }}>
                  {job.config?.slurm_partition || 'default'}
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* 主内容区 */}
      <Row gutter={16}>
        {/* 左侧：任务配置和时间线 */}
        <Col span={8}>
          {/* 分子信息 */}
          <Card
            title={<Space><ExperimentOutlined />分子信息</Space>}
            size="small"
            style={{ marginBottom: 16 }}
          >
            <Descriptions column={1} size="small">
              <Descriptions.Item label="分子名称">
                <Text strong>{job.molecule_name}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="SMILES">
                <Text copyable code style={{ fontSize: 12, wordBreak: 'break-all' }}>
                  {job.smiles}
                </Text>
              </Descriptions.Item>
              <Descriptions.Item label="分子类型">
                <Tag>{job.molecule_type}</Tag>
              </Descriptions.Item>
              <Descriptions.Item label="电荷">
                {job.charge}
              </Descriptions.Item>
              <Descriptions.Item label="自旋多重度">
                {job.spin_multiplicity}
              </Descriptions.Item>
            </Descriptions>
          </Card>

          {/* 计算参数 */}
          <Card
            title={<Space><SettingOutlined />计算参数</Space>}
            size="small"
            style={{ marginBottom: 16 }}
          >
            <Descriptions column={1} size="small">
              <Descriptions.Item label="泛函">
                <Tag color="blue">{job.functional}</Tag>
              </Descriptions.Item>
              <Descriptions.Item label="基组">
                <Tag color="green">{job.basis_set}</Tag>
              </Descriptions.Item>
              <Descriptions.Item label="溶剂模型">
                {(() => {
                  // 优先从顶层的 solvent_config 读取，如果没有则从 config.solvent_config 读取
                  const solventConfig = (job as any).solvent_config || job.config?.solvent_config;
                  if (!solventConfig) {
                    return <Tag color="default">气相</Tag>;
                  }
                  const model = solventConfig.model || 'gas';
                  const solventName = solventConfig.solvent_name;

                  if (model === 'gas') {
                    return <Tag color="default">气相</Tag>;
                  } else if (model === 'pcm') {
                    return <Tag color="cyan">PCM - {solventName || 'Water'}</Tag>;
                  } else if (model === 'smd') {
                    return <Tag color="blue">SMD - {solventName || 'Water'}</Tag>;
                  } else if (model === 'custom') {
                    return <Tag color="purple">自定义 - {solventName || '自定义溶剂'}</Tag>;
                  }
                  return <Tag color="default">{model}</Tag>;
                })()}
              </Descriptions.Item>
              <Descriptions.Item label="计算分区">
                {job.config?.slurm_partition || '-'}
              </Descriptions.Item>
              <Descriptions.Item label="CPU核数">
                {job.config?.slurm_cpus || 16} 核
              </Descriptions.Item>
            </Descriptions>
          </Card>

          {/* 任务时间线 */}
          <Card
            title={<Space><ClockCircleOutlined />任务进度</Space>}
            size="small"
            style={{ marginBottom: 16 }}
          >
            <Timeline items={getStatusTimeline()} />
            {['QUEUED', 'RUNNING'].includes(job.status) && (
              <div style={{ textAlign: 'center', color: '#999', marginTop: 8 }}>
                <HourglassOutlined spin /> 任务进行中...
              </div>
            )}
          </Card>

          {/* 关联的MD任务 */}
          {job.md_job_id && (
            <Card
              title="关联任务"
              size="small"
              style={{ marginBottom: 16 }}
            >
              <Button
                type="link"
                style={{ padding: 0 }}
                onClick={() => navigate(`/workspace/jobs/${job.md_job_id}/detail`)}
              >
                查看关联的MD任务 #{job.md_job_id}
              </Button>
            </Card>
          )}
        </Col>

        {/* 右侧：计算结果 */}
        <Col span={16}>
          {job.status === 'COMPLETED' && results.length > 0 ? (
            <QCResultsPanel results={results} job={job} />
          ) : (
            <Card
              title={<Space><ThunderboltOutlined />计算结果</Space>}
              style={{ minHeight: 400 }}
            >
              <div style={{
                textAlign: 'center',
                padding: '80px 40px',
                background: token.colorBgContainer,
                borderRadius: 8
              }}>
                {job.status === 'COMPLETED' ? (
                  <>
                    <ExperimentOutlined style={{ fontSize: 48, color: '#faad14' }} />
                    <Paragraph type="warning" style={{ marginTop: 16, fontSize: 16 }}>
                      计算已完成，但暂无结果数据
                    </Paragraph>
                    <Paragraph type="secondary" style={{ fontSize: 14, marginTop: 8 }}>
                      可能原因：
                    </Paragraph>
                    <ul style={{ textAlign: 'left', display: 'inline-block', color: '#8c8c8c' }}>
                      <li>Gaussian 计算未正常结束（检查日志文件）</li>
                      <li>Worker 解析结果失败（联系管理员查看 Worker 日志）</li>
                      <li>结果上传到云端失败（网络或认证问题）</li>
                    </ul>
                    {job.error_message && (
                      <Alert
                        type="error"
                        message="详细错误信息"
                        description={job.error_message}
                        style={{ marginTop: 16, textAlign: 'left' }}
                        showIcon
                      />
                    )}
                  </>
                ) : job.status === 'FAILED' ? (
                  <>
                    <CloseCircleOutlined style={{ fontSize: 48, color: '#ff4d4f' }} />
                    <Paragraph type="danger" style={{ marginTop: 16, fontSize: 16 }}>
                      计算失败，请检查错误信息
                    </Paragraph>
                  </>
                ) : (
                  <>
                    {job.status === 'CREATED' ? (
                      <PlayCircleOutlined style={{ fontSize: 48, color: '#d9d9d9' }} />
                    ) : (
                      <SyncOutlined spin style={{ fontSize: 48, color: '#1890ff' }} />
                    )}
                    <Paragraph type="secondary" style={{ marginTop: 16, fontSize: 16 }}>
                      {job.status === 'CREATED' ? '任务尚未提交，点击上方"提交计算"按钮开始' :
                       job.status === 'QUEUED' ? '任务正在排队等待计算资源...' :
                       job.status === 'RUNNING' ? '正在进行量子化学计算，请耐心等待...' :
                       job.status === 'POSTPROCESSING' ? '计算完成，正在处理结果数据...' :
                       '等待计算完成'}
                    </Paragraph>
                    {job.status !== 'CREATED' && (
                      <Text type="secondary">页面将自动刷新显示最新状态</Text>
                    )}
                  </>
                )}
              </div>
            </Card>
          )}
        </Col>
      </Row>
    </div>
  );
}

