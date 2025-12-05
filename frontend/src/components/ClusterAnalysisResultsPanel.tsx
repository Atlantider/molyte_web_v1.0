/**
 * Cluster 高级计算结果展示面板
 *
 * 使用 Ant Design 组件库，与项目其他组件保持一致
 */
import React, { useState, useEffect, useCallback } from 'react';
import {
  Card,
  Table,
  Tag,
  Progress,
  Alert,
  Button,
  Collapse,
  Tooltip,
  Row,
  Col,
  Statistic,
  Space,
  Spin,
  Typography,
  theme,
} from 'antd';
import {
  ReloadOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  ArrowLeftOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import {
  getClusterAnalysisJob,
  getClusterAnalysisResults,
  getClusterAnalysisQCStatus,
  CALC_TYPE_INFO,
  AdvancedClusterJob,
  ClusterAnalysisResults,
  QCStatus,
  ClusterCalcType,
} from '../api/clusterAnalysis';

const { Text, Title } = Typography;
const { Panel } = Collapse;

interface Props {
  jobId: number;
  onBack?: () => void;
}

const STATUS_CONFIG: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
  CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
  SUBMITTED: { color: 'blue', icon: <ClockCircleOutlined />, text: '已提交' },
  RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
  WAITING_QC: { color: 'orange', icon: <ClockCircleOutlined />, text: '等待 QC' },
  CALCULATING: { color: 'processing', icon: <SyncOutlined spin />, text: '计算中' },
  COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
  FAILED: { color: 'error', icon: <CloseCircleOutlined />, text: '失败' },
  CANCELLED: { color: 'default', icon: <CloseCircleOutlined />, text: '已取消' },
};

export default function ClusterAnalysisResultsPanel({ jobId, onBack }: Props) {
  const [job, setJob] = useState<AdvancedClusterJob | null>(null);
  const [results, setResults] = useState<ClusterAnalysisResults | null>(null);
  const [qcStatus, setQcStatus] = useState<QCStatus | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({});

  const fetchData = useCallback(async () => {
    try {
      setLoading(true);
      const [jobData, resultsData, qcData] = await Promise.all([
        getClusterAnalysisJob(jobId),
        getClusterAnalysisResults(jobId),
        getClusterAnalysisQCStatus(jobId),
      ]);
      setJob(jobData);
      setResults(resultsData);
      setQcStatus(qcData);
      setError(null);
    } catch (err) {
      setError((err as Error).message || '获取数据失败');
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  useEffect(() => {
    fetchData();
    // 如果任务还在进行中，定时刷新
    const interval = setInterval(() => {
      if (job && !['COMPLETED', 'FAILED', 'CANCELLED'].includes(job.status)) {
        fetchData();
      }
    }, 10000);
    return () => clearInterval(interval);
  }, [fetchData, job]);

  const toggleSection = (section: string) => {
    setExpandedSections((prev) => ({ ...prev, [section]: !prev[section] }));
  };

  const { token } = theme.useToken();

  if (loading && !job) {
    return (
      <div style={{ padding: 24, textAlign: 'center' }}>
        <Spin size="large" />
        <div style={{ marginTop: 16 }}>
          <Text type="secondary">加载中...</Text>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <Alert
        type="error"
        message="加载失败"
        description={error}
        action={<Button onClick={fetchData}>重试</Button>}
        style={{ margin: 16 }}
      />
    );
  }

  if (!job) return null;

  const statusConfig = STATUS_CONFIG[job.status] || STATUS_CONFIG.CREATED;

  return (
    <div style={{ padding: 16 }}>
      {/* 顶部信息 */}
      <Card style={{ marginBottom: 16 }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <Space align="center" style={{ marginBottom: 8 }}>
              <ExperimentOutlined style={{ fontSize: 24, color: token.colorPrimary }} />
              <Title level={4} style={{ margin: 0 }}>
                Cluster 高级计算结果 #{jobId}
              </Title>
            </Space>
            <Space>
              <Tag color={statusConfig.color} icon={statusConfig.icon}>
                {statusConfig.text}
              </Tag>
              <Text type="secondary">进度: {job.progress.toFixed(0)}%</Text>
            </Space>
          </div>
          <Space>
            <Tooltip title="刷新">
              <Button
                icon={<ReloadOutlined />}
                onClick={fetchData}
                loading={loading}
              />
            </Tooltip>
            {onBack && (
              <Button icon={<ArrowLeftOutlined />} onClick={onBack}>
                返回
              </Button>
            )}
          </Space>
        </div>

        {/* 进度条 */}
        {job.status !== 'COMPLETED' && job.status !== 'FAILED' && (
          <Progress
            percent={job.progress}
            status={job.status === 'RUNNING' || job.status === 'CALCULATING' ? 'active' : 'normal'}
            style={{ marginTop: 16 }}
          />
        )}
      </Card>

      {/* QC 任务状态 */}
      {qcStatus && qcStatus.total_qc_jobs > 0 && (
        <Card title="QC 任务状态" style={{ marginBottom: 16 }}>
          <Row gutter={16}>
            <Col span={6}>
              <Statistic
                title="已完成"
                value={qcStatus.completed}
                valueStyle={{ color: token.colorSuccess }}
                prefix={<CheckCircleOutlined />}
              />
            </Col>
            <Col span={6}>
              <Statistic
                title="运行中"
                value={qcStatus.running}
                valueStyle={{ color: token.colorPrimary }}
                prefix={<SyncOutlined spin={qcStatus.running > 0} />}
              />
            </Col>
            <Col span={6}>
              <Statistic
                title="等待中"
                value={qcStatus.pending}
                valueStyle={{ color: token.colorWarning }}
                prefix={<ClockCircleOutlined />}
              />
            </Col>
            <Col span={6}>
              <Statistic
                title="失败"
                value={qcStatus.failed}
                valueStyle={{ color: qcStatus.failed > 0 ? token.colorError : undefined }}
                prefix={<CloseCircleOutlined />}
              />
            </Col>
          </Row>
        </Card>
      )}

      {/* 计算类型和结果 */}
      <Collapse
        defaultActiveKey={job.calc_types}
        style={{ marginBottom: 16 }}
      >
        {job.calc_types.map((calcType) => {
          const info = CALC_TYPE_INFO[calcType as ClusterCalcType];
          const calcResult = results?.results?.[calcType] as Record<string, unknown>;

          const getResultStatus = () => {
            if (calcResult?.error) {
              return <Tag color="error" icon={<CloseCircleOutlined />}>失败</Tag>;
            }
            if (calcResult && Object.keys(calcResult).length > 0) {
              return <Tag color="success" icon={<CheckCircleOutlined />}>完成</Tag>;
            }
            return <Tag icon={<ClockCircleOutlined />}>等待</Tag>;
          };

          return (
            <Panel
              key={calcType}
              header={
                <Space>
                  <span>{info?.icon} {info?.label || calcType}</span>
                  {getResultStatus()}
                </Space>
              }
            >
              <div style={{ marginBottom: 8 }}>
                <Text type="secondary">{info?.description}</Text>
              </div>
              <div style={{ marginBottom: 16 }}>
                <Text code style={{ fontSize: 12 }}>{info?.formula}</Text>
              </div>

              {/* 结果展示 */}
              {calcResult && !calcResult.error && (
                <div>{renderCalcTypeResult(calcType as ClusterCalcType, calcResult)}</div>
              )}

              {calcResult?.error && (
                <Alert type="error" message={calcResult.error as string} />
              )}
            </Panel>
          );
        })}
      </Collapse>

      {/* 错误信息 */}
      {job.error_message && (
        <Alert type="error" message="任务错误" description={job.error_message} />
      )}
    </div>
  );
}

// ============================================================================
// 结果渲染函数
// ============================================================================

function renderCalcTypeResult(calcType: ClusterCalcType, result: Record<string, unknown>): React.ReactNode {
  switch (calcType) {
    case 'BINDING_TOTAL':
    case 'DESOLVATION_FULL':
      return renderBindingTotalResult(result);
    case 'BINDING_PAIRWISE':
      return renderBindingPairwiseResult(result);
    case 'DESOLVATION_STEPWISE':
      return renderDesolvationStepwiseResult(result);
    case 'REDOX':
      return renderRedoxResult(result);
    case 'REORGANIZATION':
      return renderReorganizationResult(result);
    default:
      return <pre style={{ background: '#f5f5f5', padding: 8, borderRadius: 4 }}>{JSON.stringify(result, null, 2)}</pre>;
  }
}

function renderBindingTotalResult(result: Record<string, unknown>) {
  return (
    <div>
      <Row gutter={16}>
        <Col span={8}>
          <Card size="small" style={{ textAlign: 'center' }}>
            <Statistic
              title="kcal/mol"
              value={((result.e_bind_kcal_mol as number) || 0).toFixed(2)}
              valueStyle={{ color: '#1890ff' }}
            />
          </Card>
        </Col>
        <Col span={8}>
          <Card size="small" style={{ textAlign: 'center' }}>
            <Statistic
              title="eV"
              value={((result.e_bind_ev as number) || 0).toFixed(4)}
              valueStyle={{ color: '#722ed1' }}
            />
          </Card>
        </Col>
        <Col span={8}>
          <Card size="small" style={{ textAlign: 'center' }}>
            <Statistic
              title="Hartree"
              value={((result.e_bind_au as number) || 0).toFixed(6)}
            />
          </Card>
        </Col>
      </Row>
      <div style={{ marginTop: 8 }}>
        <Text type="secondary" style={{ fontSize: 12 }}>
          E(cluster) = {((result.e_cluster_au as number) || 0).toFixed(6)} a.u.,
          E(ion) = {((result.e_ion_au as number) || 0).toFixed(6)} a.u.
        </Text>
      </div>
    </div>
  );
}

function renderBindingPairwiseResult(result: Record<string, unknown>) {
  const bindings = (result.pairwise_bindings as Array<Record<string, unknown>>) || [];

  const columns: ColumnsType<Record<string, unknown>> = [
    { title: '配体', dataIndex: 'ligand', key: 'ligand' },
    {
      title: 'E(dimer) / a.u.',
      dataIndex: 'e_dimer_au',
      key: 'e_dimer_au',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(6),
    },
    {
      title: 'E(ligand) / a.u.',
      dataIndex: 'e_ligand_au',
      key: 'e_ligand_au',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(6),
    },
    {
      title: 'E_bind / kcal·mol⁻¹',
      dataIndex: 'e_bind_kcal_mol',
      key: 'e_bind_kcal_mol',
      align: 'right',
      render: (v: number) => <Text strong>{(v || 0).toFixed(2)}</Text>,
    },
    {
      title: 'E_bind / eV',
      dataIndex: 'e_bind_ev',
      key: 'e_bind_ev',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(4),
    },
  ];

  return (
    <Table
      size="small"
      dataSource={bindings.map((b, i) => ({ ...b, key: i }))}
      columns={columns}
      pagination={false}
    />
  );
}

function renderDesolvationStepwiseResult(result: Record<string, unknown>) {
  const steps = (result.stepwise_desolvation as Array<Record<string, unknown>>) || [];

  const columns: ColumnsType<Record<string, unknown>> = [
    { title: '移除配体', dataIndex: 'ligand', key: 'ligand' },
    {
      title: 'E(cluster) / a.u.',
      dataIndex: 'e_cluster_au',
      key: 'e_cluster_au',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(6),
    },
    {
      title: 'E(minus) / a.u.',
      dataIndex: 'e_minus_au',
      key: 'e_minus_au',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(6),
    },
    {
      title: 'ΔE / kcal·mol⁻¹',
      dataIndex: 'delta_e_kcal_mol',
      key: 'delta_e_kcal_mol',
      align: 'right',
      render: (v: number) => (
        <Text strong style={{ color: v > 0 ? '#ff4d4f' : '#52c41a' }}>
          {(v || 0).toFixed(2)}
        </Text>
      ),
    },
    {
      title: 'ΔE / eV',
      dataIndex: 'delta_e_ev',
      key: 'delta_e_ev',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(4),
    },
  ];

  return (
    <Table
      size="small"
      dataSource={steps.map((s, i) => ({ ...s, key: i }))}
      columns={columns}
      pagination={false}
    />
  );
}

function renderRedoxResult(result: Record<string, unknown>) {
  const potentials = (result.redox_potentials as Array<Record<string, unknown>>) || [];

  const columns: ColumnsType<Record<string, unknown>> = [
    {
      title: '物种',
      dataIndex: 'smiles',
      key: 'smiles',
      render: (v: string) => <Text code>{v}</Text>,
    },
    {
      title: 'E(neutral,gas) / a.u.',
      dataIndex: 'e_neutral_gas_au',
      key: 'e_neutral_gas_au',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(6),
    },
    {
      title: 'E(charged,sol) / a.u.',
      dataIndex: 'e_charged_sol_au',
      key: 'e_charged_sol_au',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(6),
    },
    {
      title: 'ΔG(sol) / eV',
      dataIndex: 'delta_g_sol_ev',
      key: 'delta_g_sol_ev',
      align: 'right',
      render: (v: number) => (v || 0).toFixed(4),
    },
    {
      title: 'E° (vs SHE) / V',
      dataIndex: 'oxidation_potential_v',
      key: 'oxidation_potential_v',
      align: 'right',
      render: (v: number) => <Text strong style={{ color: '#1890ff' }}>{(v || 0).toFixed(3)}</Text>,
    },
  ];

  return (
    <Table
      size="small"
      dataSource={potentials.map((p, i) => ({ ...p, key: i }))}
      columns={columns}
      pagination={false}
    />
  );
}

function renderReorganizationResult(result: Record<string, unknown>) {
  if (result.status === 'not_implemented') {
    return <Alert type="info" message={result.message as string} />;
  }
  return <pre style={{ background: '#f5f5f5', padding: 8, borderRadius: 4 }}>{JSON.stringify(result, null, 2)}</pre>;
}

