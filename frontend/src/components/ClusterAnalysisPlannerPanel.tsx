/**
 * Cluster 高级计算规划面板
 * 
 * 功能：
 * 1. 结构筛选和选择
 * 2. 计算类型多选（Binding/Desolvation/Redox/Reorg）
 * 3. QC 任务复用预览
 * 4. 提交和追踪任务
 */
import React, { useState, useEffect, useCallback } from 'react';
import {
  Card,
  Table,
  Button,
  Space,
  Tag,
  Progress,
  Select,
  Checkbox,
  Row,
  Col,
  Typography,
  message,
  Tooltip,
  Spin,
  Empty,
  Divider,
  Alert,
  Modal,
  Statistic,
  Badge,
  Collapse,
} from 'antd';
import {
  ThunderboltOutlined,
  ReloadOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  ExclamationCircleOutlined,
  SendOutlined,
  InfoCircleOutlined,
  RocketOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import { autoSelectSolvationStructures, type AutoSelectedStructure } from '../api/jobs';
import {
  planClusterAnalysis,
  submitClusterAnalysis,
  listClusterAnalysisJobs,
  getClusterAnalysisResults,
  getClusterAnalysisQCStatus,
  CALC_TYPE_INFO,
  type ClusterCalcType,
  type ClusterAnalysisPlanResponse,
  type AdvancedClusterJob,
  type CalcTypeRequirements,
  type ClusterAnalysisResults,
  type QCStatus,
} from '../api/clusterAnalysis';

const { Text, Title, Paragraph } = Typography;
const { Panel } = Collapse;

interface Props {
  mdJobId: number;
}

// 计算类型选项
const CALC_TYPE_OPTIONS: { value: ClusterCalcType; label: string; risk: string }[] = [
  { value: 'BINDING_TOTAL', label: '🔗 总 Binding Energy', risk: 'low' },
  { value: 'BINDING_PAIRWISE', label: '⚛️ 分子-Li Binding', risk: 'low' },
  { value: 'DESOLVATION_STEPWISE', label: '📉 逐级去溶剂化', risk: 'medium' },
  { value: 'DESOLVATION_FULL', label: '🎯 完全去溶剂化', risk: 'low' },
  { value: 'REDOX', label: '⚡ 氧化还原电位', risk: 'high' },
  { value: 'REORGANIZATION', label: '🔄 Marcus 重组能', risk: 'high' },
];

export default function ClusterAnalysisPlannerPanel({ mdJobId }: Props) {
  // 状态
  const [loading, setLoading] = useState(false);
  const [structures, setStructures] = useState<AutoSelectedStructure[]>([]);
  const [selectedStructureIds, setSelectedStructureIds] = useState<number[]>([]);
  const [selectedCalcTypes, setSelectedCalcTypes] = useState<ClusterCalcType[]>([]);
  const [planResult, setPlanResult] = useState<ClusterAnalysisPlanResponse | null>(null);
  const [planLoading, setPlanLoading] = useState(false);
  const [submitting, setSubmitting] = useState(false);
  const [existingJobs, setExistingJobs] = useState<AdvancedClusterJob[]>([]);
  
  // QC 配置
  const [qcConfig, setQcConfig] = useState({
    functional: 'B3LYP',
    basis_set: '6-31G*',
    use_dispersion: true,
    charge_ion: 1,
  });

  // 加载溶剂化结构
  const loadStructures = useCallback(async () => {
    setLoading(true);
    try {
      const result = await autoSelectSolvationStructures(mdJobId, { max_per_composition: 100 });
      setStructures(result.selected);
    } catch (error) {
      message.error('加载溶剂化结构失败');
      console.error(error);
    } finally {
      setLoading(false);
    }
  }, [mdJobId]);

  // 加载已有任务
  const loadExistingJobs = useCallback(async () => {
    try {
      const jobs = await listClusterAnalysisJobs(mdJobId);
      setExistingJobs(jobs);
    } catch (error) {
      console.error('加载已有任务失败:', error);
    }
  }, [mdJobId]);

  useEffect(() => {
    loadStructures();
    loadExistingJobs();
  }, [loadStructures, loadExistingJobs]);

  // 规划计算
  const handlePlan = async () => {
    if (selectedStructureIds.length === 0) {
      message.warning('请先选择溶剂化结构');
      return;
    }
    if (selectedCalcTypes.length === 0) {
      message.warning('请选择至少一种计算类型');
      return;
    }

    setPlanLoading(true);
    try {
      const result = await planClusterAnalysis({
        md_job_id: mdJobId,
        solvation_structure_ids: selectedStructureIds,
        calc_types: selectedCalcTypes,
        qc_config: qcConfig,
      });
      setPlanResult(result);
      message.success('规划完成');
    } catch (error) {
      message.error('规划失败');
      console.error(error);
    } finally {
      setPlanLoading(false);
    }
  };

  // 提交任务
  const handleSubmit = async () => {
    if (!planResult) return;

    Modal.confirm({
      title: '确认提交计算任务',
      content: (
        <div>
          <p>将提交以下计算：</p>
          <ul>
            {selectedCalcTypes.map(ct => (
              <li key={ct}>{CALC_TYPE_INFO[ct].label}</li>
            ))}
          </ul>
          <p>
            <strong>新建 QC 任务：</strong> {planResult.total_new_qc_tasks} 个
          </p>
          <p>
            <strong>复用已有任务：</strong> {planResult.total_reused_qc_tasks} 个
          </p>
          <p>
            <strong>预估时间：</strong> {planResult.estimated_compute_hours.toFixed(1)} 小时
          </p>
        </div>
      ),
      okText: '提交',
      cancelText: '取消',
      onOk: async () => {
        setSubmitting(true);
        try {
          await submitClusterAnalysis({
            md_job_id: mdJobId,
            solvation_structure_ids: selectedStructureIds,
            calc_types: selectedCalcTypes,
            qc_config: qcConfig,
          });
          message.success('任务已提交');
          setPlanResult(null);
          loadExistingJobs();
        } catch (error) {
          message.error('提交失败');
          console.error(error);
        } finally {
          setSubmitting(false);
        }
      },
    });
  };

  // 结构表格列定义
  const structureColumns: ColumnsType<AutoSelectedStructure> = [
    {
      title: '结构 ID',
      dataIndex: 'structure_id',
      key: 'structure_id',
      width: 80,
    },
    {
      title: '配位数',
      dataIndex: 'coordination_number',
      key: 'coordination_number',
      width: 80,
      render: (cn: number) => <Tag color="blue">{cn}</Tag>,
    },
    {
      title: '组成',
      dataIndex: 'composition',
      key: 'composition',
      render: (comp: Record<string, number>) => (
        <Space size="small" wrap>
          {Object.entries(comp || {}).map(([mol, count]) => (
            count > 0 && <Tag key={mol}>{mol}: {count}</Tag>
          ))}
        </Space>
      ),
    },
    {
      title: '帧号',
      dataIndex: 'frame',
      key: 'frame',
      width: 80,
    },
  ];

  // 渲染计算类型复选框
  const renderCalcTypeCheckboxes = () => (
    <div style={{ marginBottom: 16 }}>
      <Text strong>选择计算类型：</Text>
      <Row gutter={[16, 8]} style={{ marginTop: 8 }}>
        {CALC_TYPE_OPTIONS.map(opt => (
          <Col key={opt.value} span={8}>
            <Checkbox
              checked={selectedCalcTypes.includes(opt.value)}
              onChange={e => {
                if (e.target.checked) {
                  setSelectedCalcTypes([...selectedCalcTypes, opt.value]);
                } else {
                  setSelectedCalcTypes(selectedCalcTypes.filter(t => t !== opt.value));
                }
                setPlanResult(null); // 清除之前的规划结果
              }}
            >
              <Space>
                {opt.label}
                {opt.risk === 'high' && (
                  <Tag color="red" style={{ marginLeft: 4 }}>⚠️ 高风险</Tag>
                )}
                {opt.risk === 'medium' && (
                  <Tag color="orange" style={{ marginLeft: 4 }}>中等</Tag>
                )}
              </Space>
            </Checkbox>
          </Col>
        ))}
      </Row>
    </div>
  );

  // 渲染规划结果
  const renderPlanResult = () => {
    if (!planResult) return null;

    return (
      <Card
        title={<><RocketOutlined /> QC 任务规划预览</>}
        style={{ marginTop: 16 }}
        extra={
          <Button
            type="primary"
            icon={<SendOutlined />}
            loading={submitting}
            onClick={handleSubmit}
          >
            提交计算
          </Button>
        }
      >
        {/* 汇总统计 */}
        <Row gutter={16} style={{ marginBottom: 16 }}>
          <Col span={6}>
            <Statistic
              title="选中结构"
              value={planResult.selected_structures_count}
              suffix="个"
            />
          </Col>
          <Col span={6}>
            <Statistic
              title="新建 QC 任务"
              value={planResult.total_new_qc_tasks}
              suffix="个"
              valueStyle={{ color: '#1890ff' }}
            />
          </Col>
          <Col span={6}>
            <Statistic
              title="复用已有任务"
              value={planResult.total_reused_qc_tasks}
              suffix="个"
              valueStyle={{ color: '#52c41a' }}
            />
          </Col>
          <Col span={6}>
            <Statistic
              title="预估时间"
              value={planResult.estimated_compute_hours.toFixed(1)}
              suffix="小时"
            />
          </Col>
        </Row>

        {/* 警告 */}
        {planResult.warnings.length > 0 && (
          <Alert
            type="warning"
            message="注意事项"
            description={
              <ul style={{ margin: 0, paddingLeft: 20 }}>
                {planResult.warnings.map((w, i) => <li key={i}>{w}</li>)}
              </ul>
            }
            style={{ marginBottom: 16 }}
          />
        )}

        {/* 各计算类型详情 */}
        <Collapse>
          {planResult.calc_requirements.map(req => (
            <Panel
              key={req.calc_type}
              header={
                <Space>
                  {CALC_TYPE_INFO[req.calc_type].icon} {CALC_TYPE_INFO[req.calc_type].label}
                  <Tag color="blue">新建 {req.new_tasks_count}</Tag>
                  <Tag color="green">复用 {req.reused_tasks_count}</Tag>
                </Space>
              }
            >
              <Paragraph type="secondary">
                公式：<code>{CALC_TYPE_INFO[req.calc_type].formula}</code>
              </Paragraph>
              <Table
                size="small"
                dataSource={req.required_qc_tasks}
                rowKey={(_, i) => `${req.calc_type}-${i}`}
                pagination={false}
                columns={[
                  { title: '类型', dataIndex: 'task_type', width: 100 },
                  { title: '描述', dataIndex: 'description' },
                  {
                    title: '状态',
                    dataIndex: 'status',
                    width: 100,
                    render: (status: string) => (
                      status === 'reused'
                        ? <Tag color="green"><CheckCircleOutlined /> 复用</Tag>
                        : <Tag color="blue"><ClockCircleOutlined /> 新建</Tag>
                    )
                  },
                ]}
              />
            </Panel>
          ))}
        </Collapse>
      </Card>
    );
  };

  // 查看结果的任务 ID
  const [viewingJobId, setViewingJobId] = useState<number | null>(null);

  // 渲染已有任务列表
  const renderExistingJobs = () => {
    if (existingJobs.length === 0) return null;

    return (
      <Card
        title="已有计算任务"
        style={{ marginTop: 16 }}
        size="small"
      >
        <Table
          size="small"
          dataSource={existingJobs}
          rowKey="id"
          pagination={false}
          columns={[
            { title: 'ID', dataIndex: 'id', width: 60 },
            {
              title: '计算类型',
              dataIndex: 'calc_types',
              render: (types: string[]) => (
                <Space size="small" wrap>
                  {types.map(t => (
                    <Tag key={t}>{CALC_TYPE_INFO[t as ClusterCalcType]?.icon} {t}</Tag>
                  ))}
                </Space>
              )
            },
            {
              title: '状态',
              dataIndex: 'status',
              width: 120,
              render: (status: string) => {
                const colors: Record<string, string> = {
                  COMPLETED: 'green',
                  RUNNING: 'blue',
                  WAITING_QC: 'orange',
                  FAILED: 'red',
                  SUBMITTED: 'cyan',
                };
                return <Tag color={colors[status] || 'default'}>{status}</Tag>;
              }
            },
            {
              title: '进度',
              dataIndex: 'progress',
              width: 100,
              render: (p: number) => <Progress percent={Math.round(p)} size="small" />
            },
            {
              title: '创建时间',
              dataIndex: 'created_at',
              width: 150,
              render: (t: string) => new Date(t).toLocaleString()
            },
            {
              title: '操作',
              key: 'action',
              width: 100,
              render: (_: unknown, record: AdvancedClusterJob) => (
                <Button
                  type="link"
                  size="small"
                  onClick={() => setViewingJobId(record.id)}
                >
                  查看结果
                </Button>
              )
            },
          ]}
        />

        {/* 结果查看模态框 */}
        <Modal
          title={`计算结果 #${viewingJobId}`}
          open={viewingJobId !== null}
          onCancel={() => setViewingJobId(null)}
          footer={null}
          width={900}
          destroyOnClose
        >
          {viewingJobId && (
            <ClusterAnalysisResultsView
              jobId={viewingJobId}
              onClose={() => setViewingJobId(null)}
            />
          )}
        </Modal>
      </Card>
    );
  };

  return (
    <Card
      title={
        <Space>
          <ExperimentOutlined />
          Cluster 高级计算规划
        </Space>
      }
      extra={
        <Button icon={<ReloadOutlined />} onClick={() => { loadStructures(); loadExistingJobs(); }}>
          刷新
        </Button>
      }
    >
      <Spin spinning={loading}>
        {structures.length === 0 ? (
          <Empty description="暂无溶剂化结构，请先完成 MD 计算" />
        ) : (
          <>
            {/* 步骤 1: 选择结构 */}
            <Card type="inner" title="步骤 1: 选择溶剂化结构" style={{ marginBottom: 16 }}>
              <Table
                size="small"
                dataSource={structures}
                columns={structureColumns}
                rowKey="structure_id"
                rowSelection={{
                  selectedRowKeys: selectedStructureIds,
                  onChange: keys => {
                    setSelectedStructureIds(keys as number[]);
                    setPlanResult(null);
                  },
                }}
                pagination={{ pageSize: 10 }}
              />
              <div style={{ marginTop: 8 }}>
                <Text type="secondary">
                  已选择 {selectedStructureIds.length} / {structures.length} 个结构
                </Text>
                <Button
                  type="link"
                  onClick={() => setSelectedStructureIds(structures.map(s => s.structure_id))}
                >
                  全选
                </Button>
                <Button
                  type="link"
                  onClick={() => setSelectedStructureIds([])}
                >
                  清空
                </Button>
              </div>
            </Card>

            {/* 步骤 2: 选择计算类型 */}
            <Card type="inner" title="步骤 2: 选择计算类型" style={{ marginBottom: 16 }}>
              {renderCalcTypeCheckboxes()}

              {/* 选中的计算类型说明 */}
              {selectedCalcTypes.length > 0 && (
                <Alert
                  type="info"
                  message={`已选择 ${selectedCalcTypes.length} 种计算`}
                  description={
                    <ul style={{ margin: 0, paddingLeft: 20 }}>
                      {selectedCalcTypes.map(ct => (
                        <li key={ct}>
                          <strong>{CALC_TYPE_INFO[ct].label}</strong>：{CALC_TYPE_INFO[ct].description}
                        </li>
                      ))}
                    </ul>
                  }
                />
              )}
            </Card>

            {/* 步骤 3: 规划预览 */}
            <Card type="inner" title="步骤 3: 规划与提交" style={{ marginBottom: 16 }}>
              <Space>
                <Button
                  type="primary"
                  icon={<ThunderboltOutlined />}
                  loading={planLoading}
                  onClick={handlePlan}
                  disabled={selectedStructureIds.length === 0 || selectedCalcTypes.length === 0}
                >
                  生成规划预览
                </Button>
                <Text type="secondary">
                  点击查看需要的 QC 任务和可复用的已有结果
                </Text>
              </Space>

              {renderPlanResult()}
            </Card>

            {/* 已有任务 */}
            {renderExistingJobs()}
          </>
        )}
      </Spin>
    </Card>
  );
}

// ============================================================================
// 内联结果查看组件
// ============================================================================

interface ResultsViewProps {
  jobId: number;
  onClose: () => void;
}

function ClusterAnalysisResultsView({ jobId, onClose }: ResultsViewProps) {
  const [loading, setLoading] = useState(true);
  const [results, setResults] = useState<ClusterAnalysisResults | null>(null);
  const [qcStatus, setQcStatus] = useState<QCStatus | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const fetchData = async () => {
      try {
        setLoading(true);
        const [resultsData, qcData] = await Promise.all([
          getClusterAnalysisResults(jobId),
          getClusterAnalysisQCStatus(jobId),
        ]);
        setResults(resultsData);
        setQcStatus(qcData);
      } catch (err) {
        setError((err as Error).message || '获取结果失败');
      } finally {
        setLoading(false);
      }
    };
    fetchData();
  }, [jobId]);

  if (loading) {
    return <Spin tip="加载中..." style={{ display: 'block', textAlign: 'center', padding: 40 }} />;
  }

  if (error) {
    return <Alert type="error" message={error} />;
  }

  if (!results) {
    return <Empty description="暂无结果" />;
  }

  return (
    <div>
      {/* QC 任务进度 */}
      {qcStatus && qcStatus.total_qc_jobs > 0 && (
        <Card size="small" style={{ marginBottom: 16 }}>
          <Row gutter={16}>
            <Col span={6}>
              <Statistic title="已完成" value={qcStatus.completed} valueStyle={{ color: '#52c41a' }} />
            </Col>
            <Col span={6}>
              <Statistic title="运行中" value={qcStatus.running} valueStyle={{ color: '#1890ff' }} />
            </Col>
            <Col span={6}>
              <Statistic title="等待中" value={qcStatus.pending} valueStyle={{ color: '#faad14' }} />
            </Col>
            <Col span={6}>
              <Statistic title="失败" value={qcStatus.failed} valueStyle={{ color: qcStatus.failed > 0 ? '#ff4d4f' : undefined }} />
            </Col>
          </Row>
          <Progress percent={Math.round((qcStatus.completed / qcStatus.total_qc_jobs) * 100)} style={{ marginTop: 16 }} />
        </Card>
      )}

      {/* 各类型结果 */}
      {results.calc_types.map((calcType) => {
        const info = CALC_TYPE_INFO[calcType as ClusterCalcType];
        const calcResult = results.results?.[calcType] as Record<string, unknown>;

        return (
          <Card
            key={calcType}
            size="small"
            title={<span>{info?.icon} {info?.label || calcType}</span>}
            style={{ marginBottom: 16 }}
            extra={
              calcResult?.error ? (
                <Tag color="red">失败</Tag>
              ) : calcResult && Object.keys(calcResult).length > 0 ? (
                <Tag color="green">完成</Tag>
              ) : (
                <Tag>等待</Tag>
              )
            }
          >
            <Text type="secondary">{info?.description}</Text>
            <div style={{ marginTop: 8 }}>
              <Text code>{info?.formula}</Text>
            </div>

            {calcResult?.error && (
              <Alert type="error" message={calcResult.error as string} style={{ marginTop: 8 }} />
            )}

            {calcResult && !calcResult.error && Object.keys(calcResult).length > 0 && (
              <div style={{ marginTop: 16 }}>
                {renderResultContent(calcType as ClusterCalcType, calcResult)}
              </div>
            )}
          </Card>
        );
      })}
    </div>
  );
}

function renderResultContent(calcType: ClusterCalcType, result: Record<string, unknown>): React.ReactNode {
  switch (calcType) {
    case 'BINDING_TOTAL':
    case 'DESOLVATION_FULL':
      return (
        <Row gutter={16}>
          <Col span={8}>
            <Statistic
              title="Binding Energy"
              value={(result.e_bind_kcal_mol as number)?.toFixed(2)}
              suffix="kcal/mol"
              precision={2}
            />
          </Col>
          <Col span={8}>
            <Statistic
              title="eV"
              value={(result.e_bind_ev as number)?.toFixed(4)}
              precision={4}
            />
          </Col>
          <Col span={8}>
            <Statistic
              title="Hartree"
              value={(result.e_bind_au as number)?.toFixed(6)}
              precision={6}
            />
          </Col>
        </Row>
      );

    case 'BINDING_PAIRWISE':
      const pairBindings = (result.pairwise_bindings as Array<Record<string, unknown>>) || [];
      return (
        <Table
          size="small"
          dataSource={pairBindings}
          rowKey={(_, i) => i?.toString() || '0'}
          pagination={false}
          columns={[
            { title: '配体', dataIndex: 'ligand' },
            { title: 'E_bind (kcal/mol)', dataIndex: 'e_bind_kcal_mol', render: (v: number) => v?.toFixed(2) },
            { title: 'E_bind (eV)', dataIndex: 'e_bind_ev', render: (v: number) => v?.toFixed(4) },
          ]}
        />
      );

    case 'DESOLVATION_STEPWISE':
      const steps = (result.stepwise_desolvation as Array<Record<string, unknown>>) || [];
      return (
        <Table
          size="small"
          dataSource={steps}
          rowKey={(_, i) => i?.toString() || '0'}
          pagination={false}
          columns={[
            { title: '移除配体', dataIndex: 'ligand' },
            { title: 'ΔE (kcal/mol)', dataIndex: 'delta_e_kcal_mol', render: (v: number) => v?.toFixed(2) },
            { title: 'ΔE (eV)', dataIndex: 'delta_e_ev', render: (v: number) => v?.toFixed(4) },
          ]}
        />
      );

    case 'REDOX':
      const potentials = (result.redox_potentials as Array<Record<string, unknown>>) || [];
      return (
        <Table
          size="small"
          dataSource={potentials}
          rowKey={(_, i) => i?.toString() || '0'}
          pagination={false}
          columns={[
            { title: 'SMILES', dataIndex: 'smiles', render: (s: string) => <Text code>{s}</Text> },
            { title: 'ΔG (eV)', dataIndex: 'delta_g_sol_ev', render: (v: number) => v?.toFixed(4) },
            { title: 'E° (V vs SHE)', dataIndex: 'oxidation_potential_v', render: (v: number) => v?.toFixed(3) },
          ]}
        />
      );

    case 'REORGANIZATION':
      if (result.status === 'not_implemented') {
        return <Alert type="info" message={result.message as string} />;
      }
      return <pre style={{ fontSize: 12 }}>{JSON.stringify(result, null, 2)}</pre>;

    default:
      return <pre style={{ fontSize: 12 }}>{JSON.stringify(result, null, 2)}</pre>;
  }
}
