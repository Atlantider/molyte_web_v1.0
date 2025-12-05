/**
 * 后处理详情页面
 * - 创建模式：选择 MD Job → 选择结构 → 选择计算类型 → 提交
 * - 查看模式：显示任务状态、QC 进度、计算结果
 */
import { useState, useEffect, useCallback, useMemo } from 'react';
import { useParams, useNavigate, useSearchParams } from 'react-router-dom';
import {
  Card,
  Steps,
  Button,
  Space,
  Select,
  Table,
  Checkbox,
  Tag,
  Typography,
  Row,
  Col,
  Statistic,
  Alert,
  Spin,
  Empty,
  message,
  Modal,
  Descriptions,
  Progress,
  Tooltip,
  Divider,
  theme,
} from 'antd';
import {
  ArrowLeftOutlined,
  RocketOutlined,
  ExperimentOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  CloseCircleOutlined,
  SendOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import {
  listClusterAnalysisJobs,
  getClusterAnalysisJob,
  planClusterAnalysis,
  submitClusterAnalysis,
  getClusterAnalysisResults,
  CALC_TYPE_INFO,
  type AdvancedClusterJob,
  type ClusterCalcType,
  type ClusterAnalysisPlanResponse,
} from '../api/clusterAnalysis';
import { getMDJobs, getMDJob, getSolvationStructures, type SolvationStructure } from '../api/jobs';
import type { MDJob } from '../types';
import { JobStatus } from '../types';
import ClusterAnalysisResultsPanel from '../components/ClusterAnalysisResultsPanel';
import dayjs from 'dayjs';

const { Title, Text, Paragraph } = Typography;

// 状态配置
const STATUS_CONFIG: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
  CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
  SUBMITTED: { color: 'processing', icon: <SyncOutlined spin />, text: '已提交' },
  RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
  WAITING_QC: { color: 'warning', icon: <ThunderboltOutlined />, text: '等待 QC' },
  CALCULATING: { color: 'processing', icon: <SyncOutlined spin />, text: '计算中' },
  COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
  FAILED: { color: 'error', icon: <CloseCircleOutlined />, text: '失败' },
  CANCELLED: { color: 'default', icon: <CloseCircleOutlined />, text: '已取消' },
};

// 计算类型选项
const CALC_TYPE_OPTIONS: { value: ClusterCalcType; label: string; description: string; riskLevel: string }[] = [
  { value: 'BINDING_TOTAL', label: '总 Binding Energy', description: '整个溶剂化簇的脱溶剂化能', riskLevel: 'low' },
  { value: 'BINDING_PAIRWISE', label: '分子-Li Binding', description: '单分子与 Li+ 的结合能', riskLevel: 'low' },
  { value: 'DESOLVATION_STEPWISE', label: '逐级去溶剂化', description: '逐个移除配体的能量变化', riskLevel: 'medium' },
  { value: 'DESOLVATION_FULL', label: '完全去溶剂化', description: '完整簇到裸离子的总能量', riskLevel: 'low' },
  { value: 'REDOX', label: '氧化还原电位', description: '热力学循环法计算', riskLevel: 'high' },
  { value: 'REORGANIZATION', label: 'Marcus 重组能', description: 'Marcus 理论 4 点方案', riskLevel: 'high' },
];

export default function PostProcessDetail() {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const { token } = theme.useToken();

  const isCreateMode = !id || id === 'create';

  // 通用状态
  const [loading, setLoading] = useState(false);
  const [job, setJob] = useState<AdvancedClusterJob | null>(null);

  // 创建模式状态
  const [currentStep, setCurrentStep] = useState(0);
  const [mdJobs, setMdJobs] = useState<MDJob[]>([]);
  const [selectedMdJobId, setSelectedMdJobId] = useState<number | undefined>(
    searchParams.get('md_job_id') ? Number(searchParams.get('md_job_id')) : undefined
  );
  const [selectedMdJob, setSelectedMdJob] = useState<MDJob | null>(null);
  const [structures, setStructures] = useState<SolvationStructure[]>([]);
  const [selectedStructureIds, setSelectedStructureIds] = useState<number[]>([]);
  const [selectedCalcTypes, setSelectedCalcTypes] = useState<ClusterCalcType[]>([]);
  const [planResult, setPlanResult] = useState<ClusterAnalysisPlanResponse | null>(null);
  const [planLoading, setPlanLoading] = useState(false);
  const [submitLoading, setSubmitLoading] = useState(false);

  // 筛选状态
  const [filterCoordNums, setFilterCoordNums] = useState<number[]>([]);
  const [filterAnions, setFilterAnions] = useState<string[]>([]);
  const [filterSolvents, setFilterSolvents] = useState<string[]>([]);

  // 从 composition 生成簇名称，如 Na1MP2PF6
  const generateClusterName = (centerIon: string, composition: Record<string, number>): string => {
    if (!composition || Object.keys(composition).length === 0) return centerIon || '-';
    const parts = Object.entries(composition)
      .filter(([_, count]) => count > 0)
      .sort(([a], [b]) => a.localeCompare(b))
      .map(([mol, count]) => `${count}${mol}`);
    return `${centerIon || ''}${parts.join('')}`;
  };

  // 提取所有可用的筛选选项
  const filterOptions = useMemo(() => {
    const coordNums = new Set<number>();
    const anions = new Set<string>();
    const solvents = new Set<string>();

    // 常见阴离子列表
    const anionList = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'DFOB', 'BOB', 'Cl', 'Br', 'I', 'NO3'];

    structures.forEach(s => {
      if (s.coordination_num) coordNums.add(s.coordination_num);
      if (s.composition) {
        Object.keys(s.composition).forEach(mol => {
          if (anionList.some(a => mol.toUpperCase().includes(a.toUpperCase()))) {
            anions.add(mol);
          } else {
            solvents.add(mol);
          }
        });
      }
    });

    return {
      coordNums: Array.from(coordNums).sort((a, b) => a - b),
      anions: Array.from(anions).sort(),
      solvents: Array.from(solvents).sort(),
    };
  }, [structures]);

  // 筛选后的结构
  const filteredStructures = useMemo(() => {
    return structures.filter(s => {
      // 配位数筛选
      if (filterCoordNums.length > 0 && !filterCoordNums.includes(s.coordination_num)) {
        return false;
      }
      // 阴离子筛选
      if (filterAnions.length > 0) {
        const hasAnion = filterAnions.some(anion =>
          s.composition && s.composition[anion] && s.composition[anion] > 0
        );
        if (!hasAnion) return false;
      }
      // 溶剂筛选
      if (filterSolvents.length > 0) {
        const hasSolvent = filterSolvents.some(solvent =>
          s.composition && s.composition[solvent] && s.composition[solvent] > 0
        );
        if (!hasSolvent) return false;
      }
      return true;
    });
  }, [structures, filterCoordNums, filterAnions, filterSolvents]);

  // 重置筛选
  const resetFilters = () => {
    setFilterCoordNums([]);
    setFilterAnions([]);
    setFilterSolvents([]);
  };

  // 加载已完成的 MD Jobs
  const loadMdJobs = useCallback(async () => {
    try {
      const data = await getMDJobs();
      setMdJobs(data.filter(j => j.status === JobStatus.COMPLETED));
    } catch (err) {
      console.error('Failed to load MD jobs:', err);
    }
  }, []);

  // 加载单个 MD Job 详情
  const loadMdJobDetail = useCallback(async (mdJobId: number) => {
    try {
      const data = await getMDJob(mdJobId);
      setSelectedMdJob(data);
    } catch (err) {
      console.error('Failed to load MD job detail:', err);
    }
  }, []);

  // 加载溶剂化结构
  const loadStructures = useCallback(async (mdJobId: number) => {
    setLoading(true);
    try {
      const data = await getSolvationStructures(mdJobId);
      setStructures(data);
    } catch (err) {
      console.error('Failed to load structures:', err);
      message.error('加载溶剂化结构失败');
    } finally {
      setLoading(false);
    }
  }, []);

  // 加载现有任务详情
  const loadJob = useCallback(async () => {
    if (!id || isCreateMode) return;
    setLoading(true);
    try {
      const data = await getClusterAnalysisJob(Number(id));
      setJob(data);
    } catch (err) {
      console.error('Failed to load job:', err);
      message.error('加载任务详情失败');
    } finally {
      setLoading(false);
    }
  }, [id, isCreateMode]);

  useEffect(() => {
    if (isCreateMode) {
      loadMdJobs();
    } else {
      loadJob();
    }
  }, [isCreateMode, loadMdJobs, loadJob]);

  useEffect(() => {
    if (selectedMdJobId && isCreateMode) {
      loadMdJobDetail(selectedMdJobId);
      loadStructures(selectedMdJobId);
      setCurrentStep(1);
    }
  }, [selectedMdJobId, isCreateMode, loadMdJobDetail, loadStructures]);

  // 生成规划预览
  const handlePlan = async () => {
    if (!selectedMdJobId || selectedStructureIds.length === 0 || selectedCalcTypes.length === 0) {
      message.warning('请选择结构和计算类型');
      return;
    }
    setPlanLoading(true);
    try {
      const result = await planClusterAnalysis({
        md_job_id: selectedMdJobId,
        solvation_structure_ids: selectedStructureIds,
        calc_types: selectedCalcTypes,
      });
      setPlanResult(result);
      setCurrentStep(3);
    } catch (err) {
      console.error('Failed to plan:', err);
      message.error('生成规划失败');
    } finally {
      setPlanLoading(false);
    }
  };

  // 提交任务
  const handleSubmit = async () => {
    if (!selectedMdJobId || selectedStructureIds.length === 0 || selectedCalcTypes.length === 0) {
      return;
    }
    setSubmitLoading(true);
    try {
      const job = await submitClusterAnalysis({
        md_job_id: selectedMdJobId,
        solvation_structure_ids: selectedStructureIds,
        calc_types: selectedCalcTypes,
      });
      message.success('分析任务已提交');
      navigate(`/workspace/liquid-electrolyte/analysis/${job.id}`);
    } catch (err) {
      console.error('Failed to submit:', err);
      message.error('提交任务失败');
    } finally {
      setSubmitLoading(false);
    }
  };

  // 结构选择表格列 - 精简设计
  const structureColumns: ColumnsType<SolvationStructure> = [
    {
      title: '',
      width: 32,
      render: (_, record) => (
        <Checkbox
          checked={selectedStructureIds.includes(record.id)}
          onChange={(e) => {
            if (e.target.checked) {
              setSelectedStructureIds([...selectedStructureIds, record.id]);
            } else {
              setSelectedStructureIds(selectedStructureIds.filter(id => id !== record.id));
            }
          }}
        />
      ),
    },
    {
      title: '簇名称',
      key: 'cluster_name',
      width: 160,
      render: (_, record) => (
        <Space size={4}>
          <Text type="secondary" style={{ fontSize: 11 }}>#{record.id}</Text>
          <Text strong style={{ fontFamily: 'monospace', color: token.colorPrimary }}>
            {generateClusterName(record.center_ion, record.composition)}
          </Text>
        </Space>
      ),
    },
    {
      title: '组成',
      dataIndex: 'composition',
      render: (comp: Record<string, number>) => {
        if (!comp || Object.keys(comp).length === 0) return '-';
        // 紧凑显示：溶剂用绿色，阴离子用橙色
        const anionList = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'DFOB', 'BOB'];
        return (
          <span style={{ fontSize: 12 }}>
            {Object.entries(comp)
              .filter(([_, count]) => count > 0)
              .map(([mol, count], idx) => {
                const isAnion = anionList.some(a => mol.toUpperCase().includes(a));
                return (
                  <Tag
                    key={mol}
                    color={isAnion ? 'orange' : 'green'}
                    style={{ margin: '1px 2px', padding: '0 4px', fontSize: 11 }}
                  >
                    {mol}×{count}
                  </Tag>
                );
              })}
          </span>
        );
      },
    },
    {
      title: 'CN',
      dataIndex: 'coordination_num',
      width: 45,
      align: 'center' as const,
      render: (cn: number) => <Text strong>{cn ?? '-'}</Text>,
    },
  ];

  // 渲染创建模式
  const renderCreateMode = () => (
    <>
      {/* 步骤条 */}
      <Card style={{ marginBottom: 24 }}>
        <Steps
          current={currentStep}
          items={[
            { title: '选择数据源', description: '选择 MD Job' },
            { title: '选择结构', description: '筛选溶剂化结构' },
            { title: '选择计算', description: '选择计算类型' },
            { title: '确认提交', description: '预览并提交' },
          ]}
        />
      </Card>

      {/* Step 0: 选择 MD Job */}
      {currentStep === 0 && (
        <Card title="选择数据来源">
          <Space direction="vertical" style={{ width: '100%' }}>
            <Text>选择一个已完成的 MD 模拟任务作为分析的数据来源：</Text>
            <Select
              style={{ width: 400 }}
              placeholder="请选择 MD Job"
              value={selectedMdJobId}
              onChange={(v) => setSelectedMdJobId(v)}
              options={mdJobs.map(j => ({
                value: j.id,
                label: `#${j.id} - ${j.config?.job_name || 'MD Job'} (${dayjs(j.created_at).format('YYYY-MM-DD')})`,
              }))}
              showSearch
              filterOption={(input, option) =>
                (option?.label ?? '').toLowerCase().includes(input.toLowerCase())
              }
            />
            {mdJobs.length === 0 && (
              <Alert
                type="info"
                message="暂无已完成的 MD 任务"
                description={
                  <Space>
                    <Text>请先完成 MD 模拟任务</Text>
                    <Button type="link" onClick={() => navigate('/workspace/liquid-electrolyte/md')}>
                      前往 MD 模拟
                    </Button>
                  </Space>
                }
              />
            )}
          </Space>
        </Card>
      )}

      {/* Step 1: 选择结构 */}
      {currentStep >= 1 && selectedMdJobId && (
        <Card
          size="small"
          style={{ marginBottom: 24 }}
          title={
            <Row justify="space-between" align="middle" style={{ width: '100%' }}>
              <Col>
                <Space size="small">
                  <Text strong>选择溶剂化结构</Text>
                  <Tag color="blue">#{selectedMdJobId}</Tag>
                  {selectedMdJob?.config?.job_name && (
                    <Text type="secondary" style={{ fontSize: 12 }}>{selectedMdJob.config.job_name}</Text>
                  )}
                </Space>
              </Col>
              <Col>
                <Space size="small">
                  <Text type="secondary" style={{ fontSize: 12 }}>
                    已选 <Text strong style={{ color: token.colorPrimary }}>{selectedStructureIds.length}</Text> / {filteredStructures.length} 个
                  </Text>
                  <Button size="small" type="link" onClick={() => setSelectedStructureIds(filteredStructures.map(s => s.id))}>
                    全选
                  </Button>
                  <Button size="small" type="link" onClick={() => setSelectedStructureIds([])}>
                    清空
                  </Button>
                </Space>
              </Col>
            </Row>
          }
        >
          <Spin spinning={loading}>
            {structures.length === 0 ? (
              <Empty description="未找到溶剂化结构，请先在 MD 详情页提取结构" />
            ) : (
              <>
                {/* 筛选区域 - 紧凑内联 */}
                <div style={{ marginBottom: 12, display: 'flex', alignItems: 'center', gap: 16, flexWrap: 'wrap' }}>
                  <Space size={4}>
                    <Text type="secondary" style={{ fontSize: 12 }}>配位数:</Text>
                    <Select
                      mode="multiple"
                      size="small"
                      style={{ minWidth: 100 }}
                      placeholder="全部"
                      value={filterCoordNums}
                      onChange={setFilterCoordNums}
                      options={filterOptions.coordNums.map(n => ({ value: n, label: n }))}
                      allowClear
                      maxTagCount={1}
                    />
                  </Space>
                  <Space size={4}>
                    <Text type="secondary" style={{ fontSize: 12 }}>阴离子:</Text>
                    <Select
                      mode="multiple"
                      size="small"
                      style={{ minWidth: 100 }}
                      placeholder="全部"
                      value={filterAnions}
                      onChange={setFilterAnions}
                      options={filterOptions.anions.map(a => ({ value: a, label: a }))}
                      allowClear
                      maxTagCount={1}
                    />
                  </Space>
                  <Space size={4}>
                    <Text type="secondary" style={{ fontSize: 12 }}>溶剂:</Text>
                    <Select
                      mode="multiple"
                      size="small"
                      style={{ minWidth: 100 }}
                      placeholder="全部"
                      value={filterSolvents}
                      onChange={setFilterSolvents}
                      options={filterOptions.solvents.map(s => ({ value: s, label: s }))}
                      allowClear
                      maxTagCount={1}
                    />
                  </Space>
                  {(filterCoordNums.length > 0 || filterAnions.length > 0 || filterSolvents.length > 0) && (
                    <Button size="small" type="link" onClick={resetFilters} style={{ padding: 0 }}>
                      重置
                    </Button>
                  )}
                </div>

                <Table
                  columns={structureColumns}
                  dataSource={filteredStructures}
                  rowKey="id"
                  size="small"
                  pagination={{
                    size: 'small',
                    pageSize: 15,
                    showSizeChanger: false,
                    showTotal: (t) => `${t} 条`,
                    style: { marginBottom: 0 }
                  }}
                />
                <div style={{ marginTop: 12, textAlign: 'right' }}>
                  <Button
                    type="primary"
                    disabled={selectedStructureIds.length === 0}
                    onClick={() => setCurrentStep(2)}
                  >
                    下一步：选择计算类型
                  </Button>
                </div>
              </>
            )}
          </Spin>
        </Card>
      )}

      {/* Step 2: 选择计算类型 */}
      {currentStep >= 2 && (
        <Card title="选择计算类型" style={{ marginBottom: 24 }}>
          <Row gutter={[16, 16]}>
            {CALC_TYPE_OPTIONS.map(opt => {
              const isSelected = selectedCalcTypes.includes(opt.value);
              const info = CALC_TYPE_INFO[opt.value];
              return (
                <Col key={opt.value} xs={24} sm={12} md={8}>
                  <div
                    style={{
                      padding: 16,
                      border: `2px solid ${isSelected ? token.colorPrimary : token.colorBorder}`,
                      borderRadius: 8,
                      background: isSelected ? token.colorPrimaryBg : token.colorBgContainer,
                      cursor: 'pointer',
                      transition: 'all 0.2s',
                    }}
                    onClick={() => {
                      if (isSelected) {
                        setSelectedCalcTypes(selectedCalcTypes.filter(t => t !== opt.value));
                      } else {
                        setSelectedCalcTypes([...selectedCalcTypes, opt.value]);
                      }
                    }}
                  >
                    <Space direction="vertical" style={{ width: '100%' }}>
                      <Space>
                        <Checkbox checked={isSelected} />
                        <Text strong>{info.icon} {opt.label}</Text>
                        <Tag color={opt.riskLevel === 'high' ? 'red' : opt.riskLevel === 'medium' ? 'orange' : 'green'}>
                          {opt.riskLevel === 'high' ? '高风险' : opt.riskLevel === 'medium' ? '中风险' : '低风险'}
                        </Tag>
                      </Space>
                      <Text type="secondary" style={{ fontSize: 12 }}>{opt.description}</Text>
                    </Space>
                  </div>
                </Col>
              );
            })}
          </Row>
          {selectedCalcTypes.some(t => ['REDOX', 'REORGANIZATION'].includes(t)) && (
            <Alert
              type="warning"
              message="高风险计算提示"
              description="Redox 和 Reorganization 计算需要更多的 QC 任务，计算时间较长，结果误差可能较大。"
              style={{ marginTop: 16 }}
              showIcon
            />
          )}
          <div style={{ marginTop: 16, textAlign: 'right' }}>
            <Space>
              <Button onClick={() => setCurrentStep(1)}>上一步</Button>
              <Button
                type="primary"
                disabled={selectedCalcTypes.length === 0}
                loading={planLoading}
                onClick={handlePlan}
              >
                生成规划预览
              </Button>
            </Space>
          </div>
        </Card>
      )}

      {/* Step 3: 确认提交 */}
      {currentStep >= 3 && planResult && (
        <Card title="确认提交" style={{ marginBottom: 24 }}>
          <Descriptions bordered column={2} size="small">
            <Descriptions.Item label="MD Job">#{selectedMdJobId}</Descriptions.Item>
            <Descriptions.Item label="选中结构数">{planResult.selected_structures_count}</Descriptions.Item>
            <Descriptions.Item label="新建 QC 任务">
              <Text type="warning" strong>{planResult.total_new_qc_tasks}</Text>
            </Descriptions.Item>
            <Descriptions.Item label="复用 QC 任务">
              <Text type="success" strong>{planResult.total_reused_qc_tasks}</Text>
            </Descriptions.Item>
            <Descriptions.Item label="预估计算时间" span={2}>
              约 {planResult.estimated_compute_hours.toFixed(1)} 核时
            </Descriptions.Item>
          </Descriptions>

          <Divider>计算类型详情</Divider>

          {planResult.calc_requirements.map(req => (
            <Card key={req.calc_type} size="small" style={{ marginBottom: 8 }}>
              <Row justify="space-between" align="middle">
                <Col>
                  <Space>
                    <Text strong>{CALC_TYPE_INFO[req.calc_type as ClusterCalcType]?.icon}</Text>
                    <Text strong>{CALC_TYPE_INFO[req.calc_type as ClusterCalcType]?.label}</Text>
                  </Space>
                </Col>
                <Col>
                  <Space>
                    <Tag color="blue">新建 {req.new_tasks_count}</Tag>
                    <Tag color="green">复用 {req.reused_tasks_count}</Tag>
                  </Space>
                </Col>
              </Row>
            </Card>
          ))}

          {planResult.warnings.length > 0 && (
            <Alert
              type="warning"
              message="注意事项"
              description={
                <ul style={{ margin: 0, paddingLeft: 20 }}>
                  {planResult.warnings.map((w, i) => <li key={i}>{w}</li>)}
                </ul>
              }
              style={{ marginTop: 16 }}
            />
          )}

          <div style={{ marginTop: 24, textAlign: 'right' }}>
            <Space>
              <Button onClick={() => setCurrentStep(2)}>上一步</Button>
              <Button
                type="primary"
                icon={<SendOutlined />}
                loading={submitLoading}
                onClick={handleSubmit}
              >
                提交分析任务
              </Button>
            </Space>
          </div>
        </Card>
      )}
    </>
  );

  // 渲染查看模式
  const renderViewMode = () => {
    if (!job) return <Spin spinning={loading}><Empty description="加载中..." /></Spin>;

    const statusCfg = STATUS_CONFIG[job.status] || STATUS_CONFIG.CREATED;

    return (
      <>
        {/* 任务状态卡片 */}
        <Card style={{ marginBottom: 24 }}>
          <Row gutter={24} align="middle">
            <Col flex="auto">
              <Space direction="vertical">
                <Space>
                  <Title level={4} style={{ margin: 0 }}>分析任务 #{job.id}</Title>
                  <Tag icon={statusCfg.icon} color={statusCfg.color}>{statusCfg.text}</Tag>
                </Space>
                <Space>
                  <Text type="secondary">MD Job:</Text>
                  <Button
                    type="link"
                    size="small"
                    icon={<RocketOutlined />}
                    onClick={() => navigate(`/workspace/liquid-electrolyte/md/${job.md_job_id}`)}
                  >
                    #{job.md_job_id}
                  </Button>
                  <Divider type="vertical" />
                  <Text type="secondary">创建时间: {dayjs(job.created_at).format('YYYY-MM-DD HH:mm')}</Text>
                </Space>
              </Space>
            </Col>
            <Col>
              <Row gutter={16}>
                <Col>
                  <Statistic title="选中结构" value={job.selected_structures?.count || 0} />
                </Col>
                <Col>
                  <Statistic
                    title="QC 进度"
                    value={job.qc_task_plan?.completed_qc_tasks || 0}
                    suffix={`/ ${job.qc_task_plan?.total_qc_tasks || 0}`}
                  />
                </Col>
                <Col>
                  <Statistic
                    title="总进度"
                    value={job.progress}
                    suffix="%"
                    valueStyle={{ color: job.progress === 100 ? '#52c41a' : undefined }}
                  />
                </Col>
              </Row>
            </Col>
          </Row>

          {job.error_message && (
            <Alert
              type="error"
              message="错误信息"
              description={job.error_message}
              style={{ marginTop: 16 }}
            />
          )}
        </Card>

        {/* 计算类型 */}
        <Card title="计算类型" size="small" style={{ marginBottom: 24 }}>
          <Space wrap>
            {job.calc_types.map(ct => {
              const info = CALC_TYPE_INFO[ct as ClusterCalcType];
              return (
                <Tag key={ct} color={info?.riskLevel === 'high' ? 'red' : info?.riskLevel === 'medium' ? 'orange' : 'blue'}>
                  {info?.icon} {info?.label || ct}
                </Tag>
              );
            })}
          </Space>
        </Card>

        {/* 结果展示 */}
        {job.status === 'COMPLETED' && (
          <ClusterAnalysisResultsPanel jobId={job.id} />
        )}

        {/* 等待中提示 */}
        {['SUBMITTED', 'RUNNING', 'WAITING_QC', 'CALCULATING'].includes(job.status) && (
          <Card>
            <div style={{ textAlign: 'center', padding: 40 }}>
              <Spin size="large" />
              <Title level={4} style={{ marginTop: 16 }}>任务运行中...</Title>
              <Text type="secondary">QC 任务进度: {job.qc_task_plan?.completed_qc_tasks || 0} / {job.qc_task_plan?.total_qc_tasks || 0}</Text>
              <div style={{ marginTop: 16 }}>
                <Progress
                  percent={job.progress}
                  status="active"
                  strokeColor={{ '0%': '#108ee9', '100%': '#87d068' }}
                />
              </div>
            </div>
          </Card>
        )}
      </>
    );
  };

  return (
    <div style={{ padding: 24 }}>
      {/* 返回按钮 */}
      <div style={{ marginBottom: 16 }}>
        <Button
          icon={<ArrowLeftOutlined />}
          onClick={() => navigate('/workspace/liquid-electrolyte/analysis')}
        >
          返回列表
        </Button>
      </div>

      {/* 页面标题 */}
      <Title level={4} style={{ marginBottom: 24 }}>
        {isCreateMode ? (
          <>
            <ExperimentOutlined style={{ marginRight: 8 }} />
            新建后处理分析
          </>
        ) : (
          <>
            <ExperimentOutlined style={{ marginRight: 8 }} />
            分析任务详情
          </>
        )}
      </Title>

      {isCreateMode ? renderCreateMode() : renderViewMode()}
    </div>
  );
}