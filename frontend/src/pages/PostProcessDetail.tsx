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
  BulbOutlined,
  AppstoreOutlined,
  UnorderedListOutlined,
  CalculatorOutlined,
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
import { getMDJobs, getMDJob, getSolvationStructures, autoSelectSolvationStructures, type SolvationStructure, type AutoSelectResponse } from '../api/jobs';
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

  // 视图模式: list = 列表视图, grouped = 分组视图
  const [viewMode, setViewMode] = useState<'list' | 'grouped'>('list');
  // 智能选择加载状态
  const [autoSelectLoading, setAutoSelectLoading] = useState(false);
  // 展开的分组
  const [expandedGroups, setExpandedGroups] = useState<string[]>([]);

  // 数字转下标
  const toSubscript = (num: number): string => {
    const subscripts = '₀₁₂₃₄₅₆₇₈₉';
    return String(num).split('').map(d => subscripts[parseInt(d)]).join('');
  };

  // 从 composition 生成簇名称
  // 旧格式: Na1MPN4PF6 → 新格式: Na⁺·MPN·(PF6)₄
  const generateClusterName = (centerIon: string, composition: Record<string, number>): string => {
    if (!composition || Object.keys(composition).length === 0) return centerIon ? `${centerIon}⁺` : '-';

    // 阴离子列表（需要括号包裹）
    const anionList = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'DFOB', 'BOB', 'NO3'];

    const parts = Object.entries(composition)
      .filter(([_, count]) => count > 0)
      .sort(([a], [b]) => a.localeCompare(b))
      .map(([mol, count]) => {
        const isAnion = anionList.some(a => mol.toUpperCase().includes(a.toUpperCase()));
        if (count === 1) {
          return mol;
        } else if (isAnion) {
          // 阴离子用括号: (PF6)₄
          return `(${mol})${toSubscript(count)}`;
        } else {
          // 溶剂直接加下标: MPN₂
          return `${mol}${toSubscript(count)}`;
        }
      });

    return `${centerIon || ''}⁺·${parts.join('·')}`;
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

  // 按组成分组的结构
  const groupedStructures = useMemo(() => {
    const groups: Record<string, { structures: SolvationStructure[]; count: number; percentage: number }> = {};
    const total = structures.length;

    structures.forEach(s => {
      const key = generateClusterName(s.center_ion, s.composition);
      if (!groups[key]) {
        groups[key] = { structures: [], count: 0, percentage: 0 };
      }
      groups[key].structures.push(s);
      groups[key].count++;
    });

    // 计算百分比
    Object.keys(groups).forEach(key => {
      groups[key].percentage = total > 0 ? (groups[key].count / total) * 100 : 0;
    });

    return groups;
  }, [structures]);

  // 分组列表（按数量排序）
  const sortedGroupKeys = useMemo(() => {
    return Object.keys(groupedStructures).sort(
      (a, b) => groupedStructures[b].count - groupedStructures[a].count
    );
  }, [groupedStructures]);

  // 获取结构的组成键
  const getStructureGroupKey = useCallback((s: SolvationStructure) => {
    return generateClusterName(s.center_ion, s.composition);
  }, []);

  // 智能选择：每种组成选1个
  const handleAutoSelect = useCallback(async () => {
    if (!selectedMdJobId) return;
    setAutoSelectLoading(true);
    try {
      const result = await autoSelectSolvationStructures(selectedMdJobId);
      const ids = result.selected_structures.map(s => s.id);
      setSelectedStructureIds(ids);
      message.success(`已智能选择 ${ids.length} 个结构（覆盖 ${result.unique_compositions} 种组成）`);
    } catch (err) {
      console.error('Auto select failed:', err);
      message.error('智能选择失败');
    } finally {
      setAutoSelectLoading(false);
    }
  }, [selectedMdJobId]);

  // 每种组成选N个
  const selectNPerGroup = useCallback((n: number) => {
    const ids: number[] = [];
    Object.values(groupedStructures).forEach(group => {
      const selected = group.structures.slice(0, n).map(s => s.id);
      ids.push(...selected);
    });
    setSelectedStructureIds(ids);
    message.success(`已选择 ${ids.length} 个结构（每种组成最多 ${n} 个）`);
  }, [groupedStructures]);

  // 选择某个分组的全部结构
  const selectGroup = useCallback((groupKey: string, select: boolean) => {
    const groupStructures = groupedStructures[groupKey]?.structures || [];
    const groupIds = groupStructures.map(s => s.id);
    if (select) {
      setSelectedStructureIds(prev => [...new Set([...prev, ...groupIds])]);
    } else {
      setSelectedStructureIds(prev => prev.filter(id => !groupIds.includes(id)));
    }
  }, [groupedStructures]);

  // 检查分组是否全选
  const isGroupAllSelected = useCallback((groupKey: string) => {
    const groupStructures = groupedStructures[groupKey]?.structures || [];
    return groupStructures.length > 0 && groupStructures.every(s => selectedStructureIds.includes(s.id));
  }, [groupedStructures, selectedStructureIds]);

  // 检查分组是否部分选中
  const isGroupPartiallySelected = useCallback((groupKey: string) => {
    const groupStructures = groupedStructures[groupKey]?.structures || [];
    const selectedCount = groupStructures.filter(s => selectedStructureIds.includes(s.id)).length;
    return selectedCount > 0 && selectedCount < groupStructures.length;
  }, [groupedStructures, selectedStructureIds]);

  // 计算预估 QC 任务数
  const estimatedQCTasks = useMemo(() => {
    const numStructures = selectedStructureIds.length;
    if (numStructures === 0 || selectedCalcTypes.length === 0) {
      return { total: 0, details: {} as Record<string, number> };
    }

    const details: Record<string, number> = {};
    let total = 0;

    // 估算每种计算类型需要的 QC 任务数
    selectedCalcTypes.forEach(calcType => {
      let count = 0;
      switch (calcType) {
        case 'BINDING_TOTAL':
        case 'DESOLVATION_FULL':
          // cluster + ion + 每种配体 ≈ 2-5 个任务/结构
          count = numStructures * 3;
          break;
        case 'BINDING_PAIRWISE':
          // 每个配体一个 dimer + ligand ≈ 2*配体数/结构
          count = numStructures * 6;
          break;
        case 'DESOLVATION_STEPWISE':
          // cluster + 每个配体的 (cluster-i + ligand) ≈ 配位数*2+1
          count = numStructures * 12;
          break;
        case 'REDOX':
          // gas优化 + freq + solvent = 3 * 2状态 = 6/结构
          count = numStructures * 6;
          break;
        case 'REORGANIZATION':
          // 2个几何 * 4个能量 = 8/结构
          count = numStructures * 8;
          break;
      }
      details[calcType] = count;
      total += count;
    });

    return { total, details };
  }, [selectedStructureIds, selectedCalcTypes]);

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

  // 结构选择表格列
  const structureColumns: ColumnsType<SolvationStructure> = [
    {
      title: '',
      width: 40,
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
      title: '#',
      dataIndex: 'id',
      width: 60,
      render: (id: number) => <Text type="secondary">#{id}</Text>,
    },
    {
      title: '簇名称',
      key: 'cluster_name',
      width: 150,
      render: (_, record) => (
        <Text strong style={{ fontFamily: 'monospace' }}>
          {generateClusterName(record.center_ion, record.composition)}
        </Text>
      ),
    },
    {
      title: '组成',
      dataIndex: 'composition',
      width: 200,
      render: (comp: Record<string, number>) => {
        if (!comp || Object.keys(comp).length === 0) return '-';
        return (
          <Space size={4} wrap>
            {Object.entries(comp)
              .filter(([_, count]) => count > 0)
              .map(([mol, count]) => (
                <Tag key={mol} color="cyan" style={{ margin: 0 }}>{mol}×{count}</Tag>
              ))}
          </Space>
        );
      },
    },
    {
      title: 'CN',
      dataIndex: 'coordination_num',
      width: 55,
      align: 'center' as const,
      render: (cn: number) => <Tag color="blue">{cn ?? '-'}</Tag>,
    },
    {
      title: '帧号',
      dataIndex: 'snapshot_frame',
      width: 70,
      align: 'center' as const,
      render: (frame: number) => <Text type="secondary">{frame ?? '-'}</Text>,
    },
    {
      title: '占比',
      key: 'percentage',
      width: 80,
      align: 'center' as const,
      render: (_, record) => {
        const groupKey = getStructureGroupKey(record);
        const group = groupedStructures[groupKey];
        const pct = group?.percentage || 0;
        return (
          <Tooltip title={`该组成共 ${group?.count || 0} 个`}>
            <Text type="secondary">{pct.toFixed(1)}%</Text>
          </Tooltip>
        );
      },
    },
    {
      title: '组内',
      key: 'group_count',
      width: 60,
      align: 'center' as const,
      render: (_, record) => {
        const groupKey = getStructureGroupKey(record);
        const group = groupedStructures[groupKey];
        return <Tag>{group?.count || 0}</Tag>;
      },
    },
  ];

  // 计算当前步骤的有效性
  const canProceedToStep2 = selectedMdJobId && selectedStructureIds.length > 0 && selectedCalcTypes.length > 0;
  const canSubmit = canProceedToStep2 && planResult;

  // 计算类型的详细说明和示意图
  const CALC_TYPE_DETAILS: Record<ClusterCalcType, { diagram: string; explanation: string }> = {
    BINDING_TOTAL: {
      diagram: 'Li⁺ + 溶剂₁ + 溶剂₂ + ... → [Li·溶剂化簇]',
      explanation: '计算整个溶剂化簇的形成能，反映离子与所有配体的总结合强度',
    },
    BINDING_PAIRWISE: {
      diagram: 'Li⁺ + 单个分子 → [Li-分子]',
      explanation: '分别计算每个配体与离子的结合能，比较不同分子的配位能力',
    },
    DESOLVATION_STEPWISE: {
      diagram: '[Li·ABCD] → [Li·ABC] + D → [Li·AB] + C → ...',
      explanation: '逐步移除配体，计算每一步的能量变化，分析脱溶剂化路径',
    },
    DESOLVATION_FULL: {
      diagram: '[Li·溶剂化簇] → Li⁺(裸离子) + 所有配体',
      explanation: '计算完全脱去溶剂化壳层需要的总能量',
    },
    REDOX: {
      diagram: 'Li⁺/Li⁰ 或 分子⁺/分子⁰ 氧化还原电位',
      explanation: '使用热力学循环计算电化学窗口和稳定性',
    },
    REORGANIZATION: {
      diagram: 'λ = E(R₁,Q₂) + E(R₂,Q₁) - E(R₁,Q₁) - E(R₂,Q₂)',
      explanation: 'Marcus理论4点法计算电子转移重组能',
    },
  };

  // 渲染创建模式
  const renderCreateMode = () => {
    // Step 0: 选择数据源 - 简化版
    if (currentStep === 0) {
      return (
        <Card title="选择数据来源" style={{ marginBottom: 24 }}>
          <Row gutter={16} align="middle">
            <Col flex="auto">
              <Select
                style={{ width: '100%', maxWidth: 500 }}
                placeholder="请选择一个已完成的 MD 模拟任务"
                value={selectedMdJobId}
                onChange={(v) => {
                  setSelectedMdJobId(v);
                  setSelectedStructureIds([]);
                  setSelectedCalcTypes([]);
                  setPlanResult(null);
                }}
                options={mdJobs.map(j => ({
                  value: j.id,
                  label: `#${j.id} - ${j.config?.job_name || 'MD Job'} (${dayjs(j.created_at).format('YYYY-MM-DD')})`,
                }))}
                showSearch
                filterOption={(input, option) =>
                  (option?.label ?? '').toLowerCase().includes(input.toLowerCase())
                }
              />
            </Col>
            <Col>
              <Button
                type="primary"
                disabled={!selectedMdJobId}
                onClick={() => setCurrentStep(1)}
              >
                下一步 →
              </Button>
            </Col>
          </Row>
          {mdJobs.length === 0 && (
            <Alert
              type="info"
              message="暂无已完成的 MD 任务，请先完成 MD 模拟"
              style={{ marginTop: 16 }}
              action={
                <Button size="small" onClick={() => navigate('/workspace/liquid-electrolyte/md')}>
                  前往 MD 模拟
                </Button>
              }
            />
          )}
        </Card>
      );
    }

    // Step 1: 配置分析（结构 + 计算类型 + 预览，三栏布局）
    if (currentStep === 1) {
      const coveredCompositions = Object.keys(groupedStructures).filter(k =>
        groupedStructures[k].structures.some(s => selectedStructureIds.includes(s.id))
      ).length;

      return (
        <div style={{ display: 'flex', flexDirection: 'column', height: 'calc(100vh - 180px)' }}>
          {/* 顶部进度条 */}
          <div style={{
            background: token.colorBgContainer,
            borderRadius: 12,
            padding: '16px 24px',
            marginBottom: 16,
            boxShadow: `0 2px 8px ${token.colorBgSpotlight}`,
          }}>
            <Row align="middle" justify="space-between">
              <Col>
                <Space size="large">
                  <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                    <div style={{
                      width: 32, height: 32, borderRadius: '50%',
                      background: token.colorPrimary, color: '#fff',
                      display: 'flex', alignItems: 'center', justifyContent: 'center',
                      fontWeight: 'bold',
                    }}>1</div>
                    <Text strong>配置分析</Text>
                  </div>
                  <div style={{ width: 60, height: 2, background: selectedStructureIds.length > 0 && selectedCalcTypes.length > 0 ? token.colorPrimary : token.colorBorder }} />
                  <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                    <div style={{
                      width: 32, height: 32, borderRadius: '50%',
                      background: canProceedToStep2 ? token.colorPrimary : token.colorBorder,
                      color: canProceedToStep2 ? '#fff' : token.colorTextDisabled,
                      display: 'flex', alignItems: 'center', justifyContent: 'center',
                      fontWeight: 'bold',
                    }}>2</div>
                    <Text type={canProceedToStep2 ? undefined : 'secondary'}>确认提交</Text>
                  </div>
                </Space>
              </Col>
              <Col>
                <Space>
                  <Tag color="blue">{selectedMdJob?.config?.job_name || `MD #${selectedMdJobId}`}</Tag>
                  <Button size="small" onClick={() => setCurrentStep(0)}>更换数据源</Button>
                </Space>
              </Col>
            </Row>
          </div>

          {/* 主内容区 - 三栏布局 */}
          <Row gutter={16} style={{ flex: 1, minHeight: 0 }}>
            {/* 左栏：结构选择 */}
            <Col span={10} style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
              <Card
                size="small"
                title={<Space><ExperimentOutlined /> 溶剂化结构</Space>}
                extra={
                  <Space.Compact size="small">
                    <Button
                      icon={<UnorderedListOutlined />}
                      type={viewMode === 'list' ? 'primary' : 'default'}
                      onClick={() => setViewMode('list')}
                    />
                    <Button
                      icon={<AppstoreOutlined />}
                      type={viewMode === 'grouped' ? 'primary' : 'default'}
                      onClick={() => setViewMode('grouped')}
                    />
                  </Space.Compact>
                }
                bodyStyle={{ padding: 12, display: 'flex', flexDirection: 'column', flex: 1, minHeight: 0, overflow: 'hidden' }}
                style={{ height: '100%', display: 'flex', flexDirection: 'column' }}
              >
                <Spin spinning={loading} style={{ flex: 1, display: 'flex', flexDirection: 'column' }}>
                  {/* 快速选择 */}
                  <Space wrap style={{ marginBottom: 12 }}>
                    <Button size="small" icon={<BulbOutlined />} loading={autoSelectLoading} onClick={handleAutoSelect}>
                      智能选择
                    </Button>
                    <Button size="small" onClick={() => selectNPerGroup(3)}>每种3个</Button>
                    <Button size="small" type={filteredStructures.length !== structures.length ? 'primary' : 'default'}
                      onClick={() => setSelectedStructureIds(filteredStructures.map(s => s.id))}>
                      全选 ({filteredStructures.length})
                    </Button>
                    <Button size="small" danger onClick={() => setSelectedStructureIds([])}>清空</Button>
                  </Space>

                  {/* 筛选器 */}
                  <Space wrap style={{ marginBottom: 12 }}>
                    <Select size="small" mode="multiple" style={{ minWidth: 80 }} placeholder="CN"
                      value={filterCoordNums} onChange={setFilterCoordNums} allowClear maxTagCount={1}
                      options={filterOptions.coordNums.map(n => ({ value: n, label: `CN=${n}` }))} />
                    <Select size="small" mode="multiple" style={{ minWidth: 100 }} placeholder="阴离子"
                      value={filterAnions} onChange={setFilterAnions} allowClear maxTagCount={1}
                      options={filterOptions.anions.map(a => ({ value: a, label: a }))} />
                    <Select size="small" mode="multiple" style={{ minWidth: 100 }} placeholder="溶剂"
                      value={filterSolvents} onChange={setFilterSolvents} allowClear maxTagCount={1}
                      options={filterOptions.solvents.map(s => ({ value: s, label: s }))} />
                    {(filterCoordNums.length > 0 || filterAnions.length > 0 || filterSolvents.length > 0) && (
                      <Button size="small" type="link" onClick={resetFilters}>重置</Button>
                    )}
                  </Space>

                  {/* 结构列表 - 带滚动 */}
                  <div style={{ flex: 1, overflow: 'auto', minHeight: 200 }}>
                    {viewMode === 'grouped' ? (
                      <div style={{ paddingRight: 4 }}>
                        {sortedGroupKeys.filter(key => groupedStructures[key].structures.some(s => filteredStructures.includes(s))).map(groupKey => {
                          const group = groupedStructures[groupKey];
                          const isAllSelected = isGroupAllSelected(groupKey);
                          const isPartial = isGroupPartiallySelected(groupKey);
                          const selectedInGroup = group.structures.filter(s => selectedStructureIds.includes(s.id)).length;
                          return (
                            <div key={groupKey} style={{
                              padding: '8px 12px', marginBottom: 4, borderRadius: 8,
                              background: isAllSelected ? token.colorPrimaryBg : (isPartial ? token.colorWarningBg : token.colorBgLayout),
                              border: `1px solid ${isAllSelected ? token.colorPrimary : token.colorBorder}`,
                              cursor: 'pointer',
                            }}
                              onClick={() => selectGroup(groupKey, !isAllSelected)}
                            >
                              <Row justify="space-between" align="middle">
                                <Col>
                                  <Space>
                                    <Checkbox checked={isAllSelected} indeterminate={isPartial} />
                                    <Text strong style={{ fontFamily: 'monospace', fontSize: 13 }}>{groupKey}</Text>
                                  </Space>
                                </Col>
                                <Col>
                                  <Space size={4}>
                                    <Tag color="blue" style={{ margin: 0 }}>{group.count}</Tag>
                                    <Text type="secondary" style={{ fontSize: 12 }}>{group.percentage.toFixed(1)}%</Text>
                                    {selectedInGroup > 0 && <Tag color="green" style={{ margin: 0 }}>{selectedInGroup}✓</Tag>}
                                  </Space>
                                </Col>
                              </Row>
                            </div>
                          );
                        })}
                      </div>
                    ) : (
                      <Table
                        columns={structureColumns}
                        dataSource={filteredStructures}
                        rowKey="id"
                        size="small"
                        pagination={{ pageSize: 8, size: 'small', showTotal: (t) => `共${t}条`, showSizeChanger: true, pageSizeOptions: ['8', '15', '30', '50'] }}
                        scroll={{ x: 600, y: 'calc(100vh - 500px)' }}
                      />
                    )}
                  </div>
                </Spin>
              </Card>
            </Col>

            {/* 中栏：计算类型选择 */}
            <Col span={8} style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
              <Card
                size="small"
                title={<Space><ThunderboltOutlined /> 计算类型 <Text type="secondary" style={{ fontSize: 12 }}>点击查看详情</Text></Space>}
                bodyStyle={{ padding: 12, overflow: 'auto' }}
                style={{ height: '100%' }}
              >
                <Space direction="vertical" style={{ width: '100%' }} size={8}>
                  {CALC_TYPE_OPTIONS.map(opt => {
                    const isSelected = selectedCalcTypes.includes(opt.value);
                    const info = CALC_TYPE_INFO[opt.value];
                    const details = CALC_TYPE_DETAILS[opt.value];
                    return (
                      <div
                        key={opt.value}
                        style={{
                          padding: '10px 14px',
                          borderRadius: 8,
                          border: `2px solid ${isSelected ? token.colorPrimary : token.colorBorder}`,
                          background: isSelected ? token.colorPrimaryBg : 'transparent',
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
                        <Row justify="space-between" align="middle">
                          <Col>
                            <Space>
                              <Checkbox checked={isSelected} />
                              <Text strong>{info.icon} {opt.label}</Text>
                            </Space>
                          </Col>
                          <Col>
                            <Tag color={opt.riskLevel === 'high' ? 'red' : opt.riskLevel === 'medium' ? 'orange' : 'green'} style={{ margin: 0 }}>
                              {opt.riskLevel === 'high' ? '高' : opt.riskLevel === 'medium' ? '中' : '低'}
                            </Tag>
                          </Col>
                        </Row>
                        {/* 示意图 */}
                        <div style={{
                          marginTop: 6,
                          marginLeft: 24,
                          padding: '4px 8px',
                          background: token.colorBgLayout,
                          borderRadius: 4,
                          fontFamily: 'monospace',
                          fontSize: 11,
                          color: token.colorTextSecondary,
                        }}>
                          {details.diagram}
                        </div>
                        {/* 解释 */}
                        <Text type="secondary" style={{ fontSize: 11, marginLeft: 24, display: 'block', marginTop: 4 }}>
                          {details.explanation}
                        </Text>
                      </div>
                    );
                  })}
                </Space>

                {selectedCalcTypes.some(t => ['REDOX', 'REORGANIZATION'].includes(t)) && (
                  <Alert
                    type="warning"
                    message="高风险计算需要更多 QC 任务"
                    style={{ marginTop: 12 }}
                    showIcon
                  />
                )}
              </Card>
            </Col>

            {/* 右栏：实时预览 */}
            <Col span={6} style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
              <Card
                size="small"
                title={<Space><CalculatorOutlined /> 分析预览</Space>}
                bodyStyle={{ padding: 16 }}
                style={{ height: '100%' }}
              >
                {/* 选择统计 */}
                <div style={{
                  background: token.colorPrimaryBg,
                  borderRadius: 8,
                  padding: 16,
                  marginBottom: 16,
                  textAlign: 'center',
                }}>
                  <Statistic
                    title="已选结构"
                    value={selectedStructureIds.length}
                    suffix={<Text type="secondary">/ {structures.length}</Text>}
                  />
                  <Progress
                    percent={structures.length > 0 ? Math.round((selectedStructureIds.length / structures.length) * 100) : 0}
                    size="small"
                    style={{ marginTop: 8 }}
                  />
                </div>

                <div style={{
                  background: token.colorBgLayout,
                  borderRadius: 8,
                  padding: 16,
                  marginBottom: 16,
                  textAlign: 'center',
                }}>
                  <Statistic
                    title="覆盖组成"
                    value={coveredCompositions}
                    suffix={<Text type="secondary">/ {sortedGroupKeys.length} 种</Text>}
                  />
                </div>

                <Divider style={{ margin: '12px 0' }}>预估 QC 任务</Divider>

                {selectedCalcTypes.length > 0 && selectedStructureIds.length > 0 ? (
                  <>
                    {Object.entries(estimatedQCTasks.details).map(([calcType, count]) => (
                      <Row key={calcType} justify="space-between" style={{ marginBottom: 8 }}>
                        <Col><Text>{CALC_TYPE_INFO[calcType as ClusterCalcType]?.icon} {CALC_TYPE_INFO[calcType as ClusterCalcType]?.label}</Text></Col>
                        <Col><Text strong>~{count}</Text></Col>
                      </Row>
                    ))}
                    <Divider style={{ margin: '8px 0' }} />
                    <Row justify="space-between">
                      <Col><Text strong>总计</Text></Col>
                      <Col><Text strong style={{ color: token.colorPrimary, fontSize: 18 }}>~{estimatedQCTasks.total}</Text></Col>
                    </Row>
                  </>
                ) : (
                  <Empty description="选择结构和计算类型后显示" image={Empty.PRESENTED_IMAGE_SIMPLE} />
                )}
              </Card>
            </Col>
          </Row>

          {/* 底部固定操作栏 */}
          <div style={{
            background: token.colorBgContainer,
            borderRadius: 12,
            padding: '16px 24px',
            marginTop: 16,
            boxShadow: `0 -2px 8px ${token.colorBgSpotlight}`,
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
          }}>
            <Space>
              <Text type="secondary">
                {selectedStructureIds.length > 0 ? `已选 ${selectedStructureIds.length} 个结构` : '请选择结构'}
                {selectedCalcTypes.length > 0 ? ` · ${selectedCalcTypes.length} 种计算` : ''}
              </Text>
            </Space>
            <Space>
              <Button onClick={() => setCurrentStep(0)}>← 返回</Button>
              <Button
                type="primary"
                size="large"
                disabled={!canProceedToStep2}
                loading={planLoading}
                onClick={handlePlan}
              >
                生成规划 · 下一步 →
              </Button>
            </Space>
          </div>
        </div>
      );
    }

    // Step 2: 确认提交
    if (currentStep === 2 || currentStep === 3) {
      return (
        <div style={{
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          padding: 24,
        }}>
          {/* 进度条 */}
          <div style={{
            background: token.colorBgContainer,
            borderRadius: 12,
            padding: '16px 24px',
            marginBottom: 24,
            width: '100%',
            maxWidth: 800,
            boxShadow: `0 2px 8px ${token.colorBgSpotlight}`,
          }}>
            <Row align="middle" justify="center">
              <Space size="large">
                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                  <div style={{
                    width: 32, height: 32, borderRadius: '50%',
                    background: token.colorSuccess, color: '#fff',
                    display: 'flex', alignItems: 'center', justifyContent: 'center',
                  }}><CheckCircleOutlined /></div>
                  <Text type="secondary">配置分析</Text>
                </div>
                <div style={{ width: 60, height: 2, background: token.colorPrimary }} />
                <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                  <div style={{
                    width: 32, height: 32, borderRadius: '50%',
                    background: token.colorPrimary, color: '#fff',
                    display: 'flex', alignItems: 'center', justifyContent: 'center',
                    fontWeight: 'bold',
                  }}>2</div>
                  <Text strong>确认提交</Text>
                </div>
              </Space>
            </Row>
          </div>

          {/* 规划结果卡片 */}
          <Card
            style={{
              width: '100%',
              maxWidth: 800,
              borderRadius: 16,
              boxShadow: `0 8px 32px ${token.colorBgSpotlight}`,
            }}
          >
            {planResult ? (
              <>
                <Row gutter={24} style={{ marginBottom: 24 }}>
                  <Col span={6}>
                    <Statistic title="选中结构" value={planResult.selected_structures_count} suffix="个" />
                  </Col>
                  <Col span={6}>
                    <Statistic title="新建 QC 任务" value={planResult.total_new_qc_tasks} valueStyle={{ color: token.colorWarning }} />
                  </Col>
                  <Col span={6}>
                    <Statistic title="复用 QC 任务" value={planResult.total_reused_qc_tasks} valueStyle={{ color: token.colorSuccess }} />
                  </Col>
                  <Col span={6}>
                    <Statistic title="预估计算时间" value={planResult.estimated_compute_hours.toFixed(1)} suffix="核时" />
                  </Col>
                </Row>

                <Divider>计算类型详情</Divider>

                <Row gutter={[16, 16]}>
                  {planResult.calc_requirements.map(req => (
                    <Col key={req.calc_type} span={8}>
                      <Card size="small" style={{ textAlign: 'center' }}>
                        <Text strong style={{ fontSize: 16 }}>
                          {CALC_TYPE_INFO[req.calc_type as ClusterCalcType]?.icon}
                          {' '}{CALC_TYPE_INFO[req.calc_type as ClusterCalcType]?.label}
                        </Text>
                        <div style={{ marginTop: 8 }}>
                          <Tag color="blue">新建 {req.new_tasks_count}</Tag>
                          <Tag color="green">复用 {req.reused_tasks_count}</Tag>
                        </div>
                      </Card>
                    </Col>
                  ))}
                </Row>

                {planResult.warnings.length > 0 && (
                  <Alert
                    type="warning"
                    message="注意事项"
                    description={
                      <ul style={{ margin: 0, paddingLeft: 20 }}>
                        {planResult.warnings.map((w, i) => <li key={i}>{w}</li>)}
                      </ul>
                    }
                    style={{ marginTop: 24 }}
                  />
                )}

                <div style={{ marginTop: 32, textAlign: 'center' }}>
                  <Space size="large">
                    <Button size="large" onClick={() => setCurrentStep(1)}>← 返回修改</Button>
                    <Button
                      type="primary"
                      size="large"
                      icon={<SendOutlined />}
                      loading={submitLoading}
                      onClick={handleSubmit}
                      style={{ paddingLeft: 32, paddingRight: 32, height: 48 }}
                    >
                      提交分析任务
                    </Button>
                  </Space>
                </div>
              </>
            ) : (
              <div style={{ textAlign: 'center', padding: 48 }}>
                <Spin size="large" />
                <div style={{ marginTop: 16 }}>
                  <Text type="secondary">正在生成规划预览...</Text>
                </div>
              </div>
            )}
          </Card>
        </div>
      );
    }

    return null;
  };

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