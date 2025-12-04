/**
 * 批量去溶剂化能计算面板
 *
 * 功能：
 * 1. 自动挑选不同配位组成的溶剂化结构
 * 2. 多维度筛选：配位数、阴离子数、溶剂类型
 * 3. 批量提交计算任务
 * 4. 显示任务进度和结果
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
  Collapse,
  Row,
  Col,
  InputNumber,
  Typography,
  message,
  Tooltip,
  Badge,
  Spin,
  Empty,
  Divider,
  theme,
  Alert,
} from 'antd';
import {
  ThunderboltOutlined,
  ReloadOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  ExclamationCircleOutlined,
  BulbOutlined,
  FilterOutlined,
  PlusOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import { autoSelectSolvationStructures, type AutoSelectedStructure } from '../api/jobs';
import {
  batchCreateDesolvationJobs,
  getDesolvationOverview,
  getDesolvationQCTasks,
  type QCTaskInfo,
  type DesolvationQCTasksResponse,
} from '../api/desolvation';
import type {
  DesolvationJobResponse,
  DesolvationOverviewResponse,
  SolventModel,
  SolventConfig,
} from '../types/desolvation';
import { useThemeStore } from '../stores/themeStore';
import DesolvationResultView from './DesolvationResultView';

const { Text } = Typography;

// 常用阴离子模式
const ANION_PATTERNS = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'NO3', 'OTf', 'BOB', 'Cl', 'Br', 'I'];

// 隐式溶剂选项（按介电常数分组）
const SOLVENT_OPTIONS = [
  {
    label: '📌 电池电解液常用',
    options: [
      { value: 'DiMethylCarbonate', label: 'DMC 碳酸二甲酯 ε=3.1' },
      { value: 'EthyleneCarbonate', label: 'EC 碳酸乙烯酯 ε=89.8' },
      { value: 'PropyleneCarbonate', label: 'PC 碳酸丙烯酯 ε=64.9' },
      { value: 'TetraHydroFuran', label: 'THF 四氢呋喃 ε=7.4' },
    ],
  },
  {
    label: '📌 高介电常数 (ε>40)',
    options: [
      { value: 'Water', label: '水 (Water) ε=78.4' },
      { value: 'DiMethylSulfoxide', label: 'DMSO ε=46.8' },
      { value: '1,2-EthaneDiol', label: '乙二醇 ε=40.2' },
    ],
  },
  {
    label: '📌 中等介电常数 (ε=15-40)',
    options: [
      { value: 'Acetonitrile', label: '乙腈 ε=35.7' },
      { value: 'Methanol', label: '甲醇 ε=32.6' },
      { value: 'Ethanol', label: '乙醇 ε=24.9' },
      { value: 'Acetone', label: '丙酮 ε=20.5' },
    ],
  },
  {
    label: '📌 低介电常数 (ε<15)',
    options: [
      { value: 'DiethylEther', label: '乙醚 ε=4.2' },
      { value: 'Benzene', label: '苯 ε=2.3' },
      { value: 'Toluene', label: '甲苯 ε=2.4' },
      { value: 'CycloHexane', label: '环己烷 ε=2.0' },
    ],
  },
  {
    label: '📌 自定义',
    options: [
      { value: 'custom', label: '自定义溶剂 (手动输入介电常数)' },
    ],
  },
];

interface DesolvationBatchPanelProps {
  jobId: number;  // MD Job ID
  onStructureSelect?: (structureId: number) => void;  // 选中结构时的回调
}

interface SelectedStructure extends AutoSelectedStructure {
  selected: boolean;
}

// 待提交任务项
interface PendingTask {
  structureId: number;
  structureName: string;  // 显示名称
  compositionKey: string;
  methodLevel: 'fast' | 'standard' | 'accurate';
  desolvationMode: 'stepwise' | 'full';
  solventModel: SolventModel;
  solventName?: string;
  customEps?: number;
  slurmPartition: string;
  slurmCpus: number;
  slurmTime: number;
}

export default function DesolvationBatchPanel({ jobId, onStructureSelect }: DesolvationBatchPanelProps) {
  const { token } = theme.useToken();
  const { isDark } = useThemeStore();

  // 状态
  const [loading, setLoading] = useState(false);
  const [structures, setStructures] = useState<SelectedStructure[]>([]);
  const [selectedKeys, setSelectedKeys] = useState<number[]>([]);
  const [submitting, setSubmitting] = useState(false);
  const [overview, setOverview] = useState<DesolvationOverviewResponse | null>(null);
  const [expandedJobId, setExpandedJobId] = useState<number | null>(null);
  const [expandedRowKeys, setExpandedRowKeys] = useState<number[]>([]);
  const [qcTasksCache, setQcTasksCache] = useState<Record<number, DesolvationQCTasksResponse>>({});

  // 待提交队列
  const [pendingTasks, setPendingTasks] = useState<PendingTask[]>([]);

  // 多维度筛选条件
  const [cnFilter, setCnFilter] = useState<number[]>([]);  // 配位数筛选
  const [anionCountFilter, setAnionCountFilter] = useState<number[]>([]);  // 阴离子数量筛选
  const [solventTypeFilter, setSolventTypeFilter] = useState<string[]>([]);  // 溶剂类型筛选

  // 当 jobId 变化时，重置所有状态
  useEffect(() => {
    setStructures([]);
    setSelectedKeys([]);
    setOverview(null);
    setExpandedJobId(null);
    setExpandedRowKeys([]);
    setQcTasksCache({});
    setCnFilter([]);
    setAnionCountFilter([]);
    setSolventTypeFilter([]);
    setPendingTasks([]);  // 清空待提交队列
  }, [jobId]);

  // 计算参数
  const [desolvationMode, setDesolvationMode] = useState<'stepwise' | 'full'>('stepwise');
  const [methodLevel, setMethodLevel] = useState<'fast' | 'standard' | 'accurate'>('standard');
  const [solventModel, setSolventModel] = useState<SolventModel>('gas');
  const [solventName, setSolventName] = useState<string>('Water');
  const [customEps, setCustomEps] = useState<number>(80.0);  // 自定义介电常数

  // Slurm 资源配置
  const [slurmPartition, setSlurmPartition] = useState<string>('cpu');
  const [slurmCpus, setSlurmCpus] = useState<number>(16);
  const [slurmTime, setSlurmTime] = useState<number>(7200);  // 分钟

  // 辅助函数：计算结构中的阴离子数量（只计算 count > 0 的）
  const getAnionCount = (composition: Record<string, number>): number => {
    let count = 0;
    Object.entries(composition).forEach(([mol, num]) => {
      if (num > 0 && ANION_PATTERNS.some(anion => mol.toUpperCase().includes(anion.toUpperCase()))) {
        count += num;
      }
    });
    return count;
  };

  // 辅助函数：获取结构中的溶剂类型列表（只返回数量 > 0 的非阴离子）
  const getSolventTypes = (composition: Record<string, number>): string[] => {
    return Object.entries(composition)
      .filter(([mol, count]) =>
        count > 0 && !ANION_PATTERNS.some(anion => mol.toUpperCase().includes(anion.toUpperCase()))
      )
      .map(([mol]) => mol);
  };

  // 获取所有可用的配位数选项
  const availableCNs = React.useMemo(() => {
    const cnSet = new Set<number>();
    structures.forEach(s => cnSet.add(s.coordination_num));
    return Array.from(cnSet).sort((a, b) => a - b);
  }, [structures]);

  // 获取所有可用的阴离子数量选项
  const availableAnionCounts = React.useMemo(() => {
    const countSet = new Set<number>();
    structures.forEach(s => countSet.add(getAnionCount(s.composition)));
    return Array.from(countSet).sort((a, b) => a - b);
  }, [structures]);

  // 获取所有可用的溶剂类型选项
  const availableSolventTypes = React.useMemo(() => {
    const typeSet = new Set<string>();
    structures.forEach(s => {
      getSolventTypes(s.composition).forEach(type => typeSet.add(type));
    });
    return Array.from(typeSet).sort();
  }, [structures]);

  // 根据筛选条件过滤后的结构
  const filteredStructures = React.useMemo(() => {
    return structures.filter(s => {
      // 配位数筛选
      if (cnFilter.length > 0 && !cnFilter.includes(s.coordination_num)) {
        return false;
      }
      // 阴离子数量筛选
      if (anionCountFilter.length > 0 && !anionCountFilter.includes(getAnionCount(s.composition))) {
        return false;
      }
      // 溶剂类型筛选 - 结构必须包含所有选中的溶剂类型
      if (solventTypeFilter.length > 0) {
        const solvents = getSolventTypes(s.composition);
        // 结构的溶剂必须包含所有筛选的溶剂类型
        if (!solventTypeFilter.every(type => solvents.includes(type))) {
          return false;
        }
      }
      return true;
    });
  }, [structures, cnFilter, anionCountFilter, solventTypeFilter]);

  // 当筛选条件变化时，更新选中的 keys
  useEffect(() => {
    const hasFilter = cnFilter.length > 0 || anionCountFilter.length > 0 || solventTypeFilter.length > 0;
    if (hasFilter) {
      const filteredIds = filteredStructures.map(s => s.id);
      setSelectedKeys(prev => prev.filter(id => filteredIds.includes(id)));
    }
  }, [cnFilter, anionCountFilter, solventTypeFilter, filteredStructures]);

  // 检测是否有阴离子
  const hasAnion = structures.some(s => getAnionCount(s.composition) > 0);

  // 加载自动挑选的结构
  const loadAutoSelectedStructures = useCallback(async () => {
    setLoading(true);
    try {
      const result = await autoSelectSolvationStructures(jobId);
      const selected = result.selected_structures.map(s => ({
        ...s,
        selected: true,
      }));
      setStructures(selected);
      setSelectedKeys(selected.map(s => s.id));
      message.success(`已自动挑选 ${result.unique_compositions} 种不同配位组成`);
    } catch (error) {
      message.error('加载溶剂化结构失败');
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  // 加载任务总览
  const loadOverview = useCallback(async () => {
    try {
      const data = await getDesolvationOverview(jobId);
      setOverview(data);
    } catch (error) {
      // 可能没有任务，忽略错误
    }
  }, [jobId]);

  useEffect(() => {
    loadOverview();
  }, [loadOverview]);

  // 添加到待提交队列
  const handleAddToQueue = () => {
    if (selectedKeys.length === 0) {
      message.warning('请选择要计算的溶剂化结构');
      return;
    }

    // 获取选中结构的信息
    const newTasks: PendingTask[] = selectedKeys.map(structureId => {
      const structure = structures.find(s => s.id === structureId);
      return {
        structureId,
        structureName: structure?.composition_key || `结构 #${structureId}`,
        compositionKey: structure?.composition_key || '',
        methodLevel,
        desolvationMode,
        solventModel,
        solventName: solventModel !== 'gas' ? solventName : undefined,
        customEps: solventName === 'custom' ? customEps : undefined,
        slurmPartition,
        slurmCpus,
        slurmTime,
      };
    });

    // 过滤掉已存在的任务
    const existingIds = new Set(pendingTasks.map(t => t.structureId));
    const uniqueNewTasks = newTasks.filter(t => !existingIds.has(t.structureId));

    if (uniqueNewTasks.length === 0) {
      message.warning('所选结构已在待提交队列中');
      return;
    }

    setPendingTasks(prev => [...prev, ...uniqueNewTasks]);
    setSelectedKeys([]);  // 清空选择
    message.success(`已添加 ${uniqueNewTasks.length} 个任务到待提交队列`);
  };

  // 从队列中移除任务
  const handleRemoveFromQueue = (structureId: number) => {
    setPendingTasks(prev => prev.filter(t => t.structureId !== structureId));
  };

  // 提交单个任务
  const handleSubmitSingle = async (task: PendingTask) => {
    setSubmitting(true);
    try {
      let solventConfig: SolventConfig | undefined;
      if (task.solventModel !== 'gas') {
        if (task.solventName === 'custom') {
          solventConfig = { model: 'custom', eps: task.customEps };
        } else {
          solventConfig = { model: task.solventModel, solvent_name: task.solventName };
        }
      }

      await batchCreateDesolvationJobs({
        md_job_id: jobId,
        structure_ids: [task.structureId],
        method_level: task.methodLevel,
        desolvation_mode: task.desolvationMode,
        solvent_config: solventConfig,
        slurm_partition: task.slurmPartition,
        slurm_cpus: task.slurmCpus,
        slurm_time: task.slurmTime,
      });

      // 从队列中移除
      setPendingTasks(prev => prev.filter(t => t.structureId !== task.structureId));
      message.success(`任务 ${task.structureName} 已提交`);
      loadOverview();
    } catch (error: any) {
      message.error(`提交失败: ${error.message || '未知错误'}`);
    } finally {
      setSubmitting(false);
    }
  };

  // 批量提交所有待提交任务
  const handleSubmitAll = async () => {
    if (pendingTasks.length === 0) {
      message.warning('待提交队列为空');
      return;
    }

    setSubmitting(true);
    let successCount = 0;
    let failCount = 0;

    for (const task of pendingTasks) {
      try {
        let solventConfig: SolventConfig | undefined;
        if (task.solventModel !== 'gas') {
          if (task.solventName === 'custom') {
            solventConfig = { model: 'custom', eps: task.customEps };
          } else {
            solventConfig = { model: task.solventModel, solvent_name: task.solventName };
          }
        }

        await batchCreateDesolvationJobs({
          md_job_id: jobId,
          structure_ids: [task.structureId],
          method_level: task.methodLevel,
          desolvation_mode: task.desolvationMode,
          solvent_config: solventConfig,
          slurm_partition: task.slurmPartition,
          slurm_cpus: task.slurmCpus,
          slurm_time: task.slurmTime,
        });
        successCount++;
      } catch {
        failCount++;
      }
    }

    setPendingTasks([]);  // 清空队列
    message.success(`已提交 ${successCount} 个任务${failCount > 0 ? `，${failCount} 个失败` : ''}`);
    loadOverview();
    setSubmitting(false);
  };

  // 加载 QC 子任务
  const loadQCTasks = useCallback(async (jobId: number) => {
    if (qcTasksCache[jobId]) return; // 已缓存
    try {
      const data = await getDesolvationQCTasks(jobId);
      setQcTasksCache(prev => ({ ...prev, [jobId]: data }));
    } catch (error) {
      console.error('加载QC子任务失败:', error);
    }
  }, [qcTasksCache]);

  // 展开行渲染
  const expandedRowRender = (record: DesolvationJobResponse) => {
    const qcData = qcTasksCache[record.job_id];

    if (!qcData) {
      return (
        <div style={{ padding: 16, textAlign: 'center' }}>
          <Spin tip="加载子任务..." />
        </div>
      );
    }

    const qcTasks = qcData.qc_tasks;

    // 如果没有子任务，根据状态显示不同提示
    if (qcTasks.length === 0) {
      let message = '暂无子任务';
      let icon = <ClockCircleOutlined />;
      if (record.status === 'SUBMITTED' || record.status === 'QUEUED') {
        message = '任务排队中，等待创建 QC 计算子任务...';
      } else if (record.status === 'FAILED') {
        message = record.error_message || '任务失败，未创建子任务';
        icon = <ExclamationCircleOutlined style={{ color: '#ff4d4f' }} />;
      }
      return (
        <div style={{ padding: 16, textAlign: 'center', color: token.colorTextSecondary }}>
          {icon} <span style={{ marginLeft: 8 }}>{message}</span>
        </div>
      );
    }

    // 按类型分组
    const clusterTask = qcTasks.find(t => t.task_type === 'cluster');
    const clusterMinusTasks = qcTasks.filter(t => t.task_type === 'cluster_minus');
    const ligandTasks = qcTasks.filter(t => t.task_type === 'ligand');

    const getStatusTag = (status: string, isReused?: boolean) => {
      if (isReused) {
        return <Tag color="cyan" icon={<CheckCircleOutlined />}>复用</Tag>;
      }
      const statusConfig: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
        CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
        SUBMITTED: { color: 'blue', icon: <ClockCircleOutlined />, text: '已提交' },
        QUEUED: { color: 'cyan', icon: <ClockCircleOutlined />, text: '排队中' },
        RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
        COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '完成' },
        FAILED: { color: 'error', icon: <ExclamationCircleOutlined />, text: '失败' },
      };
      const config = statusConfig[status] || { color: 'default', icon: null, text: status };
      return <Tag color={config.color} icon={config.icon}>{config.text}</Tag>;
    };

    return (
      <div style={{
        padding: '12px 16px',
        background: isDark ? 'rgba(0,0,0,0.2)' : '#fafafa',
        borderRadius: 4,
      }}>
        {/* 统计信息 */}
        <Row gutter={16} style={{ marginBottom: 12 }}>
          <Col><Text type="secondary">子任务: {qcData.total}</Text></Col>
          <Col><Text style={{ color: '#52c41a' }}>✓ {qcData.completed}</Text></Col>
          <Col><Text style={{ color: '#1890ff' }}>⟳ {qcData.running}</Text></Col>
          <Col><Text style={{ color: '#ff4d4f' }}>✗ {qcData.failed}</Text></Col>
          {qcData.reused > 0 && <Col><Text style={{ color: '#13c2c2' }}>♻ 复用 {qcData.reused}</Text></Col>}
        </Row>

        {/* Cluster 完整结构 */}
        {clusterTask && (
          <div style={{ marginBottom: 8 }}>
            <Text strong style={{ fontSize: 12 }}>完整 Cluster:</Text>
            <div style={{ marginLeft: 16, marginTop: 4 }}>
              <Space size={8}>
                <Text style={{ fontSize: 11 }}>{clusterTask.molecule_name}</Text>
                {getStatusTag(clusterTask.status, clusterTask.is_reused)}
                <Text type="secondary" style={{ fontSize: 10 }}>
                  {clusterTask.functional}/{clusterTask.basis_set}
                </Text>
              </Space>
            </div>
          </div>
        )}

        {/* Cluster-minus 结构 */}
        {clusterMinusTasks.length > 0 && (
          <div style={{ marginBottom: 8 }}>
            <Text strong style={{ fontSize: 12 }}>去配体 Cluster ({clusterMinusTasks.length}):</Text>
            <div style={{ marginLeft: 16, marginTop: 4, maxHeight: 150, overflowY: 'auto' }}>
              {clusterMinusTasks.map(task => (
                <div key={task.id} style={{ marginBottom: 6, display: 'flex', alignItems: 'center', gap: 8 }}>
                  <Text style={{ fontSize: 11, fontFamily: 'monospace', minWidth: 200 }}>
                    {task.molecule_name}
                  </Text>
                  <Text type="secondary" style={{ fontSize: 10 }}>
                    {task.functional}/{task.basis_set}
                  </Text>
                  {getStatusTag(task.status, task.is_reused)}
                  {task.error_message && (
                    <Tooltip title={task.error_message}>
                      <ExclamationCircleOutlined style={{ color: '#ff4d4f', fontSize: 11 }} />
                    </Tooltip>
                  )}
                </div>
              ))}
            </div>
          </div>
        )}

        {/* 配体分子 */}
        {ligandTasks.length > 0 && (
          <div>
            <Text strong style={{ fontSize: 12 }}>配体分子 ({ligandTasks.length}):</Text>
            <div style={{ marginLeft: 16, marginTop: 4 }}>
              <Space wrap size={4}>
                {ligandTasks.map(task => (
                  <Tooltip
                    key={task.id}
                    title={`${task.functional}/${task.basis_set} | charge=${task.charge}${task.is_reused ? ' (复用)' : ''}`}
                  >
                    <Tag
                      color={task.status === 'COMPLETED' ? 'success' : task.status === 'FAILED' ? 'error' : 'processing'}
                      style={{ fontSize: 11 }}
                    >
                      {task.molecule_name}
                      {task.is_reused && ' ♻'}
                    </Tag>
                  </Tooltip>
                ))}
              </Space>
            </div>
          </div>
        )}
      </div>
    );
  };

  // 结构表格列
  const structureColumns: ColumnsType<SelectedStructure> = [
    {
      title: '配位组成',
      dataIndex: 'composition_key',
      key: 'composition_key',
      render: (key: string, record) => (
        <Space direction="vertical" size={0}>
          <Text strong style={{ fontSize: 13 }}>{key}</Text>
          <Text type="secondary" style={{ fontSize: 11 }}>
            {record.center_ion}⁺ CN={record.coordination_num}
          </Text>
        </Space>
      ),
    },
    {
      title: '分子组成',
      dataIndex: 'composition',
      key: 'composition',
      render: (composition: Record<string, number>) => (
        <Space size={4} wrap>
          {Object.entries(composition)
            .filter(([_, count]) => count > 0)
            .map(([mol, count]) => (
              <Tag key={mol} style={{ margin: 0, fontSize: 11 }}>
                {mol}: {count}
              </Tag>
            ))}
        </Space>
      ),
    },
    {
      title: '帧号',
      dataIndex: 'frame_index',
      key: 'frame_index',
      width: 80,
      render: (frame: number) => <Text type="secondary">#{frame}</Text>,
    },
  ];

  // 任务表格列
  const jobColumns: ColumnsType<DesolvationJobResponse> = [
    {
      title: '任务名称',
      key: 'structure',
      width: 240,
      render: (_, record) => (
        <Space direction="vertical" size={0}>
          <Text strong style={{ fontSize: 12, fontFamily: 'monospace' }}>
            {record.electrolyte_name || record.composition_key || `结构 #${record.solvation_structure_id}`}
          </Text>
          {record.composition_key && record.electrolyte_name && (
            <Text type="secondary" style={{ fontSize: 10 }}>
              配位: {record.composition_key}
            </Text>
          )}
        </Space>
      ),
    },
    {
      title: '计算方法',
      key: 'method',
      width: 200,
      render: (_, record) => {
        const methodConfig: Record<string, { functional: string; basis: string; color: string }> = {
          fast: { functional: 'B3LYP', basis: '6-31G(d)', color: 'green' },
          standard: { functional: 'B3LYP', basis: '6-31++G(d,p)', color: 'blue' },
          accurate: { functional: 'ωB97XD', basis: '6-311++G(2d,2p)', color: 'purple' },
        };
        const m = methodConfig[record.method_level] || { functional: '?', basis: '?', color: 'default' };

        // 溶剂模型显示
        let solventDisplay = '气相';
        if (record.solvent_config) {
          const model = record.solvent_config.model?.toUpperCase() || '';
          if (record.solvent_config.model === 'custom') {
            solventDisplay = `${model} (ε=${record.solvent_config.eps || '?'})`;
          } else if (record.solvent_config.solvent_name) {
            solventDisplay = `${model}: ${record.solvent_config.solvent_name}`;
          } else {
            solventDisplay = model;
          }
        }

        return (
          <Space direction="vertical" size={0}>
            <Text style={{ fontSize: 12 }}>{m.functional}/{m.basis}</Text>
            <Text type="secondary" style={{ fontSize: 10 }}>
              溶剂: {solventDisplay}
            </Text>
          </Space>
        );
      },
    },
    {
      title: '状态',
      key: 'status',
      width: 150,
      render: (_, record) => {
        const statusConfig: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
          CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
          SUBMITTED: { color: 'blue', icon: <ClockCircleOutlined />, text: '已提交' },
          QUEUED: { color: 'cyan', icon: <ClockCircleOutlined />, text: '排队中' },
          RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
          COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
          FAILED: { color: 'error', icon: <ExclamationCircleOutlined />, text: '失败' },
        };
        const config = statusConfig[record.status] || { color: 'default', icon: null, text: record.status };
        
        return (
          <Space direction="vertical" size={0}>
            <Tag color={config.color} icon={config.icon}>{config.text}</Tag>
            {record.qc_progress && (
              <Progress
                percent={record.qc_progress.progress_percent}
                size="small"
                style={{ width: 100 }}
                format={() => `${record.qc_progress?.completed}/${record.qc_progress?.total}`}
              />
            )}
          </Space>
        );
      },
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 140,
      render: (time: string) => new Date(time).toLocaleString('zh-CN', {
        month: '2-digit',
        day: '2-digit',
        hour: '2-digit',
        minute: '2-digit',
      }),
    },
    {
      title: '操作',
      key: 'action',
      width: 120,
      render: (_, record) => (
        <Space size={4}>
          <Button
            type="link"
            size="small"
            disabled={record.status !== 'COMPLETED'}
            onClick={(e) => {
              e.stopPropagation();
              setExpandedJobId(expandedJobId === record.job_id ? null : record.job_id);
              if (expandedJobId !== record.job_id) {
                setExpandedRowKeys([record.job_id]);
              }
            }}
          >
            {expandedJobId === record.job_id ? '收起' : '结果'}
          </Button>
        </Space>
      ),
    },
  ];

  return (
    <Card
      title={
        <Space>
          <ThunderboltOutlined style={{ color: '#1890ff' }} />
          <span>去溶剂化能计算</span>
          {overview && overview.total_jobs > 0 && (
            <Badge
              count={overview.status_summary['RUNNING'] || 0}
              style={{ backgroundColor: '#1890ff' }}
              title="运行中的任务"
            />
          )}
        </Space>
      }
      extra={
        <Button
          icon={<ReloadOutlined />}
          onClick={loadOverview}
          size="small"
        >
          刷新
        </Button>
      }
      style={{
        background: isDark ? token.colorBgContainer : undefined,
        borderColor: token.colorBorder,
      }}
    >
      {/* 第一步：挑选结构 */}
      <Collapse
        defaultActiveKey={structures.length === 0 ? ['select'] : []}
        items={[{
          key: 'select',
          label: (
            <Space>
              <span>第一步：挑选溶剂化结构</span>
              {structures.length > 0 && (
                <Tag color="blue">{selectedKeys.length} 个已选</Tag>
              )}
            </Space>
          ),
          children: (
            <div>
              <div style={{ marginBottom: 16 }}>
                <Space wrap>
                  <Button
                    type="primary"
                    icon={<BulbOutlined />}
                    onClick={loadAutoSelectedStructures}
                    loading={loading}
                  >
                    自动挑选不同配位组成
                  </Button>
                  <Text type="secondary" style={{ fontSize: 12 }}>
                    系统会自动从所有溶剂化结构中挑选出不同配位组成的代表性结构
                  </Text>
                </Space>
              </div>

              {/* 多维度筛选器 */}
              {structures.length > 0 && (
                <div style={{
                  marginBottom: 16,
                  padding: '12px 16px',
                  background: isDark ? 'rgba(24, 144, 255, 0.05)' : '#f0f5ff',
                  border: `1px solid ${isDark ? 'rgba(24, 144, 255, 0.2)' : '#adc6ff'}`,
                  borderRadius: 8,
                }}>
                  <Space size={4} style={{ marginBottom: 12 }}>
                    <FilterOutlined style={{ color: token.colorPrimary }} />
                    <Text strong style={{ fontSize: 13 }}>筛选条件</Text>
                    {(cnFilter.length > 0 || anionCountFilter.length > 0 || solventTypeFilter.length > 0) && (
                      <Tag color="blue">筛选后 {filteredStructures.length} 个结构</Tag>
                    )}
                  </Space>

                  <Row gutter={[12, 12]}>
                    {/* 配位数筛选 */}
                    {availableCNs.length > 0 && (
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          配位数 (CN)
                        </Text>
                        <Select
                          mode="multiple"
                          placeholder="全部"
                          value={cnFilter}
                          onChange={setCnFilter}
                          style={{ width: '100%' }}
                          size="small"
                          allowClear
                          maxTagCount={2}
                          options={availableCNs.map(cn => ({
                            label: `CN=${cn}`,
                            value: cn,
                          }))}
                        />
                      </Col>
                    )}

                    {/* 阴离子数量筛选 */}
                    {availableAnionCounts.length > 0 && (
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          阴离子数量
                        </Text>
                        <Select
                          mode="multiple"
                          placeholder="全部"
                          value={anionCountFilter}
                          onChange={setAnionCountFilter}
                          style={{ width: '100%' }}
                          size="small"
                          allowClear
                          maxTagCount={2}
                          options={availableAnionCounts.map(count => ({
                            label: count === 0 ? '无阴离子' : `${count}个阴离子`,
                            value: count,
                          }))}
                        />
                      </Col>
                    )}

                    {/* 溶剂类型筛选 */}
                    {availableSolventTypes.length > 0 && (
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          溶剂类型
                        </Text>
                        <Select
                          mode="multiple"
                          placeholder="全部"
                          value={solventTypeFilter}
                          onChange={setSolventTypeFilter}
                          style={{ width: '100%' }}
                          size="small"
                          allowClear
                          maxTagCount={2}
                          options={availableSolventTypes.map(type => ({
                            label: type,
                            value: type,
                          }))}
                        />
                      </Col>
                    )}
                  </Row>

                  <div style={{ marginTop: 12 }}>
                    <Space size={8}>
                      <Button
                        size="small"
                        type="primary"
                        ghost
                        onClick={() => setSelectedKeys(filteredStructures.map(s => s.id))}
                      >
                        全选当前 ({filteredStructures.length})
                      </Button>
                      <Button
                        size="small"
                        onClick={() => setSelectedKeys([])}
                      >
                        清空选择
                      </Button>
                      <Button
                        size="small"
                        onClick={() => {
                          setCnFilter([]);
                          setAnionCountFilter([]);
                          setSolventTypeFilter([]);
                          setSelectedKeys(structures.map(s => s.id));
                        }}
                      >
                        重置筛选
                      </Button>
                    </Space>
                  </div>
                </div>
              )}

              {structures.length > 0 && (
                <Table
                  dataSource={filteredStructures}
                  columns={structureColumns}
                  rowKey="id"
                  size="small"
                  rowSelection={{
                    selectedRowKeys: selectedKeys,
                    onChange: (keys) => setSelectedKeys(keys as number[]),
                  }}
                  pagination={false}
                  scroll={{ y: 200 }}
                  onRow={(record) => ({
                    onClick: () => onStructureSelect?.(record.id),
                    style: { cursor: 'pointer' },
                  })}
                />
              )}
            </div>
          ),
        }]}
      />

      {/* 第二步：设置参数并提交 */}
      {structures.length > 0 && (
        <Collapse
          style={{ marginTop: 16 }}
          defaultActiveKey={['params']}
          items={[{
            key: 'params',
            label: '第二步：设置计算参数并提交',
            children: (
              <div>
                {/* 智能推荐 */}
                {hasAnion && (
                  <Alert
                    message={
                      <Space size={4}>
                        <BulbOutlined />
                        <span><strong>智能推荐：</strong>检测到阴离子，建议选择带弥散函数的基组（标准或精确）</span>
                      </Space>
                    }
                    type="warning"
                    showIcon={false}
                    style={{ marginBottom: 16 }}
                  />
                )}

                <Row gutter={[16, 16]}>
                  {/* 1. 计算模式 */}
                  <Col span={8}>
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                      1. 计算模式
                    </Text>
                    <Select
                      value={desolvationMode}
                      onChange={setDesolvationMode}
                      style={{ width: '100%' }}
                      options={[
                        { label: '逐级去溶剂 (推荐)', value: 'stepwise' },
                        { label: '全部去溶剂', value: 'full' },
                      ]}
                    />
                    <Text type="secondary" style={{ fontSize: 11, marginTop: 4, display: 'block' }}>
                      {desolvationMode === 'stepwise' ? '依次移除每个配体计算能量' : '一次性移除所有配体'}
                    </Text>
                  </Col>

                  {/* 2. 计算方法 */}
                  <Col span={8}>
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                      2. 计算方法
                    </Text>
                    <Select
                      value={methodLevel}
                      onChange={setMethodLevel}
                      style={{ width: '100%' }}
                      options={[
                        { label: '快速 (B3LYP/6-31G(d))', value: 'fast' },
                        { label: '标准 (B3LYP/6-31++G(d,p))', value: 'standard' },
                        { label: '精确 (ωB97XD/6-311++G(2d,2p))', value: 'accurate' },
                      ]}
                    />
                    <Text type="secondary" style={{ fontSize: 11, marginTop: 4, display: 'block' }}>
                      {methodLevel === 'fast' ? '适合快速预筛选' : methodLevel === 'standard' ? '平衡精度与速度' : '高精度计算'}
                    </Text>
                  </Col>

                  {/* 3. 溶剂模型 */}
                  <Col span={8}>
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                      3. 溶剂模型
                    </Text>
                    <Select
                      value={solventModel}
                      onChange={(value) => {
                        setSolventModel(value);
                        if (value !== 'gas' && !solventName) {
                          setSolventName('Water');
                        }
                      }}
                      style={{ width: '100%' }}
                      options={[
                        { label: '气相 (无溶剂)', value: 'gas' },
                        { label: 'PCM (极化连续介质)', value: 'pcm' },
                        { label: 'SMD (溶剂密度模型)', value: 'smd' },
                      ]}
                    />
                    <Text type="secondary" style={{ fontSize: 11, marginTop: 4, display: 'block' }}>
                      {solventModel === 'gas' ? '真空环境计算' : solventModel === 'pcm' ? '通过介电常数模拟溶剂' : '更精确的隐式溶剂'}
                    </Text>
                  </Col>
                </Row>

                {/* 隐式溶剂选择 */}
                {solventModel !== 'gas' && (
                  <Row gutter={[16, 16]} style={{ marginTop: 16 }}>
                    <Col span={8}>
                      <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                        隐式溶剂
                      </Text>
                      <Select
                        value={solventName}
                        onChange={setSolventName}
                        style={{ width: '100%' }}
                        placeholder="选择溶剂"
                        showSearch
                        optionFilterProp="label"
                        options={SOLVENT_OPTIONS}
                      />
                    </Col>
                    {solventName === 'custom' && (
                      <Col span={8}>
                        <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          介电常数 ε
                        </Text>
                        <InputNumber
                          value={customEps}
                          onChange={(v) => setCustomEps(v || 80)}
                          min={1}
                          max={200}
                          step={0.1}
                          style={{ width: '100%' }}
                          placeholder="输入介电常数"
                        />
                      </Col>
                    )}
                    <Col span={solventName === 'custom' ? 8 : 16}>
                      <div style={{
                        padding: '8px 12px',
                        background: isDark ? 'rgba(24, 144, 255, 0.1)' : '#e6f4ff',
                        borderRadius: 6,
                        marginTop: 20,
                      }}>
                        <Text style={{ fontSize: 11 }}>
                          💡 <strong>提示：</strong>电池电解液建议选 EC (ε=89.8) 或 PC (ε=64.9)
                        </Text>
                      </div>
                    </Col>
                  </Row>
                )}

                {/* 计算资源配置 */}
                <Card
                  size="small"
                  title="计算资源配置"
                  style={{ marginTop: 16, background: isDark ? 'rgba(0,0,0,0.2)' : '#fafafa' }}
                >
                  <Row gutter={16}>
                    <Col span={8}>
                      <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                        队列/分区
                      </Text>
                      <Select
                        value={slurmPartition}
                        onChange={setSlurmPartition}
                        style={{ width: '100%' }}
                        options={[
                          { label: 'cpu', value: 'cpu' },
                          { label: 'gpu', value: 'gpu' },
                          { label: 'debug', value: 'debug' },
                        ]}
                      />
                    </Col>
                    <Col span={8}>
                      <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                        CPU 核心数
                      </Text>
                      <InputNumber
                        value={slurmCpus}
                        onChange={(v) => setSlurmCpus(v || 16)}
                        min={1}
                        max={64}
                        style={{ width: '100%' }}
                      />
                    </Col>
                    <Col span={8}>
                      <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                        最大时间 (分钟)
                      </Text>
                      <InputNumber
                        value={slurmTime}
                        onChange={(v) => setSlurmTime(v || 7200)}
                        min={10}
                        max={43200}
                        style={{ width: '100%' }}
                      />
                    </Col>
                  </Row>
                  <Text type="secondary" style={{ fontSize: 11, marginTop: 8, display: 'block' }}>
                    💡 QC计算通常使用 16 核，标准基组需要 30分钟~数小时，大基组可能需要更长时间
                  </Text>
                </Card>

                <div style={{ marginTop: 20 }}>
                  <Button
                    type="primary"
                    icon={<PlusOutlined />}
                    onClick={handleAddToQueue}
                    disabled={selectedKeys.length === 0}
                    size="large"
                  >
                    添加到待提交队列 ({selectedKeys.length} 个)
                  </Button>
                </div>
              </div>
            ),
          }]}
        />
      )}

      {/* 第三步：待提交队列 */}
      {pendingTasks.length > 0 && (
        <Card
          size="small"
          title={
            <Space>
              <span>待提交队列</span>
              <Tag color="orange">{pendingTasks.length} 个任务</Tag>
            </Space>
          }
          style={{ marginTop: 16 }}
          extra={
            <Space>
              <Button
                danger
                size="small"
                onClick={() => setPendingTasks([])}
              >
                清空队列
              </Button>
              <Button
                type="primary"
                icon={<ThunderboltOutlined />}
                onClick={handleSubmitAll}
                loading={submitting}
              >
                全部提交
              </Button>
            </Space>
          }
        >
          <Table
            dataSource={pendingTasks}
            rowKey="structureId"
            size="small"
            pagination={false}
            columns={[
              {
                title: '结构',
                key: 'structure',
                render: (_, task) => (
                  <Space direction="vertical" size={0}>
                    <Text strong style={{ fontSize: 12 }}>{task.structureName}</Text>
                    <Text type="secondary" style={{ fontSize: 10 }}>{task.compositionKey}</Text>
                  </Space>
                ),
              },
              {
                title: '计算方法',
                key: 'method',
                width: 180,
                render: (_, task) => {
                  const methodConfig: Record<string, string> = {
                    fast: 'B3LYP/6-31G(d)',
                    standard: 'B3LYP/6-31++G(d,p)',
                    accurate: 'ωB97XD/6-311++G(2d,2p)',
                  };
                  return <Text style={{ fontSize: 12 }}>{methodConfig[task.methodLevel]}</Text>;
                },
              },
              {
                title: '溶剂模型',
                key: 'solvent',
                width: 120,
                render: (_, task) => {
                  if (task.solventModel === 'gas') return <Text style={{ fontSize: 12 }}>气相</Text>;
                  if (task.solventName === 'custom') {
                    return <Text style={{ fontSize: 12 }}>{task.solventModel.toUpperCase()} (ε={task.customEps})</Text>;
                  }
                  return <Text style={{ fontSize: 12 }}>{task.solventModel.toUpperCase()}: {task.solventName}</Text>;
                },
              },
              {
                title: '资源',
                key: 'resource',
                width: 140,
                render: (_, task) => (
                  <Text type="secondary" style={{ fontSize: 11 }}>
                    {task.slurmPartition} / {task.slurmCpus}核 / {task.slurmTime}分钟
                  </Text>
                ),
              },
              {
                title: '操作',
                key: 'action',
                width: 120,
                render: (_, task) => (
                  <Space size={4}>
                    <Button
                      type="primary"
                      size="small"
                      onClick={() => handleSubmitSingle(task)}
                      loading={submitting}
                    >
                      提交
                    </Button>
                    <Button
                      danger
                      size="small"
                      onClick={() => handleRemoveFromQueue(task.structureId)}
                    >
                      移除
                    </Button>
                  </Space>
                ),
              },
            ]}
          />
        </Card>
      )}

      {/* 第四步：任务监控 */}
      {overview && overview.total_jobs > 0 && (
        <div style={{ marginTop: 16 }}>
          <Divider orientation="left">
            <Space>
              任务监控
              <Tag color="blue">{overview.total_jobs} 个任务</Tag>
              {overview.status_summary['COMPLETED'] > 0 && (
                <Tag color="success">{overview.status_summary['COMPLETED']} 完成</Tag>
              )}
              {(overview.status_summary['RUNNING'] || 0) + (overview.status_summary['QUEUED'] || 0) > 0 && (
                <Tag color="processing">
                  {(overview.status_summary['RUNNING'] || 0) + (overview.status_summary['QUEUED'] || 0)} 进行中
                </Tag>
              )}
              {overview.status_summary['FAILED'] > 0 && (
                <Tag color="error">{overview.status_summary['FAILED']} 失败</Tag>
              )}
            </Space>
          </Divider>

          <Table
            dataSource={overview.jobs}
            columns={jobColumns}
            rowKey="job_id"
            size="small"
            pagination={{ pageSize: 5, size: 'small' }}
            expandable={{
              expandedRowKeys: expandedRowKeys,
              onExpand: (expanded, record) => {
                if (expanded) {
                  setExpandedRowKeys([record.job_id]);
                  loadQCTasks(record.job_id);
                } else {
                  setExpandedRowKeys([]);
                }
              },
              expandedRowRender: (record) => {
                // 如果有结果且要查看结果
                if (expandedJobId === record.job_id && record.result) {
                  return <DesolvationResultView result={record.result} />;
                }
                // 否则展示子任务
                return expandedRowRender(record);
              },
              rowExpandable: () => true,
            }}
          />
        </div>
      )}
    </Card>
  );
}

