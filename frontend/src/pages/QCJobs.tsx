/**
 * QC量子化学计算任务管理页面
 */
import { useState, useEffect, useCallback, useRef } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Button,
  Space,
  message,
  Modal,
  Form,
  Input,
  Select,
  InputNumber,
  Spin,
  Empty,
  Tabs,
  Typography,
  Card,
  Table,
  Tag,
  Tooltip,
  Popconfirm,
  Row,
  Col,
  Divider,
  Switch,
  Alert,
  Collapse,
  AutoComplete,
  Descriptions,
  Checkbox,
  Steps,
  Result,
} from 'antd';
import {
  PlusOutlined,
  ReloadOutlined,
  ExperimentOutlined,
  DeleteOutlined,
  PlayCircleOutlined,
  EyeOutlined,
  EditOutlined,
  ThunderboltOutlined,
  SettingOutlined,
  BulbOutlined,
  RedoOutlined,
  CopyOutlined,
  CheckSquareOutlined,
  UploadOutlined,
  DownloadOutlined,
  FileExcelOutlined,
  InboxOutlined,
  StopOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import {
  getQCJobs,
  createQCJob,
  updateQCJob,
  deleteQCJob,
  submitQCJob,
  getBasisSets,
  getFunctionals,
  getAccuracyLevels,
  getSolventModels,
  getSolvents,
  getCommonMolecules,
  calculateSpinMultiplicity,
  batchDeleteQCJobs,
  batchSubmitQCJobs,
  batchCancelQCJobs,
  checkDuplicateCalculations,
} from '../api/qc';
import { getPartitions, getSlurmSuggestion, PartitionInfo } from '../api/slurm';
import { downloadTemplate, batchImportUpload, BatchImportResult } from '../api/batchImport';
import type { QCJob, QCJobCreate, SolventConfig } from '../types/qc';

const { Title, Text } = Typography;
const { TextArea } = Input;
const { Panel } = Collapse;
const { Step } = Steps;

// QC任务状态映射
const statusMap: Record<string, { color: string; text: string }> = {
  CREATED: { color: 'default', text: '已创建' },
  QUEUED: { color: 'processing', text: '排队中' },
  RUNNING: { color: 'processing', text: '运行中' },
  POSTPROCESSING: { color: 'processing', text: '后处理中' },
  COMPLETED: { color: 'success', text: '已完成' },
  FAILED: { color: 'error', text: '失败' },
  CANCELLED: { color: 'warning', text: '已取消' },
};

// 分子类型映射
const moleculeTypeMap: Record<string, string> = {
  solvent: '溶剂',
  cation: '阳离子',
  anion: '阴离子',
  custom: '自定义',
};

// 精度等级类型
interface AccuracyLevel {
  value: string;
  label: string;
  functional: string | null;
  basis_set: string | null;
  description: string;
  estimated_time: string;
}

// 溶剂模型类型
interface SolventModelOption {
  value: string;
  label: string;
  description: string;
}

// 溶剂选项类型
interface SolventOption {
  value: string;
  label: string;
  eps: number;
  description: string;
}

// 常用分子类型
interface CommonMolecule {
  name: string;       // 英文名称，用于生成Gaussian文件
  label?: string;     // 中文显示标签
  smiles: string;
  charge: number;
}

interface MoleculeCategory {
  name: string;
  molecules: CommonMolecule[];
}

export default function QCJobs() {
  const navigate = useNavigate();
  const [jobs, setJobs] = useState<QCJob[]>([]);
  const [total, setTotal] = useState(0);
  const [loading, setLoading] = useState(false);
  const [modalVisible, setModalVisible] = useState(false);
  const [submitting, setSubmitting] = useState(false);
  const [activeTab, setActiveTab] = useState('all');
  const [basisSets, setBasisSets] = useState<{ value: string; label: string; category?: string; description?: string }[]>([]);
  const [functionals, setFunctionals] = useState<{ value: string; label: string; category?: string; description?: string }[]>([]);
  const [accuracyLevels, setAccuracyLevels] = useState<AccuracyLevel[]>([]);
  const [solventModels, setSolventModels] = useState<SolventModelOption[]>([]);
  const [solvents, setSolvents] = useState<SolventOption[]>([]);
  const [commonMolecules, setCommonMolecules] = useState<MoleculeCategory[]>([]);
  const [selectedAccuracy, setSelectedAccuracy] = useState('standard');
  const [selectedSolventModel, setSelectedSolventModel] = useState('gas');
  const [autoSpin, setAutoSpin] = useState(true);
  const [calculatedSpin, setCalculatedSpin] = useState<number | null>(null);
  const [smilesOptions, setSmilesOptions] = useState<{ value: string; label: string; charge?: number }[]>([]);
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);
  const [useRecommendedParams, setUseRecommendedParams] = useState(true);
  const [recommendationReason, setRecommendationReason] = useState<string>('');
  const [moleculeCollapseKey, setMoleculeCollapseKey] = useState<string[]>([]);

  // 批量选择相关状态
  const [selectedRowKeys, setSelectedRowKeys] = useState<React.Key[]>([]);
  const [editingJob, setEditingJob] = useState<QCJob | null>(null);
  const [form] = Form.useForm();
  const [pagination, setPagination] = useState({ current: 1, pageSize: 20 });
  const pollingRef = useRef<ReturnType<typeof setInterval> | null>(null);

  // 批量导入相关状态
  const [batchImportVisible, setBatchImportVisible] = useState(false);
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [importing, setImporting] = useState(false);
  const [importStep, setImportStep] = useState(0);
  const [importResult, setImportResult] = useState<BatchImportResult | null>(null);

  // 加载QC任务列表
  const loadJobs = useCallback(async (status?: string) => {
    try {
      const skip = (pagination.current - 1) * pagination.pageSize;
      const data = await getQCJobs({
        status: status === 'all' ? undefined : status,
        skip,
        limit: pagination.pageSize,
      });
      setJobs(data.jobs);
      setTotal(data.total);
    } catch (error: any) {
      console.error('加载QC任务列表失败:', error);
    }
  }, [pagination]);

  // 加载所有配置选项
  const loadOptions = async () => {
    try {
      const [basisData, funcData, accuracyData, solventModelData, solventData, moleculesData, partitionsData] = await Promise.all([
        getBasisSets(),
        getFunctionals(),
        getAccuracyLevels(),
        getSolventModels(),
        getSolvents(),
        getCommonMolecules(),
        getPartitions().catch(() => []),  // 分区加载失败时使用空数组
      ]);
      setBasisSets(basisData.basis_sets);
      setFunctionals(funcData.functionals);
      setAccuracyLevels(accuracyData.levels);
      setSolventModels(solventModelData.models);
      setSolvents(solventData.solvents);
      setCommonMolecules(moleculesData.categories);
      setPartitions(partitionsData);

      // 构建SMILES下拉选项（包含charge信息用于自动判断分子类型）
      const options: { value: string; label: string; charge?: number }[] = [];
      moleculesData.categories.forEach((cat: MoleculeCategory) => {
        cat.molecules.forEach((mol: CommonMolecule) => {
          options.push({
            value: mol.smiles,
            label: `${mol.label || mol.name} (${mol.smiles})`,
            charge: mol.charge,
          });
        });
      });
      setSmilesOptions(options);
    } catch (error) {
      console.error('加载选项失败:', error);
    }
  };

  // 检查是否有活跃任务
  const hasActiveJobs = useCallback(() => {
    return jobs.some(job =>
      job.status === 'QUEUED' || job.status === 'RUNNING' || job.status === 'POSTPROCESSING'
    );
  }, [jobs]);

  // 获取默认分区
  const getDefaultPartition = () => {
    if (partitions.length > 0) {
      const upPartition = partitions.find(p => p.state === 'up');
      return upPartition?.name || partitions[0].name;
    }
    return 'cpu';
  };

  // 获取推荐配置
  const handleGetSuggestion = async () => {
    try {
      const suggestion = await getSlurmSuggestion({ job_type: 'qc' });
      form.setFieldsValue({
        slurm_partition: suggestion.partition,
        slurm_cpus: suggestion.cpus_per_task,
      });
      message.success(`已应用推荐配置: ${suggestion.reason}`);
    } catch (error: any) {
      message.error('获取推荐配置失败');
    }
  };

  // 计算自旋多重度
  const handleCalculateSpin = async () => {
    const smiles = form.getFieldValue('smiles');
    const charge = form.getFieldValue('charge') || 0;
    if (!smiles) {
      message.warning('请先输入SMILES');
      return;
    }
    try {
      const result = await calculateSpinMultiplicity(smiles, charge);
      setCalculatedSpin(result.spin_multiplicity);
      form.setFieldsValue({ spin_multiplicity: result.spin_multiplicity });
      message.success(result.description);
    } catch (error: any) {
      message.error('计算自旋多重度失败: ' + (error.response?.data?.detail || error.message));
    }
  };

  // 从分子名称中提取英文名（如果存在）
  // 例如: "碳酸乙烯酯(EC)" -> "EC", "1-乙基-3-甲基咪唑(EMIm+)" -> "EMIm+"
  // 如果没有括号中的英文名，则返回原名称
  const extractEnglishName = (name: string): string => {
    // 匹配括号中的英文名，支持中文括号和英文括号
    const match = name.match(/[（(]([A-Za-z0-9+\-_]+)[）)]/);
    if (match) {
      return match[1];
    }
    // 如果名称本身就是纯英文/数字，直接返回
    if (/^[A-Za-z0-9+\-_\s]+$/.test(name)) {
      return name;
    }
    // 否则返回原名称
    return name;
  };

  // 根据精度等级获取基础参数
  const getBaseParams = (accuracy: string) => {
    switch (accuracy) {
      case 'fast':
        return { functional: 'HF', basis_set: 'STO-3G' };
      case 'standard':
        return { functional: 'B3LYP', basis_set: '6-31G(d)' };
      case 'accurate':
        return { functional: 'B3LYP', basis_set: '6-311++G(d,p)' };
      default:
        return { functional: form.getFieldValue('functional') || 'B3LYP', basis_set: form.getFieldValue('basis_set') || '6-31G(d)' };
    }
  };

  // 根据分子类型调整参数（智能推荐核心逻辑 - 改进版）
  const adjustParamsForMoleculeType = (molType: string, baseParams: { functional: string; basis_set: string }, currentSolventModel: string) => {
    let { functional, basis_set } = baseParams;
    let solventModel = currentSolventModel;
    let reason = '';
    let changes: string[] = [];

    if (molType === 'anion') {
      // 阴离子需要弥散函数 - 在原基组基础上添加
      if (!basis_set.includes('+')) {
        // 升级基组：添加弥散函数
        if (basis_set === 'STO-3G') {
          basis_set = '6-31++G(d,p)'; // STO-3G 太小，直接升级
          changes.push('基组升级为 6-31++G(d,p)（包含弥散函数）');
        } else if (basis_set === '6-31G(d)') {
          basis_set = '6-31++G(d,p)';
          changes.push('基组升级为 6-31++G(d,p)（添加弥散函数）');
        } else if (basis_set === '6-31G(d,p)') {
          basis_set = '6-31++G(d,p)';
          changes.push('基组升级为 6-31++G(d,p)（添加弥散函数）');
        } else if (basis_set === '6-311G(d,p)') {
          basis_set = '6-311++G(d,p)';
          changes.push('基组升级为 6-311++G(d,p)（添加弥散函数）');
        }
      }

      // 阴离子在气相可能不稳定
      if (solventModel === 'gas') {
        solventModel = 'pcm';
        changes.push('溶剂模型切换为 PCM（阴离子在气相不稳定）');
      }

      reason = `阴离子具有扩展的电子密度分布。${changes.join('；')}`;
    } else if (molType === 'cation') {
      // 阳离子确保有极化函数
      if (!basis_set.includes('(d')) {
        if (basis_set === 'STO-3G') {
          basis_set = '6-31G(d,p)';
          changes.push('基组升级为 6-31G(d,p)（包含极化函数）');
        } else if (basis_set === '6-31G(d)') {
          basis_set = '6-31G(d,p)';
          changes.push('基组升级为 6-31G(d,p)（添加极化函数）');
        }
      }

      if (solventModel === 'gas') {
        solventModel = 'pcm';
        changes.push('溶剂模型切换为 PCM（阳离子在气相不稳定）');
      }

      reason = `阳离子电子结构紧凑。${changes.join('；')}`;
    } else if (molType === 'solvent') {
      reason = '中性分子使用标准计算参数';
    } else {
      reason = '自定义分子类型，使用指定参数';
    }

    return { functional, basis_set, solventModel, reason };
  };

  // 应用智能参数推荐
  const applyRecommendedParams = (molType: string, accuracy?: string) => {
    const currentAccuracy = accuracy || selectedAccuracy;
    const baseParams = getBaseParams(currentAccuracy);

    if (!useRecommendedParams) {
      // 不启用智能推荐时，只使用基础参数
      form.setFieldsValue({
        functional: baseParams.functional,
        basis_set: baseParams.basis_set,
      });
      setRecommendationReason('');
      return;
    }

    // 启用智能推荐：根据分子类型调整参数
    const adjusted = adjustParamsForMoleculeType(molType, baseParams, selectedSolventModel);

    // 更新表单
    form.setFieldsValue({
      functional: adjusted.functional,
      basis_set: adjusted.basis_set,
    });

    // 更新溶剂模型（仅当从gas切换到pcm/smd时才设置默认溶剂）
    if (adjusted.solventModel !== selectedSolventModel) {
      const wasGas = selectedSolventModel === 'gas';
      setSelectedSolventModel(adjusted.solventModel);
      // 只有从气相切换到隐式溶剂且当前没有设置溶剂名称时，才设置默认值
      if (wasGas && (adjusted.solventModel === 'pcm' || adjusted.solventModel === 'smd')) {
        const currentSolvent = form.getFieldValue('solvent_name');
        if (!currentSolvent) {
          form.setFieldsValue({ solvent_name: 'Water' });
        }
      }
    }

    setRecommendationReason(adjusted.reason);
  };

  // 选择常用分子
  const handleSelectMolecule = (molecule: CommonMolecule) => {
    // 根据电荷自动判断分子类型
    let moleculeType = 'custom';
    if (molecule.charge > 0.5) {
      moleculeType = 'cation';
    } else if (molecule.charge < -0.5) {
      moleculeType = 'anion';
    } else {
      moleculeType = 'solvent';
    }

    // 使用英文名称 (name字段)
    form.setFieldsValue({
      molecule_name: molecule.name,
      smiles: molecule.smiles,
      charge: molecule.charge,
      molecule_type: moleculeType,
    });
    // 自动计算自旋多重度
    if (autoSpin) {
      setTimeout(() => handleCalculateSpin(), 100);
    }
    // 应用智能参数推荐
    applyRecommendedParams(moleculeType);
    // 关闭 Collapse 面板
    setMoleculeCollapseKey([]);
  };

  // 分子类型变化时应用推荐参数
  const handleMoleculeTypeChange = (value: string) => {
    applyRecommendedParams(value);
  };

  // 精度等级变化
  const handleAccuracyChange = (value: string) => {
    setSelectedAccuracy(value);
    // 获取分子类型，重新应用推荐参数
    const molType = form.getFieldValue('molecule_type') || 'solvent';
    applyRecommendedParams(molType, value);
  };

  // 智能推荐开关变化
  const handleRecommendedParamsChange = (enabled: boolean) => {
    setUseRecommendedParams(enabled);
    const molType = form.getFieldValue('molecule_type') || 'solvent';
    if (enabled) {
      // 重新计算推荐参数
      const baseParams = getBaseParams(selectedAccuracy);
      const adjusted = adjustParamsForMoleculeType(molType, baseParams, selectedSolventModel);
      form.setFieldsValue({
        functional: adjusted.functional,
        basis_set: adjusted.basis_set,
      });
      if (adjusted.solventModel !== selectedSolventModel) {
        const wasGas = selectedSolventModel === 'gas';
        setSelectedSolventModel(adjusted.solventModel);
        // 只有从气相切换到隐式溶剂且当前没有设置溶剂名称时，才设置默认值
        if (wasGas && (adjusted.solventModel === 'pcm' || adjusted.solventModel === 'smd')) {
          const currentSolvent = form.getFieldValue('solvent_name');
          if (!currentSolvent) {
            form.setFieldsValue({ solvent_name: 'Water' });
          }
        }
      }
      setRecommendationReason(adjusted.reason);
    } else {
      // 恢复基础参数
      const baseParams = getBaseParams(selectedAccuracy);
      form.setFieldsValue({
        functional: baseParams.functional,
        basis_set: baseParams.basis_set,
      });
      setRecommendationReason('');
    }
  };

  useEffect(() => {
    loadJobs(activeTab);
    loadOptions();
  }, [activeTab, pagination, loadJobs]);

  // 轮询活跃任务
  useEffect(() => {
    if (hasActiveJobs()) {
      pollingRef.current = setInterval(() => {
        loadJobs(activeTab);
      }, 10000);
    }
    return () => {
      if (pollingRef.current) {
        clearInterval(pollingRef.current);
      }
    };
  }, [hasActiveJobs, loadJobs, activeTab]);

  // 创建或更新QC任务
  const handleCreateOrUpdate = async (values: any) => {
    setSubmitting(true);
    try {
      // 构建溶剂配置
      let solventConfig: SolventConfig | undefined;
      if (selectedSolventModel !== 'gas') {
        // 检查是否选择了自定义溶剂参数（在PCM/SMD模型下）
        const isCustomSolvent = values.solvent_name === 'custom';

        if (isCustomSolvent) {
          // 使用自定义溶剂参数
          solventConfig = {
            model: selectedSolventModel,
            solvent_name: values.custom_solvent_name || 'CustomSolvent',
            eps: values.custom_eps,
            eps_inf: values.custom_eps_inf,
          };
        } else if (selectedSolventModel === 'custom') {
          // 旧的自定义溶剂模型（保留兼容性）
          solventConfig = {
            model: selectedSolventModel,
            solvent_name: values.solvent_name,
            eps: values.eps,
            eps_inf: values.eps_inf,
            hbond_acidity: values.hbond_acidity,
            hbond_basicity: values.hbond_basicity,
            surface_tension: values.surface_tension,
            carbon_aromaticity: values.carbon_aromaticity,
            halogenicity: values.halogenicity,
          };
        } else {
          // 预设溶剂
          solventConfig = {
            model: selectedSolventModel,
            solvent_name: values.solvent_name,
          };
        }
      }

      // 获取精度等级对应的泛函和基组
      const accuracyLevel = accuracyLevels.find(l => l.value === selectedAccuracy);
      const functional = selectedAccuracy === 'custom' ? values.functional : (accuracyLevel?.functional || 'B3LYP');
      const basisSet = selectedAccuracy === 'custom' ? values.basis_set : (accuracyLevel?.basis_set || '6-31G(d)');

      const jobData: QCJobCreate = {
        molecule_name: values.molecule_name,
        smiles: values.smiles,
        molecule_type: values.molecule_type || 'solvent',
        basis_set: basisSet,
        functional: functional,
        charge: values.charge || 0,
        spin_multiplicity: values.spin_multiplicity || 1,
        accuracy_level: selectedAccuracy,
        solvent_config: solventConfig,
        auto_spin: autoSpin,
        config: {
          accuracy_level: selectedAccuracy,
          solvent_config: solventConfig,
          auto_spin: autoSpin,
        },
        // Slurm 资源配置
        slurm_partition: values.slurm_partition || 'cpu',
        slurm_cpus: values.slurm_cpus || 16,
        slurm_time: values.slurm_time || 7200,
      };

      if (editingJob) {
        // 更新现有任务
        await updateQCJob(editingJob.id, jobData);
        message.success('QC任务更新成功');
        setModalVisible(false);
        setEditingJob(null);
        form.resetFields();
        setSelectedAccuracy('standard');
        setSelectedSolventModel('gas');
        loadJobs(activeTab);
      } else {
        // 直接创建任务，后端会进行查重检查
        await createQCJob(jobData);
        message.success('QC任务创建成功');
        setModalVisible(false);
        setEditingJob(null);
        form.resetFields();
        setSelectedAccuracy('standard');
        setSelectedSolventModel('gas');
        loadJobs(activeTab);
      }
    } catch (error: any) {
      // 处理409冲突错误（重复计算检测）
      if (error.response?.status === 409) {
        const errorDetail = error.response?.data?.detail;
        const existingJobId = errorDetail?.existing_job_id;
        const errorMessage = typeof errorDetail === 'string' ? errorDetail : errorDetail?.message;

        Modal.warning({
          title: '检测到重复计算',
          width: 600,
          content: (
            <div>
              <Alert
                message="已存在完全相同参数的QC计算"
                description={errorMessage || '为避免浪费计算资源，请直接查看已有结果。'}
                type="warning"
                showIcon
                style={{ marginBottom: 16 }}
              />
              {existingJobId && (
                <p style={{ marginTop: 12 }}>
                  <strong>已有任务ID:</strong> {existingJobId}
                </p>
              )}
            </div>
          ),
          okText: '查看已有结果',
          onOk: () => {
            if (existingJobId) {
              navigate(`/workspace/qc-jobs/${existingJobId}`);
              setModalVisible(false);
              setEditingJob(null);
              form.resetFields();
            }
          },
        });
      } else {
        message.error(error.response?.data?.detail || (editingJob ? '更新失败' : '创建失败'));
      }
    } finally {
      setSubmitting(false);
    }
  };

  // 打开编辑弹窗
  const handleEdit = (job: QCJob) => {
    setEditingJob(job);
    // 设置溶剂模型
    const solventModel = job.solvent_config?.model || 'gas';
    setSelectedSolventModel(solventModel);
    // 设置精度等级
    setSelectedAccuracy(job.accuracy_level || 'standard');
    // 填充表单
    form.setFieldsValue({
      molecule_name: job.molecule_name,
      smiles: job.smiles,
      molecule_type: job.molecule_type || 'solvent',
      charge: job.charge,
      spin_multiplicity: job.spin_multiplicity,
      functional: job.functional,
      basis_set: job.basis_set,
      solvent_model: solventModel,
      solvent_name: job.solvent_config?.solvent_name,
      slurm_partition: job.slurm_partition || 'cpu',
      slurm_cpus: job.slurm_cpus || 16,
      slurm_time: job.slurm_time || 7200,
    });
    setModalVisible(true);
  };

  // 提交任务到集群
  const handleSubmit = async (jobId: number) => {
    try {
      await submitQCJob(jobId);
      message.success('任务已提交到计算集群');
      loadJobs(activeTab);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '提交失败');
    }
  };

  // 删除/取消任务
  const handleDelete = async (jobId: number) => {
    try {
      await deleteQCJob(jobId);
      message.success('任务已删除');
      loadJobs(activeTab);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  // 批量删除
  const handleBatchDelete = async () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择要删除的任务');
      return;
    }

    Modal.confirm({
      title: '批量删除确认',
      content: `确定要删除/取消选中的 ${selectedRowKeys.length} 个任务吗？运行中的任务将被取消。`,
      okText: '确定删除',
      okButtonProps: { danger: true },
      cancelText: '取消',
      onOk: async () => {
        try {
          const result = await batchDeleteQCJobs(selectedRowKeys as number[]);
          message.success(result.message);
          setSelectedRowKeys([]);
          loadJobs(activeTab);
        } catch (error: any) {
          message.error(error.response?.data?.detail || '批量删除失败');
        }
      },
    });
  };

  // 批量提交
  const handleBatchSubmit = async () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择要提交的任务');
      return;
    }

    Modal.confirm({
      title: '批量提交确认',
      content: `确定要提交选中的 ${selectedRowKeys.length} 个任务到计算集群吗？`,
      okText: '确定提交',
      cancelText: '取消',
      onOk: async () => {
        try {
          const result = await batchSubmitQCJobs(selectedRowKeys as number[]);
          if (result.success_count > 0) {
            message.success(result.message);
          }
          if (result.failed_count > 0 && result.errors.length > 0) {
            Modal.warning({
              title: '部分任务提交失败',
              content: (
                <div>
                  {result.errors.map((err, i) => (
                    <div key={i}>任务ID {err.job_id}: {err.error}</div>
                  ))}
                </div>
              ),
            });
          }
          setSelectedRowKeys([]);
          loadJobs(activeTab);
        } catch (error: any) {
          message.error(error.response?.data?.detail || '批量提交失败');
        }
      },
    });
  };

  // 批量取消
  const handleBatchCancel = async () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择要取消的任务');
      return;
    }

    Modal.confirm({
      title: '批量取消确认',
      content: `确定要取消选中的 ${selectedRowKeys.length} 个任务吗？`,
      okText: '确定取消',
      okButtonProps: { danger: true },
      cancelText: '返回',
      onOk: async () => {
        try {
          const result = await batchCancelQCJobs(selectedRowKeys as number[]);
          if (result.success_count > 0) {
            message.success(result.message);
          }
          if (result.failed_count > 0 && result.errors.length > 0) {
            Modal.warning({
              title: '部分任务取消失败',
              content: (
                <div>
                  {result.errors.map((err, i) => (
                    <div key={i}>任务ID {err.job_id}: {err.error}</div>
                  ))}
                </div>
              ),
            });
          }
          setSelectedRowKeys([]);
          loadJobs(activeTab);
        } catch (error: any) {
          message.error(error.response?.data?.detail || '批量取消失败');
        }
      },
    });
  };

  // 克隆任务（基于已有任务创建新任务，用于重新计算或修改参数）
  const handleCloneJob = (job: QCJob) => {
    // 设置溶剂模型
    const solventConfig = job.solvent_config || job.config?.solvent_config;
    const solventModel = solventConfig?.model || 'gas';
    setSelectedSolventModel(solventModel);
    // 设置精度等级
    setSelectedAccuracy(job.accuracy_level || 'custom');
    // 填充表单，但不设置editingJob，这样会创建新任务而不是更新
    form.setFieldsValue({
      molecule_name: job.molecule_name,
      smiles: job.smiles,
      molecule_type: job.molecule_type || 'solvent',
      charge: job.charge,
      spin_multiplicity: job.spin_multiplicity,
      functional: job.functional,
      basis_set: job.basis_set,
      solvent_model: solventModel,
      solvent_name: solventConfig?.solvent_name,
      slurm_partition: job.slurm_partition || 'cpu',
      slurm_cpus: job.slurm_cpus || 16,
      slurm_time: job.slurm_time || 7200,
    });
    setEditingJob(null);  // 确保是创建新任务
    setModalVisible(true);
    message.info('已加载任务参数，可修改后创建新任务');
  };

  // 批量导入处理函数
  const handleDownloadTemplate = async () => {
    try {
      const blob = await downloadTemplate('qc', true);
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'QC计算批量导入模板.xlsx';
      document.body.appendChild(a);
      a.click();
      window.URL.revokeObjectURL(url);
      document.body.removeChild(a);
      message.success('模板下载成功');
    } catch (error) {
      message.error('模板下载失败');
    }
  };

  const handleBatchImport = async () => {
    if (!uploadFile) {
      message.error('请先选择文件');
      return;
    }

    setImporting(true);
    try {
      const result = await batchImportUpload(
        uploadFile,
        undefined,
        undefined,
        undefined,
        'QC计算'
      );
      setImportResult(result);
      setImportStep(1); // 跳转到结果页面
      loadJobs(activeTab);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '批量导入失败');
    } finally {
      setImporting(false);
    }
  };

  const resetBatchImport = () => {
    setImportStep(0);
    setUploadFile(null);
    setImportResult(null);
    setBatchImportVisible(false);
  };

  // 表格列定义
  const columns: ColumnsType<QCJob> = [
    {
      title: '分子名称',
      dataIndex: 'molecule_name',
      key: 'molecule_name',
      width: 150,
      ellipsis: true,
    },
    {
      title: 'SMILES',
      dataIndex: 'smiles',
      key: 'smiles',
      width: 200,
      ellipsis: true,
      render: (smiles: string) => (
        <Tooltip title={smiles}>
          <Text copyable={{ text: smiles }} style={{ maxWidth: 180 }}>
            {smiles.length > 25 ? `${smiles.slice(0, 25)}...` : smiles}
          </Text>
        </Tooltip>
      ),
    },
    {
      title: '类型',
      dataIndex: 'molecule_type',
      key: 'molecule_type',
      width: 80,
      render: (type: string) => moleculeTypeMap[type] || type,
    },
    {
      title: '基组',
      dataIndex: 'basis_set',
      key: 'basis_set',
      width: 120,
    },
    {
      title: '泛函',
      dataIndex: 'functional',
      key: 'functional',
      width: 80,
    },
    {
      title: '溶剂模型',
      key: 'solvent_model',
      width: 140,
      render: (_: unknown, record: QCJob) => {
        const solventConfig = record.solvent_config || record.config?.solvent_config;
        if (!solventConfig || solventConfig.model === 'gas') {
          return <Tag>气相</Tag>;
        }
        const model = solventConfig.model?.toUpperCase() || 'N/A';
        const solventName = solventConfig.solvent_name || '';
        return (
          <Space size={4}>
            <Tag color="blue">{model}</Tag>
            {solventName && <span style={{ fontSize: '12px', color: '#666' }}>{solventName}</span>}
          </Space>
        );
      },
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      width: 100,
      render: (status: string) => {
        const { color, text } = statusMap[status] || { color: 'default', text: status };
        return <Tag color={color}>{text}</Tag>;
      },
    },
    {
      title: '进度',
      dataIndex: 'progress',
      key: 'progress',
      width: 80,
      render: (progress: number) => `${progress.toFixed(0)}%`,
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 160,
      render: (time: string) => new Date(time).toLocaleString('zh-CN'),
    },
    {
      title: '操作',
      key: 'actions',
      width: 220,
      fixed: 'right',
      render: (_, record) => (
        <Space size="small">
          {record.status === 'CREATED' && (
            <>
              <Tooltip title="编辑">
                <Button
                  type="link"
                  size="small"
                  icon={<EditOutlined />}
                  onClick={() => handleEdit(record)}
                />
              </Tooltip>
              <Tooltip title="提交任务">
                <Button
                  type="link"
                  size="small"
                  icon={<PlayCircleOutlined />}
                  onClick={() => handleSubmit(record.id)}
                />
              </Tooltip>
            </>
          )}
          {(record.status === 'COMPLETED' || record.status === 'FAILED') && (
            <Tooltip title="基于此任务创建新任务（可修改参数重新计算）">
              <Button
                type="link"
                size="small"
                icon={<CopyOutlined />}
                onClick={() => handleCloneJob(record)}
              />
            </Tooltip>
          )}
          <Tooltip title="查看详情">
            <Button
              type="link"
              size="small"
              icon={<EyeOutlined />}
              onClick={() => navigate(`/workspace/qc-jobs/${record.id}`)}
            />
          </Tooltip>
          <Popconfirm
            title="确定要删除这个任务吗？"
            onConfirm={() => handleDelete(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Tooltip title="删除">
              <Button type="link" size="small" danger icon={<DeleteOutlined />} />
            </Tooltip>
          </Popconfirm>
        </Space>
      ),
    },
  ];

  // Tab项
  const tabItems = [
    { key: 'all', label: '全部' },
    { key: 'CREATED', label: '待提交' },
    { key: 'QUEUED', label: '排队中' },
    { key: 'RUNNING', label: '运行中' },
    { key: 'COMPLETED', label: '已完成' },
    { key: 'FAILED', label: '失败' },
  ];

  return (
    <div style={{ padding: 24, position: 'relative', background: '#f5f7fb', minHeight: 'calc(100vh - 64px)' }}>
      {/* 浮动批量操作栏 */}
      {selectedRowKeys.length > 0 && (
        <div style={{
          position: 'fixed',
          bottom: 32,
          left: '50%',
          transform: 'translateX(-50%)',
          zIndex: 1000,
          background: '#fff',
          padding: '16px 24px',
          borderRadius: 16,
          boxShadow: '0 8px 32px rgba(0, 0, 0, 0.12), 0 2px 8px rgba(0, 0, 0, 0.08)',
          border: '1px solid #e8e8e8',
          display: 'flex',
          alignItems: 'center',
          gap: 16,
          animation: 'slideUp 0.3s ease-out',
        }}>
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: 8,
            paddingRight: 16,
            borderRight: '1px solid #e8e8e8',
          }}>
            <CheckSquareOutlined style={{ fontSize: 18, color: '#1677ff' }} />
            <Text strong style={{ fontSize: 14 }}>
              已选择 <span style={{ color: '#1677ff', fontSize: 16 }}>{selectedRowKeys.length}</span> 个任务
            </Text>
          </div>
          <Space size={12}>
            <Button
              type="primary"
              icon={<PlayCircleOutlined />}
              onClick={handleBatchSubmit}
              style={{ borderRadius: 8 }}
            >
              批量提交
            </Button>
            <Button
              icon={<StopOutlined />}
              onClick={handleBatchCancel}
              style={{ borderRadius: 8 }}
            >
              批量取消
            </Button>
            <Button
              danger
              icon={<DeleteOutlined />}
              onClick={handleBatchDelete}
              style={{ borderRadius: 8 }}
            >
              批量删除
            </Button>
            <Button
              onClick={() => setSelectedRowKeys([])}
              style={{ borderRadius: 8 }}
            >
              取消选择
            </Button>
          </Space>
        </div>
      )}
      <style>{`
        @keyframes slideUp {
          from {
            opacity: 0;
            transform: translateX(-50%) translateY(20px);
          }
          to {
            opacity: 1;
            transform: translateX(-50%) translateY(0);
          }
        }
      `}</style>

      {/* 页面标题区域 */}
      <div style={{ marginBottom: 24 }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
              <ExperimentOutlined style={{ marginRight: 12, color: '#1890ff' }} />
              量子化学计算
            </Title>
            <Text type="secondary">管理QC量子化学计算任务，进行精确能量和性质计算</Text>
          </div>
          <Space>
            <Button icon={<ReloadOutlined />} onClick={() => loadJobs(activeTab)} style={{ borderRadius: 8 }}>
              刷新
            </Button>
            <Button icon={<UploadOutlined />} onClick={() => setBatchImportVisible(true)} style={{ borderRadius: 8 }}>
              批量导入
            </Button>
            <Button type="primary" icon={<PlusOutlined />} onClick={() => setModalVisible(true)} style={{ borderRadius: 8 }}>
              新建QC任务
            </Button>
          </Space>
        </div>
      </div>

      <Card>

        <Tabs
          activeKey={activeTab}
          onChange={setActiveTab}
          items={tabItems}
          style={{ marginBottom: 16 }}
        />

        <Spin spinning={loading}>
          {jobs.length === 0 ? (
            <Empty description="暂无QC计算任务" />
          ) : (
            <Table
              columns={columns}
              dataSource={jobs}
              rowKey="id"
              rowSelection={{
                selectedRowKeys,
                onChange: (keys) => setSelectedRowKeys(keys),
              }}
              pagination={{
                current: pagination.current,
                pageSize: pagination.pageSize,
                total,
                showSizeChanger: true,
                showTotal: (total) => `共 ${total} 条`,
                onChange: (page, pageSize) => setPagination({ current: page, pageSize }),
              }}
              scroll={{ x: 1200 }}
            />
          )}
        </Spin>
      </Card>

      {/* 创建/编辑QC任务弹窗 */}
      <Modal
        title={
          <Space>
            <ExperimentOutlined style={{ color: '#722ed1' }} />
            <span style={{ fontWeight: 600 }}>{editingJob ? '编辑QC计算任务' : '新建QC计算任务'}</span>
          </Space>
        }
        open={modalVisible}
        onCancel={() => {
          setModalVisible(false);
          setEditingJob(null);
          form.resetFields();
          setSelectedAccuracy('standard');
          setSelectedSolventModel('gas');
        }}
        footer={null}
        width={800}
        centered
        styles={{
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' }
        }}
      >
        <Form
          form={form}
          layout="vertical"
          onFinish={handleCreateOrUpdate}
          initialValues={{
            basis_set: '6-31G(d)',
            functional: 'B3LYP',
            charge: 0,
            spin_multiplicity: 1,
            molecule_type: 'solvent',
          }}
        >
          {/* 分子信息 */}
          <Divider orientation="left">
            <Space>
              <ExperimentOutlined style={{ color: '#722ed1' }} />
              分子信息
            </Space>
          </Divider>

          {/* 常用分子选择 */}
          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#d3adf7',
              background: 'linear-gradient(135deg, #f9f0ff 0%, #fff 100%)'
            }}
          >
            <Collapse
              ghost
              style={{ margin: -12 }}
              activeKey={moleculeCollapseKey}
              onChange={(keys) => setMoleculeCollapseKey(keys as string[])}
            >
              <Panel
                header={
                  <Text style={{ color: '#722ed1' }}>
                    <ThunderboltOutlined style={{ marginRight: 8 }} />
                    快速选择常用分子
                  </Text>
                }
                key="molecules"
              >
                <Row gutter={[8, 12]}>
                  {commonMolecules.map(category => (
                    <Col span={24} key={category.name}>
                      <Text strong style={{ marginBottom: 8, display: 'block', color: '#595959' }}>
                        {category.name}
                      </Text>
                      <Space wrap>
                        {category.molecules.map(mol => (
                          <Tooltip key={mol.smiles} title={`SMILES: ${mol.smiles}`}>
                            <Tag
                              color="purple"
                              style={{ cursor: 'pointer', borderRadius: 4 }}
                              onClick={() => handleSelectMolecule(mol)}
                            >
                              {mol.label || mol.name}
                            </Tag>
                          </Tooltip>
                        ))}
                      </Space>
                    </Col>
                  ))}
                </Row>
              </Panel>
            </Collapse>
          </Card>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                name="molecule_name"
                label="分子名称"
                rules={[
                  { required: true, message: '请输入分子名称' },
                  {
                    pattern: /^[A-Za-z0-9+\-_\s]+$/,
                    message: '分子名称仅支持英文、数字和符号，不支持中文'
                  }
                ]}
                tooltip="请使用英文名称，例如: EC, DMC, LiPF6"
                extra="仅支持英文、数字和符号（+、-、_），不支持中文"
              >
                <Input placeholder="例如: EC, DMC, LiPF6" />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item name="molecule_type" label="分子类型">
                <Select onChange={handleMoleculeTypeChange}>
                  <Select.Option value="solvent">溶剂</Select.Option>
                  <Select.Option value="cation">阳离子</Select.Option>
                  <Select.Option value="anion">阴离子</Select.Option>
                  <Select.Option value="custom">自定义</Select.Option>
                </Select>
              </Form.Item>
            </Col>
          </Row>

          <Form.Item
            name="smiles"
            label="SMILES表达式"
            rules={[{ required: true, message: '请输入SMILES' }]}
            extra="可输入自定义SMILES或从下拉列表选择常用分子"
          >
            <AutoComplete
              options={smilesOptions}
              placeholder="输入或选择SMILES，例如: C1COC(=O)O1 (EC)"
              filterOption={(inputValue, option) =>
                option!.label.toLowerCase().indexOf(inputValue.toLowerCase()) !== -1 ||
                option!.value.toLowerCase().indexOf(inputValue.toLowerCase()) !== -1
              }
              onSelect={(value: string) => {
                // 选择时自动填充分子名称（提取英文名）
                const mol = smilesOptions.find(opt => opt.value === value);
                if (mol) {
                  const name = mol.label.split(' (')[0];
                  const englishName = extractEnglishName(name);

                  // 根据电荷自动判断分子类型
                  let moleculeType = 'custom';
                  if (mol.charge !== undefined) {
                    if (mol.charge > 0.5) {
                      moleculeType = 'cation';
                    } else if (mol.charge < -0.5) {
                      moleculeType = 'anion';
                    } else {
                      moleculeType = 'solvent';
                    }
                  }

                  form.setFieldsValue({
                    molecule_name: englishName,
                    molecule_type: moleculeType,
                    charge: mol.charge,
                  });
                  // 自动计算自旋
                  if (autoSpin) {
                    setTimeout(() => handleCalculateSpin(), 100);
                  }
                  // 应用智能参数推荐
                  applyRecommendedParams(moleculeType);
                }
              }}
            />
          </Form.Item>

          <Divider orientation="left">
            <Space>
              <ThunderboltOutlined style={{ color: '#fa8c16' }} />
              计算精度
            </Space>
          </Divider>

          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#ffd591',
              background: 'linear-gradient(135deg, #fff7e6 0%, #fff 100%)',
              boxShadow: '0 2px 8px rgba(0,0,0,0.08)'
            }}
          >
            <Form.Item label="精度等级" style={{ marginBottom: 16 }}>
              <Select
                value={selectedAccuracy}
                onChange={handleAccuracyChange}
                style={{ width: '100%' }}
              >
                {accuracyLevels.map(level => (
                  <Select.Option key={level.value} value={level.value}>
                    <Space>
                      <span>{level.label}</span>
                      <Text type="secondary" style={{ fontSize: 12 }}>
                        {level.description}
                      </Text>
                    </Space>
                  </Select.Option>
                ))}
              </Select>
            </Form.Item>

            {/* 精度等级信息卡片 */}
            {selectedAccuracy !== 'custom' && (
              <div style={{ padding: '12px', background: '#fff', borderRadius: 6, border: '1px solid #ffe7ba', marginBottom: 12 }}>
                <Row gutter={16}>
                  <Col span={8}>
                    <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>预计耗时</div>
                    <Tag color="orange" style={{ fontSize: 12, padding: '4px 12px' }}>
                      {accuracyLevels.find(l => l.value === selectedAccuracy)?.estimated_time || '-'}
                    </Tag>
                  </Col>
                  <Col span={8}>
                    <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>推荐泛函</div>
                    <Tag color="blue" style={{ fontSize: 12, padding: '4px 12px' }}>
                      {selectedAccuracy === 'fast' ? 'HF' : 'B3LYP'}
                    </Tag>
                  </Col>
                  <Col span={8}>
                    <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>推荐基组</div>
                    <Tag color="green" style={{ fontSize: 12, padding: '4px 12px' }}>
                      {selectedAccuracy === 'fast' ? 'STO-3G' : selectedAccuracy === 'accurate' ? '6-311++G(d,p)' : '6-31G(d)'}
                    </Tag>
                  </Col>
                </Row>
              </div>
            )}

            {/* 自定义精度等级 */}
            {selectedAccuracy === 'custom' && (
              <div style={{ padding: '12px', background: '#fff', borderRadius: 6, border: '1px solid #ffe7ba' }}>
                <div style={{ fontSize: 12, color: '#666', marginBottom: 12, fontWeight: 500 }}>自定义计算参数</div>
                <Row gutter={16}>
                  <Col span={12}>
                    <Form.Item name="functional" label="泛函" rules={[{ required: true, message: '请选择泛函' }]} style={{ marginBottom: 0 }}>
                      <Select placeholder="选择泛函">
                        {functionals.map(f => (
                          <Select.Option key={f.value} value={f.value}>
                            <Tooltip title={f.description}>{f.label}</Tooltip>
                          </Select.Option>
                        ))}
                      </Select>
                    </Form.Item>
                  </Col>
                  <Col span={12}>
                    <Form.Item name="basis_set" label="基组" rules={[{ required: true, message: '请选择基组' }]} style={{ marginBottom: 0 }}>
                      <Select placeholder="选择基组">
                        {basisSets.map(bs => (
                          <Select.Option key={bs.value} value={bs.value}>
                            <Tooltip title={bs.description}>{bs.label}</Tooltip>
                          </Select.Option>
                        ))}
                      </Select>
                    </Form.Item>
                  </Col>
                </Row>
              </div>
            )}
          </Card>

          <Divider orientation="left">
            <Space>
              <SettingOutlined style={{ color: '#13c2c2' }} />
              溶剂环境
            </Space>
          </Divider>

          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#87e8de',
              background: 'linear-gradient(135deg, #e6fffb 0%, #fff 100%)',
              boxShadow: '0 2px 8px rgba(0,0,0,0.08)'
            }}
          >
            <Form.Item label="溶剂模型" style={{ marginBottom: 16 }}>
              <Select
                value={selectedSolventModel}
                onChange={setSelectedSolventModel}
                style={{ width: '100%' }}
              >
                {solventModels.map(model => (
                  <Select.Option key={model.value} value={model.value}>
                    <Space>
                      <span>{model.label}</span>
                      <Text type="secondary" style={{ fontSize: 12 }}>{model.description}</Text>
                    </Space>
                  </Select.Option>
                ))}
              </Select>
            </Form.Item>

            {(selectedSolventModel === 'pcm' || selectedSolventModel === 'smd') && (
              <Form.Item name="solvent_name" label="溶剂" style={{ marginBottom: 12 }}>
                <Select
                  showSearch
                  placeholder="选择溶剂"
                  optionFilterProp="children"
                >
                  <Select.OptGroup label="常用溶剂">
                    <Select.Option value="Water">水 (Water, ε=78.4)</Select.Option>
                    <Select.Option value="Acetonitrile">乙腈 (Acetonitrile, ε=35.7)</Select.Option>
                    <Select.Option value="DMSO">二甲亚砜 (DMSO, ε=46.8)</Select.Option>
                    <Select.Option value="Methanol">甲醇 (Methanol, ε=32.6)</Select.Option>
                    <Select.Option value="Ethanol">乙醇 (Ethanol, ε=24.9)</Select.Option>
                  </Select.OptGroup>
                  <Select.OptGroup label="电解液常用">
                    <Select.Option value="Acetone">丙酮 (Acetone, ε=20.5)</Select.Option>
                    <Select.Option value="Dichloromethane">二氯甲烷 (DCM, ε=8.9)</Select.Option>
                    <Select.Option value="THF">四氢呋喃 (THF, ε=7.4)</Select.Option>
                    <Select.Option value="DiethylEther">乙醚 (ε=4.2)</Select.Option>
                    <Select.Option value="PropyleneCarbonate">碳酸丙烯酯 (PC, ε=64.9)</Select.Option>
                    <Select.Option value="DimethylCarbonate">碳酸二甲酯 (DMC, ε=3.1)</Select.Option>
                    <Select.Option value="EthylMethylCarbonate">碳酸甲乙酯 (EMC, ε=2.9)</Select.Option>
                    <Select.Option value="EthyleneCarbonate">碳酸乙烯酯 (EC, ε=90.5)</Select.Option>
                    <Select.Option value="DiethylCarbonate">碳酸二乙酯 (DEC, ε=2.8)</Select.Option>
                    <Select.Option value="GammaButyrolactone">γ-丁内酯 (GBL, ε=42.0)</Select.Option>
                    <Select.Option value="12Dimethoxyethane">1,2-二甲氧基乙烷 (DME, ε=7.2)</Select.Option>
                  </Select.OptGroup>
                  <Select.OptGroup label="自定义">
                    <Select.Option value="custom">自定义溶剂参数...</Select.Option>
                  </Select.OptGroup>
                </Select>
              </Form.Item>
            )}

            {/* PCM/SMD模型下的自定义溶剂参数 */}
            {(selectedSolventModel === 'pcm' || selectedSolventModel === 'smd') && (
              <Form.Item
                noStyle
                shouldUpdate={(prevValues, currentValues) => prevValues.solvent_name !== currentValues.solvent_name}
              >
                {({ getFieldValue }) =>
                  getFieldValue('solvent_name') === 'custom' ? (
                    <Card size="small" style={{ background: '#fff', marginBottom: 0 }}>
                      <Text type="secondary" style={{ display: 'block', marginBottom: 12 }}>
                        自定义溶剂参数
                      </Text>
                      <Row gutter={[8, 8]}>
                        <Col span={8}>
                          <Form.Item
                            name="custom_eps"
                            label="介电常数 ε"
                            rules={[{ required: true, message: '请输入介电常数' }]}
                            style={{ marginBottom: 8 }}
                            tooltip="溶剂的静态介电常数，如水为78.4"
                          >
                            <InputNumber min={1} max={200} step={0.1} style={{ width: '100%' }} placeholder="例如: 78.4" />
                          </Form.Item>
                        </Col>
                        <Col span={8}>
                          <Form.Item
                            name="custom_eps_inf"
                            label="光学介电常数 ε∞"
                            style={{ marginBottom: 8 }}
                            tooltip="溶剂的光学介电常数（可选），默认为1.0"
                          >
                            <InputNumber min={1} max={10} step={0.01} style={{ width: '100%' }} placeholder="例如: 1.78" />
                          </Form.Item>
                        </Col>
                        <Col span={8}>
                          <Form.Item
                            name="custom_solvent_name"
                            label="溶剂名称"
                            style={{ marginBottom: 8 }}
                            tooltip="自定义溶剂的名称（用于记录）"
                          >
                            <Input placeholder="例如: 高浓LiTFSI" />
                          </Form.Item>
                        </Col>
                      </Row>
                    </Card>
                  ) : null
                }
              </Form.Item>
            )}

            {selectedSolventModel === 'custom' && (
              <Card size="small" style={{ background: '#fff' }}>
                <Text type="secondary" style={{ display: 'block', marginBottom: 12 }}>
                  自定义溶剂参数
                </Text>
                <Row gutter={[8, 8]}>
                  <Col span={8}>
                    <Form.Item name="eps" label="介电常数 ε" style={{ marginBottom: 8 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="78.35" />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="eps_inf" label="光学介电常数 n²" style={{ marginBottom: 8 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="1.778" />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="hbond_acidity" label="氢键酸度 α" style={{ marginBottom: 8 }}>
                      <InputNumber style={{ width: '100%' }} min={0} max={1} step={0.01} placeholder="0.82" />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="hbond_basicity" label="氢键碱度 β" style={{ marginBottom: 8 }}>
                      <InputNumber style={{ width: '100%' }} min={0} max={1} step={0.01} placeholder="0.35" />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="surface_tension" label="表面张力 γ" style={{ marginBottom: 8 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="71.99" />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="carbon_aromaticity" label="芳香碳比例 φ" style={{ marginBottom: 8 }}>
                      <InputNumber style={{ width: '100%' }} min={0} max={1} step={0.01} placeholder="0.0" />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="halogenicity" label="卤素比例 ψ" style={{ marginBottom: 0 }}>
                      <InputNumber style={{ width: '100%' }} min={0} max={1} step={0.01} placeholder="0.0" />
                    </Form.Item>
                  </Col>
                </Row>
              </Card>
            )}
          </Card>

          {/* 智能参数推荐 - 重构版 */}
          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#b7eb8f',
              background: 'linear-gradient(135deg, #f6ffed 0%, #fff 100%)',
              boxShadow: '0 2px 8px rgba(0,0,0,0.08)'
            }}
          >
            <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: 16 }}>
              <Space>
                <BulbOutlined style={{ color: '#52c41a', fontSize: 18 }} />
                <Text strong style={{ fontSize: 14 }}>智能参数推荐</Text>
              </Space>
              <Checkbox
                checked={useRecommendedParams}
                onChange={(e) => handleRecommendedParamsChange(e.target.checked)}
              >
                启用智能推荐
              </Checkbox>
            </div>

            {/* 显示推荐参数详情 */}
            <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
              prevValues.molecule_type !== currentValues.molecule_type ||
              prevValues.accuracy_level !== currentValues.accuracy_level ||
              prevValues.solvent_model !== currentValues.solvent_model ||
              prevValues.functional !== currentValues.functional ||
              prevValues.basis_set !== currentValues.basis_set
            }>
              {({ getFieldValue }) => {
                const moleculeType = getFieldValue('molecule_type') || 'solvent';
                const accuracyLevel = getFieldValue('accuracy_level') || 'standard';
                const solventModel = getFieldValue('solvent_model') || 'gas';
                const currentFunctional = getFieldValue('functional') || 'B3LYP';
                const currentBasisSet = getFieldValue('basis_set') || '6-31G(d)';

                // 获取基础参数
                const getBaseParamsForDisplay = (level: string) => {
                  switch (level) {
                    case 'fast': return { basis_set: 'STO-3G', functional: 'HF' };
                    case 'standard': return { basis_set: '6-31G(d)', functional: 'B3LYP' };
                    case 'accurate': return { basis_set: '6-311++G(d,p)', functional: 'B3LYP' };
                    default: return {
                      basis_set: currentBasisSet,
                      functional: currentFunctional
                    };
                  }
                };

                // 获取推荐参数
                const getRecommendedDisplay = () => {
                  const baseParams = getBaseParamsForDisplay(accuracyLevel);
                  let params = {
                    functional: baseParams.functional,
                    basis_set: baseParams.basis_set,
                    solvent_model: solventModel,
                    reason: ''
                  };

                  if (!useRecommendedParams) {
                    return params;
                  }

                  if (moleculeType === 'anion') {
                    // 阴离子需要弥散函数
                    if (!params.basis_set.includes('+')) {
                      params.basis_set = accuracyLevel === 'accurate' ? '6-311++G(d,p)' : '6-31++G(d,p)';
                    }
                    params.reason = '阴离子具有扩展的电子密度分布，需要使用带弥散函数(++)的基组以准确描述。';
                    if (params.solvent_model === 'gas') {
                      params.solvent_model = 'pcm';
                      params.reason += '阴离子在气相中可能不稳定，已自动切换到PCM隐式溶剂模型。';
                    }
                  } else if (moleculeType === 'cation') {
                    // 阳离子确保有极化函数
                    if (!params.basis_set.includes('(d')) {
                      params.basis_set = '6-31G(d,p)';
                    }
                    params.reason = '阳离子电子结构较为紧凑，使用带极化函数(d,p)的基组即可获得良好精度。';
                    if (params.solvent_model === 'gas') {
                      params.solvent_model = 'pcm';
                      params.reason += '已自动切换到PCM隐式溶剂模型。';
                    }
                  } else {
                    params.reason = '中性分子使用标准计算参数。';
                  }

                  return params;
                };

                const recommended = getRecommendedDisplay();
                const moleculeTypeLabel = moleculeType === 'solvent' ? '溶剂' : moleculeType === 'cation' ? '阳离子' : '阴离子';
                const moleculeTypeColor = moleculeType === 'solvent' ? 'blue' : moleculeType === 'cation' ? 'green' : 'orange';

                // 判断是否有变化
                const hasChanges = useRecommendedParams && (
                  recommended.basis_set !== currentBasisSet ||
                  recommended.solvent_model !== solventModel
                );

                return (
                  <div>
                    {/* 分子类型和状态 */}
                    <div style={{ marginBottom: 16, padding: '12px', background: '#fff', borderRadius: 6, border: '1px solid #e8f5e9' }}>
                      <Row gutter={16}>
                        <Col span={12}>
                          <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>分子类型</div>
                          <Tag color={moleculeTypeColor} style={{ fontSize: 12, padding: '4px 12px' }}>
                            {moleculeTypeLabel}
                          </Tag>
                        </Col>
                        <Col span={12}>
                          <div style={{ fontSize: 12, color: '#666', marginBottom: 4 }}>推荐状态</div>
                          {useRecommendedParams ? (
                            <Tag color={hasChanges ? 'green' : 'default'} icon={<BulbOutlined />} style={{ fontSize: 12, padding: '4px 12px' }}>
                              {hasChanges ? '已优化' : '无需调整'}
                            </Tag>
                          ) : (
                            <Tag color="default" style={{ fontSize: 12, padding: '4px 12px' }}>已禁用</Tag>
                          )}
                        </Col>
                      </Row>
                    </div>

                    {/* 参数对比 */}
                    {useRecommendedParams && (
                      <div style={{ marginBottom: 16 }}>
                        <div style={{ fontSize: 12, color: '#666', marginBottom: 8, fontWeight: 500 }}>参数对比</div>
                        <Row gutter={12}>
                          {/* 泛函 */}
                          <Col span={8}>
                            <div style={{ padding: '8px', background: '#f5f5f5', borderRadius: 4, border: '1px solid #e8e8e8' }}>
                              <div style={{ fontSize: 11, color: '#999', marginBottom: 4 }}>泛函</div>
                              <div style={{ fontSize: 13, fontWeight: 500, color: '#1890ff' }}>{recommended.functional}</div>
                              {recommended.functional !== currentFunctional && (
                                <div style={{ fontSize: 10, color: '#999', marginTop: 2 }}>原: {currentFunctional}</div>
                              )}
                            </div>
                          </Col>

                          {/* 基组 */}
                          <Col span={8}>
                            <div style={{ padding: '8px', background: '#f5f5f5', borderRadius: 4, border: '1px solid #e8e8e8' }}>
                              <div style={{ fontSize: 11, color: '#999', marginBottom: 4 }}>基组</div>
                              <div style={{ fontSize: 13, fontWeight: 500, color: '#52c41a' }}>{recommended.basis_set}</div>
                              {recommended.basis_set !== currentBasisSet && (
                                <div style={{ fontSize: 10, color: '#999', marginTop: 2 }}>原: {currentBasisSet}</div>
                              )}
                            </div>
                          </Col>

                          {/* 溶剂模型 */}
                          <Col span={8}>
                            <div style={{ padding: '8px', background: '#f5f5f5', borderRadius: 4, border: '1px solid #e8e8e8' }}>
                              <div style={{ fontSize: 11, color: '#999', marginBottom: 4 }}>溶剂模型</div>
                              <div style={{ fontSize: 13, fontWeight: 500, color: '#fa8c16' }}>
                                {recommended.solvent_model === 'gas' ? '气相' : recommended.solvent_model.toUpperCase()}
                              </div>
                              {recommended.solvent_model !== solventModel && (
                                <div style={{ fontSize: 10, color: '#999', marginTop: 2 }}>
                                  原: {solventModel === 'gas' ? '气相' : solventModel.toUpperCase()}
                                </div>
                              )}
                            </div>
                          </Col>
                        </Row>
                      </div>
                    )}

                    {/* 推荐理由 */}
                    {useRecommendedParams && recommended.reason && (
                      <Alert
                        message={recommended.reason}
                        type={hasChanges ? 'success' : 'info'}
                        icon={<BulbOutlined />}
                        showIcon
                        style={{ marginTop: 8 }}
                      />
                    )}
                  </div>
                );
              }}
            </Form.Item>
          </Card>

          <Divider orientation="left">
            <Space>
              <ThunderboltOutlined style={{ color: '#eb2f96' }} />
              电荷与自旋
            </Space>
          </Divider>

          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#ffadd2',
              background: 'linear-gradient(135deg, #fff0f6 0%, #fff 100%)',
              boxShadow: '0 2px 8px rgba(0,0,0,0.08)'
            }}
          >
            <Row gutter={16}>
              <Col span={8}>
                <Form.Item name="charge" label="电荷" style={{ marginBottom: 0 }}>
                  <InputNumber
                    style={{ width: '100%' }}
                    onChange={() => {
                      // 电荷改变时，如果开启自动计算，则重新计算自旋多重度
                      if (autoSpin) {
                        setTimeout(() => handleCalculateSpin(), 100);
                      }
                    }}
                  />
                </Form.Item>
              </Col>
              <Col span={8}>
                <Form.Item name="spin_multiplicity" label="自旋多重度" style={{ marginBottom: 0 }}>
                  <InputNumber min={1} style={{ width: '100%' }} />
                </Form.Item>
              </Col>
              <Col span={8}>
                <Form.Item label="自动计算" style={{ marginBottom: 0 }}>
                  <Space>
                    <Switch
                      checked={autoSpin}
                      onChange={(checked) => {
                        setAutoSpin(checked);
                        // 开启自动计算时立即计算
                        if (checked) {
                          setTimeout(() => handleCalculateSpin(), 100);
                        }
                      }}
                      checkedChildren="开"
                      unCheckedChildren="关"
                    />
                    <Button size="small" onClick={handleCalculateSpin} type="dashed">
                      计算
                    </Button>
                  </Space>
                </Form.Item>
              </Col>
            </Row>
          </Card>

          <Divider orientation="left">
            <Space>
              <SettingOutlined style={{ color: '#1677ff' }} />
              计算资源配置
              <Button type="link" size="small" onClick={handleGetSuggestion} icon={<ThunderboltOutlined />}>
                获取推荐配置
              </Button>
            </Space>
          </Divider>

          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#91caff',
              background: 'linear-gradient(135deg, #e6f4ff 0%, #fff 100%)'
            }}
          >
            <Row gutter={16}>
              <Col span={8}>
                <Form.Item
                  name="slurm_partition"
                  label="队列/分区"
                  tooltip="显示管理员分配给您的可用队列，队列状态实时从集群获取"
                  initialValue={getDefaultPartition()}
                  style={{ marginBottom: 0 }}
                >
                  <Select
                    placeholder={partitions.length > 0 ? "选择队列" : "暂无可用队列"}
                    disabled={partitions.length === 0}
                  >
                    {partitions.length > 0 ? (
                      partitions.map(p => (
                        <Select.Option
                          key={p.name}
                          value={p.name}
                          disabled={p.state !== 'up'}
                        >
                          <span style={{ color: p.state === 'up' ? 'inherit' : '#999' }}>
                            {p.name} {p.state === 'up'
                              ? `(可用 ${p.available_cpus}/${p.total_cpus} CPUs)`
                              : '(不可用)'}
                          </span>
                        </Select.Option>
                      ))
                    ) : (
                      <>
                        <Select.Option value="cpu">cpu</Select.Option>
                        <Select.Option value="gpu">gpu</Select.Option>
                        <Select.Option value="debug">debug</Select.Option>
                      </>
                    )}
                  </Select>
                </Form.Item>
              </Col>
              <Col span={8}>
                <Form.Item name="slurm_cpus" label="CPU核心数" initialValue={16} style={{ marginBottom: 0 }}>
                  <InputNumber min={1} max={64} style={{ width: '100%' }} />
                </Form.Item>
              </Col>
              <Col span={8}>
                <Form.Item name="slurm_time" label="最大运行时间(分钟)" initialValue={7200} style={{ marginBottom: 0 }}>
                  <InputNumber min={10} max={43200} style={{ width: '100%' }} />
                </Form.Item>
              </Col>
            </Row>
          </Card>

          {partitions.length === 0 && (
            <Alert
              message="暂无可用队列"
              description="请联系管理员分配队列权限，或等待集群信息加载"
              type="warning"
              showIcon
              style={{ marginBottom: 16 }}
            />
          )}

          <Alert
            message="QC计算资源说明"
            description="Gaussian计算为单节点多核任务。标准QC计算(6-31G(d)基组)通常需要16核、30分钟到数小时。较大基组可能需要更多资源。"
            type="info"
            showIcon
            style={{ marginBottom: 16 }}
          />

          <div style={{
            borderTop: '1px solid #f0f0f0',
            paddingTop: 16,
            marginTop: 8,
            textAlign: 'right'
          }}>
            <Space size="middle">
              <Button onClick={() => {
                setModalVisible(false);
                setEditingJob(null);
              }} size="large">
                取消
              </Button>
              <Button
                type="primary"
                htmlType="submit"
                loading={submitting}
                size="large"
                icon={<ExperimentOutlined />}
                style={{
                  background: 'linear-gradient(135deg, #722ed1 0%, #9254de 100%)',
                  border: 'none',
                  boxShadow: '0 4px 12px rgba(114, 46, 209, 0.3)'
                }}
              >
                {editingJob ? '保存修改' : '创建任务'}
              </Button>
            </Space>
          </div>
        </Form>
      </Modal>

      {/* 批量导入对话框 - 步骤式流程 */}
      <Modal
        title={
          <Space>
            <UploadOutlined style={{ color: '#1677ff' }} />
            <span style={{ fontWeight: 600 }}>批量导入QC计算任务</span>
          </Space>
        }
        open={batchImportVisible}
        onCancel={resetBatchImport}
        footer={importStep === 1 ? [
          <Button key="close" type="primary" onClick={resetBatchImport}>
            完成
          </Button>
        ] : [
          <Button key="cancel" onClick={resetBatchImport}>
            取消
          </Button>,
          <Button key="import" type="primary" onClick={handleBatchImport} loading={importing} disabled={!uploadFile}>
            开始导入
          </Button>,
        ]}
        width={700}
        centered
        destroyOnClose
      >
        <div style={{ padding: '24px 0' }}>
          <Steps current={importStep} style={{ marginBottom: 32 }}>
            <Step title="选择文件" />
            <Step title="导入结果" />
          </Steps>

          {/* 步骤1: 选择文件 */}
          {importStep === 0 && (
            <Space direction="vertical" size="large" style={{ width: '100%' }}>
              <Alert
                message="批量导入说明"
                description="下载模板后，在Excel中填写分子信息，包括名称、SMILES、电荷、泛函、基组等参数。支持 .xlsx、.xls、.csv 格式。"
                type="info"
                showIcon
              />
              <Card size="small" title="下载模板">
                <Button
                  icon={<DownloadOutlined />}
                  onClick={handleDownloadTemplate}
                  type="primary"
                  ghost
                >
                  下载QC计算导入模板 (含示例)
                </Button>
              </Card>
              <Card
                size="small"
                title="上传文件"
                style={{
                  border: uploadFile ? '1px solid #52c41a' : '1px dashed #d9d9d9',
                  backgroundColor: uploadFile ? '#f6ffed' : undefined
                }}
              >
                <div
                  style={{
                    padding: 24,
                    textAlign: 'center',
                    cursor: 'pointer',
                    borderRadius: 8,
                  }}
                  onClick={() => document.getElementById('qc-batch-file-input')?.click()}
                >
                  <input
                    id="qc-batch-file-input"
                    type="file"
                    accept=".xlsx,.xls,.csv"
                    style={{ display: 'none' }}
                    onChange={(e) => {
                      const file = e.target.files?.[0];
                      if (file) {
                        setUploadFile(file);
                      }
                    }}
                  />
                  {uploadFile ? (
                    <Space direction="vertical">
                      <FileExcelOutlined style={{ fontSize: 48, color: '#52c41a' }} />
                      <Text strong>{uploadFile.name}</Text>
                      <Text type="secondary">点击重新选择</Text>
                    </Space>
                  ) : (
                    <Space direction="vertical">
                      <InboxOutlined style={{ fontSize: 48, color: '#1677ff' }} />
                      <Text>点击选择文件</Text>
                      <Text type="secondary">支持 .xlsx, .xls, .csv</Text>
                    </Space>
                  )}
                </div>
              </Card>
            </Space>
          )}

          {/* 步骤2: 导入结果 */}
          {importStep === 1 && importResult && (
            <Space direction="vertical" size="large" style={{ width: '100%' }}>
              <Result
                status={importResult.failed_qc_jobs === 0 ? 'success' : 'warning'}
                title={importResult.failed_qc_jobs === 0 ? '导入成功' : '导入完成（部分失败）'}
                subTitle={`共 ${importResult.total_qc_jobs} 个任务`}
              />

              <Row gutter={16}>
                <Col span={12}>
                  <Card size="small">
                    <div style={{ textAlign: 'center' }}>
                      <div style={{ fontSize: 24, fontWeight: 600, color: '#52c41a' }}>
                        {importResult.success_qc_jobs}
                      </div>
                      <div style={{ color: '#666' }}>成功创建</div>
                    </div>
                  </Card>
                </Col>
                <Col span={12}>
                  <Card size="small">
                    <div style={{ textAlign: 'center' }}>
                      <div style={{ fontSize: 24, fontWeight: 600, color: '#ff4d4f' }}>
                        {importResult.failed_qc_jobs}
                      </div>
                      <div style={{ color: '#666' }}>创建失败</div>
                    </div>
                  </Card>
                </Col>
              </Row>

              {/* 错误详情 */}
              {importResult.errors.length > 0 && (
                <Card size="small" title={<span style={{ color: '#ff4d4f' }}>错误详情</span>}>
                  <Table
                    size="small"
                    dataSource={importResult.errors.map((e, i) => ({ ...e, key: i }))}
                    columns={[
                      { title: '行号', dataIndex: 'row', width: 60 },
                      { title: '类型', dataIndex: 'type', width: 120 },
                      {
                        title: '错误信息',
                        dataIndex: 'message',
                        render: (text: string) => (
                          <Text type="danger" style={{ fontSize: 12 }}>{text}</Text>
                        )
                      },
                    ]}
                    pagination={false}
                    scroll={{ y: 150 }}
                  />
                </Card>
              )}

              {/* 成功详情 */}
              {importResult.qc_job_results.length > 0 && (
                <Card size="small" title={<span style={{ color: '#52c41a' }}>成功导入的任务</span>}>
                  <div style={{ maxHeight: 150, overflow: 'auto' }}>
                    {importResult.qc_job_results.map((r, i) => (
                      <Tag key={i} color="green" style={{ margin: 4 }}>
                        {r.molecule_name} (ID: {r.job_id})
                      </Tag>
                    ))}
                  </div>
                </Card>
              )}
            </Space>
          )}
        </div>
      </Modal>
    </div>
  );
}

