/**
 * 计算任务管理页面
 */
import { useState, useEffect, useRef, useCallback } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import {
  Button,
  Space,
  message,
  Modal,
  Form,
  Row,
  Col,
  Spin,
  Empty,
  Select,
  InputNumber,
  Input,
  Tabs,
  Divider,
  Alert,
  Tooltip,
  Typography,
  Card,
  Checkbox,
  Tag,
  DatePicker,
  Statistic,
  Table,
  Switch,
} from 'antd';
import {
  PlusOutlined,
  ReloadOutlined,
  ThunderboltOutlined,
  RocketOutlined,
  ExperimentOutlined,
  FolderAddOutlined,
  SearchOutlined,
  FilterOutlined,
  AppstoreOutlined,
  UnorderedListOutlined,
  ClockCircleOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
} from '@ant-design/icons';
import JobCard from '../components/JobCard';
import { getMDJobs, createMDJob, cancelMDJob, deleteMDJob, resubmitMDJob, updateMDJobConfig } from '../api/jobs';
import { getElectrolytes, createElectrolyteNew } from '../api/electrolytes';
import { getPartitions, getSlurmSuggestion, type PartitionInfo } from '../api/slurm';
import { getProjects } from '../api/projects';
import type { MDJob, MDJobCreate, ElectrolyteSystem, Project } from '../types';
import { JobStatus } from '../types';
import ElectrolyteFormOptimized from '../components/ElectrolyteFormOptimized';
import AccuracyLevelSelector from '../components/AccuracyLevelSelector';
import dayjs, { Dayjs } from 'dayjs';

const { Title, Text } = Typography;
const { RangePicker } = DatePicker;

export default function Jobs() {
  const location = useLocation();
  const navigate = useNavigate();
  const [jobs, setJobs] = useState<MDJob[]>([]);
  const [electrolytes, setElectrolytes] = useState<ElectrolyteSystem[]>([]);
  const [projects, setProjects] = useState<Project[]>([]);
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);
  const [loading, setLoading] = useState(false);
  const [modalVisible, setModalVisible] = useState(false);
  const [resubmitModalVisible, setResubmitModalVisible] = useState(false);
  const [activeTab, setActiveTab] = useState('all');
  const [form] = Form.useForm();
  const [resubmitForm] = Form.useForm();
  const [resubmittingJob, setResubmittingJob] = useState<MDJob | null>(null);
  const [lastRefresh, setLastRefresh] = useState<Date>(new Date());
  const pollingRef = useRef<ReturnType<typeof setInterval> | null>(null);

  // 新建配方相关状态
  const [electrolyteModalVisible, setElectrolyteModalVisible] = useState(false);
  const [electrolyteForm] = Form.useForm();
  const [selectedCations, setSelectedCations] = useState<any[]>([]);
  const [selectedAnions, setSelectedAnions] = useState<any[]>([]);

  // 精度等级相关状态
  const [selectedAccuracyLevel, setSelectedAccuracyLevel] = useState<string>('standard');
  const [accuracyDefaults, setAccuracyDefaults] = useState<any>(null);

  // 筛选和视图状态
  const [searchText, setSearchText] = useState('');
  const [projectFilter, setProjectFilter] = useState<number | undefined>(undefined);
  const [electrolyteFilter, setElectrolyteFilter] = useState<number | undefined>(undefined);
  const [partitionFilter, setPartitionFilter] = useState<string | undefined>(undefined);
  const [dateRange, setDateRange] = useState<[Dayjs | null, Dayjs | null] | null>(null);
  const [viewMode, setViewMode] = useState<'card' | 'table'>(() => {
    // 从localStorage读取视图模式
    const saved = localStorage.getItem('md-jobs-view-mode');
    return (saved === 'card' || saved === 'table') ? saved : 'card';
  });
  const [sortBy, setSortBy] = useState<'created_at' | 'updated_at' | 'id'>('created_at');
  const [sortOrder, setSortOrder] = useState<'asc' | 'desc'>('desc');

  // 保存视图模式到localStorage
  const handleViewModeChange = (mode: 'card' | 'table') => {
    setViewMode(mode);
    localStorage.setItem('md-jobs-view-mode', mode);
  };

  // 加载任务列表
  const loadJobs = useCallback(async () => {
    try {
      const data = await getMDJobs();
      setJobs(data);
      setLastRefresh(new Date());
    } catch (error: any) {
      console.error('加载任务列表失败:', error);
    }
  }, []);

  // 加载电解质配方
  const loadElectrolytes = async () => {
    try {
      const data = await getElectrolytes();
      setElectrolytes(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载电解质配方列表失败');
    }
  };

  // 加载项目列表
  const loadProjects = async () => {
    try {
      const data = await getProjects();
      setProjects(data);
    } catch (error: any) {
      console.error('加载项目列表失败:', error);
    }
  };

  // 加载 Slurm 分区信息
  const loadPartitions = async () => {
    try {
      const data = await getPartitions();
      setPartitions(data);
    } catch (error: any) {
      console.error('加载分区信息失败:', error);
      // 使用默认分区
      setPartitions([
        { name: 'cpu', state: 'up', total_nodes: 0, available_nodes: 0, total_cpus: 0, available_cpus: 0 },
        { name: 'gpu', state: 'up', total_nodes: 0, available_nodes: 0, total_cpus: 0, available_cpus: 0 },
      ]);
    }
  };

  const loadData = async () => {
    setLoading(true);
    try {
      await Promise.all([loadJobs(), loadElectrolytes(), loadProjects(), loadPartitions()]);
    } finally {
      setLoading(false);
    }
  };

  // 检查是否有活跃任务（需要轮询）
  const hasActiveJobs = useCallback(() => {
    return jobs.some(job =>
      job.status === JobStatus.QUEUED ||
      job.status === JobStatus.RUNNING ||
      job.status === JobStatus.POSTPROCESSING
    );
  }, [jobs]);

  useEffect(() => {
    loadData();
  }, []);

  // 加载精度等级配置
  useEffect(() => {
    const loadAccuracyLevels = async () => {
      try {
        const token = localStorage.getItem('access_token');
        const response = await fetch('/api/v1/jobs/accuracy-levels', {
          headers: { Authorization: `Bearer ${token}` }
        });
        const data = await response.json();
        setAccuracyDefaults(data);
      } catch (error) {
        console.error('加载精度等级配置失败:', error);
      }
    };
    loadAccuracyLevels();
  }, []);

  // 智能轮询：只有在有活跃任务时才轮询
  useEffect(() => {
    // 清除之前的轮询
    if (pollingRef.current) {
      clearInterval(pollingRef.current);
      pollingRef.current = null;
    }

    // 如果有活跃任务，启动轮询（每 10 秒刷新一次）
    if (hasActiveJobs()) {
      pollingRef.current = setInterval(() => {
        loadJobs();
      }, 10000);
    }

    // 清理轮询
    return () => {
      if (pollingRef.current) {
        clearInterval(pollingRef.current);
      }
    };
  }, [hasActiveJobs, loadJobs]);

  // 获取默认分区
  const getDefaultPartition = () => {
    if (partitions.length > 0) {
      const upPartition = partitions.find(p => p.state === 'up');
      return upPartition?.name || partitions[0].name;
    }
    return 'cpu';
  };

  // 检查是否需要自动打开创建对话框
  useEffect(() => {
    if (location.state?.openCreateModal) {
      // 直接设置 modal 可见并初始化表单
      form.resetFields();
      form.setFieldsValue({
        nsteps_npt: 100000,
        nsteps_nvt: 500000,
        timestep: 1.0,
        slurm_partition: getDefaultPartition(),
        slurm_nodes: 1,
        slurm_ntasks: 8,
        slurm_cpus_per_task: 8,
        slurm_time: 7200,
      });
      setModalVisible(true);
      // 清除 state，避免刷新时重复打开
      window.history.replaceState({}, document.title);
    }
  }, [location, partitions]);

  // 打开创建对话框
  const handleOpenModal = () => {
    form.resetFields();
    setSelectedAccuracyLevel('standard');
    form.setFieldsValue({
      job_name: '',
      accuracy_level: 'standard',
      // 模拟参数留空，使用精度等级默认值
      nsteps_npt: undefined,
      nsteps_nvt: undefined,
      timestep: undefined,
      temperature: undefined,
      pressure: undefined,
      freq_trj_npt: undefined,
      freq_trj_nvt: undefined,
      thermo_freq: undefined,
      slurm_partition: getDefaultPartition(),
      slurm_nodes: 1,
      slurm_ntasks: 8,
      slurm_cpus_per_task: 8,
      slurm_time: 7200,
      // QC计算选项 - 多选模式
      qc_enabled: false,
      qc_functionals: ['B3LYP'],
      qc_basis_sets: ['6-31++g(d,p)'],
      qc_solvent_models: ['pcm'],
      qc_solvents: ['Water'],
    });
    setModalVisible(true);
  };

  // 获取推荐配置
  const handleGetSuggestion = async (formInstance: typeof form) => {
    try {
      const suggestion = await getSlurmSuggestion({ job_type: 'md' });
      formInstance.setFieldsValue({
        slurm_partition: suggestion.partition,
        slurm_ntasks: suggestion.ntasks,
        slurm_cpus_per_task: suggestion.cpus_per_task,
      });
      message.success(`已应用推荐配置: ${suggestion.reason}`);
    } catch (error: any) {
      message.error('获取推荐配置失败');
    }
  };

  // 关闭对话框
  const handleCloseModal = () => {
    setModalVisible(false);
    form.resetFields();
  };

  // 打开新建配方对话框
  const handleOpenElectrolyteModal = () => {
    setElectrolyteModalVisible(true);
    electrolyteForm.resetFields();
    setSelectedCations([]);
    setSelectedAnions([]);
  };

  // 关闭新建配方对话框
  const handleCloseElectrolyteModal = () => {
    setElectrolyteModalVisible(false);
    electrolyteForm.resetFields();
    setSelectedCations([]);
    setSelectedAnions([]);
  };

  // 创建配方
  const handleCreateElectrolyte = async () => {
    try {
      const values = await electrolyteForm.validateFields();

      // 获取盒子尺寸
      const boxSize = values.box_size || 40;
      const box = {
        type: 'cubic' as const,
        dimensions: [boxSize],
      };

      // 构建请求数据（新格式 - 使用浓度）
      const electrolyteData = {
        project_id: values.project_id,
        name: values.name,
        description: values.description,
        temperature: values.temperature || 298.15,
        pressure: 1.0,
        nsteps_npt: 5000000,
        nsteps_nvt: 10000000,
        timestep: 1.0,
        force_field: 'OPLS',
        solvents: values.solvents || [],
        box: box,
        // 使用 charge 和 concentration，而不是 smiles 和 count
        cations: selectedCations.map(cat => ({
          name: cat.name,
          charge: cat.charge,
          concentration: cat.concentration,
        })),
        anions: selectedAnions.map(an => ({
          name: an.name,
          charge: an.charge,
          concentration: an.concentration,
        })),
      };

      console.log('=== Jobs.tsx 创建电解质请求数据 ===');
      console.log('electrolyteData:', JSON.stringify(electrolyteData, null, 2));

      const newElectrolyte = await createElectrolyteNew(electrolyteData);
      message.success('配方创建成功');

      // 重新加载配方列表
      await loadElectrolytes();

      // 自动选择新创建的配方
      form.setFieldsValue({ electrolyte_id: newElectrolyte.id });

      handleCloseElectrolyteModal();
    } catch (error: any) {
      console.error('=== Jobs.tsx 创建电解质失败 ===');
      console.error('error:', error);
      console.error('error.response:', error.response);
      console.error('error.response.data:', error.response?.data);
      if (error.response) {
        const detail = error.response?.data?.detail;
        if (Array.isArray(detail)) {
          // Pydantic validation errors
          const errorMessages = detail.map((err: any) =>
            `${err.loc.join('.')}: ${err.msg}`
          ).join('; ');
          message.error(`验证失败: ${errorMessages}`);
        } else {
          message.error(detail || '创建配方失败');
        }
      }
    }
  };

  // 提交表单
  const handleSubmit = async () => {
    try {
      const values = await form.validateFields();
      const data: MDJobCreate = {
        system_id: values.electrolyte_id,
        job_name: values.job_name || undefined,
        accuracy_level: values.accuracy_level || 'standard',
        nsteps_npt: values.nsteps_npt || undefined,
        nsteps_nvt: values.nsteps_nvt || undefined,
        timestep: values.timestep,
        temperature: values.temperature,
        pressure: values.pressure,
        freq_trj_npt: values.freq_trj_npt || undefined,
        freq_trj_nvt: values.freq_trj_nvt || undefined,
        thermo_freq: values.thermo_freq || undefined,
        // Slurm 资源配置
        slurm_partition: values.slurm_partition || 'cpu',
        slurm_nodes: values.slurm_nodes || 1,
        slurm_ntasks: values.slurm_ntasks || 8,
        slurm_cpus_per_task: values.slurm_cpus_per_task || 8,
        slurm_time: values.slurm_time || 7200,
      };
      // QC计算选项 - 支持多选
      if (values.qc_enabled) {
        data.qc_options = {
          enabled: true,
          // 使用复数形式的数组字段（后端已支持）
          functionals: values.qc_functionals || ['B3LYP'],
          basis_sets: values.qc_basis_sets || ['6-31++g(d,p)'],
          solvent_models: values.qc_solvent_models || ['pcm'],
          solvents: values.qc_solvents || ['Water'],
          molecules: [], // 将由后端从电解质配方中提取
          // 兼容旧版字段（取第一个值）
          functional: values.qc_functionals?.[0] || 'B3LYP',
          basis_set: values.qc_basis_sets?.[0] || '6-31++g(d,p)',
          solvent_model: values.qc_solvent_models?.[0] || 'pcm',
          solvent_name: values.qc_solvents?.[0] || 'Water',
        } as any;
      }
      await createMDJob(data);
      message.success('任务创建成功');
      handleCloseModal();
      loadJobs();
    } catch (error: any) {
      if (error.response) {
        message.error(error.response?.data?.detail || '创建失败');
      }
    }
  };

  // 取消任务
  const handleCancel = async (id: number) => {
    try {
      await cancelMDJob(id);
      message.success('任务已取消');
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '取消失败');
    }
  };

  // 删除任务
  const handleDelete = async (id: number) => {
    try {
      await deleteMDJob(id);
      message.success('任务已删除');
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  // 打开重新提交对话框
  const handleOpenResubmitModal = (job: MDJob) => {
    setResubmittingJob(job);

    // 从任务配置中读取参数并填充表单
    const config = job.config || {};
    resubmitForm.setFieldsValue({
      nsteps_npt: config.nsteps_npt || 100000,
      nsteps_nvt: config.nsteps_nvt || 500000,
      timestep: config.timestep || 1.0,
      slurm_partition: config.slurm_partition || 'cpu',
      slurm_nodes: config.slurm_nodes || 1,
      slurm_ntasks: config.slurm_ntasks || 8,
      slurm_cpus_per_task: config.slurm_cpus_per_task || 8,
      slurm_time: config.slurm_time || 7200,
    });

    setResubmitModalVisible(true);
  };

  // 关闭重新提交对话框
  const handleCloseResubmitModal = () => {
    setResubmitModalVisible(false);
    setResubmittingJob(null);
    resubmitForm.resetFields();
  };

  // 提交重新提交表单
  const handleResubmitSubmit = async () => {
    if (!resubmittingJob) return;

    try {
      const values = await resubmitForm.validateFields();

      // 更新任务配置
      const updatedConfig = {
        ...resubmittingJob.config,
        nsteps_npt: values.nsteps_npt,
        nsteps_nvt: values.nsteps_nvt,
        timestep: values.timestep,
        slurm_partition: values.slurm_partition,
        slurm_nodes: values.slurm_nodes,
        slurm_ntasks: values.slurm_ntasks,
        slurm_cpus_per_task: values.slurm_cpus_per_task,
        slurm_time: values.slurm_time,
      };

      // 调用 API 更新配置并重新提交
      await updateMDJobConfig(resubmittingJob.id, updatedConfig);
      await resubmitMDJob(resubmittingJob.id);

      message.success('任务配置已更新并重新提交到集群');
      handleCloseResubmitModal();
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '重新提交失败');
    }
  };

  // 过滤任务
  const getFilteredJobs = () => {
    let filtered = [...jobs];

    // 按状态标签页筛选
    if (activeTab === 'created') {
      filtered = filtered.filter((j) => j.status === JobStatus.CREATED);
    } else if (activeTab === 'running') {
      filtered = filtered.filter((j) =>
        j.status === JobStatus.QUEUED ||
        j.status === JobStatus.RUNNING ||
        j.status === JobStatus.POSTPROCESSING
      );
    } else if (activeTab === 'completed') {
      filtered = filtered.filter((j) => j.status === JobStatus.COMPLETED);
    } else if (activeTab === 'failed') {
      filtered = filtered.filter((j) => j.status === JobStatus.FAILED || j.status === JobStatus.CANCELLED);
    }

    // 搜索筛选（任务名称、Slurm Job ID）
    if (searchText) {
      const search = searchText.toLowerCase();
      filtered = filtered.filter((job) =>
        (job.config?.job_name && job.config.job_name.toLowerCase().includes(search)) ||
        (job.slurm_job_id && job.slurm_job_id.toString().includes(search)) ||
        (job.id && job.id.toString().includes(search))
      );
    }

    // 项目筛选
    if (projectFilter !== undefined) {
      filtered = filtered.filter((job) => {
        const electrolyte = electrolytes.find(e => e.id === job.system_id);
        return electrolyte?.project_id === projectFilter;
      });
    }

    // 配方筛选
    if (electrolyteFilter !== undefined) {
      filtered = filtered.filter((job) => job.system_id === electrolyteFilter);
    }

    // 分区筛选
    if (partitionFilter) {
      filtered = filtered.filter((job) => job.config?.slurm_partition === partitionFilter);
    }

    // 时间范围筛选
    if (dateRange && dateRange[0] && dateRange[1]) {
      const startDate = dateRange[0].startOf('day');
      const endDate = dateRange[1].endOf('day');
      filtered = filtered.filter((job) => {
        const jobDate = dayjs(job.created_at);
        return jobDate.isAfter(startDate) && jobDate.isBefore(endDate);
      });
    }

    // 排序
    filtered.sort((a, b) => {
      let aValue: any, bValue: any;

      if (sortBy === 'created_at') {
        aValue = new Date(a.created_at).getTime();
        bValue = new Date(b.created_at).getTime();
      } else if (sortBy === 'updated_at') {
        aValue = new Date(a.updated_at).getTime();
        bValue = new Date(b.updated_at).getTime();
      } else {
        aValue = a.id;
        bValue = b.id;
      }

      return sortOrder === 'asc' ? aValue - bValue : bValue - aValue;
    });

    return filtered;
  };

  const filteredJobs = getFilteredJobs();

  // 重置筛选
  const handleResetFilters = () => {
    setSearchText('');
    setProjectFilter(undefined);
    setElectrolyteFilter(undefined);
    setPartitionFilter(undefined);
    setDateRange(null);
    setSortBy('created_at');
    setSortOrder('desc');
  };

  // 计算各状态任务数量
  const createdCount = jobs.filter((j) => j.status === JobStatus.CREATED).length;
  const runningCount = jobs.filter((j) =>
    j.status === JobStatus.QUEUED ||
    j.status === JobStatus.RUNNING ||
    j.status === JobStatus.POSTPROCESSING
  ).length;
  const completedCount = jobs.filter((j) => j.status === JobStatus.COMPLETED).length;
  const failedCount = jobs.filter((j) => j.status === JobStatus.FAILED || j.status === JobStatus.CANCELLED).length;

  return (
    <div style={{
      padding: '24px',
      background: '#f5f7fb',
      minHeight: 'calc(100vh - 64px)'
    }}>
      {/* 页面标题区域 */}
      <div style={{ marginBottom: 24 }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
          <div>
            <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
              <RocketOutlined style={{ marginRight: 12, color: '#1677ff' }} />
              计算任务管理
            </Title>
            <Space>
              <Text type="secondary">管理分子动力学模拟任务，监控计算进度</Text>
              <Text type="secondary">|</Text>
              <Text type="secondary" style={{ fontSize: 12 }}>
                最后更新: {lastRefresh.toLocaleTimeString()}
                {hasActiveJobs() && <Text type="success" style={{ marginLeft: 8 }}>(自动刷新中)</Text>}
              </Text>
            </Space>
          </div>
          <Space>
            <Tooltip title="刷新任务列表">
              <Button
                icon={<ReloadOutlined />}
                onClick={loadJobs}
                style={{ borderRadius: 8 }}
              >
                刷新
              </Button>
            </Tooltip>
            <Button
              type="primary"
              icon={<PlusOutlined />}
              onClick={handleOpenModal}
              size="large"
              style={{
                borderRadius: 8,
                boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)',
              }}
            >
              创建新任务
            </Button>
          </Space>
        </div>
      </div>

      {/* 统计卡片 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
        }}
      >
        <Row gutter={24} align="middle" justify="space-around">
          <Col>
            <div style={{ textAlign: 'center' }}>
              <div style={{
                fontSize: 28,
                fontWeight: 700,
                color: '#1677ff',
                lineHeight: 1.2
              }}>
                {jobs.length}
              </div>
              <Text type="secondary" style={{ fontSize: 12 }}>全部任务</Text>
            </div>
          </Col>
          <Col>
            <div style={{ textAlign: 'center' }}>
              <div style={{
                fontSize: 28,
                fontWeight: 700,
                color: '#faad14',
                lineHeight: 1.2
              }}>
                {createdCount}
              </div>
              <Text type="secondary" style={{ fontSize: 12 }}>待配置</Text>
            </div>
          </Col>
          <Col>
            <div style={{ textAlign: 'center' }}>
              <div style={{
                fontSize: 28,
                fontWeight: 700,
                color: '#52c41a',
                lineHeight: 1.2
              }}>
                {runningCount}
              </div>
              <Text type="secondary" style={{ fontSize: 12 }}>运行中</Text>
            </div>
          </Col>
          <Col>
            <div style={{ textAlign: 'center' }}>
              <div style={{
                fontSize: 28,
                fontWeight: 700,
                color: '#722ed1',
                lineHeight: 1.2
              }}>
                {completedCount}
              </div>
              <Text type="secondary" style={{ fontSize: 12 }}>已完成</Text>
            </div>
          </Col>
          <Col>
            <div style={{ textAlign: 'center' }}>
              <div style={{
                fontSize: 28,
                fontWeight: 700,
                color: '#ff4d4f',
                lineHeight: 1.2
              }}>
                {failedCount}
              </div>
              <Text type="secondary" style={{ fontSize: 12 }}>失败/取消</Text>
            </div>
          </Col>
        </Row>
      </Card>

      {/* 任务分类标签 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
        }}
        styles={{ body: { padding: '12px 24px' } }}
      >
        <Tabs
          activeKey={activeTab}
          onChange={setActiveTab}
          items={[
            { key: 'all', label: `全部 (${jobs.length})` },
            { key: 'created', label: `待配置 (${createdCount})` },
            { key: 'running', label: `运行中 (${runningCount})` },
            { key: 'completed', label: `已完成 (${completedCount})` },
            { key: 'failed', label: `失败/取消 (${failedCount})` },
          ]}
        />
      </Card>

      {/* 筛选栏 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
          background: '#fafafa',
        }}
      >
        <Row gutter={[16, 16]}>
          <Col xs={24} sm={12} md={8} lg={6}>
            <Input
              placeholder="搜索任务名称、Job ID"
              prefix={<SearchOutlined />}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              allowClear
            />
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Select
              placeholder="项目"
              value={projectFilter}
              onChange={setProjectFilter}
              allowClear
              style={{ width: '100%' }}
            >
              {projects.map((p) => (
                <Select.Option key={p.id} value={p.id}>
                  {p.name}
                </Select.Option>
              ))}
            </Select>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Select
              placeholder="配方"
              value={electrolyteFilter}
              onChange={setElectrolyteFilter}
              allowClear
              style={{ width: '100%' }}
            >
              {electrolytes.map((e) => (
                <Select.Option key={e.id} value={e.id}>
                  {e.name}
                </Select.Option>
              ))}
            </Select>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Select
              placeholder="分区"
              value={partitionFilter}
              onChange={setPartitionFilter}
              allowClear
              style={{ width: '100%' }}
            >
              {partitions.map((p) => (
                <Select.Option key={p.name} value={p.name}>
                  {p.name}
                </Select.Option>
              ))}
            </Select>
          </Col>
          <Col xs={24} sm={12} md={8} lg={6}>
            <RangePicker
              value={dateRange}
              onChange={setDateRange}
              style={{ width: '100%' }}
              placeholder={['开始日期', '结束日期']}
            />
          </Col>
        </Row>
        <Row gutter={[16, 16]} style={{ marginTop: 16 }}>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Select
              placeholder="排序方式"
              value={sortBy}
              onChange={setSortBy}
              style={{ width: '100%' }}
            >
              <Select.Option value="created_at">创建时间</Select.Option>
              <Select.Option value="updated_at">更新时间</Select.Option>
              <Select.Option value="id">任务ID</Select.Option>
            </Select>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Select
              placeholder="排序顺序"
              value={sortOrder}
              onChange={setSortOrder}
              style={{ width: '100%' }}
            >
              <Select.Option value="desc">降序</Select.Option>
              <Select.Option value="asc">升序</Select.Option>
            </Select>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Space>
              <Button onClick={handleResetFilters} icon={<ReloadOutlined />}>
                重置筛选
              </Button>
              <Button
                type={viewMode === 'card' ? 'primary' : 'default'}
                icon={<AppstoreOutlined />}
                onClick={() => handleViewModeChange('card')}
              />
              <Button
                type={viewMode === 'table' ? 'primary' : 'default'}
                icon={<UnorderedListOutlined />}
                onClick={() => handleViewModeChange('table')}
              />
            </Space>
          </Col>
          <Col xs={24} sm={12} md={8} lg={12} style={{ textAlign: 'right' }}>
            <Text type="secondary">
              显示 {filteredJobs.length} / {jobs.length} 个任务
            </Text>
          </Col>
        </Row>
      </Card>

      {/* 任务列表 */}
      <Spin spinning={loading}>
        {filteredJobs.length === 0 ? (
          <Card
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
            }}
          >
            <Empty
              image={<RocketOutlined style={{ fontSize: 64, color: '#d9d9d9' }} />}
              description={
                <Space direction="vertical" size={8}>
                  <Text type="secondary" style={{ fontSize: 16 }}>
                    {activeTab === 'all' ? '还没有任务' : '没有符合条件的任务'}
                  </Text>
                  {activeTab === 'all' && (
                    <Text type="secondary">点击上方按钮创建第一个任务</Text>
                  )}
                </Space>
              }
              style={{ padding: '60px 0' }}
            >
              {activeTab === 'all' && (
                <Button
                  type="primary"
                  icon={<PlusOutlined />}
                  onClick={handleOpenModal}
                >
                  创建新任务
                </Button>
              )}
            </Empty>
          </Card>
        ) : viewMode === 'card' ? (
          <Row gutter={[16, 16]}>
            {filteredJobs.map((job) => {
              const electrolyte = electrolytes.find(e => e.id === job.system_id);
              return (
                <Col xs={24} sm={24} md={12} lg={8} key={job.id}>
                  <JobCard
                    job={job}
                    electrolyte={electrolyte}
                    onCancel={handleCancel}
                    onResubmit={handleOpenResubmitModal}
                    onDelete={handleDelete}
                  />
                </Col>
              );
            })}
          </Row>
        ) : (
          <Card
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
            }}
          >
            <Table
              dataSource={filteredJobs}
              rowKey="id"
              pagination={{
                pageSize: 20,
                showSizeChanger: true,
                showTotal: (total) => `共 ${total} 个任务`,
                pageSizeOptions: ['10', '20', '50', '100'],
              }}
              scroll={{ x: 1200 }}
              columns={[
                {
                  title: 'ID',
                  dataIndex: 'id',
                  key: 'id',
                  width: 80,
                  fixed: 'left',
                },
                {
                  title: '任务名称',
                  key: 'job_name',
                  width: 200,
                  fixed: 'left',
                  render: (_: any, record: MDJob) => (
                    <a onClick={() => navigate(`/workspace/jobs/${record.id}/detail`)}>
                      {record.config?.job_name || `任务 #${record.id}`}
                    </a>
                  ),
                },
                {
                  title: '状态',
                  dataIndex: 'status',
                  key: 'status',
                  width: 120,
                  render: (status: JobStatus) => {
                    const statusConfig: Record<JobStatus, { color: string; text: string; icon: any }> = {
                      [JobStatus.CREATED]: { color: 'default', text: '待配置', icon: <ClockCircleOutlined /> },
                      [JobStatus.QUEUED]: { color: 'processing', text: '排队中', icon: <ClockCircleOutlined /> },
                      [JobStatus.RUNNING]: { color: 'processing', text: '运行中', icon: <RocketOutlined /> },
                      [JobStatus.POSTPROCESSING]: { color: 'processing', text: '后处理', icon: <RocketOutlined /> },
                      [JobStatus.COMPLETED]: { color: 'success', text: '已完成', icon: <CheckCircleOutlined /> },
                      [JobStatus.FAILED]: { color: 'error', text: '失败', icon: <CloseCircleOutlined /> },
                      [JobStatus.CANCELLED]: { color: 'default', text: '已取消', icon: <CloseCircleOutlined /> },
                    };
                    const config = statusConfig[status] || statusConfig[JobStatus.CREATED];
                    return (
                      <Tag color={config.color} icon={config.icon}>
                        {config.text}
                      </Tag>
                    );
                  },
                },
                {
                  title: '配方',
                  dataIndex: 'system_id',
                  key: 'system_id',
                  width: 200,
                  render: (systemId: number) => {
                    const electrolyte = electrolytes.find(e => e.id === systemId);
                    return electrolyte?.name || '-';
                  },
                },
                {
                  title: 'Slurm Job ID',
                  dataIndex: 'slurm_job_id',
                  key: 'slurm_job_id',
                  width: 120,
                  render: (id: number | null) => id || '-',
                },
                {
                  title: '分区',
                  key: 'slurm_partition',
                  width: 100,
                  render: (_: any, record: MDJob) => {
                    const partition = record.config?.slurm_partition;
                    return partition ? <Tag color="blue">{partition}</Tag> : '-';
                  },
                },
                {
                  title: '创建时间',
                  dataIndex: 'created_at',
                  key: 'created_at',
                  width: 180,
                  render: (time: string) => new Date(time).toLocaleString('zh-CN'),
                },
                {
                  title: '操作',
                  key: 'actions',
                  width: 200,
                  fixed: 'right',
                  render: (_: any, record: MDJob) => (
                    <Space size="small">
                      <Button
                        type="link"
                        size="small"
                        onClick={() => navigate(`/workspace/jobs/${record.id}/detail`)}
                      >
                        详情
                      </Button>
                      {(record.status === JobStatus.QUEUED || record.status === JobStatus.RUNNING) && (
                        <Button
                          type="link"
                          size="small"
                          danger
                          onClick={() => handleCancel(record.id)}
                        >
                          取消
                        </Button>
                      )}
                      {(record.status === JobStatus.FAILED || record.status === JobStatus.CANCELLED) && (
                        <Button
                          type="link"
                          size="small"
                          onClick={() => handleOpenResubmitModal(record)}
                        >
                          重新提交
                        </Button>
                      )}
                    </Space>
                  ),
                },
              ]}
            />
          </Card>
        )}
      </Spin>

      {/* 创建任务对话框 */}
      <Modal
        title={
          <Space>
            <RocketOutlined style={{ color: '#1677ff' }} />
            创建新计算任务
          </Space>
        }
        open={modalVisible}
        onOk={handleSubmit}
        onCancel={handleCloseModal}
        okText="创建"
        cancelText="取消"
        width={800}
        centered
        styles={{
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' },
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
        }}
      >
        <Form form={form} layout="vertical" style={{ marginTop: 24 }}>
          <Form.Item
            name="electrolyte_id"
            label="选择电解质配方"
            rules={[{ required: true, message: '请选择电解质配方' }]}
          >
            <Select
              placeholder="选择要计算的电解质配方"
              notFoundContent={
                electrolytes.length === 0 ? (
                  <div style={{ textAlign: 'center', padding: '20px 0' }}>
                    <Empty
                      image={Empty.PRESENTED_IMAGE_SIMPLE}
                      description="暂无配方"
                      style={{ marginBottom: 12 }}
                    />
                    <Button
                      type="primary"
                      icon={<PlusOutlined />}
                      onClick={handleOpenElectrolyteModal}
                      size="small"
                    >
                      新建配方
                    </Button>
                  </div>
                ) : undefined
              }
              dropdownRender={(menu) => {
                const hasElectrolytes = electrolytes && electrolytes.length > 0;
                return (
                  <>
                    {menu}
                    {hasElectrolytes && (
                      <>
                        <Divider style={{ margin: '8px 0' }} />
                        <div style={{ padding: '4px 8px' }}>
                          <Button
                            type="link"
                            icon={<PlusOutlined />}
                            onClick={handleOpenElectrolyteModal}
                            style={{ width: '100%', textAlign: 'left' }}
                          >
                            新建配方
                          </Button>
                        </div>
                      </>
                    )}
                  </>
                );
              }}
            >
              {electrolytes.map((e) => (
                <Select.Option key={e.id} value={e.id}>
                  {e.name} ({e.temperature} K)
                </Select.Option>
              ))}
            </Select>
          </Form.Item>

          <Form.Item
            label="自定义名称（可选）"
            name="job_name"
            tooltip="可选的自定义名称后缀，将添加到自动生成的任务名称后面"
          >
            <Input placeholder="留空或输入自定义名称后缀（如：高温测试）" allowClear />
          </Form.Item>

          <Divider orientation="left">精度等级</Divider>

          <Form.Item name="accuracy_level">
            <AccuracyLevelSelector
              value={selectedAccuracyLevel}
              onChange={(value) => {
                setSelectedAccuracyLevel(value);
                form.setFieldsValue({ accuracy_level: value });

                // 如果切换到自定义模式，自动填充参考值（标准模式的参数）
                if (value === 'custom' && accuracyDefaults?.custom) {
                  const customDefaults = accuracyDefaults.custom;
                  form.setFieldsValue({
                    nsteps_npt: customDefaults.nsteps_npt,
                    nsteps_nvt: customDefaults.nsteps_nvt,
                    timestep: customDefaults.timestep,
                    temperature: customDefaults.temperature,
                    pressure: customDefaults.pressure,
                    freq_trj_npt: customDefaults.freq_trj_npt,
                    freq_trj_nvt: customDefaults.freq_trj_nvt,
                    thermo_freq: customDefaults.thermo_freq,
                  });
                } else {
                  // 切换到其他模式时，清空这些字段，让后端使用默认值
                  form.setFieldsValue({
                    nsteps_npt: undefined,
                    nsteps_nvt: undefined,
                    timestep: undefined,
                    temperature: undefined,
                    pressure: undefined,
                    freq_trj_npt: undefined,
                    freq_trj_nvt: undefined,
                    thermo_freq: undefined,
                  });
                }
              }}
            />
          </Form.Item>

          {selectedAccuracyLevel !== 'custom' && accuracyDefaults && (
            <Alert
              message="提示"
              description="选择精度等级后，系统会自动设置模拟参数。如果需要自定义参数，请选择「自定义」模式。"
              type="info"
              showIcon
              style={{ marginBottom: 16 }}
            />
          )}

          {selectedAccuracyLevel === 'custom' && (
            <>
              <Alert
                message="自定义模式"
                description="您选择了自定义模式。下方已自动填充标准模式的参数作为参考，您可以根据需要修改。"
                type="warning"
                showIcon
                style={{ marginBottom: 16 }}
              />

              <Divider orientation="left">模拟参数设置（必填）</Divider>

              <Row gutter={16}>
                <Col span={12}>
                  <Form.Item
                    name="nsteps_npt"
                    label="NPT 步数"
                    rules={[{ required: true, message: '请输入 NPT 步数' }]}
                    tooltip="等压等温系综的模拟步数"
                  >
                    <InputNumber min={1000} max={100000000} step={100000} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
                <Col span={12}>
                  <Form.Item
                    name="nsteps_nvt"
                    label="NVT 步数"
                    rules={[{ required: true, message: '请输入 NVT 步数' }]}
                    tooltip="等容等温系综的模拟步数"
                  >
                    <InputNumber min={1000} max={100000000} step={100000} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
              </Row>

              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item
                    name="timestep"
                    label="时间步长 (fs)"
                    rules={[{ required: true, message: '请输入时间步长' }]}
                  >
                    <InputNumber min={0.1} max={10} step={0.1} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="temperature"
                    label="温度 (K)"
                    rules={[{ required: true, message: '请输入温度' }]}
                  >
                    <InputNumber min={200} max={500} step={1} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="pressure"
                    label="压力 (atm)"
                    rules={[{ required: true, message: '请输入压力' }]}
                  >
                    <InputNumber min={0.1} max={100} step={0.1} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
              </Row>

              <Divider orientation="left">输出频率设置</Divider>

              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item
                    name="freq_trj_npt"
                    label="NPT 轨迹输出频率"
                    rules={[{ required: true, message: '请输入频率' }]}
                  >
                    <InputNumber min={100} max={10000000} step={100} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="freq_trj_nvt"
                    label="NVT 轨迹输出频率"
                    rules={[{ required: true, message: '请输入频率' }]}
                  >
                    <InputNumber min={100} max={10000000} step={100} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="thermo_freq"
                    label="热力学输出频率"
                    rules={[{ required: true, message: '请输入频率' }]}
                  >
                    <InputNumber min={100} max={10000000} step={100} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
              </Row>
            </>
          )}

          {/* QC量子化学计算选项 - 放在资源配置前面 */}
          <Divider orientation="left">
            <Space>
              <ExperimentOutlined style={{ color: '#722ed1' }} />
              量子化学计算 (可选)
            </Space>
          </Divider>

          <Card
            size="small"
            style={{
              marginBottom: 16,
              borderColor: '#d3adf7',
              background: 'linear-gradient(135deg, #f9f0ff 0%, #fff 100%)'
            }}
          >
            <Form.Item
              name="qc_enabled"
              valuePropName="checked"
              initialValue={false}
              style={{ marginBottom: 8 }}
            >
              <Checkbox>
                <Space>
                  <ExperimentOutlined style={{ color: '#722ed1' }} />
                  <Text strong>启用QC计算</Text>
                </Space>
              </Checkbox>
            </Form.Item>
            <Text type="secondary" style={{ fontSize: 12 }}>
              勾选后将对电解质中的溶剂分子进行量子化学计算，获取HOMO、LUMO、ESP等性质
            </Text>
          </Card>

          <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
            prevValues.qc_enabled !== currentValues.qc_enabled ||
            prevValues.electrolyte_id !== currentValues.electrolyte_id ||
            prevValues.qc_functionals !== currentValues.qc_functionals ||
            prevValues.qc_basis_sets !== currentValues.qc_basis_sets ||
            prevValues.qc_solvent_models !== currentValues.qc_solvent_models ||
            prevValues.qc_solvents !== currentValues.qc_solvents
          }>
            {({ getFieldValue }) => {
              const qcEnabled = getFieldValue('qc_enabled');
              if (!qcEnabled) return null;

              const electrolyteId = getFieldValue('electrolyte_id');
              const selectedElectrolyte = electrolytes.find(e => e.id === electrolyteId);

              // 收集将要计算的分子列表
              const moleculesToCalc: Array<{name: string, smiles: string, type: string, charge: number}> = [];
              if (selectedElectrolyte) {
                // 溶剂分子
                selectedElectrolyte.solvents?.forEach((sol: any) => {
                  if (sol.smiles && !moleculesToCalc.find(m => m.smiles === sol.smiles)) {
                    moleculesToCalc.push({ name: sol.name, smiles: sol.smiles, type: 'solvent', charge: 0 });
                  }
                });
                // 阳离子
                selectedElectrolyte.cations?.forEach((cat: any) => {
                  if (cat.smiles && !moleculesToCalc.find(m => m.smiles === cat.smiles)) {
                    moleculesToCalc.push({ name: cat.name, smiles: cat.smiles, type: 'cation', charge: 1 });
                  }
                });
                // 阴离子
                selectedElectrolyte.anions?.forEach((an: any) => {
                  if (an.smiles && !moleculesToCalc.find(m => m.smiles === an.smiles)) {
                    moleculesToCalc.push({ name: an.name, smiles: an.smiles, type: 'anion', charge: -1 });
                  }
                });
              }

              return (
                <Card size="small" style={{ marginBottom: 16 }}>
                  <Row gutter={16}>
                    <Col span={12}>
                      <Form.Item
                        name="qc_functionals"
                        label="泛函"
                        initialValue={['B3LYP']}
                        style={{ marginBottom: 8 }}
                        tooltip="可选择多个泛函进行对比计算"
                      >
                        <Select mode="multiple" placeholder="选择泛函（可多选）">
                          <Select.Option value="B3LYP">B3LYP (混合泛函)</Select.Option>
                          <Select.Option value="M062X">M06-2X (Minnesota泛函)</Select.Option>
                          <Select.Option value="wB97XD">ωB97X-D (长程校正)</Select.Option>
                          <Select.Option value="PBE0">PBE0 (混合GGA)</Select.Option>
                          <Select.Option value="CAM-B3LYP">CAM-B3LYP (长程校正)</Select.Option>
                        </Select>
                      </Form.Item>
                    </Col>
                    <Col span={12}>
                      <Form.Item
                        name="qc_basis_sets"
                        label="基组"
                        initialValue={['6-31++g(d,p)']}
                        style={{ marginBottom: 8 }}
                        tooltip="可选择多个基组进行对比计算"
                      >
                        <Select mode="multiple" placeholder="选择基组（可多选）">
                          <Select.Option value="6-31g(d,p)">6-31G(d,p) (标准)</Select.Option>
                          <Select.Option value="6-31++g(d,p)">6-31++G(d,p) (含弥散)</Select.Option>
                          <Select.Option value="6-311g(d,p)">6-311G(d,p) (三重劈裂)</Select.Option>
                          <Select.Option value="6-311++g(d,p)">6-311++G(d,p) (三重劈裂+弥散)</Select.Option>
                          <Select.Option value="Def2TZVP">Def2-TZVP (高精度)</Select.Option>
                        </Select>
                      </Form.Item>
                    </Col>
                  </Row>

                  <Row gutter={16}>
                    <Col span={12}>
                      <Form.Item
                        name="qc_solvent_models"
                        label="溶剂环境"
                        initialValue={['pcm']}
                        style={{ marginBottom: 8 }}
                        tooltip={
                          <div>
                            <p><strong>气相 (Gas)</strong>: 真空环境，无溶剂效应</p>
                            <p><strong>PCM</strong>: 极化连续介质模型，使用介电常数描述溶剂</p>
                            <p><strong>SMD</strong>: 溶剂密度模型，更精确但计算量更大</p>
                            <p>可多选进行对比计算</p>
                          </div>
                        }
                      >
                        <Select mode="multiple" placeholder="选择溶剂环境（可多选）">
                          <Select.Option value="gas">气相 (Gas Phase) - 无溶剂效应</Select.Option>
                          <Select.Option value="pcm">PCM - 极化连续介质模型</Select.Option>
                          <Select.Option value="smd">SMD - 溶剂密度模型（更精确）</Select.Option>
                        </Select>
                      </Form.Item>
                    </Col>
                    <Col span={12}>
                      <Form.Item
                        name="qc_solvents"
                        label="隐式溶剂"
                        initialValue={['Water']}
                        style={{ marginBottom: 8 }}
                        tooltip={
                          <div>
                            <p><strong>选择原则</strong>：选择介电常数(ε)接近您电解液的溶剂</p>
                            <hr style={{ margin: '4px 0', borderColor: 'rgba(255,255,255,0.3)' }} />
                            <p>• <strong>水系电解液</strong>: 选择 Water (ε=78.4)</p>
                            <p>• <strong>高浓电解液</strong>: 选择 Acetone (ε=20.5)</p>
                            <p>• <strong>EC基电解液</strong>: 选择 Water 或 PC (ε≈65-90)</p>
                            <p>• <strong>DMC/EMC/DEC电解液</strong>: 选择 Chloroform (ε≈3-5)</p>
                            <p>• <strong>离子液体</strong>: 选择 DMSO (ε=46.8)</p>
                          </div>
                        }
                      >
                        <Select mode="multiple" placeholder="选择隐式溶剂（可多选）" showSearch>
                          <Select.OptGroup label="📌 水系电解液 (ε>50)">
                            <Select.Option value="Water">水 (Water) ε=78.4</Select.Option>
                          </Select.OptGroup>
                          <Select.OptGroup label="📌 高介电常数碳酸酯 (ε=40-90)">
                            <Select.Option value="DiMethylSulfoxide">DMSO ε=46.8 (离子液体参考)</Select.Option>
                            <Select.Option value="1,2-EthaneDiol">乙二醇 ε=40.2</Select.Option>
                          </Select.OptGroup>
                          <Select.OptGroup label="📌 中等介电常数 (ε=15-40)">
                            <Select.Option value="Acetonitrile">乙腈 ε=35.7</Select.Option>
                            <Select.Option value="Methanol">甲醇 ε=32.6</Select.Option>
                            <Select.Option value="Ethanol">乙醇 ε=24.9</Select.Option>
                            <Select.Option value="Acetone">丙酮 ε=20.5 (高浓电解液)</Select.Option>
                            <Select.Option value="1-Propanol">正丙醇 ε=20.5</Select.Option>
                          </Select.OptGroup>
                          <Select.OptGroup label="📌 低介电常数 (ε<15) - DMC/EMC/DEC体系">
                            <Select.Option value="DiChloroEthane">二氯乙烷 ε=10.1</Select.Option>
                            <Select.Option value="Dichloromethane">二氯甲烷 ε=8.9</Select.Option>
                            <Select.Option value="TetraHydroFuran">四氢呋喃 (THF) ε=7.4</Select.Option>
                            <Select.Option value="Chloroform">氯仿 ε=4.7 (线性碳酸酯参考)</Select.Option>
                            <Select.Option value="DiethylEther">乙醚 ε=4.2</Select.Option>
                            <Select.Option value="CarbonTetraChloride">四氯化碳 ε=2.2</Select.Option>
                            <Select.Option value="Toluene">甲苯 ε=2.4</Select.Option>
                            <Select.Option value="Benzene">苯 ε=2.3</Select.Option>
                          </Select.OptGroup>
                        </Select>
                      </Form.Item>
                    </Col>
                  </Row>

                  {/* 溶剂选择提示 */}
                  <Alert
                    type="info"
                    showIcon
                    style={{ marginBottom: 12 }}
                    message={
                      <Text style={{ fontSize: 12 }}>
                        <strong>隐式溶剂选择提示：</strong>选择介电常数(ε)接近您电解液的溶剂。
                        例如：EC体系选Water(ε≈78)，DMC/EMC体系选Chloroform(ε≈4.7)，高浓电解液选Acetone(ε≈20)。
                      </Text>
                    }
                  />

                  {/* 显示将要计算的分子列表和任务数量 */}
                  <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
                    prevValues.qc_functionals !== currentValues.qc_functionals ||
                    prevValues.qc_basis_sets !== currentValues.qc_basis_sets ||
                    prevValues.qc_solvent_models !== currentValues.qc_solvent_models ||
                    prevValues.qc_solvents !== currentValues.qc_solvents
                  }>
                    {({ getFieldValue: getFieldValueInner }) => {
                      const functionals = getFieldValueInner('qc_functionals') || ['B3LYP'];
                      const basisSets = getFieldValueInner('qc_basis_sets') || ['6-31++g(d,p)'];
                      const solventModels = getFieldValueInner('qc_solvent_models') || ['pcm'];
                      const solvents = getFieldValueInner('qc_solvents') || ['Water'];

                      // 计算溶剂组合数
                      let solventCombinations = 0;
                      if (solventModels.includes('gas')) {
                        solventCombinations += 1;
                      }
                      const nonGasModels = solventModels.filter((m: string) => m !== 'gas');
                      solventCombinations += nonGasModels.length * solvents.length;

                      const totalJobs = moleculesToCalc.length * functionals.length * basisSets.length * solventCombinations;

                      return moleculesToCalc.length > 0 ? (
                        <Alert
                          type="info"
                          showIcon
                          style={{ marginTop: 8 }}
                          message={
                            <div>
                              <strong>将创建 {totalJobs} 个 QC 任务</strong>
                              <Text type="secondary" style={{ fontSize: 12, marginLeft: 8 }}>
                                ({moleculesToCalc.length} 分子 × {functionals.length} 泛函 × {basisSets.length} 基组 × {solventCombinations} 溶剂组合)
                              </Text>
                            </div>
                          }
                          description={
                            <div style={{ marginTop: 8 }}>
                              <div style={{ marginBottom: 8 }}>
                                <Text strong style={{ fontSize: 12 }}>分子列表：</Text>
                              </div>
                              {moleculesToCalc.map((mol, index) => (
                                <div key={index} style={{
                                  display: 'inline-block',
                                  marginRight: 8,
                                  marginBottom: 4,
                                  padding: '2px 8px',
                                  background: mol.type === 'solvent' ? '#f6ffed' :
                                             mol.type === 'cation' ? '#fff2f0' : '#f0f5ff',
                                  borderRadius: 4,
                                  border: `1px solid ${mol.type === 'solvent' ? '#b7eb8f' :
                                                      mol.type === 'cation' ? '#ffccc7' : '#adc6ff'}`
                                }}>
                                  <Text style={{ fontSize: 12 }}>
                                    {mol.name}
                                    <Text type="secondary" style={{ fontSize: 11, marginLeft: 4 }}>
                                      ({mol.type === 'solvent' ? '溶剂' :
                                        mol.type === 'cation' ? '阳离子' : '阴离子'})
                                    </Text>
                                  </Text>
                                </div>
                              ))}
                            </div>
                          }
                        />
                      ) : null;
                    }}
                  </Form.Item>
                </Card>
              );
            }}
          </Form.Item>

          <Divider orientation="left">
            计算资源配置
            <Button
              type="link"
              size="small"
              icon={<ThunderboltOutlined />}
              onClick={() => handleGetSuggestion(form)}
              style={{ marginLeft: 8 }}
            >
              获取推荐配置
            </Button>
          </Divider>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                name="slurm_partition"
                label="队列/分区"
                tooltip="Slurm 队列名称"
              >
                <Select>
                  {partitions.length > 0 ? (
                    partitions.map(p => (
                      <Select.Option key={p.name} value={p.name}>
                        {p.name} ({p.state === 'up' ? `可用 ${p.available_cpus} CPUs` : '不可用'})
                      </Select.Option>
                    ))
                  ) : (
                    <>
                      <Select.Option value="cpu">cpu</Select.Option>
                      <Select.Option value="gpu">gpu</Select.Option>
                    </>
                  )}
                </Select>
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item
                name="slurm_nodes"
                label="节点数"
                initialValue={1}
                tooltip="使用的计算节点数量"
              >
                <InputNumber min={1} max={10} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
          </Row>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                name="slurm_ntasks"
                label="任务数"
                initialValue={8}
                tooltip="Slurm 任务数（通常对应 MPI 进程数的一部分）"
              >
                <InputNumber min={1} max={128} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item
                name="slurm_cpus_per_task"
                label="每任务 CPU 数"
                initialValue={8}
                tooltip="每个任务使用的 CPU 核心数"
              >
                <InputNumber min={1} max={64} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
          </Row>

          <Form.Item
            name="slurm_time"
            label="最大运行时间 (分钟)"
            initialValue={7200}
            tooltip="任务的最大运行时间，超时将被终止"
          >
            <InputNumber min={60} max={43200} step={60} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
            prevValues.slurm_ntasks !== currentValues.slurm_ntasks ||
            prevValues.slurm_cpus_per_task !== currentValues.slurm_cpus_per_task
          }>
            {({ getFieldValue }) => {
              const ntasks = getFieldValue('slurm_ntasks') || 8;
              const cpusPerTask = getFieldValue('slurm_cpus_per_task') || 8;
              const totalProcesses = ntasks * cpusPerTask;

              return (
                <Alert
                  message="总 MPI 进程数 = 任务数 × 每任务 CPU 数"
                  description={`当前配置将使用 ${ntasks} × ${cpusPerTask} = ${totalProcesses} 个 MPI 进程`}
                  type="info"
                  showIcon
                  style={{ marginTop: 16 }}
                />
              );
            }}
          </Form.Item>
        </Form>
      </Modal>

      {/* 重新提交任务对话框 */}
      <Modal
        title={
          <Space>
            <ReloadOutlined style={{ color: '#1677ff' }} />
            {`重新提交任务 - ${resubmittingJob?.config?.job_name || ''}`}
          </Space>
        }
        open={resubmitModalVisible}
        onOk={handleResubmitSubmit}
        onCancel={handleCloseResubmitModal}
        okText="重新提交"
        cancelText="取消"
        width={800}
        centered
        styles={{
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' },
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
        }}
      >
        <Form form={resubmitForm} layout="vertical">
          <Divider orientation="left">模拟参数</Divider>

          <Form.Item
            name="nsteps_npt"
            label="NPT 步数"
            rules={[{ required: true, message: '请输入 NPT 步数' }]}
          >
            <InputNumber min={1000} max={100000000} step={1000} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="nsteps_nvt"
            label="NVT 步数"
            rules={[{ required: true, message: '请输入 NVT 步数' }]}
          >
            <InputNumber min={1000} max={100000000} step={1000} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="timestep"
            label="时间步长 (fs)"
            rules={[{ required: true, message: '请输入时间步长' }]}
          >
            <InputNumber min={0.1} max={10} step={0.1} style={{ width: '100%' }} />
          </Form.Item>

          <Divider orientation="left">
            计算资源配置
            <Button
              type="link"
              size="small"
              icon={<ThunderboltOutlined />}
              onClick={() => handleGetSuggestion(resubmitForm)}
              style={{ marginLeft: 8 }}
            >
              获取推荐配置
            </Button>
          </Divider>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                name="slurm_partition"
                label="队列/分区"
                tooltip="Slurm 队列名称"
              >
                <Select>
                  {partitions.length > 0 ? (
                    partitions.map(p => (
                      <Select.Option key={p.name} value={p.name}>
                        {p.name} ({p.state === 'up' ? `可用 ${p.available_cpus} CPUs` : '不可用'})
                      </Select.Option>
                    ))
                  ) : (
                    <>
                      <Select.Option value="cpu">cpu</Select.Option>
                      <Select.Option value="gpu">gpu</Select.Option>
                    </>
                  )}
                </Select>
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item
                name="slurm_nodes"
                label="节点数"
                tooltip="使用的计算节点数量"
              >
                <InputNumber min={1} max={10} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
          </Row>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                name="slurm_ntasks"
                label="任务数"
                tooltip="Slurm 任务数（通常对应 MPI 进程数的一部分）"
              >
                <InputNumber min={1} max={128} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item
                name="slurm_cpus_per_task"
                label="每任务 CPU 数"
                tooltip="每个任务使用的 CPU 核心数"
              >
                <InputNumber min={1} max={64} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
          </Row>

          <Form.Item
            name="slurm_time"
            label="最大运行时间 (分钟)"
            tooltip="任务的最大运行时间，超时将被终止"
          >
            <InputNumber min={60} max={43200} step={60} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
            prevValues.slurm_ntasks !== currentValues.slurm_ntasks ||
            prevValues.slurm_cpus_per_task !== currentValues.slurm_cpus_per_task
          }>
            {({ getFieldValue }) => {
              const ntasks = getFieldValue('slurm_ntasks') || 8;
              const cpusPerTask = getFieldValue('slurm_cpus_per_task') || 8;
              const totalProcesses = ntasks * cpusPerTask;

              return (
                <Alert
                  message="总 MPI 进程数 = 任务数 × 每任务 CPU 数"
                  description={`当前配置将使用 ${ntasks} × ${cpusPerTask} = ${totalProcesses} 个 MPI 进程`}
                  type="info"
                  showIcon
                  style={{ marginTop: 16 }}
                />
              );
            }}
          </Form.Item>
        </Form>
      </Modal>

      {/* 新建配方对话框 */}
      <Modal
        title={
          <Space>
            <ExperimentOutlined style={{ color: '#1677ff' }} />
            新建电解质配方
          </Space>
        }
        open={electrolyteModalVisible}
        onOk={handleCreateElectrolyte}
        onCancel={handleCloseElectrolyteModal}
        okText="创建"
        cancelText="取消"
        width={900}
        centered
        styles={{
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' },
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
        }}
      >
        <ElectrolyteFormOptimized
          form={electrolyteForm}
          projects={projects}
          initialCations={selectedCations}
          initialAnions={selectedAnions}
          onIonsChange={(cations, anions) => {
            setSelectedCations(cations);
            setSelectedAnions(anions);
          }}
        />
      </Modal>
    </div>
  );
}


