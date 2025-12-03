/**
 * 电解质体系管理页面
 */
import { useState, useEffect } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import {
  Button,
  Input,
  Space,
  message,
  Modal,
  Form,
  Row,
  Col,
  Spin,
  Empty,
  Tabs,
  Typography,
  Card,
  Select,
  Dropdown,
  Steps,
  Table,
  Tag,
  Alert,
  Radio,
  Divider,
  Result,
  InputNumber,
  Checkbox,
  DatePicker,
  Statistic,
  Switch,
  theme,
} from 'antd';
import type { MenuProps } from 'antd';
import {
  PlusOutlined,
  SearchOutlined,
  ReloadOutlined,
  ExperimentOutlined,
  DeleteOutlined,
  CheckSquareOutlined,
  CloseSquareOutlined,
  FolderOutlined,
  MoreOutlined,
  UploadOutlined,
  DownloadOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  FileExcelOutlined,
  InboxOutlined,
  FilterOutlined,
  AppstoreOutlined,
  UnorderedListOutlined,
} from '@ant-design/icons';
import ElectrolyteCard from '../components/ElectrolyteCard';
import ElectrolyteFormOptimized from '../components/ElectrolyteFormOptimized';
import {
  getElectrolytes,
  createElectrolyteNew,
  updateElectrolyteNew,
  getElectrolyteEditable,
  deleteElectrolyte,
  batchDeleteElectrolytes,
  batchUpdateProject,
} from '../api/electrolytes';
import { getMDJobs, batchCreateMDJobs } from '../api/jobs';
import { getProjects, createProject } from '../api/projects';
import { downloadTemplate, batchImportUpload, BatchImportResult } from '../api/batchImport';
import type { ElectrolyteSystem, Project, MDJob, ProjectCreate } from '../types';
import { JobStatus } from '../types';
import { useThemeStore } from '../stores/themeStore';
import dayjs, { Dayjs } from 'dayjs';

const { Title, Text } = Typography;
const { Step } = Steps;
const { RangePicker } = DatePicker;

export default function Electrolytes() {
  const location = useLocation();
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const [electrolytes, setElectrolytes] = useState<ElectrolyteSystem[]>([]);
  const [filteredElectrolytes, setFilteredElectrolytes] = useState<ElectrolyteSystem[]>([]);
  const [projects, setProjects] = useState<Project[]>([]);
  const [jobs, setJobs] = useState<MDJob[]>([]);
  const [loading, setLoading] = useState(false);
  const [searchText, setSearchText] = useState('');
  const [activeTab, setActiveTab] = useState('all');
  const [modalVisible, setModalVisible] = useState(false);
  const [editingElectrolyte, setEditingElectrolyte] = useState<ElectrolyteSystem | null>(null);
  const [copyingElectrolyte, setCopyingElectrolyte] = useState<ElectrolyteSystem | null>(null);
  const [form] = Form.useForm();
  const [selectedCations, setSelectedCations] = useState<any[]>([]);
  const [selectedAnions, setSelectedAnions] = useState<any[]>([]);

  // 筛选和视图状态
  const [projectFilter, setProjectFilter] = useState<number | undefined>(undefined);
  const [solventFilter, setSolventFilter] = useState<string | undefined>(undefined);
  const [ionFilter, setIonFilter] = useState<string | undefined>(undefined);
  const [dateRange, setDateRange] = useState<[Dayjs | null, Dayjs | null] | null>(null);
  const [viewMode, setViewMode] = useState<'card' | 'table'>(() => {
    // 从localStorage读取视图模式
    const saved = localStorage.getItem('electrolytes-view-mode');
    return (saved === 'card' || saved === 'table') ? saved : 'card';
  });
  const [sortBy, setSortBy] = useState<'created_at' | 'name'>('created_at');
  const [sortOrder, setSortOrder] = useState<'asc' | 'desc'>('desc');

  // 保存视图模式到localStorage
  const handleViewModeChange = (mode: 'card' | 'table') => {
    setViewMode(mode);
    localStorage.setItem('electrolytes-view-mode', mode);
  };

  // 批量选择相关状态
  const [selectMode, setSelectMode] = useState(false);
  const [selectedIds, setSelectedIds] = useState<number[]>([]);

  // 批量导入相关状态
  const [batchImportVisible, setBatchImportVisible] = useState(false);
  const [uploadFile, setUploadFile] = useState<File | null>(null);
  const [importing, setImporting] = useState(false);
  const [importStep, setImportStep] = useState(0);
  const [importResult, setImportResult] = useState<BatchImportResult | null>(null);
  const [projectMode, setProjectMode] = useState<'existing' | 'new'>('new');
  const [selectedProjectId, setSelectedProjectId] = useState<number | undefined>(undefined);
  const [newProjectName, setNewProjectName] = useState('');
  const [newProjectDesc, setNewProjectDesc] = useState('');

  // 批量创建MD任务相关状态
  const [batchMDModalVisible, setBatchMDModalVisible] = useState(false);
  const [batchMDForm] = Form.useForm();
  const [creatingBatchMD, setCreatingBatchMD] = useState(false);
  const [batchMDAccuracyLevel, setBatchMDAccuracyLevel] = useState<string>('standard');
  const [batchMDPartitions, setBatchMDPartitions] = useState<any[]>([]);
  const [batchMDQuota, setBatchMDQuota] = useState<any>(null);
  const [createdMDJobIds, setCreatedMDJobIds] = useState<number[]>([]);
  const [batchMDSubmitToCluster, setBatchMDSubmitToCluster] = useState<boolean>(false);

  // 新建项目相关状态
  const [projectModalVisible, setProjectModalVisible] = useState(false);
  const [projectForm] = Form.useForm();

  // 加载数据
  const loadData = async () => {
    setLoading(true);
    try {
      const [electrolyteData, projectData, jobData] = await Promise.all([
        getElectrolytes(),
        getProjects(),
        getMDJobs(),
      ]);
      setElectrolytes(electrolyteData);
      setFilteredElectrolytes(electrolyteData);
      setProjects(projectData);
      setJobs(jobData);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadData();
  }, []);

  // 获取配方的状态分类
  const getElectrolyteCategory = (electrolyte: ElectrolyteSystem) => {
    const relatedJobs = jobs.filter(job => job.system_id === electrolyte.id);

    if (relatedJobs.length === 0) {
      return 'draft'; // 草稿：没有任务
    }

    const hasRunning = relatedJobs.some(job =>
      [JobStatus.CREATED, JobStatus.QUEUED, JobStatus.RUNNING, JobStatus.POSTPROCESSING].includes(job.status)
    );

    if (hasRunning) {
      return 'running'; // 进行中：有任务在运行
    }

    return 'completed'; // 已完成：所有任务都已完成
  };

  // 获取配方中的所有溶剂
  const getSolvents = (electrolyte: ElectrolyteSystem): string[] => {
    if (!electrolyte.solvents || electrolyte.solvents.length === 0) {
      return [];
    }
    return electrolyte.solvents.map(s => s.name);
  };

  // 获取配方中的所有离子
  const getIons = (electrolyte: ElectrolyteSystem): string[] => {
    const ions: string[] = [];
    if (electrolyte.cations) {
      ions.push(...electrolyte.cations.map(c => c.name));
    }
    if (electrolyte.anions) {
      ions.push(...electrolyte.anions.map(a => a.name));
    }
    return ions;
  };

  // 搜索和分类过滤
  useEffect(() => {
    let filtered = [...electrolytes];

    // 按状态标签页筛选
    if (activeTab !== 'all') {
      filtered = filtered.filter((e) => getElectrolyteCategory(e) === activeTab);
    }

    // 搜索过滤（配方名称、组成成分）
    if (searchText) {
      const search = searchText.toLowerCase();
      filtered = filtered.filter((e) => {
        const nameMatch = e.name.toLowerCase().includes(search);
        const solvents = getSolvents(e);
        const ions = getIons(e);
        const solventMatch = solvents.some(s => s.toLowerCase().includes(search));
        const ionMatch = ions.some(i => i.toLowerCase().includes(search));
        return nameMatch || solventMatch || ionMatch;
      });
    }

    // 项目筛选
    if (projectFilter !== undefined) {
      filtered = filtered.filter((e) => e.project_id === projectFilter);
    }

    // 溶剂筛选
    if (solventFilter) {
      filtered = filtered.filter((e) => {
        const solvents = getSolvents(e);
        return solvents.some(s => s === solventFilter);
      });
    }

    // 离子筛选
    if (ionFilter) {
      filtered = filtered.filter((e) => {
        const ions = getIons(e);
        return ions.some(i => i === ionFilter);
      });
    }

    // 时间范围筛选
    if (dateRange && dateRange[0] && dateRange[1]) {
      const startDate = dateRange[0].startOf('day');
      const endDate = dateRange[1].endOf('day');
      filtered = filtered.filter((e) => {
        const createdDate = dayjs(e.created_at);
        return createdDate.isAfter(startDate) && createdDate.isBefore(endDate);
      });
    }

    // 排序
    filtered.sort((a, b) => {
      let aValue: any, bValue: any;

      if (sortBy === 'created_at') {
        aValue = new Date(a.created_at).getTime();
        bValue = new Date(b.created_at).getTime();
        return sortOrder === 'asc' ? aValue - bValue : bValue - aValue;
      } else {
        // name
        aValue = a.name.toLowerCase();
        bValue = b.name.toLowerCase();
        return sortOrder === 'asc'
          ? aValue.localeCompare(bValue)
          : bValue.localeCompare(aValue);
      }
    });

    setFilteredElectrolytes(filtered);
  }, [searchText, electrolytes, jobs, activeTab, projectFilter, solventFilter, ionFilter, dateRange, sortBy, sortOrder]);

  // 重置筛选
  const handleResetFilters = () => {
    setSearchText('');
    setProjectFilter(undefined);
    setSolventFilter(undefined);
    setIonFilter(undefined);
    setDateRange(null);
    setSortBy('created_at');
    setSortOrder('desc');
  };

  // 获取所有唯一的溶剂
  const getAllSolvents = (): string[] => {
    const solventsSet = new Set<string>();
    electrolytes.forEach(e => {
      getSolvents(e).forEach(s => solventsSet.add(s));
    });
    return Array.from(solventsSet).sort();
  };

  // 获取所有唯一的离子
  const getAllIons = (): string[] => {
    const ionsSet = new Set<string>();
    electrolytes.forEach(e => {
      getIons(e).forEach(i => ionsSet.add(i));
    });
    return Array.from(ionsSet).sort();
  };

  // 检查是否需要自动打开创建对话框
  useEffect(() => {
    if (location.state?.openCreateModal) {
      // 如果有指定项目ID，设置默认项目
      if (location.state?.projectId) {
        form.setFieldsValue({ project_id: location.state.projectId });
      }
      setModalVisible(true);
      // 清除 state，避免刷新时重复打开
      window.history.replaceState({}, document.title);
    }
  }, [location]);

  // 打开创建/编辑/复制对话框
  const handleOpenModal = async (electrolyte?: ElectrolyteSystem, isCopy: boolean = false) => {
    if (electrolyte) {
      if (isCopy) {
        setCopyingElectrolyte(electrolyte);
      } else {
        setEditingElectrolyte(electrolyte);
      }
      setModalVisible(true);

      // 获取可编辑格式的数据
      try {
        const editableData = await getElectrolyteEditable(electrolyte.id);

        // 设置表单值
        form.setFieldsValue({
          project_id: editableData.project_id,
          name: isCopy ? `${editableData.name} (副本)` : editableData.name,
          description: editableData.description || '',
          temperature: editableData.temperature,
          box_type: editableData.box.type,
          box_size: editableData.box.type === 'cubic' ? editableData.box.dimensions[0] : undefined,
          box_dimensions: editableData.box.type === 'rectangular' ? editableData.box.dimensions : undefined,
          solvents: editableData.solvents || [],
        });

        // 设置选中的离子
        setSelectedCations(editableData.cations || []);
        setSelectedAnions(editableData.anions || []);
      } catch (error: any) {
        message.error('加载配方数据失败: ' + (error.response?.data?.detail || error.message));
        setModalVisible(false);
        setEditingElectrolyte(null);
      }
    } else {
      setEditingElectrolyte(null);
      form.resetFields();
      setSelectedCations([]);
      setSelectedAnions([]);
      // 设置默认值
      form.setFieldsValue({
        temperature: 298.15,
        box_type: 'cubic',
        box_size: 40,
        solvents: [],
      });
      setModalVisible(true);
    }
  };

  // 关闭对话框
  const handleCloseModal = () => {
    setModalVisible(false);
    setEditingElectrolyte(null);
    setCopyingElectrolyte(null);
    form.resetFields();
  };

  // 打开新建项目对话框
  const handleOpenProjectModal = () => {
    projectForm.resetFields();
    setProjectModalVisible(true);
  };

  // 关闭新建项目对话框
  const handleCloseProjectModal = () => {
    setProjectModalVisible(false);
    projectForm.resetFields();
  };

  // 创建新项目
  const handleCreateProject = async () => {
    try {
      const values = await projectForm.validateFields();
      const newProject = await createProject(values as ProjectCreate);
      message.success('项目创建成功');

      // 重新加载项目列表
      const projectData = await getProjects();
      setProjects(projectData);

      // 自动选择新创建的项目
      form.setFieldsValue({ project_id: newProject.id });

      handleCloseProjectModal();
    } catch (error: any) {
      if (error.response) {
        message.error(error.response?.data?.detail || '创建项目失败');
      }
    }
  };

  // 提交表单
  const handleSubmit = async () => {
    try {
      // 验证至少选择了一个阳离子和一个阴离子
      if (selectedCations.length === 0) {
        message.error('请至少选择一种阳离子');
        return;
      }
      if (selectedAnions.length === 0) {
        message.error('请至少选择一种阴离子');
        return;
      }

      // 获取当前表单值（不验证）
      const allValues = form.getFieldsValue();
      const boxType = allValues.box_type || 'cubic';

      // 根据 box_type 清理不需要的字段
      if (boxType === 'cubic') {
        form.setFieldsValue({ box_dimensions: undefined });
      } else {
        form.setFieldsValue({ box_size: undefined });
      }

      // 现在验证表单
      let values;
      try {
        values = await form.validateFields();
      } catch (error: any) {
        console.error('表单验证失败:', error);
        if (error.errorFields && error.errorFields.length > 0) {
          console.error('失败的字段:', error.errorFields.map((f: any) => ({
            name: f.name,
            errors: f.errors
          })));
          message.error(`表单验证失败: ${error.errorFields[0].errors[0]}`);
        }
        return;
      }

      // 构建盒子配置
      const box = {
        type: boxType,
        dimensions: boxType === 'cubic'
          ? [values.box_size || 40]
          : (values.box_dimensions || [40, 40, 40]),
      };

      // 构建新格式数据
      const data: any = {
        project_id: values.project_id,
        name: values.name,
        description: values.description,
        cations: selectedCations,
        anions: selectedAnions,
        solvents: values.solvents || [],
        box: box,
        temperature: values.temperature || 298.15,
        pressure: 1.0,
        nsteps_npt: 5000000,
        nsteps_nvt: 10000000,
        timestep: 1.0,
        force_field: 'OPLS',
      };

      console.log('=== 创建电解质请求数据 ===');
      console.log('data:', JSON.stringify(data, null, 2));
      console.log('selectedCations:', selectedCations);
      console.log('selectedAnions:', selectedAnions);
      console.log('values.solvents:', values.solvents);

      if (editingElectrolyte && !copyingElectrolyte) {
        await updateElectrolyteNew(editingElectrolyte.id, data);
        message.success('电解质配方更新成功');
      } else {
        await createElectrolyteNew(data);
        message.success(copyingElectrolyte ? '电解质配方复制成功' : '电解质配方创建成功');
      }
      handleCloseModal();
      loadData();
    } catch (error: any) {
      console.error('=== 创建电解质失败 ===');
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
          message.error(detail || '操作失败');
        }
      } else {
        console.error('Submit error:', error);
      }
    }
  };

  // 删除电解质配方
  const handleDelete = async (id: number) => {
    try {
      await deleteElectrolyte(id);
      message.success('电解质配方删除成功');
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  // 创建任务 - 跳转到配置页面
  const handleCreateJob = (electrolyte: ElectrolyteSystem) => {
    // 跳转到任务配置页面，传递电解质系统信息
    navigate(`/workspace/jobs/create/${electrolyte.id}`, {
      state: { electrolyte }
    });
  };

  // 批量删除
  const handleBatchDelete = async () => {
    if (selectedIds.length === 0) {
      message.warning('请先选择要删除的配方');
      return;
    }

    Modal.confirm({
      title: '批量删除确认',
      content: `确定要删除选中的 ${selectedIds.length} 个配方吗？此操作不可恢复。`,
      okText: '确定删除',
      okButtonProps: { danger: true },
      cancelText: '取消',
      onOk: async () => {
        try {
          const result = await batchDeleteElectrolytes(selectedIds);

          if (result.deleted_count > 0) {
            message.success(result.message);
          }

          // 如果有失败的，显示详细原因
          if (result.failed_ids && result.failed_ids.length > 0) {
            const failedReasons = result.failed_reasons || {};
            const reasonMessages = Object.entries(failedReasons)
              .map(([id, reason]) => `配方ID ${id}: ${reason}`)
              .join('\n');

            if (result.deleted_count === 0) {
              message.error(`删除失败：\n${reasonMessages}`);
            } else {
              message.warning(`部分删除失败：\n${reasonMessages}`);
            }
          }

          setSelectedIds([]);
          setSelectMode(false);
          loadData();
        } catch (error: any) {
          message.error(error.response?.data?.detail || '批量删除失败');
        }
      },
    });
  };

  // 切换选择状态
  const toggleSelect = (id: number) => {
    setSelectedIds(prev =>
      prev.includes(id)
        ? prev.filter(i => i !== id)
        : [...prev, id]
    );
  };

  // 全选/取消全选
  const toggleSelectAll = () => {
    if (selectedIds.length === filteredElectrolytes.length) {
      setSelectedIds([]);
    } else {
      setSelectedIds(filteredElectrolytes.map(e => e.id));
    }
  };

  // 批量更改项目归属
  const handleBatchChangeProject = async (projectId: number) => {
    if (selectedIds.length === 0) {
      message.warning('请先选择要移动的配方');
      return;
    }

    const targetProject = projects.find(p => p.id === projectId);
    Modal.confirm({
      title: '批量移动确认',
      content: `确定要将选中的 ${selectedIds.length} 个配方移动到项目 "${targetProject?.name}" 吗？`,
      okText: '确定移动',
      cancelText: '取消',
      onOk: async () => {
        try {
          const result = await batchUpdateProject(selectedIds, projectId);
          message.success(result.message);
          setSelectedIds([]);
          setSelectMode(false);
          loadData();
        } catch (error: any) {
          message.error(error.response?.data?.detail || '批量移动失败');
        }
      },
    });
  };

  // 批量导入处理函数
  const handleDownloadTemplate = async () => {
    try {
      const blob = await downloadTemplate('electrolyte', true);
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = '配方批量导入模板.xlsx';
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

    // 验证项目设置
    if (projectMode === 'existing' && !selectedProjectId) {
      message.error('请选择一个项目');
      return;
    }
    if (projectMode === 'new' && !newProjectName.trim()) {
      message.error('请输入项目名称');
      return;
    }

    setImporting(true);
    try {
      const result = await batchImportUpload(
        uploadFile,
        projectMode === 'existing' ? selectedProjectId : undefined,
        projectMode === 'new' ? newProjectName : undefined,
        projectMode === 'new' ? newProjectDesc : undefined,
        '配方'
      );
      setImportResult(result);
      setImportStep(2); // 跳转到结果页面
      loadData();
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
    setProjectMode('new');
    setSelectedProjectId(undefined);
    setNewProjectName('');
    setNewProjectDesc('');
    setBatchImportVisible(false);
  };

  // 批量操作菜单项
  const batchMenuItems: MenuProps['items'] = [
    {
      key: 'delete',
      label: '批量删除',
      icon: <DeleteOutlined />,
      danger: true,
      onClick: handleBatchDelete,
    },
    {
      key: 'project',
      label: '移动到项目',
      icon: <FolderOutlined />,
      children: projects.map(p => ({
        key: `project-${p.id}`,
        label: p.name,
        onClick: () => handleBatchChangeProject(p.id),
      })),
    },
  ];

  // 计算各分类的数量
  const draftCount = electrolytes.filter(e => getElectrolyteCategory(e) === 'draft').length;
  const runningCount = electrolytes.filter(e => getElectrolyteCategory(e) === 'running').length;
  const completedCount = electrolytes.filter(e => getElectrolyteCategory(e) === 'completed').length;

  return (
    <div style={{
      padding: '24px',
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      position: 'relative',
      transition: 'background 0.3s',
    }}>
      {/* 浮动批量操作栏 */}
      {selectMode && selectedIds.length > 0 && (
        <div style={{
          position: 'fixed',
          bottom: 32,
          left: '50%',
          transform: 'translateX(-50%)',
          zIndex: 1000,
          background: token.colorBgContainer,
          padding: '16px 24px',
          borderRadius: 16,
          boxShadow: isDark ? '0 8px 32px rgba(0, 0, 0, 0.4)' : '0 8px 32px rgba(0, 0, 0, 0.12), 0 2px 8px rgba(0, 0, 0, 0.08)',
          border: `1px solid ${token.colorBorder}`,
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
              已选择 <span style={{ color: '#1677ff', fontSize: 16 }}>{selectedIds.length}</span> 项
            </Text>
          </div>
          <Space size={12}>
            <Button
              icon={selectedIds.length === filteredElectrolytes.length ? <CloseSquareOutlined /> : <CheckSquareOutlined />}
              onClick={toggleSelectAll}
              style={{ borderRadius: 8 }}
            >
              {selectedIds.length === filteredElectrolytes.length ? '取消全选' : '全选'}
            </Button>
            <Dropdown
              menu={{ items: batchMenuItems }}
              trigger={['click']}
            >
              <Button
                type="primary"
                icon={<MoreOutlined />}
                style={{ borderRadius: 8 }}
              >
                批量操作
              </Button>
            </Dropdown>
            <Button
              onClick={() => { setSelectMode(false); setSelectedIds([]); }}
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
              <ExperimentOutlined style={{ marginRight: 12, color: '#1677ff' }} />
              配方管理
            </Title>
            <Text type="secondary">创建和管理电解质配方，配置分子组成和模拟参数</Text>
          </div>
          <Space>
            <Button
              icon={<ReloadOutlined />}
              onClick={loadData}
              style={{ borderRadius: 8 }}
            >
              刷新
            </Button>
            {selectMode ? (
              <Button
                onClick={() => { setSelectMode(false); setSelectedIds([]); }}
                style={{ borderRadius: 8 }}
              >
                退出批量管理
              </Button>
            ) : (
              <>
                <Button
                  icon={<CheckSquareOutlined />}
                  onClick={() => setSelectMode(true)}
                  style={{ borderRadius: 8 }}
                >
                  批量管理
                </Button>
                <Button
                  icon={<UploadOutlined />}
                  onClick={() => setBatchImportVisible(true)}
                  style={{ borderRadius: 8 }}
                >
                  批量导入
                </Button>
                <Button
                  type="primary"
                  icon={<PlusOutlined />}
                  onClick={() => handleOpenModal()}
                  size="large"
                  style={{
                    borderRadius: 8,
                    boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)',
                  }}
                >
                  创建新配方
                </Button>
              </>
            )}
          </Space>
        </div>
      </div>

      {/* 搜索和统计卡片 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
        }}
      >
        <Row gutter={24} align="middle">
          <Col flex="auto">
            <Input
              placeholder="搜索配方名称..."
              prefix={<SearchOutlined style={{ color: '#bfbfbf' }} />}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              style={{ maxWidth: 400, borderRadius: 8 }}
              allowClear
              size="large"
            />
          </Col>
          <Col>
            <Space size={24}>
              <div style={{ textAlign: 'center' }}>
                <div style={{
                  fontSize: 24,
                  fontWeight: 700,
                  color: '#1677ff',
                  lineHeight: 1.2
                }}>
                  {electrolytes.length}
                </div>
                <Text type="secondary" style={{ fontSize: 12 }}>全部配方</Text>
              </div>
              <div style={{ textAlign: 'center' }}>
                <div style={{
                  fontSize: 24,
                  fontWeight: 700,
                  color: '#faad14',
                  lineHeight: 1.2
                }}>
                  {draftCount}
                </div>
                <Text type="secondary" style={{ fontSize: 12 }}>草稿</Text>
              </div>
              <div style={{ textAlign: 'center' }}>
                <div style={{
                  fontSize: 24,
                  fontWeight: 700,
                  color: '#52c41a',
                  lineHeight: 1.2
                }}>
                  {runningCount}
                </div>
                <Text type="secondary" style={{ fontSize: 12 }}>进行中</Text>
              </div>
              <div style={{ textAlign: 'center' }}>
                <div style={{
                  fontSize: 24,
                  fontWeight: 700,
                  color: '#722ed1',
                  lineHeight: 1.2
                }}>
                  {completedCount}
                </div>
                <Text type="secondary" style={{ fontSize: 12 }}>已完成</Text>
              </div>
            </Space>
          </Col>
        </Row>
      </Card>

      {/* 分类 Tabs */}
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
            {
              key: 'all',
              label: `全部 (${electrolytes.length})`,
            },
            {
              key: 'draft',
              label: `草稿 (${draftCount})`,
            },
            {
              key: 'running',
              label: `进行中 (${runningCount})`,
            },
            {
              key: 'completed',
              label: `已完成 (${completedCount})`,
            },
          ]}
        />
      </Card>

      {/* 筛选栏 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
          border: `1px solid ${token.colorBorder}`,
          background: token.colorBgContainer,
        }}
      >
        <Row gutter={[16, 16]}>
          <Col xs={24} sm={12} md={8} lg={6}>
            <Input
              placeholder="搜索配方名称、成分"
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
              placeholder="溶剂"
              value={solventFilter}
              onChange={setSolventFilter}
              allowClear
              style={{ width: '100%' }}
              showSearch
            >
              {getAllSolvents().map((s) => (
                <Select.Option key={s} value={s}>
                  {s}
                </Select.Option>
              ))}
            </Select>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Select
              placeholder="离子"
              value={ionFilter}
              onChange={setIonFilter}
              allowClear
              style={{ width: '100%' }}
              showSearch
            >
              {getAllIons().map((i) => (
                <Select.Option key={i} value={i}>
                  {i}
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
              <Select.Option value="name">配方名称</Select.Option>
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
              显示 {filteredElectrolytes.length} / {electrolytes.length} 个配方
            </Text>
          </Col>
        </Row>
      </Card>

      {/* 电解质列表 */}
      <Spin spinning={loading}>
        {filteredElectrolytes.length === 0 ? (
          <Card
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
            }}
          >
            <Empty
              image={<ExperimentOutlined style={{ fontSize: 64, color: '#d9d9d9' }} />}
              description={
                <Space direction="vertical" size={8}>
                  <Text type="secondary" style={{ fontSize: 16 }}>
                    {searchText ? '没有找到匹配的配方' : '还没有电解质配方'}
                  </Text>
                  {!searchText && (
                    <Text type="secondary">点击上方按钮创建第一个配方</Text>
                  )}
                </Space>
              }
              style={{ padding: '60px 0' }}
            >
              {!searchText && (
                <Button
                  type="primary"
                  icon={<PlusOutlined />}
                  onClick={() => handleOpenModal()}
                >
                  创建新配方
                </Button>
              )}
            </Empty>
          </Card>
        ) : viewMode === 'card' ? (
          <Row gutter={[16, 16]}>
            {filteredElectrolytes.map((electrolyte) => (
              <Col xs={24} sm={24} md={12} lg={8} key={electrolyte.id}>
                <div
                  style={{
                    position: 'relative',
                    cursor: selectMode ? 'pointer' : 'default',
                  }}
                  onClick={selectMode ? () => toggleSelect(electrolyte.id) : undefined}
                >
                  {selectMode && (
                    <div style={{
                      position: 'absolute',
                      top: 8,
                      left: 8,
                      zIndex: 10,
                      width: 24,
                      height: 24,
                      borderRadius: 4,
                      backgroundColor: selectedIds.includes(electrolyte.id) ? '#1677ff' : '#fff',
                      border: selectedIds.includes(electrolyte.id) ? '2px solid #1677ff' : '2px solid #d9d9d9',
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center',
                      color: '#fff',
                      fontWeight: 'bold',
                    }}>
                      {selectedIds.includes(electrolyte.id) && '✓'}
                    </div>
                  )}
                  <div style={{
                    opacity: selectMode && !selectedIds.includes(electrolyte.id) ? 0.7 : 1,
                    border: selectedIds.includes(electrolyte.id) ? '2px solid #1677ff' : '2px solid transparent',
                    borderRadius: 14,
                    transition: 'all 0.2s',
                  }}>
                    <ElectrolyteCard
                      electrolyte={electrolyte}
                      jobs={jobs}
                      onEdit={(e) => !selectMode && handleOpenModal(e, false)}
                      onCopy={(e) => !selectMode && handleOpenModal(e, true)}
                      onDelete={!selectMode ? handleDelete : () => {}}
                      onCreateJob={!selectMode ? handleCreateJob : () => {}}
                    />
                  </div>
                </div>
              </Col>
            ))}
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
              dataSource={filteredElectrolytes}
              rowKey="id"
              pagination={{
                pageSize: 20,
                showSizeChanger: true,
                showTotal: (total) => `共 ${total} 个配方`,
                pageSizeOptions: ['10', '20', '50', '100'],
              }}
              scroll={{ x: 1400 }}
              columns={[
                {
                  title: 'ID',
                  dataIndex: 'id',
                  key: 'id',
                  width: 80,
                  fixed: 'left',
                },
                {
                  title: '配方名称',
                  dataIndex: 'name',
                  key: 'name',
                  width: 250,
                  fixed: 'left',
                  ellipsis: true,
                  render: (name: string, record: ElectrolyteSystem) => (
                    <a onClick={() => navigate(`/workspace/electrolytes/${record.id}`)}>
                      {name}
                    </a>
                  ),
                },
                {
                  title: '项目',
                  dataIndex: 'project_id',
                  key: 'project_id',
                  width: 150,
                  render: (projectId: number) => {
                    const project = projects.find(p => p.id === projectId);
                    return project?.name || '-';
                  },
                },
                {
                  title: '阳离子',
                  key: 'cations',
                  width: 150,
                  render: (_: any, record: ElectrolyteSystem) => {
                    if (!record.cations || record.cations.length === 0) return '-';
                    return (
                      <Space size={4} wrap>
                        {record.cations.map((c, idx) => (
                          <Tag key={idx} color="blue">{c.name}</Tag>
                        ))}
                      </Space>
                    );
                  },
                },
                {
                  title: '阴离子',
                  key: 'anions',
                  width: 150,
                  render: (_: any, record: ElectrolyteSystem) => {
                    if (!record.anions || record.anions.length === 0) return '-';
                    return (
                      <Space size={4} wrap>
                        {record.anions.map((a, idx) => (
                          <Tag key={idx} color="orange">{a.name}</Tag>
                        ))}
                      </Space>
                    );
                  },
                },
                {
                  title: '溶剂',
                  key: 'solvents',
                  width: 150,
                  render: (_: any, record: ElectrolyteSystem) => {
                    if (!record.solvents || record.solvents.length === 0) return '-';
                    return (
                      <Space size={4} wrap>
                        {record.solvents.map((s, idx) => (
                          <Tag key={idx} color="green">{s.name}</Tag>
                        ))}
                      </Space>
                    );
                  },
                },
                {
                  title: '状态',
                  key: 'status',
                  width: 100,
                  render: (_: any, record: ElectrolyteSystem) => {
                    const category = getElectrolyteCategory(record);
                    const statusConfig: Record<string, { color: string; text: string }> = {
                      draft: { color: 'default', text: '草稿' },
                      running: { color: 'processing', text: '进行中' },
                      completed: { color: 'success', text: '已完成' },
                    };
                    const config = statusConfig[category] || statusConfig.draft;
                    return <Tag color={config.color}>{config.text}</Tag>;
                  },
                },
                {
                  title: '任务数',
                  key: 'job_count',
                  width: 100,
                  render: (_: any, record: ElectrolyteSystem) => {
                    const relatedJobs = jobs.filter(j => j.system_id === record.id);
                    return relatedJobs.length;
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
                  width: 250,
                  fixed: 'right',
                  render: (_: any, record: ElectrolyteSystem) => (
                    <Space size="small">
                      <Button
                        type="link"
                        size="small"
                        onClick={() => navigate(`/workspace/electrolytes/${record.id}`)}
                      >
                        详情
                      </Button>
                      <Button
                        type="link"
                        size="small"
                        onClick={() => handleOpenModal(record, false)}
                      >
                        编辑
                      </Button>
                      <Button
                        type="link"
                        size="small"
                        onClick={() => handleOpenModal(record, true)}
                      >
                        复制
                      </Button>
                      <Button
                        type="link"
                        size="small"
                        onClick={() => handleCreateJob(record)}
                      >
                        创建任务
                      </Button>
                    </Space>
                  ),
                },
              ]}
            />
          </Card>
        )}
      </Spin>

      {/* 创建/编辑/复制对话框 */}
      <Modal
        title={
          <Space>
            <ExperimentOutlined style={{ color: '#1677ff' }} />
            <span style={{ fontWeight: 600 }}>
              {copyingElectrolyte
                ? '复制电解质配方'
                : editingElectrolyte
                ? '编辑电解质配方'
                : '创建新电解质配方'}
            </span>
          </Space>
        }
        open={modalVisible}
        onOk={handleSubmit}
        onCancel={handleCloseModal}
        okText="保存"
        cancelText="取消"
        width={1200}
        centered
        styles={{
          body: { maxHeight: '70vh', overflowY: 'auto', padding: '24px' },
          header: { borderBottom: '1px solid #f0f0f0', paddingBottom: 16 },
        }}
      >
        <div style={{ marginTop: 24 }}>
          <ElectrolyteFormOptimized
            form={form}
            projects={projects}
            initialCations={selectedCations}
            initialAnions={selectedAnions}
            onIonsChange={(cations, anions) => {
              setSelectedCations(cations);
              setSelectedAnions(anions);
            }}
            onCreateProject={handleOpenProjectModal}
          />
        </div>
      </Modal>

      {/* 新建项目对话框 */}
      <Modal
        title="创建新项目"
        open={projectModalVisible}
        onOk={handleCreateProject}
        onCancel={handleCloseProjectModal}
        okText="创建"
        cancelText="取消"
        width={600}
      >
        <Form form={projectForm} layout="vertical" style={{ marginTop: 24 }}>
          <Form.Item
            name="name"
            label="项目名称"
            rules={[
              { required: true, message: '请输入项目名称' },
              { max: 100, message: '名称不能超过100个字符' }
            ]}
          >
            <Input placeholder="例如：锂离子电池电解液研究" />
          </Form.Item>
          <Form.Item
            name="description"
            label="项目描述"
            rules={[{ max: 500, message: '描述不能超过500个字符' }]}
          >
            <Input.TextArea
              rows={4}
              placeholder="简要描述项目的研究目标和内容"
            />
          </Form.Item>
        </Form>
      </Modal>

      {/* 批量导入对话框 - 步骤式流程 */}
      <Modal
        title={
          <Space>
            <UploadOutlined style={{ color: '#1677ff' }} />
            <span style={{ fontWeight: 600 }}>批量导入配方</span>
          </Space>
        }
        open={batchImportVisible}
        onCancel={resetBatchImport}
        footer={importStep === 2 ? [
          <Button key="close" type="primary" onClick={resetBatchImport}>
            完成
          </Button>
        ] : [
          <Button key="cancel" onClick={resetBatchImport}>
            取消
          </Button>,
          importStep > 0 && importStep < 2 && (
            <Button key="prev" onClick={() => setImportStep(importStep - 1)}>
              上一步
            </Button>
          ),
          importStep < 1 && (
            <Button
              key="next"
              type="primary"
              onClick={() => {
                if (!uploadFile) {
                  message.error('请先选择文件');
                  return;
                }
                setImportStep(1);
              }}
            >
              下一步
            </Button>
          ),
          importStep === 1 && (
            <Button key="import" type="primary" onClick={handleBatchImport} loading={importing}>
              开始导入
            </Button>
          ),
        ]}
        width={700}
        centered
        destroyOnClose
      >
        <div style={{ padding: '24px 0' }}>
          <Steps current={importStep} style={{ marginBottom: 32 }}>
            <Step title="选择文件" />
            <Step title="配置项目" />
            <Step title="导入结果" />
          </Steps>

          {/* 步骤1: 选择文件 */}
          {importStep === 0 && (
            <Space direction="vertical" size="large" style={{ width: '100%' }}>
              <Alert
                message="批量导入说明"
                description="下载模板后，在Excel中填写配方数据，多个离子/溶剂用分号(;)分隔。支持 .xlsx、.xls、.csv 格式。"
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
                  下载配方导入模板 (含示例)
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
                  onClick={() => document.getElementById('batch-file-input')?.click()}
                >
                  <input
                    id="batch-file-input"
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

          {/* 步骤2: 配置项目 */}
          {importStep === 1 && (
            <Space direction="vertical" size="large" style={{ width: '100%' }}>
              <Alert
                message={`已选择文件: ${uploadFile?.name}`}
                type="success"
                showIcon
              />
              <Card size="small" title="项目设置">
                <Radio.Group
                  value={projectMode}
                  onChange={(e) => setProjectMode(e.target.value)}
                  style={{ width: '100%' }}
                >
                  <Space direction="vertical" style={{ width: '100%' }}>
                    <Radio value="new">创建新项目</Radio>
                    {projectMode === 'new' && (
                      <div style={{ marginLeft: 24, marginTop: 8 }}>
                        <Space direction="vertical" style={{ width: '100%' }}>
                          <Input
                            placeholder="项目名称 *"
                            value={newProjectName}
                            onChange={(e) => setNewProjectName(e.target.value)}
                            style={{ maxWidth: 400 }}
                          />
                          <Input.TextArea
                            placeholder="项目描述（可选）"
                            value={newProjectDesc}
                            onChange={(e) => setNewProjectDesc(e.target.value)}
                            rows={2}
                            style={{ maxWidth: 400 }}
                          />
                        </Space>
                      </div>
                    )}
                    <Divider style={{ margin: '12px 0' }} />
                    <Radio value="existing">添加到现有项目</Radio>
                    {projectMode === 'existing' && (
                      <div style={{ marginLeft: 24, marginTop: 8 }}>
                        <Select
                          placeholder="选择项目"
                          value={selectedProjectId}
                          onChange={(value) => setSelectedProjectId(value)}
                          style={{ width: 300 }}
                          options={projects.map(p => ({
                            label: p.name,
                            value: p.id
                          }))}
                        />
                      </div>
                    )}
                  </Space>
                </Radio.Group>
              </Card>
            </Space>
          )}

          {/* 步骤3: 导入结果 */}
          {importStep === 2 && importResult && (
            <Space direction="vertical" size="large" style={{ width: '100%' }}>
              <Result
                status={importResult.failed_electrolytes === 0 && importResult.failed_md_jobs === 0 ? 'success' : 'warning'}
                title={importResult.failed_electrolytes === 0 && importResult.failed_md_jobs === 0 ? '导入成功' : '导入完成（部分失败）'}
                subTitle={`项目: ${importResult.project_name}`}
              />

              <Row gutter={16}>
                <Col span={8}>
                  <Card size="small">
                    <div style={{ textAlign: 'center' }}>
                      <div style={{ fontSize: 24, fontWeight: 600, color: '#52c41a' }}>
                        {importResult.success_electrolytes}
                      </div>
                      <div style={{ color: '#666' }}>配方成功</div>
                    </div>
                  </Card>
                </Col>
                <Col span={8}>
                  <Card size="small">
                    <div style={{ textAlign: 'center' }}>
                      <div style={{ fontSize: 24, fontWeight: 600, color: '#ff4d4f' }}>
                        {importResult.failed_electrolytes}
                      </div>
                      <div style={{ color: '#666' }}>配方失败</div>
                    </div>
                  </Card>
                </Col>
                <Col span={8}>
                  <Card size="small">
                    <div style={{ textAlign: 'center' }}>
                      <div style={{ fontSize: 24, fontWeight: 600, color: '#1677ff' }}>
                        {importResult.success_md_jobs}
                      </div>
                      <div style={{ color: '#666' }}>MD任务已创建</div>
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
              {importResult.electrolyte_results.length > 0 && (
                <Card
                  size="small"
                  title={<span style={{ color: '#52c41a' }}>成功导入的配方</span>}
                  extra={
                    <Space>
                      <Button
                        type="primary"
                        size="small"
                        onClick={async () => {
                          // 从导入结果中提取成功的配方ID
                          const successIds = importResult.electrolyte_results
                            ?.filter(r => r.status === 'success' && r.id)
                            .map(r => r.id) || [];
                          if (successIds.length > 0) {
                            setSelectedIds(successIds);

                            // 加载分区和配额信息
                            try {
                              const [partitionsData, quotaData] = await Promise.all([
                                fetch('/api/v1/slurm/partitions', {
                                  headers: { Authorization: `Bearer ${localStorage.getItem('access_token')}` }
                                }).then(res => res.json()),
                                fetch('/api/v1/jobs/quota/check', {
                                  headers: { Authorization: `Bearer ${localStorage.getItem('access_token')}` }
                                }).then(res => res.json())
                              ]);
                              setBatchMDPartitions(partitionsData);
                              setBatchMDQuota(quotaData);

                              // 设置默认分区
                              const defaultPartition = partitionsData.find((p: any) => p.state === 'up')?.name || 'cpu';
                              batchMDForm.setFieldsValue({ slurm_partition: defaultPartition });
                            } catch (error) {
                              console.error('加载分区/配额信息失败:', error);
                            }

                            setBatchMDModalVisible(true);
                          } else {
                            message.warning('没有可用的配方创建MD任务');
                          }
                        }}
                      >
                        批量创建MD任务
                      </Button>
                      <Button
                        size="small"
                        onClick={() => {
                          setBatchImportVisible(false);
                          setImportStep(0);
                          setImportResult(null);
                          loadData();
                        }}
                      >
                        完成并刷新列表
                      </Button>
                    </Space>
                  }
                >
                  <div style={{ maxHeight: 150, overflow: 'auto' }}>
                    {importResult.electrolyte_results.map((r, i) => (
                      <Tag key={i} color="green" style={{ margin: 4 }}>
                        {r.name} (ID: {r.id})
                      </Tag>
                    ))}
                  </div>
                </Card>
              )}
            </Space>
          )}
        </div>
      </Modal>

      {/* 批量创建MD任务Modal */}
      <Modal
        title="批量创建MD任务"
        open={batchMDModalVisible}
        onCancel={() => {
          setBatchMDModalVisible(false);
          batchMDForm.resetFields();
        }}
        width={600}
        footer={[
          <Button key="cancel" onClick={() => {
            setBatchMDModalVisible(false);
            batchMDForm.resetFields();
          }}>
            取消
          </Button>,
          <Button
            key="submit"
            type="primary"
            loading={creatingBatchMD}
            onClick={async () => {
              try {
                const values = await batchMDForm.validateFields();
                setCreatingBatchMD(true);

                // 构建QC选项
                const qcOptions = values.qc_enabled ? {
                  enabled: true,
                  functionals: values.qc_functionals || ['B3LYP'],
                  basis_sets: values.qc_basis_sets || ['6-31++g(d,p)'],
                  solvent_models: values.qc_solvent_models || ['pcm'],
                  solvents: values.qc_solvents || ['Water'],
                  molecules: [], // 将由后端从电解质配方中提取
                } : undefined;

                const result = await batchCreateMDJobs(selectedIds, {
                  job_name: values.job_name || undefined,
                  accuracy_level: values.accuracy_level || 'standard',
                  nsteps_npt: values.nsteps_npt || undefined,
                  nsteps_nvt: values.nsteps_nvt || undefined,
                  timestep: values.timestep || undefined,
                  temperature: values.temperature || undefined,
                  pressure: values.pressure || undefined,
                  freq_trj_npt: values.freq_trj_npt || undefined,
                  freq_trj_nvt: values.freq_trj_nvt || undefined,
                  thermo_freq: values.thermo_freq || undefined,
                  submit_to_cluster: false,
                  slurm_partition: values.slurm_partition || 'cpu',
                  slurm_nodes: values.slurm_nodes || 1,
                  slurm_ntasks: values.slurm_ntasks || 8,
                  slurm_cpus_per_task: values.slurm_cpus_per_task || 8,
                  slurm_time: values.slurm_time || 7200,
                  qc_options: qcOptions,
                });

                if (result.quota_exceeded) {
                  Modal.warning({
                    title: '配额不足',
                    content: result.message,
                  });
                  return;
                }

                if (result.success) {
                  // 保存创建成功的任务ID
                  const jobIds = result.success_jobs?.map((j: any) => j.job_id) || [];
                  setCreatedMDJobIds(jobIds);

                  setBatchMDModalVisible(false);
                  batchMDForm.resetFields();
                  setSelectedIds([]);
                  loadData();

                  // 询问是否批量提交
                  Modal.confirm({
                    title: '批量创建成功',
                    content: `成功创建 ${result.success_count} 个MD任务，是否立即批量提交到集群？`,
                    okText: '立即提交',
                    cancelText: '稍后提交',
                    onOk: async () => {
                      try {
                        // 批量提交任务
                        const submitPromises = jobIds.map((jobId: number) =>
                          fetch(`/api/v1/jobs/${jobId}/submit`, {
                            method: 'POST',
                            headers: {
                              'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
                              'Content-Type': 'application/json'
                            }
                          }).then(res => res.json())
                        );

                        const results = await Promise.allSettled(submitPromises);
                        const successCount = results.filter(r => r.status === 'fulfilled').length;
                        const failedCount = results.filter(r => r.status === 'rejected').length;

                        if (failedCount === 0) {
                          message.success(`成功提交 ${successCount} 个MD任务到集群`);
                        } else {
                          message.warning(`提交完成：成功 ${successCount} 个，失败 ${failedCount} 个`);
                        }
                        loadData();
                      } catch (error: any) {
                        message.error('批量提交失败: ' + (error.message || '未知错误'));
                      }
                    },
                    onCancel: () => {
                      message.info('任务已创建，您可以稍后在任务列表中提交');
                    }
                  });
                } else {
                  message.warning(`创建完成：成功 ${result.success_count} 个，失败 ${result.failed_count} 个`);
                  if (result.errors && result.errors.length > 0) {
                    Modal.error({
                      title: '部分任务创建失败',
                      content: (
                        <div>
                          {result.errors.map((err: any, i: number) => (
                            <div key={i}>配方ID {err.system_id}: {err.error}</div>
                          ))}
                        </div>
                      ),
                    });
                  }

                  // 即使部分失败，也询问是否提交成功的任务
                  if (result.success_count > 0) {
                    const jobIds = result.success_jobs?.map((j: any) => j.job_id) || [];
                    setCreatedMDJobIds(jobIds);

                    Modal.confirm({
                      title: '部分任务创建成功',
                      content: `成功创建 ${result.success_count} 个MD任务，是否立即批量提交到集群？`,
                      okText: '立即提交',
                      cancelText: '稍后提交',
                      onOk: async () => {
                        try {
                          const submitPromises = jobIds.map((jobId: number) =>
                            fetch(`/api/v1/jobs/${jobId}/submit`, {
                              method: 'POST',
                              headers: {
                                'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
                                'Content-Type': 'application/json'
                              }
                            }).then(res => res.json())
                          );

                          const results = await Promise.allSettled(submitPromises);
                          const successCount = results.filter(r => r.status === 'fulfilled').length;
                          const failedCount = results.filter(r => r.status === 'rejected').length;

                          if (failedCount === 0) {
                            message.success(`成功提交 ${successCount} 个MD任务到集群`);
                          } else {
                            message.warning(`提交完成：成功 ${successCount} 个，失败 ${failedCount} 个`);
                          }
                          loadData();
                        } catch (error: any) {
                          message.error('批量提交失败: ' + (error.message || '未知错误'));
                        }
                      }
                    });
                  }
                  loadData();
                }
              } catch (error: any) {
                message.error(error.response?.data?.detail || '批量创建失败');
              } finally {
                setCreatingBatchMD(false);
              }
            }}
          >
            创建 {selectedIds.length} 个MD任务
          </Button>,
        ]}
      >
        <Form
          form={batchMDForm}
          layout="vertical"
          initialValues={{
            accuracy_level: 'standard',
            slurm_partition: 'cpu',
            slurm_nodes: 1,
            slurm_ntasks: 8,
            slurm_cpus_per_task: 8,
            slurm_time: 7200,
          }}
        >
          <Alert
            message={`将为 ${selectedIds.length} 个配方创建MD任务`}
            type="info"
            showIcon
            style={{ marginBottom: 16 }}
          />

          {/* 配额检查 */}
          {batchMDQuota && !batchMDQuota.can_create && (
            <Alert
              message="配额不足"
              description={`当前已创建 ${batchMDQuota.current_count}/${batchMDQuota.limit} 个任务，剩余配额：${batchMDQuota.remaining} 个`}
              type="warning"
              showIcon
              style={{ marginBottom: 16 }}
            />
          )}
          {batchMDQuota && batchMDQuota.can_create && selectedIds.length > batchMDQuota.remaining && (
            <Alert
              message="配额警告"
              description={`批量创建 ${selectedIds.length} 个任务将超过剩余配额（${batchMDQuota.remaining} 个）`}
              type="warning"
              showIcon
              style={{ marginBottom: 16 }}
            />
          )}

          <Form.Item
            label="自定义名称后缀（可选）"
            name="job_name"
            tooltip="可选的自定义名称后缀，将添加到自动生成的任务名称后面"
            extra="留空：MD-日期-序号-配方名 | 填写：MD-日期-序号-配方名-自定义名称"
          >
            <Input placeholder="留空或输入自定义名称后缀（如：批量测试）" allowClear />
          </Form.Item>

          <Divider />

          <Form.Item
            label="精度等级"
            name="accuracy_level"
            tooltip="选择计算精度等级，影响步数和输出频率。选择自定义模式可手动配置所有参数"
          >
            <Select onChange={(value) => setBatchMDAccuracyLevel(value)}>
              <Select.Option value="fast">快速模式 (0.6 ns)</Select.Option>
              <Select.Option value="standard">标准模式 (15.5 ns)</Select.Option>
              <Select.Option value="accurate">精确模式 (55 ns)</Select.Option>
              <Select.Option value="custom">自定义模式</Select.Option>
            </Select>
          </Form.Item>

          {batchMDAccuracyLevel === 'custom' && (
            <>
              <Alert
                message="自定义模式"
                description="您可以手动配置所有模拟参数。留空的参数将使用标准模式的默认值。"
                type="warning"
                showIcon
                style={{ marginBottom: 16 }}
              />

              <Row gutter={16}>
                <Col span={12}>
                  <Form.Item
                    label="NPT步数"
                    name="nsteps_npt"
                    tooltip="NPT系综模拟步数"
                  >
                    <InputNumber
                      min={1000}
                      max={100000000}
                      step={100000}
                      style={{ width: '100%' }}
                      placeholder="默认: 15,000,000"
                    />
                  </Form.Item>
                </Col>
                <Col span={12}>
                  <Form.Item
                    label="NVT步数"
                    name="nsteps_nvt"
                    tooltip="NVT系综模拟步数"
                  >
                    <InputNumber
                      min={1000}
                      max={100000000}
                      step={100000}
                      style={{ width: '100%' }}
                      placeholder="默认: 500,000"
                    />
                  </Form.Item>
                </Col>
              </Row>

              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item
                    label="时间步长(fs)"
                    name="timestep"
                    tooltip="模拟时间步长"
                  >
                    <InputNumber
                      min={0.1}
                      max={10}
                      step={0.1}
                      style={{ width: '100%' }}
                      placeholder="默认: 1.0"
                    />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    label="温度(K)"
                    name="temperature"
                    tooltip="模拟温度"
                  >
                    <InputNumber
                      min={0}
                      max={1000}
                      step={10}
                      style={{ width: '100%' }}
                      placeholder="默认: 298.15"
                    />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    label="压力(atm)"
                    name="pressure"
                    tooltip="模拟压力"
                  >
                    <InputNumber
                      min={0}
                      max={1000}
                      step={0.1}
                      style={{ width: '100%' }}
                      placeholder="默认: 1.0"
                    />
                  </Form.Item>
                </Col>
              </Row>

              <Row gutter={16}>
                <Col span={12}>
                  <Form.Item
                    label="NPT轨迹输出频率"
                    name="freq_trj_npt"
                    tooltip="NPT阶段轨迹输出频率（步）"
                  >
                    <InputNumber
                      min={100}
                      max={10000000}
                      step={1000}
                      style={{ width: '100%' }}
                      placeholder="默认: 10,000"
                    />
                  </Form.Item>
                </Col>
                <Col span={12}>
                  <Form.Item
                    label="NVT轨迹输出频率"
                    name="freq_trj_nvt"
                    tooltip="NVT阶段轨迹输出频率（步）"
                  >
                    <InputNumber
                      min={100}
                      max={10000000}
                      step={1000}
                      style={{ width: '100%' }}
                      placeholder="默认: 1,000"
                    />
                  </Form.Item>
                </Col>
              </Row>

              <Form.Item
                label="热力学输出频率"
                name="thermo_freq"
                tooltip="热力学数据输出频率（步）"
              >
                <InputNumber
                  min={100}
                  max={10000000}
                  step={100}
                  style={{ width: '100%' }}
                  placeholder="默认: 1,000"
                />
              </Form.Item>

              <Divider />
            </>
          )}

          <Divider orientation="left">Slurm资源配置</Divider>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                label="队列/分区"
                name="slurm_partition"
                tooltip="显示管理员分配给您的可用队列，队列状态实时从集群获取"
                rules={[{ required: true, message: '请选择队列' }]}
              >
                <Select
                  placeholder={batchMDPartitions.length > 0 ? "选择队列" : "暂无可用队列"}
                  disabled={batchMDPartitions.length === 0}
                >
                  {batchMDPartitions.map((p: any) => (
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
                  ))}
                </Select>
              </Form.Item>
              {batchMDPartitions.length === 0 && (
                <Alert
                  message="暂无可用队列"
                  description="请联系管理员分配队列权限"
                  type="warning"
                  showIcon
                  style={{ marginBottom: 16 }}
                />
              )}
            </Col>
            <Col span={12}>
              <Form.Item label="节点数" name="slurm_nodes">
                <InputNumber min={1} max={10} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
          </Row>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item label="任务数" name="slurm_ntasks">
                <InputNumber min={1} max={128} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item label="每任务CPU数" name="slurm_cpus_per_task">
                <InputNumber min={1} max={64} style={{ width: '100%' }} />
              </Form.Item>
            </Col>
          </Row>

          <Form.Item label="最大运行时间(分钟)" name="slurm_time">
            <InputNumber min={10} max={43200} style={{ width: '100%' }} />
          </Form.Item>

          {/* 总MPI进程数显示 */}
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
                  style={{ marginBottom: 16 }}
                />
              );
            }}
          </Form.Item>

          <Divider orientation="left">量子化学计算 (可选)</Divider>

          <Form.Item
            name="qc_enabled"
            valuePropName="checked"
            initialValue={false}
          >
            <Checkbox>
              <Space>
                <ExperimentOutlined style={{ color: '#722ed1' }} />
                <Text strong>启用QC计算</Text>
              </Space>
            </Checkbox>
          </Form.Item>

          <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 16 }}>
            勾选后将对电解质中的分子进行量子化学计算，获取HOMO、LUMO、ESP等电子结构性质。
            计算将在MD任务创建后自动进行。
          </Text>

          <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
            prevValues.qc_enabled !== currentValues.qc_enabled
          }>
            {({ getFieldValue }) => {
              const qcEnabled = getFieldValue('qc_enabled');
              if (!qcEnabled) return null;

              return (
                <Card size="small" style={{ marginBottom: 16 }}>
                  <Row gutter={16}>
                    <Col span={12}>
                      <Form.Item
                        name="qc_functionals"
                        label="泛函"
                        initialValue={['B3LYP']}
                        tooltip="可选择多个泛函进行对比计算"
                      >
                        <Select mode="multiple" placeholder="选择泛函（可多选）">
                          <Select.Option value="B3LYP">B3LYP</Select.Option>
                          <Select.Option value="M062X">M06-2X</Select.Option>
                          <Select.Option value="wB97XD">ωB97X-D</Select.Option>
                          <Select.Option value="PBE0">PBE0</Select.Option>
                        </Select>
                      </Form.Item>
                    </Col>
                    <Col span={12}>
                      <Form.Item
                        name="qc_basis_sets"
                        label="基组"
                        initialValue={['6-31++g(d,p)']}
                        tooltip="可选择多个基组进行对比计算"
                      >
                        <Select mode="multiple" placeholder="选择基组（可多选）">
                          <Select.Option value="6-31g(d,p)">6-31G(d,p)</Select.Option>
                          <Select.Option value="6-31++g(d,p)">6-31++G(d,p)</Select.Option>
                          <Select.Option value="6-311g(d,p)">6-311G(d,p)</Select.Option>
                          <Select.Option value="Def2TZVP">Def2-TZVP</Select.Option>
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
                    style={{ marginBottom: 8 }}
                    message={
                      <Text style={{ fontSize: 12 }}>
                        <strong>隐式溶剂选择提示：</strong>选择介电常数(ε)接近您电解液的溶剂。
                        例如：EC体系选Water(ε≈78)，DMC/EMC体系选Chloroform(ε≈4.7)，高浓电解液选Acetone(ε≈20)。
                      </Text>
                    }
                  />

                  <Alert
                    message="QC任务数量"
                    description="将根据配方中的分子数量、泛函、基组和溶剂组合自动创建QC任务"
                    type="info"
                    showIcon
                    style={{ marginTop: 8 }}
                  />
                </Card>
              );
            }}
          </Form.Item>
        </Form>
      </Modal>
    </div>
  );
}

