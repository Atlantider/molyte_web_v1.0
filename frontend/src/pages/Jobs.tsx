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
} from 'antd';
import { PlusOutlined, ReloadOutlined, ThunderboltOutlined, RocketOutlined, ExperimentOutlined, FolderAddOutlined } from '@ant-design/icons';
import JobCard from '../components/JobCard';
import { getMDJobs, createMDJob, cancelMDJob, deleteMDJob, resubmitMDJob, updateMDJobConfig } from '../api/jobs';
import { getElectrolytes, createElectrolyteNew } from '../api/electrolytes';
import { getPartitions, getSlurmSuggestion, type PartitionInfo } from '../api/slurm';
import { getProjects } from '../api/projects';
import type { MDJob, MDJobCreate, ElectrolyteSystem, Project } from '../types';
import { JobStatus } from '../types';
import ElectrolyteFormOptimized from '../components/ElectrolyteFormOptimized';

const { Title, Text } = Typography;

export default function Jobs() {
  const location = useLocation();
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
        nsteps_npt: values.nsteps_npt,
        nsteps_nvt: values.nsteps_nvt,
        timestep: values.timestep,
        // Slurm 资源配置
        slurm_partition: values.slurm_partition,
        slurm_nodes: values.slurm_nodes,
        slurm_ntasks: values.slurm_ntasks,
        slurm_cpus_per_task: values.slurm_cpus_per_task,
        slurm_time: values.slurm_time,
      };
      // QC计算选项
      if (values.qc_enabled) {
        data.qc_options = {
          enabled: true,
          accuracy_level: values.qc_accuracy_level || 'standard',
          basis_set: values.qc_basis_set || '6-31++g(d,p)',
          functional: values.qc_functional || 'B3LYP',
          solvent_model: values.qc_solvent_model || 'pcm',
          solvent_name: values.qc_solvent_name || 'water',
          use_recommended_params: values.qc_use_recommended_params !== false,
        };
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
    if (activeTab === 'all') return jobs;
    if (activeTab === 'created') return jobs.filter((j) => j.status === JobStatus.CREATED);
    if (activeTab === 'running')
      return jobs.filter((j) =>
        j.status === JobStatus.QUEUED ||
        j.status === JobStatus.RUNNING ||
        j.status === JobStatus.POSTPROCESSING
      );
    if (activeTab === 'completed') return jobs.filter((j) => j.status === JobStatus.COMPLETED);
    if (activeTab === 'failed')
      return jobs.filter((j) => j.status === JobStatus.FAILED || j.status === JobStatus.CANCELLED);
    return jobs;
  };

  const filteredJobs = getFilteredJobs();

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
        ) : (
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

          <Divider orientation="left">模拟参数</Divider>

          <Form.Item
            name="nsteps_npt"
            label="NPT 步数"
            rules={[{ required: true, message: '请输入 NPT 步数' }]}
            tooltip="等压等温系综的模拟步数"
          >
            <InputNumber min={1000} max={10000000} step={1000} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="nsteps_nvt"
            label="NVT 步数"
            rules={[{ required: true, message: '请输入 NVT 步数' }]}
            tooltip="等容等温系综的模拟步数"
          >
            <InputNumber min={1000} max={10000000} step={1000} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="timestep"
            label="时间步长 (fs)"
            rules={[{ required: true, message: '请输入时间步长' }]}
            tooltip="分子动力学模拟的时间步长，单位为飞秒"
          >
            <InputNumber min={0.1} max={10} step={0.1} style={{ width: '100%' }} />
          </Form.Item>

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
            prevValues.qc_accuracy_level !== currentValues.qc_accuracy_level ||
            prevValues.qc_solvent_model !== currentValues.qc_solvent_model ||
            prevValues.qc_solvent_name !== currentValues.qc_solvent_name ||
            prevValues.qc_use_recommended_params !== currentValues.qc_use_recommended_params ||
            prevValues.electrolyte_id !== currentValues.electrolyte_id
          }>
            {({ getFieldValue }) => {
              const qcEnabled = getFieldValue('qc_enabled');
              if (!qcEnabled) return null;

              const accuracyLevel = getFieldValue('qc_accuracy_level') || 'standard';
              const solventModel = getFieldValue('qc_solvent_model') || 'pcm';
              const useRecommendedParams = getFieldValue('qc_use_recommended_params') !== false;
              const electrolyteId = getFieldValue('electrolyte_id');
              const selectedElectrolyte = electrolytes.find(e => e.id === electrolyteId);

              // 根据精度等级获取默认参数
              const getDefaultParams = (level: string) => {
                switch (level) {
                  case 'fast': return { basis_set: 'STO-3G', functional: 'HF' };
                  case 'standard': return { basis_set: '6-31G(d)', functional: 'B3LYP' };
                  case 'accurate': return { basis_set: '6-311++G(d,p)', functional: 'B3LYP' };
                  default: return { basis_set: getFieldValue('qc_basis_set') || '6-31++G(d,p)', functional: getFieldValue('qc_functional') || 'B3LYP' };
                }
              };

              // 根据分子类型获取推荐参数
              const getRecommendedParamsForMolecule = (molType: string, baseParams: { basis_set: string; functional: string }) => {
                if (!useRecommendedParams) {
                  return { ...baseParams, solvent_model: solventModel, reason: '' };
                }

                let params = { ...baseParams, solvent_model: solventModel, reason: '' };

                if (molType === 'anion') {
                  // 阴离子需要弥散函数
                  if (!params.basis_set.includes('+')) {
                    params.basis_set = accuracyLevel === 'accurate' ? '6-311++G(d,p)' : '6-31++G(d,p)';
                  }
                  params.reason = '阴离子：使用弥散函数(++)描述扩展电子密度';
                  if (params.solvent_model === 'gas') {
                    params.solvent_model = 'pcm';
                    params.reason += '，使用PCM溶剂模型稳定电子结构';
                  }
                } else if (molType === 'cation') {
                  params.reason = '阳离子：使用极化函数描述紧凑电子结构';
                  if (params.solvent_model === 'gas') {
                    params.solvent_model = 'pcm';
                    params.reason += '，使用PCM溶剂模型';
                  }
                } else {
                  params.reason = '中性分子：使用标准参数';
                }

                return params;
              };

              const baseParams = getDefaultParams(accuracyLevel);

              // 提取所有分子
              const allMolecules: { name: string; smiles: string; type: string; params: any }[] = [];
              if (selectedElectrolyte) {
                selectedElectrolyte.solvents?.forEach((s: any) => {
                  const params = getRecommendedParamsForMolecule('solvent', baseParams);
                  allMolecules.push({ name: s.name, smiles: s.smiles, type: 'solvent', params });
                });
                selectedElectrolyte.cations?.forEach((c: any) => {
                  const params = getRecommendedParamsForMolecule('cation', baseParams);
                  allMolecules.push({ name: c.name, smiles: c.smiles, type: 'cation', params });
                });
                selectedElectrolyte.anions?.forEach((a: any) => {
                  const params = getRecommendedParamsForMolecule('anion', baseParams);
                  allMolecules.push({ name: a.name, smiles: a.smiles, type: 'anion', params });
                });
              }

              return (
                <Card size="small" style={{ marginBottom: 16 }}>
                  {/* 精度等级选择 */}
                  <Form.Item
                    name="qc_accuracy_level"
                    label="精度等级"
                    initialValue="standard"
                    style={{ marginBottom: 12 }}
                  >
                    <Select>
                      <Select.Option value="fast">
                        <Space>
                          <Tag color="green">快速</Tag>
                          <Text type="secondary">HF/STO-3G (~5分钟)</Text>
                        </Space>
                      </Select.Option>
                      <Select.Option value="standard">
                        <Space>
                          <Tag color="blue">标准</Tag>
                          <Text type="secondary">B3LYP/6-31G(d) (~30分钟)</Text>
                        </Space>
                      </Select.Option>
                      <Select.Option value="accurate">
                        <Space>
                          <Tag color="orange">精确</Tag>
                          <Text type="secondary">B3LYP/6-311++G(d,p) (~2小时)</Text>
                        </Space>
                      </Select.Option>
                      <Select.Option value="custom">
                        <Space>
                          <Tag color="purple">自定义</Tag>
                          <Text type="secondary">自定义泛函和基组</Text>
                        </Space>
                      </Select.Option>
                    </Select>
                  </Form.Item>

                  {/* 自定义泛函和基组 - 仅在自定义模式下显示 */}
                  {accuracyLevel === 'custom' && (
                    <Row gutter={16} style={{ marginBottom: 12 }}>
                      <Col span={12}>
                        <Form.Item
                          name="qc_basis_set"
                          label="基组"
                          initialValue="6-31++g(d,p)"
                          style={{ marginBottom: 0 }}
                          tooltip="阴离子建议使用带弥散函数(++)的基组"
                        >
                          <Select>
                            <Select.OptGroup label="标准基组">
                              <Select.Option value="6-31g(d)">6-31G(d) - 几何优化推荐</Select.Option>
                              <Select.Option value="6-31g(d,p)">6-31G(d,p) - 含氢体系</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="带弥散函数 (阴离子推荐)">
                              <Select.Option value="6-31++g(d,p)">6-31++G(d,p) - 阴离子/弱相互作用</Select.Option>
                              <Select.Option value="6-311++g(d,p)">6-311++G(d,p) - 高精度</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="Def2系列">
                              <Select.Option value="Def2SVP">Def2-SVP - 平衡精度效率</Select.Option>
                              <Select.Option value="Def2TZVP">Def2-TZVP - 高精度</Select.Option>
                            </Select.OptGroup>
                          </Select>
                        </Form.Item>
                      </Col>
                      <Col span={12}>
                        <Form.Item
                          name="qc_functional"
                          label="泛函"
                          initialValue="B3LYP"
                          style={{ marginBottom: 0 }}
                        >
                          <Select>
                            <Select.OptGroup label="杂化泛函 (推荐)">
                              <Select.Option value="B3LYP">B3LYP - 常用</Select.Option>
                              <Select.Option value="PBE0">PBE0 - 无经验参数</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="Minnesota泛函">
                              <Select.Option value="M062X">M06-2X - 主族化学</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="色散校正">
                              <Select.Option value="wB97XD">ωB97X-D - 弱相互作用</Select.Option>
                            </Select.OptGroup>
                          </Select>
                        </Form.Item>
                      </Col>
                    </Row>
                  )}

                  {/* 溶剂环境设置 */}
                  <Row gutter={16} style={{ marginBottom: 12 }}>
                    <Col span={12}>
                      <Form.Item
                        name="qc_solvent_model"
                        label="溶剂环境"
                        initialValue="pcm"
                        style={{ marginBottom: 0 }}
                        tooltip={
                          <div>
                            <p><strong>气相 (Gas)</strong>: 真空环境，无溶剂效应</p>
                            <p><strong>PCM</strong>: 极化连续介质模型，使用介电常数描述溶剂</p>
                            <p><strong>SMD</strong>: 溶剂密度模型，更精确但计算量更大</p>
                            <p>离子在气相中可能不稳定，建议使用PCM/SMD</p>
                          </div>
                        }
                      >
                        <Select>
                          <Select.Option value="gas">气相 (Gas Phase) - 无溶剂效应</Select.Option>
                          <Select.Option value="pcm">PCM - 极化连续介质模型 (推荐)</Select.Option>
                          <Select.Option value="smd">SMD - 溶剂密度模型 (更精确)</Select.Option>
                        </Select>
                      </Form.Item>
                    </Col>
                    <Col span={12}>
                      {(solventModel === 'pcm' || solventModel === 'smd') && (
                        <Form.Item
                          name="qc_solvent_name"
                          label="隐式溶剂"
                          initialValue="Water"
                          style={{ marginBottom: 0 }}
                          tooltip={
                            <div>
                              <p><strong>选择原则</strong>：选择介电常数(ε)接近您电解液的溶剂</p>
                              <p>• 水系电解液 → Water (ε=78.4)</p>
                              <p>• 高浓电解液 → Acetone (ε=20.5)</p>
                              <p>• DMC/EMC体系 → Chloroform (ε≈4.7)</p>
                            </div>
                          }
                        >
                          <Select showSearch optionFilterProp="children">
                            <Select.OptGroup label="📌 水系电解液 (ε>50)">
                              <Select.Option value="Water">水 (Water) ε=78.4</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="📌 高介电常数 (ε=40-90)">
                              <Select.Option value="DiMethylSulfoxide">DMSO ε=46.8 (离子液体参考)</Select.Option>
                              <Select.Option value="1,2-EthaneDiol">乙二醇 ε=40.2</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="📌 中等介电常数 (ε=15-40)">
                              <Select.Option value="Acetonitrile">乙腈 ε=35.7</Select.Option>
                              <Select.Option value="Methanol">甲醇 ε=32.6</Select.Option>
                              <Select.Option value="Ethanol">乙醇 ε=24.9</Select.Option>
                              <Select.Option value="Acetone">丙酮 ε=20.5 (高浓电解液)</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="📌 低介电常数 (ε<15) - DMC/EMC/DEC体系">
                              <Select.Option value="DiChloroEthane">二氯乙烷 ε=10.1</Select.Option>
                              <Select.Option value="Dichloromethane">二氯甲烷 ε=8.9</Select.Option>
                              <Select.Option value="TetraHydroFuran">四氢呋喃 (THF) ε=7.4</Select.Option>
                              <Select.Option value="Chloroform">氯仿 ε=4.7 (线性碳酸酯参考)</Select.Option>
                              <Select.Option value="DiethylEther">乙醚 ε=4.2</Select.Option>
                            </Select.OptGroup>
                            <Select.OptGroup label="自定义">
                              <Select.Option value="custom">自定义溶剂参数...</Select.Option>
                            </Select.OptGroup>
                          </Select>
                        </Form.Item>
                      )}
                    </Col>
                  </Row>

                  {/* 自定义溶剂参数输入 */}
                  {(solventModel === 'pcm' || solventModel === 'smd') && getFieldValue('qc_solvent_name') === 'custom' && (
                    <Row gutter={16} style={{ marginBottom: 12 }}>
                      <Col span={8}>
                        <Form.Item
                          name="qc_custom_solvent_eps"
                          label="介电常数 ε"
                          rules={[{ required: true, message: '请输入介电常数' }]}
                          style={{ marginBottom: 0 }}
                          tooltip="溶剂的静态介电常数，如水为78.4"
                        >
                          <InputNumber min={1} max={200} step={0.1} style={{ width: '100%' }} placeholder="例如: 78.4" />
                        </Form.Item>
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          name="qc_custom_solvent_epsinf"
                          label="光学介电常数 ε∞"
                          style={{ marginBottom: 0 }}
                          tooltip="溶剂的光学介电常数（可选），默认为1.0"
                        >
                          <InputNumber min={1} max={10} step={0.01} style={{ width: '100%' }} placeholder="例如: 1.78" />
                        </Form.Item>
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          name="qc_custom_solvent_name"
                          label="溶剂名称"
                          style={{ marginBottom: 0 }}
                          tooltip="自定义溶剂的名称（用于记录）"
                        >
                          <Input placeholder="例如: 高浓LiTFSI" />
                        </Form.Item>
                      </Col>
                    </Row>
                  )}

                  {/* 智能参数推荐 */}
                  <Form.Item
                    name="qc_use_recommended_params"
                    valuePropName="checked"
                    initialValue={true}
                    style={{ marginBottom: 12 }}
                  >
                    <Checkbox>
                      <Space>
                        <Text>智能参数推荐</Text>
                        <Text type="secondary" style={{ fontSize: 12 }}>
                          (自动为阴离子添加弥散函数，为离子选择合适的溶剂模型)
                        </Text>
                      </Space>
                    </Checkbox>
                  </Form.Item>

                  {/* 分子参数详情 */}
                  {allMolecules.length > 0 && (
                    <Alert
                      type="info"
                      showIcon
                      message={`将对 ${allMolecules.length} 个分子进行QC计算`}
                      description={
                        <div style={{ marginTop: 8 }}>
                          {allMolecules.map((mol, idx) => (
                            <div key={idx} style={{
                              padding: '8px 12px',
                              marginBottom: 8,
                              background: '#fafafa',
                              borderRadius: 6,
                              border: '1px solid #f0f0f0'
                            }}>
                              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 4 }}>
                                <Space>
                                  <Text strong>{mol.name}</Text>
                                  <Tag color={mol.type === 'solvent' ? 'blue' : mol.type === 'cation' ? 'green' : 'orange'}>
                                    {mol.type === 'solvent' ? '溶剂' : mol.type === 'cation' ? '阳离子' : '阴离子'}
                                  </Tag>
                                </Space>
                              </div>
                              <div style={{ fontSize: 12, color: '#666' }}>
                                <Space split={<span style={{ color: '#d9d9d9' }}>|</span>}>
                                  <span>泛函: <Text code>{mol.params.functional}</Text></span>
                                  <span>基组: <Text code>{mol.params.basis_set}</Text></span>
                                  <span>溶剂: <Text code>{mol.params.solvent_model === 'gas' ? '气相' : mol.params.solvent_model.toUpperCase()}</Text></span>
                                </Space>
                              </div>
                              {mol.params.reason && (
                                <div style={{ fontSize: 11, color: '#999', marginTop: 4 }}>
                                  💡 {mol.params.reason}
                                </div>
                              )}
                            </div>
                          ))}
                        </div>
                      }
                    />
                  )}
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


