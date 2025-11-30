/**
 * 任务配置和创建页面
 */
import { useState, useEffect } from 'react';
import { useNavigate, useParams, useLocation } from 'react-router-dom';
import {
  Card,
  Form,
  Input,
  InputNumber,
  Button,
  Space,
  Typography,
  Descriptions,
  message,
  Alert,
  Divider,
  Row,
  Col,
  Select,
  Checkbox,
} from 'antd';
import { ArrowLeftOutlined, ThunderboltOutlined, SaveOutlined, WarningOutlined, ExperimentOutlined } from '@ant-design/icons';
import type { ElectrolyteSystem, MDJobCreate } from '../types';
import { getElectrolyte } from '../api/electrolytes';
import { createMDJob, checkJobQuota } from '../api/jobs';
import { getPartitions, type PartitionInfo } from '../api/slurm';
import AccuracyLevelSelector from '../components/AccuracyLevelSelector';

const { Title, Text } = Typography;

export default function JobCreate() {
  const navigate = useNavigate();
  const location = useLocation();
  const { systemId } = useParams<{ systemId: string }>();
  const [form] = Form.useForm();
  const [loading, setLoading] = useState(false);
  const [electrolyte, setElectrolyte] = useState<ElectrolyteSystem | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [quota, setQuota] = useState<{
    can_create: boolean;
    current_count: number;
    limit: number;
    remaining: number;
  } | null>(null);
  const [accuracyDefaults, setAccuracyDefaults] = useState<any>(null);
  const [selectedAccuracyLevel, setSelectedAccuracyLevel] = useState<string>('standard');
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);

  // 从 location.state 或 API 加载电解质系统信息
  useEffect(() => {
    const loadElectrolyte = async () => {
      setLoading(true);
      try {
        // 优先使用 location.state 传递的数据
        if (location.state?.electrolyte) {
          setElectrolyte(location.state.electrolyte);
        } else if (systemId) {
          const data = await getElectrolyte(Number(systemId));
          setElectrolyte(data);
        }
      } catch (error: any) {
        message.error('加载配方信息失败: ' + (error.response?.data?.detail || error.message));
        navigate('/workspace/electrolytes');
      } finally {
        setLoading(false);
      }
    };

    loadElectrolyte();
  }, [systemId, location.state]);

  // 加载配额信息
  useEffect(() => {
    const loadQuota = async () => {
      try {
        const quotaData = await checkJobQuota();
        setQuota(quotaData);
      } catch (error: any) {
        console.error('加载配额信息失败:', error);
      }
    };

    loadQuota();
  }, []);

  // 加载集群分区信息
  useEffect(() => {
    const loadPartitions = async () => {
      try {
        const data = await getPartitions();
        setPartitions(data);
      } catch (error) {
        console.error('Failed to load partitions:', error);
        // 如果加载失败，使用默认值
        setPartitions([]);
      }
    };

    loadPartitions();
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

  // 获取默认分区
  const getDefaultPartition = () => {
    if (partitions.length > 0) {
      const upPartition = partitions.find(p => p.state === 'up');
      return upPartition?.name || partitions[0].name;
    }
    return 'cpu';
  };

  // 设置默认值
  useEffect(() => {
    if (electrolyte) {
      form.setFieldsValue({
        job_name: '', // 留空，让后端自动生成
        accuracy_level: 'standard', // 默认标准模式
        // 所有参数留空，使用精度等级的默认值
        nsteps_npt: undefined,
        nsteps_nvt: undefined,
        timestep: undefined,
        temperature: undefined,
        pressure: undefined,
        freq_trj_npt: undefined,
        freq_trj_nvt: undefined,
        thermo_freq: undefined,
        // Slurm 资源配置默认值
        slurm_partition: getDefaultPartition(),
        slurm_nodes: 1,
        slurm_ntasks: 8,
        slurm_cpus_per_task: 8,
        slurm_time: 7200,
      });
    }
  }, [electrolyte, form, partitions]);

  // 获取当前精度等级的默认值
  const getCurrentDefaults = () => {
    return accuracyDefaults?.[selectedAccuracyLevel] || {};
  };

  // 检查是否为自定义模式
  const isCustomMode = () => {
    return selectedAccuracyLevel === 'custom';
  };

  // 将步数转换为时间，自动选择合适的单位（ns、ps或fs）
  const stepsToTime = (steps: number | null | undefined, timestep: number = 1.0): string => {
    if (!steps) return '-';

    const timeFs = steps * timestep; // 总时间（飞秒）

    // 根据时间大小选择合适的单位
    if (timeFs >= 1_000_000) {
      // >= 1 ns，使用 ns
      const timeNs = timeFs / 1_000_000;
      return `${timeNs.toFixed(1)} ns`;
    } else if (timeFs >= 1_000) {
      // >= 1 ps，使用 ps
      const timePs = timeFs / 1_000;
      return `${timePs.toFixed(1)} ps`;
    } else {
      // < 1 ps，使用 fs
      return `${timeFs.toFixed(0)} fs`;
    }
  };

  // 格式化步数显示（步数 + 时间）
  const formatStepsWithTime = (steps: number | null | undefined, timestep: number = 1.0): string => {
    if (!steps) return '-';
    return `${steps.toLocaleString()} 步 (${stepsToTime(steps, timestep)})`;
  };

  // 生成字段配置的辅助函数
  const getFieldConfig = (fieldName: string, label: string, unit: string = '') => {
    const customMode = isCustomMode();
    const defaultValue = getCurrentDefaults()[fieldName];

    return {
      label: `${label}${customMode ? '（可修改）' : '（可选）'}`,
      rules: customMode ? [{ required: true, message: `请输入${label}` }] : [],
      tooltip: customMode
        ? '已填充参考值，您可以根据需要修改'
        : '留空使用精度等级的默认值',
      extra: customMode && defaultValue
        ? `参考值（标准模式）: ${typeof defaultValue === 'number' ? defaultValue.toLocaleString() : defaultValue}${unit}`
        : (!customMode && accuracyDefaults && defaultValue
            ? `默认值: ${typeof defaultValue === 'number' ? defaultValue.toLocaleString() : defaultValue}${unit}`
            : undefined),
      placeholder: customMode
        ? `参考值: ${defaultValue || ''}${unit}`
        : (defaultValue ? `默认: ${typeof defaultValue === 'number' ? defaultValue.toLocaleString() : defaultValue}${unit}` : '请选择精度等级'),
    };
  };

  // 构建任务数据
  const buildJobData = (values: any, submitToCluster: boolean = false): MDJobCreate => {
    const jobData: MDJobCreate = {
      system_id: electrolyte!.id,
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
      submit_to_cluster: submitToCluster,
      // Slurm 资源配置
      slurm_partition: values.slurm_partition || 'cpu',
      slurm_nodes: values.slurm_nodes || 1,
      slurm_ntasks: values.slurm_ntasks || 8,
      slurm_cpus_per_task: values.slurm_cpus_per_task || 8,
      slurm_time: values.slurm_time || 7200,
    };

    // QC计算选项
    if (values.qc_enabled) {
      jobData.qc_options = {
        enabled: true,
        functional: values.qc_functionals?.[0] || 'B3LYP',
        basis_set: values.qc_basis_sets?.[0] || '6-31++g(d,p)',
        solvent_model: values.qc_solvent_models?.[0] || 'pcm',
        solvent_name: values.qc_solvents?.[0] || 'Water',
        molecules: [], // 将由后端从电解质配方中提取
      } as any;
    }

    return jobData;
  };

  // 通用提交处理
  const handleJobSubmit = async (submitToCluster: boolean) => {
    try {
      const values = await form.validateFields();
      setSubmitting(true);

      const jobData = buildJobData(values, submitToCluster);
      const job = await createMDJob(jobData);

      const successMsg = submitToCluster ? '计算任务已提交到集群！' : '计算任务已保存！';
      const targetPath = submitToCluster ? `/workspace/jobs/${job.id}/detail` : `/workspace/jobs/${job.id}/submit`;

      message.success(successMsg);
      navigate(targetPath);
    } catch (error: any) {
      if (error.errorFields) {
        message.error('请检查表单填写');
      } else {
        const action = submitToCluster ? '提交' : '保存';
        message.error(`${action}任务失败: ` + (error.response?.data?.detail || error.message));
      }
    } finally {
      setSubmitting(false);
    }
  };

  // 保存任务（不提交到集群）
  const handleSave = () => handleJobSubmit(false);

  // 提交任务到集群
  const handleSubmit = () => handleJobSubmit(true);

  if (loading || !electrolyte) {
    return (
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        height: 'calc(100vh - 64px)',
        background: '#f5f7fb',
      }}>
        加载中...
      </div>
    );
  }

  return (
    <div style={{ padding: 24, background: '#f5f7fb', minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面头部 */}
      <div style={{ marginBottom: 24 }}>
        <Space style={{ marginBottom: 16 }}>
          <Button
            icon={<ArrowLeftOutlined />}
            onClick={() => navigate('/workspace/electrolytes')}
            style={{ borderRadius: 8 }}
          >
            返回
          </Button>
        </Space>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <ThunderboltOutlined style={{ marginRight: 12, color: '#1677ff' }} />
          创建计算任务
        </Title>
        <Text type="secondary">配置分子动力学模拟参数并提交计算任务</Text>
      </div>

      {/* 配额提示 */}
      {quota && (
        <Alert
          message={
            <Space>
              <span>今日任务配额</span>
              <Text strong>
                {quota.current_count} / {quota.limit}
              </Text>
              <Text type="secondary">（剩余 {quota.remaining} 个）</Text>
            </Space>
          }
          type={quota.can_create ? 'info' : 'warning'}
          showIcon
          icon={!quota.can_create ? <WarningOutlined /> : undefined}
          description={
            !quota.can_create
              ? '您今日的任务创建配额已用完，请明天再试或联系管理员。'
              : '请根据需要配置分子动力学模拟的参数，然后提交任务到计算集群'
          }
          style={{ marginBottom: 24, borderRadius: 8 }}
        />
      )}

      {/* 配方信息 */}
      <Card
        title="电解质配方信息"
        style={{
          marginBottom: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}
      >
        <Descriptions bordered column={2}>
          <Descriptions.Item label="配方名称">{electrolyte.name}</Descriptions.Item>
          <Descriptions.Item label="配方 ID">#{electrolyte.id}</Descriptions.Item>
          <Descriptions.Item label="温度">{electrolyte.temperature} K</Descriptions.Item>
          <Descriptions.Item label="压力">{electrolyte.pressure} atm</Descriptions.Item>
          <Descriptions.Item label="盒子大小">
            {electrolyte.box_size ? Number(electrolyte.box_size).toFixed(1) : '-'} Å
          </Descriptions.Item>
          <Descriptions.Item label="力场">{electrolyte.force_field || 'OPLS-AA'}</Descriptions.Item>
        </Descriptions>
      </Card>

      {/* 计算参数配置表单 */}
      <Card
        title="计算参数配置"
        style={{
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}
      >
        <Form form={form} layout="vertical">
          <Form.Item
            label="自定义名称（可选）"
            name="job_name"
            tooltip="可选的自定义名称后缀，将添加到自动生成的任务名称后面"
            extra="留空：MD-20251119-0001-配方名 | 填写：MD-20251119-0001-配方名-自定义名称"
          >
            <Input placeholder="留空或输入自定义名称后缀（如：高温测试）" allowClear />
          </Form.Item>

          <Divider />

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

          {!isCustomMode() && (
            <Alert
              message="提示"
              description="选择精度等级后，系统会自动设置下方的模拟参数。如果需要自定义参数，可以在下方手动修改。"
              type="info"
              showIcon
              style={{ marginBottom: 24 }}
            />
          )}

          {isCustomMode() && (
            <Alert
              message="自定义模式"
              description="您选择了自定义模式。下方已自动填充标准模式的参数作为参考，您可以根据需要修改任何参数。"
              type="warning"
              showIcon
              style={{ marginBottom: 24 }}
            />
          )}

          {!isCustomMode() && accuracyDefaults && (
            <>
              <Divider>模拟参数（由精度等级自动配置）</Divider>
              <Descriptions bordered column={1} size="small" style={{ marginBottom: 24 }}>
                <Descriptions.Item label="电荷产生方式">
                  {getCurrentDefaults().charge_method === 'resp' ? (
                    <span style={{ color: '#f5222d' }}>🎯 RESP（高精度电荷）</span>
                  ) : (
                    <span style={{ color: '#52c41a' }}>🚀 LigParGen（快速生成）</span>
                  )}
                </Descriptions.Item>
                <Descriptions.Item label="NPT 系综模拟">
                  {formatStepsWithTime(getCurrentDefaults().nsteps_npt, getCurrentDefaults().timestep)}
                </Descriptions.Item>
                <Descriptions.Item label="NVT 系综模拟">
                  {formatStepsWithTime(getCurrentDefaults().nsteps_nvt, getCurrentDefaults().timestep)}
                </Descriptions.Item>
                <Descriptions.Item label="时间步长">
                  {getCurrentDefaults().timestep || '-'} fs
                </Descriptions.Item>
                <Descriptions.Item label="温度">
                  {getCurrentDefaults().temperature || '-'} K
                </Descriptions.Item>
                <Descriptions.Item label="压力">
                  {getCurrentDefaults().pressure || '-'} atm
                </Descriptions.Item>
                <Descriptions.Item label="NPT 轨迹输出频率">
                  {formatStepsWithTime(getCurrentDefaults().freq_trj_npt, getCurrentDefaults().timestep)}
                </Descriptions.Item>
                <Descriptions.Item label="NVT 轨迹输出频率">
                  {formatStepsWithTime(getCurrentDefaults().freq_trj_nvt, getCurrentDefaults().timestep)}
                </Descriptions.Item>
                <Descriptions.Item label="热力学输出频率">
                  {formatStepsWithTime(getCurrentDefaults().thermo_freq, getCurrentDefaults().timestep)}
                </Descriptions.Item>
              </Descriptions>
            </>
          )}

          {isCustomMode() && (
            <>
              <Divider>电荷计算方法</Divider>

              <Form.Item
                name="charge_method"
                label="电荷产生方式"
                rules={[{ required: true, message: '请选择电荷产生方式' }]}
                tooltip="LigParGen: 快速生成，适合一般研究；RESP: 高精度电荷，适合论文发表"
                initialValue="ligpargen"
              >
                <Select>
                  <Select.Option value="ligpargen">
                    <Space>
                      <span>🚀 LigParGen</span>
                      <span style={{ color: '#888', fontSize: 12 }}>快速生成，适合一般研究</span>
                    </Space>
                  </Select.Option>
                  <Select.Option value="resp">
                    <Space>
                      <span>🎯 RESP</span>
                      <span style={{ color: '#888', fontSize: 12 }}>高精度电荷，适合论文发表（计算时间更长）</span>
                    </Space>
                  </Select.Option>
                </Select>
              </Form.Item>

              <Divider>模拟参数设置（必填）</Divider>

              <Form.Item
                name="nsteps_npt"
                label="NPT 系综模拟步数（可修改）"
                rules={[{ required: true, message: '请输入 NPT 步数' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${formatStepsWithTime(getCurrentDefaults().nsteps_npt, getCurrentDefaults().timestep)}`}
              >
                <InputNumber
                  min={1000}
                  max={100000000}
                  step={100000}
                  style={{ width: '100%' }}
                  formatter={(value) => `${value}`.replace(/\B(?=(\d{3})+(?!\d))/g, ',')}
                  addonAfter="步"
                />
              </Form.Item>

              <Form.Item
                name="nsteps_nvt"
                label="NVT 系综模拟步数（可修改）"
                rules={[{ required: true, message: '请输入 NVT 步数' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${formatStepsWithTime(getCurrentDefaults().nsteps_nvt, getCurrentDefaults().timestep)}`}
              >
                <InputNumber
                  min={1000}
                  max={100000000}
                  step={100000}
                  style={{ width: '100%' }}
                  formatter={(value) => `${value}`.replace(/\B(?=(\d{3})+(?!\d))/g, ',')}
                  addonAfter="步"
                />
              </Form.Item>

              <Form.Item
                name="timestep"
                label="时间步长（可修改）"
                rules={[{ required: true, message: '请输入时间步长' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${getCurrentDefaults().timestep} fs`}
              >
                <InputNumber
                  min={0.1}
                  max={10}
                  step={0.1}
                  style={{ width: '100%' }}
                  addonAfter="fs"
                />
              </Form.Item>

              <Form.Item
                name="temperature"
                label="温度（可修改）"
                rules={[{ required: true, message: '请输入温度' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${getCurrentDefaults().temperature} K`}
              >
                <InputNumber
                  min={200}
                  max={500}
                  step={1}
                  style={{ width: '100%' }}
                  addonAfter="K"
                />
              </Form.Item>

              <Form.Item
                name="pressure"
                label="压力（可修改）"
                rules={[{ required: true, message: '请输入压力' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${getCurrentDefaults().pressure} atm`}
              >
                <InputNumber
                  min={0.1}
                  max={100}
                  step={0.1}
                  style={{ width: '100%' }}
                  addonAfter="atm"
                />
              </Form.Item>

              <Divider>输出频率设置（必填）</Divider>

              <Form.Item
                name="freq_trj_npt"
                label="NPT 轨迹输出频率（可修改）"
                rules={[{ required: true, message: '请输入 NPT 轨迹输出频率' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${formatStepsWithTime(getCurrentDefaults().freq_trj_npt, getCurrentDefaults().timestep)}`}
              >
                <InputNumber
                  min={100}
                  max={10000000}
                  step={100}
                  style={{ width: '100%' }}
                  formatter={(value) => `${value}`.replace(/\B(?=(\d{3})+(?!\d))/g, ',')}
                  addonAfter="步"
                />
              </Form.Item>

              <Form.Item
                name="freq_trj_nvt"
                label="NVT 轨迹输出频率（可修改）"
                rules={[{ required: true, message: '请输入 NVT 轨迹输出频率' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${formatStepsWithTime(getCurrentDefaults().freq_trj_nvt, getCurrentDefaults().timestep)}`}
              >
                <InputNumber
                  min={100}
                  max={10000000}
                  step={100}
                  style={{ width: '100%' }}
                  formatter={(value) => `${value}`.replace(/\B(?=(\d{3})+(?!\d))/g, ',')}
                  addonAfter="步"
                />
              </Form.Item>

              <Form.Item
                name="thermo_freq"
                label="热力学输出频率（可修改）"
                rules={[{ required: true, message: '请输入热力学输出频率' }]}
                tooltip="已填充参考值，您可以根据需要修改"
                extra={`参考值（标准模式）: ${formatStepsWithTime(getCurrentDefaults().thermo_freq, getCurrentDefaults().timestep)}`}
              >
                <InputNumber
                  min={100}
                  max={10000000}
                  step={100}
                  style={{ width: '100%' }}
                  formatter={(value) => `${value}`.replace(/\B(?=(\d{3})+(?!\d))/g, ',')}
                  addonAfter="步"
                />
              </Form.Item>
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
              style={{ marginBottom: 12 }}
            >
              <Checkbox>
                <Space>
                  <ExperimentOutlined style={{ color: '#722ed1' }} />
                  <Text strong>启用QC计算</Text>
                </Space>
              </Checkbox>
            </Form.Item>

            <Text type="secondary" style={{ fontSize: 12 }}>
              勾选后将对电解质中的分子进行量子化学计算，获取HOMO、LUMO、ESP等电子结构性质。
              计算将在MD任务提交后自动进行。
            </Text>
          </Card>

          <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
            prevValues.qc_enabled !== currentValues.qc_enabled
          }>
            {({ getFieldValue }) => {
              const qcEnabled = getFieldValue('qc_enabled');
              if (!qcEnabled) return null;

              // 收集将要计算的分子列表
              const moleculesToCalc: Array<{name: string, smiles: string, type: string, charge: number}> = [];
              if (electrolyte) {
                // 溶剂分子
                electrolyte.solvents?.forEach(sol => {
                  if (sol.smiles && !moleculesToCalc.find(m => m.smiles === sol.smiles)) {
                    moleculesToCalc.push({ name: sol.name, smiles: sol.smiles, type: 'solvent', charge: 0 });
                  }
                });
                // 阳离子
                electrolyte.cations?.forEach(cat => {
                  if (cat.smiles && !moleculesToCalc.find(m => m.smiles === cat.smiles)) {
                    moleculesToCalc.push({ name: cat.name, smiles: cat.smiles, type: 'cation', charge: 1 });
                  }
                });
                // 阴离子
                electrolyte.anions?.forEach(an => {
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
                        label="隐式溶剂模型"
                        initialValue={['pcm']}
                        style={{ marginBottom: 8 }}
                        tooltip="可选择多个模型进行对比计算"
                      >
                        <Select mode="multiple" placeholder="选择溶剂模型（可多选）">
                          <Select.Option value="gas">气相 (Gas Phase)</Select.Option>
                          <Select.Option value="pcm">PCM (极化连续介质模型)</Select.Option>
                          <Select.Option value="smd">SMD (溶剂密度模型)</Select.Option>
                        </Select>
                      </Form.Item>
                    </Col>
                    <Col span={12}>
                      <Form.Item
                        name="qc_solvents"
                        label="溶剂"
                        initialValue={['Water']}
                        style={{ marginBottom: 8 }}
                        tooltip="可选择多个溶剂进行对比计算"
                      >
                        <Select mode="multiple" placeholder="选择溶剂（可多选）" showSearch>
                          <Select.OptGroup label="常用溶剂">
                            <Select.Option value="Water">水 (Water, ε=78.4)</Select.Option>
                            <Select.Option value="Acetonitrile">乙腈 (Acetonitrile, ε=35.7)</Select.Option>
                            <Select.Option value="DMSO">二甲亚砜 (DMSO, ε=46.8)</Select.Option>
                            <Select.Option value="Methanol">甲醇 (Methanol, ε=32.6)</Select.Option>
                            <Select.Option value="Ethanol">乙醇 (Ethanol, ε=24.9)</Select.Option>
                          </Select.OptGroup>
                          <Select.OptGroup label="电解液常用">
                            <Select.Option value="Acetone">丙酮 (Acetone, ε=20.5)</Select.Option>
                            <Select.Option value="EC">碳酸乙烯酯 (EC, ε=89.8)</Select.Option>
                            <Select.Option value="DMC">碳酸二甲酯 (DMC, ε=3.1)</Select.Option>
                            <Select.Option value="EMC">碳酸甲乙酯 (EMC, ε=2.9)</Select.Option>
                            <Select.Option value="DEC">碳酸二乙酯 (DEC, ε=2.8)</Select.Option>
                            <Select.Option value="PC">碳酸丙烯酯 (PC, ε=64.9)</Select.Option>
                          </Select.OptGroup>
                          <Select.OptGroup label="其他">
                            <Select.Option value="Chloroform">氯仿 (Chloroform, ε=4.7)</Select.Option>
                            <Select.Option value="Dichloromethane">二氯甲烷 (DCM, ε=8.9)</Select.Option>
                            <Select.Option value="THF">四氢呋喃 (THF, ε=7.4)</Select.Option>
                            <Select.Option value="Toluene">甲苯 (Toluene, ε=2.4)</Select.Option>
                          </Select.OptGroup>
                        </Select>
                      </Form.Item>
                    </Col>
                  </Row>

                  {/* 显示将要计算的分子列表和任务数量 */}
                  <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
                    prevValues.qc_functionals !== currentValues.qc_functionals ||
                    prevValues.qc_basis_sets !== currentValues.qc_basis_sets ||
                    prevValues.qc_solvent_models !== currentValues.qc_solvent_models ||
                    prevValues.qc_solvents !== currentValues.qc_solvents
                  }>
                    {({ getFieldValue }) => {
                      const functionals = getFieldValue('qc_functionals') || ['B3LYP'];
                      const basisSets = getFieldValue('qc_basis_sets') || ['6-31++g(d,p)'];
                      const solventModels = getFieldValue('qc_solvent_models') || ['pcm'];
                      const solvents = getFieldValue('qc_solvents') || ['Water'];

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

          <Divider>计算资源配置</Divider>

          <Row gutter={16}>
            <Col span={12}>
              <Form.Item
                name="slurm_partition"
                label="队列/分区"
                tooltip="显示管理员分配给您的可用队列，队列状态实时从集群获取"
                rules={[{ required: true, message: '请选择队列' }]}
              >
                <Select
                  placeholder={partitions.length > 0 ? "选择队列" : "暂无可用队列"}
                  disabled={partitions.length === 0}
                >
                  {partitions.map(p => (
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
              {partitions.length === 0 && (
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
                  style={{ marginBottom: 24 }}
                />
              );
            }}
          </Form.Item>

          <Form.Item>
            <Space>
              <Button
                icon={<SaveOutlined />}
                onClick={handleSave}
                loading={submitting}
                size="large"
                disabled={quota ? !quota.can_create : false}
              >
                保存计算任务
              </Button>
              <Button
                type="primary"
                icon={<ThunderboltOutlined />}
                onClick={handleSubmit}
                loading={submitting}
                size="large"
                disabled={quota ? !quota.can_create : false}
              >
                提交到集群
              </Button>
              <Button onClick={() => navigate('/workspace/electrolytes')} size="large">
                取消
              </Button>
            </Space>
          </Form.Item>
        </Form>
      </Card>
    </div>
  );
}


