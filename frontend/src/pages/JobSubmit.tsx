/**
 * 任务提交页面 - 检查配置并提交到集群
 */
import { useState, useEffect } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
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
  Modal,
  Tag,
  Select,
  Switch,
  Row,
  Col,
} from 'antd';
import { ArrowLeftOutlined, ThunderboltOutlined, EditOutlined, WalletOutlined, ExperimentOutlined, BulbOutlined, CheckCircleOutlined, SyncOutlined } from '@ant-design/icons';
import type { ElectrolyteSystem, MDJob, MDJobCreate } from '../types';
import { getElectrolyte } from '../api/electrolytes';
import { getMDJob, updateMDJobConfig, submitJobToCluster, createMDJob } from '../api/jobs';
import { checkCanSubmit } from '../api/billing';
import { checkDuplicateCalculations, MoleculeCheckResult, DuplicateCheckResponse } from '../api/qc';

const { Title, Text } = Typography;

export default function JobSubmit() {
  const navigate = useNavigate();
  const { jobId } = useParams<{ jobId: string }>();
  const [form] = Form.useForm();
  const [loading, setLoading] = useState(false);
  const [submitting, setSubmitting] = useState(false);
  const [job, setJob] = useState<MDJob | null>(null);
  const [electrolyte, setElectrolyte] = useState<ElectrolyteSystem | null>(null);
  const [editMode, setEditMode] = useState(false);
  const [isSubmittedJob, setIsSubmittedJob] = useState(false); // 是否是已提交的任务
  // QC 分子参数编辑状态
  const [moleculeParams, setMoleculeParams] = useState<Record<string, {
    functional: string;
    basis_set: string;
    solvent_model: string;
  }>>({});
  const [editingMolecule, setEditingMolecule] = useState<string | null>(null); // 当前正在编辑的分子key
  // QC 全局参数编辑状态
  const [editingGlobalQC, setEditingGlobalQC] = useState(false);
  const [globalQCParams, setGlobalQCParams] = useState<{
    accuracy_level: string;
    solvent_model: string;
    solvent_name: string;
    use_recommended_params: boolean;
  } | null>(null);
  // 重复计算检查状态
  const [duplicateCheckResult, setDuplicateCheckResult] = useState<DuplicateCheckResponse | null>(null);
  const [checkingDuplicates, setCheckingDuplicates] = useState(false);

  // 加载任务和配方信息
  useEffect(() => {
    const loadData = async () => {
      setLoading(true);
      try {
        const jobData = await getMDJob(Number(jobId));
        setJob(jobData);

        // 检查任务是否已提交
        // CREATED 和 CANCELLED 状态可以编辑和提交
        setIsSubmittedJob(jobData.status !== 'CREATED' && jobData.status !== 'CANCELLED');

        const electrolyteData = await getElectrolyte(jobData.system_id);
        setElectrolyte(electrolyteData);

        // 设置表单值
        if (jobData.config) {
          form.setFieldsValue(jobData.config);
        }
      } catch (error: any) {
        message.error('加载任务信息失败: ' + (error.response?.data?.detail || error.message));
        navigate('/workspace/jobs');
      } finally {
        setLoading(false);
      }
    };

    if (jobId) {
      loadData();
    }
  }, [jobId]);

  // 检查QC重复计算
  const checkQCDuplicates = async () => {
    if (!job?.config?.qc_enabled || !electrolyte) return null;

    const molecules: any[] = [];
    const config = job.config;
    const solventModel = config.qc_solvent_model || 'pcm';
    const solventName = config.qc_solvent_name || 'water';
    const functional = config.qc_functional || 'B3LYP';
    const basisSet = config.qc_basis_set || '6-31G(d)';

    // 收集所有需要计算的分子
    electrolyte.solvents?.forEach((s: any) => {
      const customParams = moleculeParams[`solvent_${electrolyte.solvents?.indexOf(s)}`];
      molecules.push({
        smiles: s.smiles,
        molecule_name: s.name,
        functional: customParams?.functional || functional,
        basis_set: customParams?.basis_set || basisSet,
        solvent_model: customParams?.solvent_model || solventModel,
        solvent_name: solventModel !== 'gas' ? solventName : undefined,
        charge: 0,
        spin_multiplicity: 1,
      });
    });

    electrolyte.cations?.forEach((c: any) => {
      const customParams = moleculeParams[`cation_${electrolyte.cations?.indexOf(c)}`];
      molecules.push({
        smiles: c.smiles,
        molecule_name: c.name,
        functional: customParams?.functional || functional,
        basis_set: customParams?.basis_set || basisSet,
        solvent_model: customParams?.solvent_model || solventModel,
        solvent_name: solventModel !== 'gas' ? solventName : undefined,
        charge: 1,
        spin_multiplicity: 1,
      });
    });

    electrolyte.anions?.forEach((a: any) => {
      const customParams = moleculeParams[`anion_${electrolyte.anions?.indexOf(a)}`];
      molecules.push({
        smiles: a.smiles,
        molecule_name: a.name,
        functional: customParams?.functional || functional,
        basis_set: customParams?.basis_set || basisSet,
        solvent_model: customParams?.solvent_model || solventModel,
        solvent_name: solventModel !== 'gas' ? solventName : undefined,
        charge: -1,
        spin_multiplicity: 1,
      });
    });

    if (molecules.length === 0) return null;

    try {
      setCheckingDuplicates(true);
      const result = await checkDuplicateCalculations(molecules);
      setDuplicateCheckResult(result);
      return result;
    } catch (error) {
      console.error('检查重复计算失败:', error);
      return null;
    } finally {
      setCheckingDuplicates(false);
    }
  };

  // 提交到集群
  const handleSubmit = async () => {
    // 先检查余额
    try {
      const canSubmitResult = await checkCanSubmit();
      if (!canSubmitResult.can_submit) {
        Modal.confirm({
          title: '余额不足',
          icon: <WalletOutlined style={{ color: '#faad14' }} />,
          content: (
            <div>
              <p>{canSubmitResult.reason}</p>
              <p>是否前往充值？</p>
            </div>
          ),
          okText: '前往充值',
          cancelText: '取消',
          onOk: () => navigate('/workspace/recharge'),
        });
        return;
      }
    } catch (error) {
      console.error('检查余额失败:', error);
      // 如果检查失败，继续提交（后端会再次验证）
    }

    // 如果启用了QC计算，检查重复
    let duplicateInfo = duplicateCheckResult;
    if (job?.config?.qc_enabled && !duplicateInfo) {
      duplicateInfo = await checkQCDuplicates();
    }

    // 构建确认消息
    let confirmContent = '确定要将此任务提交到 Slurm 集群执行吗？提交后将开始计算。';
    if (duplicateInfo && duplicateInfo.existing_count > 0) {
      confirmContent = (
        <div>
          <p>确定要将此任务提交到 Slurm 集群执行吗？</p>
          <Alert
            type="success"
            showIcon
            icon={<CheckCircleOutlined />}
            style={{ marginTop: 12 }}
            message={
              <span>
                检测到 <strong>{duplicateInfo.existing_count}</strong> 个分子已有计算结果，
                将直接复用，无需重复计算！
              </span>
            }
            description={
              duplicateInfo.new_count > 0 ? (
                <span>另外 {duplicateInfo.new_count} 个分子将执行新计算。</span>
              ) : (
                <span>所有QC计算都将复用已有结果，节省计算时间和资源！</span>
              )
            }
          />
        </div>
      ) as any;
    }

    Modal.confirm({
      title: '确认提交任务到集群',
      content: confirmContent,
      okText: '确定提交',
      cancelText: '取消',
      onOk: async () => {
        try {
          setSubmitting(true);
          await submitJobToCluster(Number(jobId));
          message.success('任务已提交到集群！');
          navigate(`/workspace/jobs/${jobId}/detail`);
        } catch (error: any) {
          const detail = error.response?.data?.detail || error.message;
          // 检查是否是余额不足错误
          if (error.response?.status === 402) {
            Modal.confirm({
              title: '余额不足',
              icon: <WalletOutlined style={{ color: '#faad14' }} />,
              content: (
                <div>
                  <p>{detail}</p>
                  <p>是否前往充值？</p>
                </div>
              ),
              okText: '前往充值',
              cancelText: '取消',
              onOk: () => navigate('/workspace/recharge'),
            });
          } else {
            message.error('提交失败: ' + detail);
          }
        } finally {
          setSubmitting(false);
        }
      },
    });
  };

  // 生成复制任务的名称（添加 -copy 后缀）
  const generateCopyName = (originalName: string) => {
    // 检查原名称是否已经有 copy 标记
    const copyMatch = originalName.match(/-copy(-(\d+))?$/);
    if (copyMatch) {
      const copyNumber = copyMatch[2] ? parseInt(copyMatch[2]) + 1 : 2;
      return originalName.replace(/-copy(-\d+)?$/, `-copy-${copyNumber}`);
    }
    return `${originalName}-copy`;
  };

  // 保存修改（已提交的任务会创建新任务）
  const handleSaveChanges = async () => {
    try {
      const values = await form.validateFields();
      setSubmitting(true);

      // 构建QC配置（如果启用）
      const qcOptions = job?.config?.qc_enabled ? {
        enabled: true,
        accuracy_level: values.qc_accuracy_level || job.config?.qc_accuracy_level || 'standard',
        basis_set: values.qc_basis_set || job.config?.qc_basis_set || '6-31++g(d,p)',
        functional: values.qc_functional || job.config?.qc_functional || 'B3LYP',
        solvent_model: values.qc_solvent_model || job.config?.qc_solvent_model || 'pcm',
        solvent_name: values.qc_solvent_name || job.config?.qc_solvent_name || 'water',
        use_recommended_params: values.qc_use_recommended_params !== undefined ? values.qc_use_recommended_params : (job.config?.qc_use_recommended_params !== false),
      } : undefined;

      if (isSubmittedJob) {
        // 已提交的任务：创建新任务
        const originalName = values.job_name || job!.config?.job_name || '';
        const newJobName = generateCopyName(originalName);
        const newJobData: MDJobCreate = {
          system_id: job!.system_id,
          job_name: newJobName,
          nsteps_npt: values.nsteps_npt,
          nsteps_nvt: values.nsteps_nvt,
          timestep: values.timestep,
          temperature: values.temperature,
          pressure: values.pressure,
          freq_trj_npt: values.freq_trj_npt,
          freq_trj_nvt: values.freq_trj_nvt,
          thermo_freq: values.thermo_freq,
          submit_to_cluster: false, // 新任务默认不提交
          qc_options: qcOptions,
        };

        const newJob = await createMDJob(newJobData);
        message.success(`已创建新任务：${newJobName}`);
        navigate(`/workspace/jobs/${newJob.id}/submit`);
      } else {
        // 未提交的任务：直接更新（包含QC配置）
        const updateData = {
          ...values,
          // 如果启用QC，确保QC配置字段被包含
          ...(job?.config?.qc_enabled && {
            qc_accuracy_level: values.qc_accuracy_level || job.config?.qc_accuracy_level,
            qc_basis_set: values.qc_basis_set || job.config?.qc_basis_set,
            qc_functional: values.qc_functional || job.config?.qc_functional,
            qc_solvent_model: values.qc_solvent_model || job.config?.qc_solvent_model,
            qc_solvent_name: values.qc_solvent_name || job.config?.qc_solvent_name,
            qc_use_recommended_params: values.qc_use_recommended_params !== undefined ? values.qc_use_recommended_params : job.config?.qc_use_recommended_params,
          }),
        };
        const updatedJob = await updateMDJobConfig(Number(jobId), updateData);
        setJob(updatedJob);
        message.success('配置已更新！');
        setEditMode(false);
      }
    } catch (error: any) {
      if (error.errorFields) {
        message.error('请检查表单填写');
      } else {
        message.error('保存失败: ' + (error.response?.data?.detail || error.message));
      }
    } finally {
      setSubmitting(false);
    }
  };

  // 处理编辑按钮点击
  const handleEditClick = () => {
    if (isSubmittedJob) {
      Modal.confirm({
        title: '创建新任务',
        content: '该任务已提交，无法直接修改。是否要基于当前配置创建一个新任务？',
        okText: '创建新任务',
        cancelText: '取消',
        onOk: () => {
          setEditMode(true);
        },
      });
    } else {
      setEditMode(true);
    }
  };

  if (loading || !job || !electrolyte) {
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
            onClick={() => navigate('/workspace/jobs')}
            style={{ borderRadius: 8 }}
          >
            返回任务列表
          </Button>
        </Space>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <ThunderboltOutlined style={{ marginRight: 12, color: '#1677ff' }} />
          {isSubmittedJob ? '查看任务配置' : '提交计算任务'}
        </Title>
        <Text type="secondary">
          {isSubmittedJob ? '查看已提交任务的配置参数' : '检查并确认计算参数后提交到集群'}
        </Text>
      </div>

      {isSubmittedJob ? (
        <Alert
          message="任务已提交"
          description="该任务已提交到集群，无法直接修改。如需修改参数，可以基于当前配置创建新任务。"
          type="warning"
          showIcon
          style={{ marginBottom: 24, borderRadius: 8 }}
        />
      ) : (
        <Alert
          message="检查计算参数"
          description="请仔细检查以下计算参数，确认无误后提交到集群执行"
          type="info"
          showIcon
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

      {/* MD计算参数 */}
      <Card
        title={
          <Space>
            <span>MD计算参数配置</span>
            {job?.config?.qc_enabled && <Tag color="blue">MD</Tag>}
          </Space>
        }
        extra={
          !editMode && (
            <Button icon={<EditOutlined />} onClick={handleEditClick} style={{ borderRadius: 8 }}>
              {isSubmittedJob ? '基于此配置创建新任务' : '修改参数'}
            </Button>
          )
        }
        style={{
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}
      >
        {editMode ? (
          <Form form={form} layout="vertical">
            <Form.Item
              label="任务名称"
              name="job_name"
              rules={[{ required: true, message: '请输入任务名称' }]}
            >
              <Input placeholder="输入任务名称" />
            </Form.Item>

            <Divider>模拟步数设置</Divider>

            <Form.Item
              label="NPT 系综步数"
              name="nsteps_npt"
              rules={[{ required: true, message: '请输入 NPT 步数' }]}
            >
              <InputNumber min={0} style={{ width: '100%' }} />
            </Form.Item>

            <Form.Item
              label="NVT 系综步数"
              name="nsteps_nvt"
              rules={[{ required: true, message: '请输入 NVT 步数' }]}
            >
              <InputNumber min={0} style={{ width: '100%' }} />
            </Form.Item>

            <Form.Item
              label="时间步长 (fs)"
              name="timestep"
              rules={[{ required: true, message: '请输入时间步长' }]}
            >
              <InputNumber min={0} step={0.1} style={{ width: '100%' }} />
            </Form.Item>

            <Divider>热力学参数</Divider>

            <Form.Item
              label="温度 (K)"
              name="temperature"
              rules={[{ required: true, message: '请输入温度' }]}
            >
              <InputNumber min={0} style={{ width: '100%' }} />
            </Form.Item>

            <Form.Item
              label="压力 (atm)"
              name="pressure"
              rules={[{ required: true, message: '请输入压力' }]}
            >
              <InputNumber min={0} step={0.1} style={{ width: '100%' }} />
            </Form.Item>

            <Divider>输出频率设置</Divider>

            <Form.Item
              label="NPT 轨迹输出频率"
              name="freq_trj_npt"
              rules={[{ required: true, message: '请输入输出频率' }]}
            >
              <InputNumber min={0} style={{ width: '100%' }} />
            </Form.Item>

            <Form.Item
              label="NVT 轨迹输出频率"
              name="freq_trj_nvt"
              rules={[{ required: true, message: '请输入输出频率' }]}
            >
              <InputNumber min={0} style={{ width: '100%' }} />
            </Form.Item>

            <Form.Item
              label="热力学输出频率"
              name="thermo_freq"
              rules={[{ required: true, message: '请输入输出频率' }]}
            >
              <InputNumber min={0} style={{ width: '100%' }} />
            </Form.Item>

            <Form.Item>
              <Space>
                <Button type="primary" onClick={handleSaveChanges} loading={submitting}>
                  {isSubmittedJob ? '创建新任务' : '保存修改'}
                </Button>
                <Button onClick={() => setEditMode(false)}>取消</Button>
              </Space>
            </Form.Item>
          </Form>
        ) : (
          <Descriptions bordered column={2}>
            <Descriptions.Item label="任务名称" span={2}>
              {job.config?.job_name || '-'}
            </Descriptions.Item>
            <Descriptions.Item label="NPT 步数">
              {job.config?.nsteps_npt?.toLocaleString() || '-'}
            </Descriptions.Item>
            <Descriptions.Item label="NVT 步数">
              {job.config?.nsteps_nvt?.toLocaleString() || '-'}
            </Descriptions.Item>
            <Descriptions.Item label="时间步长">
              {job.config?.timestep || '-'} fs
            </Descriptions.Item>
            <Descriptions.Item label="温度">{job.config?.temperature || '-'} K</Descriptions.Item>
            <Descriptions.Item label="压力">{job.config?.pressure || '-'} atm</Descriptions.Item>
            <Descriptions.Item label="NPT 轨迹输出频率">
              {job.config?.freq_trj_npt?.toLocaleString() || '-'}
            </Descriptions.Item>
            <Descriptions.Item label="NVT 轨迹输出频率">
              {job.config?.freq_trj_nvt?.toLocaleString() || '-'}
            </Descriptions.Item>
            <Descriptions.Item label="热力学输出频率">
              {job.config?.thermo_freq?.toLocaleString() || '-'}
            </Descriptions.Item>
          </Descriptions>
        )}
      </Card>

      {/* QC计算配置卡片 */}
      {job.config?.qc_enabled && (() => {
        // 根据精度等级获取默认参数
        const getDefaultParams = (level: string) => {
          switch (level) {
            case 'fast': return { basis_set: 'STO-3G', functional: 'HF' };
            case 'standard': return { basis_set: '6-31G(d)', functional: 'B3LYP' };
            case 'accurate': return { basis_set: '6-311++G(d,p)', functional: 'B3LYP' };
            default: return { basis_set: job.config?.qc_basis_set || '6-31++G(d,p)', functional: job.config?.qc_functional || 'B3LYP' };
          }
        };

        // 根据分子类型获取推荐参数
        const getRecommendedParamsForMolecule = (molType: string, baseParams: { basis_set: string; functional: string }) => {
          const useRecommendedParams = job.config?.qc_use_recommended_params !== false;
          const solventModel = job.config?.qc_solvent_model || 'pcm';

          if (!useRecommendedParams) {
            return { ...baseParams, solvent_model: solventModel, reason: '' };
          }

          let params = { ...baseParams, solvent_model: solventModel, reason: '' };

          if (molType === 'anion') {
            if (!params.basis_set.includes('+')) {
              const accuracyLevel = job.config?.qc_accuracy_level || 'standard';
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

        const baseParams = getDefaultParams(job.config?.qc_accuracy_level || 'standard');

        // 提取所有分子并获取参数
        const allMolecules: { key: string; name: string; smiles: string; type: string; params: any }[] = [];
        if (electrolyte) {
          electrolyte.solvents?.forEach((s: any, idx: number) => {
            const key = `solvent_${idx}`;
            const defaultParams = getRecommendedParamsForMolecule('solvent', baseParams);
            const customParams = moleculeParams[key];
            allMolecules.push({
              key,
              name: s.name,
              smiles: s.smiles,
              type: 'solvent',
              params: customParams ? { ...customParams, reason: '' } : defaultParams
            });
          });
          electrolyte.cations?.forEach((c: any, idx: number) => {
            const key = `cation_${idx}`;
            const defaultParams = getRecommendedParamsForMolecule('cation', baseParams);
            const customParams = moleculeParams[key];
            allMolecules.push({
              key,
              name: c.name,
              smiles: c.smiles,
              type: 'cation',
              params: customParams ? { ...customParams, reason: '' } : defaultParams
            });
          });
          electrolyte.anions?.forEach((a: any, idx: number) => {
            const key = `anion_${idx}`;
            const defaultParams = getRecommendedParamsForMolecule('anion', baseParams);
            const customParams = moleculeParams[key];
            allMolecules.push({
              key,
              name: a.name,
              smiles: a.smiles,
              type: 'anion',
              params: customParams ? { ...customParams, reason: '' } : defaultParams
            });
          });
        }

        // 更新单个分子参数
        const handleMoleculeParamChange = (key: string, field: string, value: string) => {
          const mol = allMolecules.find(m => m.key === key);
          if (mol) {
            setMoleculeParams(prev => ({
              ...prev,
              [key]: {
                ...(prev[key] || mol.params),
                [field]: value
              }
            }));
          }
        };

        // 保存分子编辑
        const handleSaveMolecule = (key: string) => {
          setEditingMolecule(null);
          message.success('分子参数已更新');
        };

        // 获取当前有效的全局参数（优先使用编辑中的值）
        const currentAccuracyLevel = globalQCParams?.accuracy_level ?? job.config?.qc_accuracy_level ?? 'standard';
        const currentSolventModel = globalQCParams?.solvent_model ?? job.config?.qc_solvent_model ?? 'pcm';
        const currentSolventName = globalQCParams?.solvent_name ?? job.config?.qc_solvent_name ?? 'water';
        const currentUseRecommended = globalQCParams?.use_recommended_params ?? job.config?.qc_use_recommended_params ?? true;

        // 初始化全局参数编辑
        const handleEditGlobalQC = () => {
          setGlobalQCParams({
            accuracy_level: job.config?.qc_accuracy_level || 'standard',
            solvent_model: job.config?.qc_solvent_model || 'pcm',
            solvent_name: job.config?.qc_solvent_name || 'water',
            use_recommended_params: job.config?.qc_use_recommended_params !== false,
          });
          setEditingGlobalQC(true);
        };

        // 保存全局参数
        const handleSaveGlobalQC = async () => {
          if (!globalQCParams) return;
          try {
            const updatedConfig = {
              ...job.config,
              qc_accuracy_level: globalQCParams.accuracy_level,
              qc_solvent_model: globalQCParams.solvent_model,
              qc_solvent_name: globalQCParams.solvent_name,
              qc_use_recommended_params: globalQCParams.use_recommended_params,
            };
            await updateMDJobConfig(Number(jobId), updatedConfig);
            setJob({ ...job, config: updatedConfig });
            setEditingGlobalQC(false);
            message.success('QC全局参数已更新');
          } catch (error: any) {
            message.error('更新失败: ' + (error.response?.data?.detail || error.message));
          }
        };

        return (
          <Card
            title={
              <Space>
                <ExperimentOutlined style={{ color: '#722ed1' }} />
                <span>QC计算参数配置</span>
                <Tag color="purple">QC</Tag>
              </Space>
            }
            extra={!isSubmittedJob && !editingGlobalQC && (
              <Button size="small" type="link" icon={<EditOutlined />} onClick={handleEditGlobalQC}>
                编辑全局参数
              </Button>
            )}
            style={{ marginTop: 24, borderLeft: '4px solid #722ed1' }}
          >
            {/* 全局配置展示/编辑 */}
            {editingGlobalQC && globalQCParams ? (
              <div style={{ background: '#f6f0ff', padding: 16, borderRadius: 8, marginBottom: 16 }}>
                <Row gutter={[16, 12]}>
                  <Col span={6}>
                    <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>精度等级</div>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={globalQCParams.accuracy_level}
                      onChange={(v) => setGlobalQCParams({ ...globalQCParams, accuracy_level: v })}
                    >
                      <Select.Option value="fast">快速</Select.Option>
                      <Select.Option value="standard">标准</Select.Option>
                      <Select.Option value="accurate">精确</Select.Option>
                      <Select.Option value="custom">自定义</Select.Option>
                    </Select>
                  </Col>
                  <Col span={6}>
                    <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>溶剂模型</div>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={globalQCParams.solvent_model}
                      onChange={(v) => setGlobalQCParams({ ...globalQCParams, solvent_model: v })}
                    >
                      <Select.Option value="gas">气相</Select.Option>
                      <Select.Option value="pcm">PCM</Select.Option>
                      <Select.Option value="smd">SMD</Select.Option>
                    </Select>
                  </Col>
                  <Col span={6}>
                    <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>隐式溶剂</div>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={globalQCParams.solvent_name}
                      onChange={(v) => setGlobalQCParams({ ...globalQCParams, solvent_name: v })}
                      disabled={globalQCParams.solvent_model === 'gas'}
                    >
                      <Select.OptGroup label="常用溶剂">
                        <Select.Option value="water">Water</Select.Option>
                        <Select.Option value="acetonitrile">Acetonitrile</Select.Option>
                        <Select.Option value="acetone">Acetone</Select.Option>
                        <Select.Option value="dmso">DMSO</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="电池电解液溶剂">
                        <Select.Option value="ec">EC</Select.Option>
                        <Select.Option value="dmc">DMC</Select.Option>
                        <Select.Option value="emc">EMC</Select.Option>
                        <Select.Option value="dec">DEC</Select.Option>
                        <Select.Option value="pc">PC</Select.Option>
                        <Select.Option value="gbl">GBL</Select.Option>
                        <Select.Option value="dme">DME</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="其他溶剂">
                        <Select.Option value="methanol">Methanol</Select.Option>
                        <Select.Option value="ethanol">Ethanol</Select.Option>
                        <Select.Option value="thf">THF</Select.Option>
                        <Select.Option value="dcm">DCM</Select.Option>
                        <Select.Option value="diethyl_ether">Diethyl Ether</Select.Option>
                      </Select.OptGroup>
                    </Select>
                  </Col>
                  <Col span={6}>
                    <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>智能推荐</div>
                    <Switch
                      checked={globalQCParams.use_recommended_params}
                      onChange={(v) => setGlobalQCParams({ ...globalQCParams, use_recommended_params: v })}
                      checkedChildren="开启"
                      unCheckedChildren="关闭"
                    />
                  </Col>
                </Row>
                <div style={{ marginTop: 12, textAlign: 'right' }}>
                  <Space>
                    <Button size="small" onClick={() => setEditingGlobalQC(false)}>取消</Button>
                    <Button size="small" type="primary" onClick={handleSaveGlobalQC}>保存</Button>
                  </Space>
                </div>
              </div>
            ) : (
              <Descriptions bordered column={2} size="small" style={{ marginBottom: 16 }}>
                <Descriptions.Item label="精度等级">
                  <Tag color={
                    currentAccuracyLevel === 'fast' ? 'green' :
                    currentAccuracyLevel === 'standard' ? 'blue' :
                    currentAccuracyLevel === 'accurate' ? 'orange' : 'purple'
                  }>
                    {currentAccuracyLevel === 'fast' ? '快速' :
                     currentAccuracyLevel === 'standard' ? '标准' :
                     currentAccuracyLevel === 'accurate' ? '精确' : '自定义'}
                  </Tag>
                </Descriptions.Item>
                <Descriptions.Item label="默认溶剂模型">
                  <Tag color={
                    currentSolventModel === 'gas' ? 'default' :
                    currentSolventModel === 'pcm' ? 'blue' : 'cyan'
                  }>
                    {currentSolventModel === 'gas' ? '气相' :
                     currentSolventModel === 'pcm' ? 'PCM' : 'SMD'}
                  </Tag>
                </Descriptions.Item>
                <Descriptions.Item label="智能推荐">
                  <Tag color={currentUseRecommended ? 'green' : 'default'}>
                    {currentUseRecommended ? '已启用' : '未启用'}
                  </Tag>
                </Descriptions.Item>
                {currentSolventModel !== 'gas' && (
                  <Descriptions.Item label="隐式溶剂">
                    <Text code>{currentSolventName}</Text>
                  </Descriptions.Item>
                )}
              </Descriptions>
            )}

            {/* 重复计算检查提示 */}
            {duplicateCheckResult && (
              <Alert
                type={duplicateCheckResult.existing_count === duplicateCheckResult.total_molecules ? 'success' :
                      duplicateCheckResult.existing_count > 0 ? 'info' : 'warning'}
                showIcon
                icon={duplicateCheckResult.existing_count > 0 ? <CheckCircleOutlined /> : <SyncOutlined />}
                style={{ marginBottom: 16 }}
                message={
                  duplicateCheckResult.existing_count === duplicateCheckResult.total_molecules
                    ? `所有 ${duplicateCheckResult.total_molecules} 个分子都已有计算结果，将直接复用！`
                    : duplicateCheckResult.existing_count > 0
                    ? `${duplicateCheckResult.existing_count} 个分子已有结果（将复用），${duplicateCheckResult.new_count} 个分子需要新计算`
                    : `${duplicateCheckResult.total_molecules} 个分子都需要新计算`
                }
                description={
                  duplicateCheckResult.existing_count > 0 && (
                    <div style={{ fontSize: 12, marginTop: 4 }}>
                      <Text type="secondary">已有结果的分子：</Text>
                      {duplicateCheckResult.results
                        .filter(r => r.has_existing_result)
                        .map((r, idx) => (
                          <Tag key={idx} color="green" style={{ margin: '2px 4px' }}>
                            {r.molecule_name || r.smiles.substring(0, 20)}
                          </Tag>
                        ))}
                    </div>
                  )
                }
              />
            )}

            {/* 检查重复按钮 */}
            {!isSubmittedJob && allMolecules.length > 0 && (
              <div style={{ marginBottom: 16 }}>
                <Button
                  size="small"
                  icon={<SyncOutlined spin={checkingDuplicates} />}
                  onClick={checkQCDuplicates}
                  loading={checkingDuplicates}
                >
                  检查是否有已完成的计算结果
                </Button>
                <Text type="secondary" style={{ marginLeft: 8, fontSize: 12 }}>
                  （全局共享的QC计算结果可直接复用，节省计算时间）
                </Text>
              </div>
            )}

            {/* 分子列表 - 每个分子可单独编辑 */}
            {allMolecules.length > 0 && (
              <div>
                <Divider orientation="left" style={{ margin: '16px 0 12px' }}>
                  <Space>
                    <BulbOutlined />
                    将对 {allMolecules.length} 个分子进行QC计算（点击编辑可调整参数）
                  </Space>
                </Divider>
                <div>
                  {allMolecules.map((mol) => (
                    <div key={mol.key} style={{
                      padding: '12px 16px',
                      marginBottom: 10,
                      background: editingMolecule === mol.key ? '#f6f0ff' : '#fafafa',
                      borderRadius: 8,
                      border: editingMolecule === mol.key ? '1px solid #722ed1' : '1px solid #f0f0f0',
                      transition: 'all 0.3s'
                    }}>
                      {/* 分子标题行 */}
                      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 8 }}>
                        <Space>
                          <Text strong style={{ fontSize: 14 }}>{mol.name}</Text>
                          <Tag color={mol.type === 'solvent' ? 'blue' : mol.type === 'cation' ? 'green' : 'orange'}>
                            {mol.type === 'solvent' ? '溶剂' : mol.type === 'cation' ? '阳离子' : '阴离子'}
                          </Tag>
                          {mol.smiles && (
                            <Text type="secondary" style={{ fontSize: 11 }}>{mol.smiles}</Text>
                          )}
                        </Space>
                        {!isSubmittedJob && (
                          editingMolecule === mol.key ? (
                            <Space size="small">
                              <Button size="small" type="primary" onClick={() => handleSaveMolecule(mol.key)}>
                                确定
                              </Button>
                              <Button size="small" onClick={() => setEditingMolecule(null)}>
                                取消
                              </Button>
                            </Space>
                          ) : (
                            <Button
                              size="small"
                              type="text"
                              icon={<EditOutlined />}
                              onClick={() => setEditingMolecule(mol.key)}
                            >
                              编辑
                            </Button>
                          )
                        )}
                      </div>

                      {/* 参数展示/编辑 */}
                      {editingMolecule === mol.key ? (
                        // 编辑模式
                        <Row gutter={12}>
                          <Col span={8}>
                            <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>泛函</div>
                            <Select
                              size="small"
                              style={{ width: '100%' }}
                              value={moleculeParams[mol.key]?.functional || mol.params.functional}
                              onChange={(v) => handleMoleculeParamChange(mol.key, 'functional', v)}
                            >
                              <Select.OptGroup label="杂化泛函">
                                <Select.Option value="B3LYP">B3LYP</Select.Option>
                                <Select.Option value="B3PW91">B3PW91</Select.Option>
                                <Select.Option value="PBE0">PBE0</Select.Option>
                              </Select.OptGroup>
                              <Select.OptGroup label="范围分离泛函">
                                <Select.Option value="wB97X-D">ωB97X-D</Select.Option>
                                <Select.Option value="CAM-B3LYP">CAM-B3LYP</Select.Option>
                              </Select.OptGroup>
                              <Select.OptGroup label="Minnesota 泛函">
                                <Select.Option value="M06-2X">M06-2X</Select.Option>
                                <Select.Option value="M06">M06</Select.Option>
                              </Select.OptGroup>
                              <Select.OptGroup label="基础方法">
                                <Select.Option value="HF">HF</Select.Option>
                              </Select.OptGroup>
                            </Select>
                          </Col>
                          <Col span={8}>
                            <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>基组</div>
                            <Select
                              size="small"
                              style={{ width: '100%' }}
                              value={moleculeParams[mol.key]?.basis_set || mol.params.basis_set}
                              onChange={(v) => handleMoleculeParamChange(mol.key, 'basis_set', v)}
                            >
                              <Select.OptGroup label="Pople 基组">
                                <Select.Option value="STO-3G">STO-3G</Select.Option>
                                <Select.Option value="6-31G(d)">6-31G(d)</Select.Option>
                                <Select.Option value="6-31+G(d,p)">6-31+G(d,p)</Select.Option>
                                <Select.Option value="6-31++G(d,p)">6-31++G(d,p)</Select.Option>
                                <Select.Option value="6-311++G(d,p)">6-311++G(d,p)</Select.Option>
                              </Select.OptGroup>
                              <Select.OptGroup label="Def2 基组">
                                <Select.Option value="Def2-SVP">Def2-SVP</Select.Option>
                                <Select.Option value="Def2-TZVP">Def2-TZVP</Select.Option>
                              </Select.OptGroup>
                            </Select>
                          </Col>
                          <Col span={8}>
                            <div style={{ marginBottom: 4, fontSize: 12, color: '#666' }}>溶剂模型</div>
                            <Select
                              size="small"
                              style={{ width: '100%' }}
                              value={moleculeParams[mol.key]?.solvent_model || mol.params.solvent_model}
                              onChange={(v) => handleMoleculeParamChange(mol.key, 'solvent_model', v)}
                            >
                              <Select.Option value="gas">气相</Select.Option>
                              <Select.Option value="pcm">PCM</Select.Option>
                              <Select.Option value="smd">SMD</Select.Option>
                            </Select>
                          </Col>
                        </Row>
                      ) : (
                        // 展示模式
                        <>
                          <div style={{ fontSize: 12, color: '#666' }}>
                            <Space split={<span style={{ color: '#d9d9d9' }}>|</span>}>
                              <span>泛函: <Text code style={{ fontSize: 11 }}>{mol.params.functional}</Text></span>
                              <span>基组: <Text code style={{ fontSize: 11 }}>{mol.params.basis_set}</Text></span>
                              <span>溶剂: <Text code style={{ fontSize: 11 }}>{mol.params.solvent_model === 'gas' ? '气相' : mol.params.solvent_model.toUpperCase()}</Text></span>
                            </Space>
                          </div>
                          {mol.params.reason && (
                            <div style={{ fontSize: 11, color: '#999', marginTop: 4 }}>
                              💡 {mol.params.reason}
                            </div>
                          )}
                        </>
                      )}
                    </div>
                  ))}
                </div>
              </div>
            )}
          </Card>
        );
      })()}

      {/* 提交按钮 - 只有未提交的任务才显示 */}
      {!editMode && !isSubmittedJob && (
        <Card style={{ marginTop: 24 }}>
          <Space>
            <Button
              type="primary"
              size="large"
              icon={<ThunderboltOutlined />}
              onClick={handleSubmit}
              loading={submitting}
            >
              提交到集群
            </Button>
            <Button size="large" onClick={() => navigate('/workspace/jobs')}>
              返回任务列表
            </Button>
          </Space>
        </Card>
      )}

      {/* 已提交任务的返回按钮 */}
      {!editMode && isSubmittedJob && (
        <Card style={{ marginTop: 24 }}>
          <Button size="large" onClick={() => navigate('/workspace/jobs')}>
            返回任务列表
          </Button>
        </Card>
      )}
    </div>
  );
}

