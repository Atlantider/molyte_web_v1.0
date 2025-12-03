/**
 * QC 重新计算对话框组件
 */
import { useState } from 'react';
import {
  Modal,
  Form,
  Select,
  message,
  Descriptions,
  Card,
  Space,
  Alert,
  Row,
  Col,
  InputNumber,
  theme,
} from 'antd';
import { ExperimentOutlined, ThunderboltOutlined } from '@ant-design/icons';
import type { QCJob } from '../types/qc';
import { recalculateQCJob } from '../api/qc';
import { getSolventModelText } from '../utils/qc';
import { useThemeStore } from '../stores/themeStore';

interface QCRecalculateModalProps {
  visible: boolean;
  job: QCJob | null;
  onClose: () => void;
  onSuccess: (newJob: QCJob) => void;
}

export default function QCRecalculateModal({
  visible,
  job,
  onClose,
  onSuccess,
}: QCRecalculateModalProps) {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [form] = Form.useForm();
  const [loading, setLoading] = useState(false);
  const [selectedSolventModel, setSelectedSolventModel] = useState<string>('gas');

  // 泛函选项
  const functionals = [
    { value: 'HF', label: 'HF', description: 'Hartree-Fock' },
    { value: 'B3LYP', label: 'B3LYP', description: '混合泛函' },
    { value: 'M062X', label: 'M06-2X', description: 'Minnesota泛函' },
    { value: 'wB97XD', label: 'ωB97X-D', description: '长程校正泛函' },
    { value: 'PBE0', label: 'PBE0', description: 'PBE混合泛函' },
  ];

  // 基组选项
  const basisSets = [
    { value: 'STO-3G', label: 'STO-3G', description: '最小基组' },
    { value: '6-31G(d)', label: '6-31G(d)', description: '标准基组' },
    { value: '6-31G(d,p)', label: '6-31G(d,p)', description: '标准基组+极化' },
    { value: '6-31++G(d,p)', label: '6-31++G(d,p)', description: '标准基组+弥散+极化' },
    { value: '6-311G(d,p)', label: '6-311G(d,p)', description: '大基组' },
    { value: '6-311++G(d,p)', label: '6-311++G(d,p)', description: '大基组+弥散+极化' },
    { value: 'Def2-TZVP', label: 'Def2-TZVP', description: 'Def2三重ζ基组' },
  ];

  // 溶剂模型选项
  const solventModels = [
    { value: 'gas', label: '气相 (Gas Phase)', description: '无溶剂效应' },
    { value: 'pcm', label: 'PCM', description: '极化连续介质模型' },
    { value: 'smd', label: 'SMD', description: '溶剂密度模型（更精确）' },
    { value: 'custom', label: '自定义', description: '手动设置介电常数等参数' },
  ];

  // 常用溶剂 - 按介电常数分组
  const solventGroups = [
    {
      label: '📌 水系电解液 (ε>50)',
      options: [
        { value: 'Water', label: '水 (Water) ε=78.4' },
      ],
    },
    {
      label: '📌 高介电常数 (ε=40-90)',
      options: [
        { value: 'DiMethylSulfoxide', label: 'DMSO ε=46.8 (离子液体参考)' },
        { value: '1,2-EthaneDiol', label: '乙二醇 ε=40.2' },
      ],
    },
    {
      label: '📌 中等介电常数 (ε=15-40)',
      options: [
        { value: 'Acetonitrile', label: '乙腈 ε=35.7' },
        { value: 'Methanol', label: '甲醇 ε=32.6' },
        { value: 'Ethanol', label: '乙醇 ε=24.9' },
        { value: 'Acetone', label: '丙酮 ε=20.5 (高浓电解液)' },
        { value: '1-Propanol', label: '正丙醇 ε=20.5' },
      ],
    },
    {
      label: '📌 低介电常数 (ε<15) - DMC/EMC/DEC体系',
      options: [
        { value: 'DiChloroEthane', label: '二氯乙烷 ε=10.1' },
        { value: 'Dichloromethane', label: '二氯甲烷 ε=8.9' },
        { value: 'TetraHydroFuran', label: '四氢呋喃 (THF) ε=7.4' },
        { value: 'Chloroform', label: '氯仿 ε=4.7 (线性碳酸酯参考)' },
        { value: 'DiethylEther', label: '乙醚 ε=4.2' },
        { value: 'CarbonTetraChloride', label: '四氯化碳 ε=2.2' },
        { value: 'Toluene', label: '甲苯 ε=2.4' },
        { value: 'Benzene', label: '苯 ε=2.3' },
      ],
    },
  ];

  const handleSubmit = async () => {
    if (!job) return;

    try {
      const values = await form.validateFields();
      setLoading(true);

      // 构建溶剂配置
      let solventConfig = undefined;
      if (values.solvent_model === 'custom') {
        // 自定义溶剂参数
        solventConfig = {
          model: 'custom',
          solvent_name: 'Custom',
          eps: values.custom_eps,
          eps_inf: values.custom_eps_inf,
          hbond_acidity: values.custom_hbond_acidity,
          hbond_basicity: values.custom_hbond_basicity,
          surface_tension: values.custom_surface_tension,
          carbon_aromaticity: values.custom_carbon_aromaticity,
          halogenicity: values.custom_halogenicity,
        };
      } else if (values.solvent_model !== 'gas') {
        solventConfig = {
          model: values.solvent_model,
          solvent_name: values.solvent_name,
        };
      }

      const newJob = await recalculateQCJob(job.id, {
        functional: values.functional,
        basis_set: values.basis_set,
        solvent_config: solventConfig,
        slurm_partition: values.slurm_partition,
        slurm_cpus: values.slurm_cpus,
        slurm_time: values.slurm_time,
      });

      message.success(`重新计算任务已创建 (ID: ${newJob.id})，请前往任务列表提交`);
      form.resetFields();
      onSuccess(newJob);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '创建重新计算任务失败');
    } finally {
      setLoading(false);
    }
  };

  const handleCancel = () => {
    form.resetFields();
    onClose();
  };

  if (!job) return null;

  return (
    <Modal
      title={
        <Space>
          <ThunderboltOutlined style={{ color: '#722ed1' }} />
          <span>重新计算 QC 任务</span>
        </Space>
      }
      open={visible}
      onOk={handleSubmit}
      onCancel={handleCancel}
      confirmLoading={loading}
      width={700}
      okText="创建重新计算任务"
      cancelText="取消"
      okButtonProps={{
        icon: <ExperimentOutlined />,
        style: {
          background: 'linear-gradient(135deg, #722ed1 0%, #9254de 100%)',
          border: 'none',
        },
      }}
    >
      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
        {/* 原任务信息 */}
        <Card
          size="small"
          title="原任务信息"
          style={{ background: token.colorBgContainer }}
        >
          <Descriptions size="small" column={2}>
            <Descriptions.Item label="任务ID">{job.id}</Descriptions.Item>
            <Descriptions.Item label="分子名称">{job.molecule_name}</Descriptions.Item>
            <Descriptions.Item label="SMILES" span={2}>
              <code style={{ fontSize: 12 }}>{job.smiles}</code>
            </Descriptions.Item>
            <Descriptions.Item label="电荷">{job.charge}</Descriptions.Item>
            <Descriptions.Item label="自旋多重度">{job.spin_multiplicity}</Descriptions.Item>
            <Descriptions.Item label="泛函">{job.functional}</Descriptions.Item>
            <Descriptions.Item label="基组">{job.basis_set}</Descriptions.Item>
            <Descriptions.Item label="溶剂模型" span={2}>
              {getSolventModelText(job.config)}
            </Descriptions.Item>
          </Descriptions>
        </Card>

        <Alert
          message="提示"
          description="重新计算将创建一个新的QC任务，复用原任务的分子信息，但使用新的计算参数。新任务创建后需要手动提交。"
          type="info"
          showIcon
        />

        {/* 新计算参数 */}
        <Card size="small" title="新计算参数">
          <Form
            form={form}
            layout="vertical"
            initialValues={{
              functional: job.functional,
              basis_set: job.basis_set,
              solvent_model: job.config?.solvent_config?.model || 'gas',
              solvent_name: job.config?.solvent_config?.solvent_name || 'Water',
              slurm_partition: 'cpu',
              slurm_cpus: 16,
              slurm_time: 7200,
            }}
          >
            <Row gutter={16}>
              <Col span={12}>
                <Form.Item
                  name="functional"
                  label="泛函"
                  rules={[{ required: true, message: '请选择泛函' }]}
                >
                  <Select placeholder="选择泛函">
                    {functionals.map(f => (
                      <Select.Option key={f.value} value={f.value}>
                        {f.label} - {f.description}
                      </Select.Option>
                    ))}
                  </Select>
                </Form.Item>
              </Col>
              <Col span={12}>
                <Form.Item
                  name="basis_set"
                  label="基组"
                  rules={[{ required: true, message: '请选择基组' }]}
                >
                  <Select placeholder="选择基组">
                    {basisSets.map(bs => (
                      <Select.Option key={bs.value} value={bs.value}>
                        {bs.label} - {bs.description}
                      </Select.Option>
                    ))}
                  </Select>
                </Form.Item>
              </Col>
            </Row>

            <Row gutter={16}>
              <Col span={12}>
                <Form.Item
                  name="solvent_model"
                  label="溶剂环境"
                  rules={[{ required: true, message: '请选择溶剂环境' }]}
                  tooltip={
                    <div>
                      <p><strong>气相 (Gas)</strong>: 真空环境，无溶剂效应</p>
                      <p><strong>PCM</strong>: 极化连续介质模型</p>
                      <p><strong>SMD</strong>: 溶剂密度模型（更精确）</p>
                      <p><strong>自定义</strong>: 手动设置介电常数等参数</p>
                    </div>
                  }
                >
                  <Select
                    placeholder="选择溶剂环境"
                    onChange={setSelectedSolventModel}
                  >
                    {solventModels.map(sm => (
                      <Select.Option key={sm.value} value={sm.value}>
                        {sm.label} - {sm.description}
                      </Select.Option>
                    ))}
                  </Select>
                </Form.Item>
              </Col>
              <Col span={12}>
                {(selectedSolventModel === 'pcm' || selectedSolventModel === 'smd') && (
                  <Form.Item
                    name="solvent_name"
                    label="隐式溶剂"
                    tooltip={
                      <div>
                        <p><strong>选择原则</strong>：选择介电常数(ε)接近您电解液的溶剂</p>
                        <p>• 水系电解液 → Water (ε=78.4)</p>
                        <p>• 高浓电解液 → Acetone (ε=20.5)</p>
                        <p>• DMC/EMC体系 → Chloroform (ε≈4.7)</p>
                      </div>
                    }
                    rules={[{ required: true, message: '请选择溶剂' }]}
                  >
                    <Select placeholder="选择隐式溶剂" showSearch>
                      {solventGroups.map(group => (
                        <Select.OptGroup key={group.label} label={group.label}>
                          {group.options.map(s => (
                            <Select.Option key={s.value} value={s.value}>
                              {s.label}
                            </Select.Option>
                          ))}
                        </Select.OptGroup>
                      ))}
                    </Select>
                  </Form.Item>
                )}
              </Col>
            </Row>

            {/* 自定义溶剂参数 */}
            {selectedSolventModel === 'custom' && (
              <Card size="small" style={{ marginBottom: 16, background: mode === 'dark' ? 'rgba(250, 173, 20, 0.15)' : '#fffbe6', borderColor: token.colorWarning }}>
                <div style={{ marginBottom: 8, fontWeight: 500 }}>🔧 自定义溶剂参数（SMD模型）</div>
                <Row gutter={[8, 8]}>
                  <Col span={8}>
                    <Form.Item name="custom_eps" label="介电常数 ε" style={{ marginBottom: 4 }} rules={[{ required: true, message: '请输入介电常数' }]}>
                      <InputNumber style={{ width: '100%' }} placeholder="如: 89.6 (EC)" step={0.1} min={1} />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="custom_eps_inf" label="光学介电常数 n²" style={{ marginBottom: 4 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="如: 2.2" step={0.01} min={1} />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="custom_hbond_acidity" label="氢键酸度 α" style={{ marginBottom: 4 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="custom_hbond_basicity" label="氢键碱度 β" style={{ marginBottom: 4 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="custom_surface_tension" label="表面张力 γ" style={{ marginBottom: 4 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="cal/mol·Å²" step={0.1} />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="custom_carbon_aromaticity" label="芳香碳比例 φ" style={{ marginBottom: 4 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                    </Form.Item>
                  </Col>
                  <Col span={8}>
                    <Form.Item name="custom_halogenicity" label="卤素比例 ψ" style={{ marginBottom: 4 }}>
                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                    </Form.Item>
                  </Col>
                </Row>
                <Alert
                  type="info"
                  showIcon
                  style={{ marginTop: 8 }}
                  message={<span style={{ fontSize: 11 }}>常用电解液介电常数：EC(ε≈89.6), PC(ε≈64.9), DMC(ε≈3.1), EMC(ε≈2.9), DEC(ε≈2.8)</span>}
                />
              </Card>
            )}

            <Card
              size="small"
              title="计算资源配置"
              style={{ marginTop: 16, background: '#f5f5f5' }}
            >
              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item name="slurm_partition" label="队列">
                    <Select>
                      <Select.Option value="cpu">cpu</Select.Option>
                      <Select.Option value="gpu">gpu</Select.Option>
                    </Select>
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item name="slurm_cpus" label="CPU核心数">
                    <InputNumber min={1} max={64} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item name="slurm_time" label="最大时间(分钟)">
                    <InputNumber min={10} max={43200} style={{ width: '100%' }} />
                  </Form.Item>
                </Col>
              </Row>
            </Card>
          </Form>
        </Card>
      </Space>
    </Modal>
  );
}

