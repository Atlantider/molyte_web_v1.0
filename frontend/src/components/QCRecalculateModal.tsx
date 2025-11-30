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
} from 'antd';
import { ExperimentOutlined, ThunderboltOutlined } from '@ant-design/icons';
import type { QCJob } from '../types/qc';
import { recalculateQCJob } from '../api/qc';
import { getSolventModelText } from '../utils/qc';

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
    { value: 'gas', label: '气相', description: '无溶剂效应' },
    { value: 'pcm', label: 'PCM', description: '极化连续介质模型' },
    { value: 'smd', label: 'SMD', description: '溶剂密度模型' },
  ];

  // 常用溶剂
  const solvents = [
    { value: 'Water', label: '水 (Water)' },
    { value: 'Acetonitrile', label: '乙腈 (Acetonitrile)' },
    { value: 'Methanol', label: '甲醇 (Methanol)' },
    { value: 'Ethanol', label: '乙醇 (Ethanol)' },
    { value: 'Acetone', label: '丙酮 (Acetone)' },
    { value: 'DiMethylSulfoxide', label: 'DMSO' },
    { value: 'Dichloromethane', label: '二氯甲烷 (DCM)' },
    { value: 'Chloroform', label: '氯仿 (Chloroform)' },
    { value: 'TetraHydroFuran', label: 'THF' },
  ];

  const handleSubmit = async () => {
    if (!job) return;

    try {
      const values = await form.validateFields();
      setLoading(true);

      // 构建溶剂配置
      let solventConfig = undefined;
      if (values.solvent_model !== 'gas') {
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
          style={{ background: '#fafafa' }}
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
                  label="溶剂模型"
                  rules={[{ required: true, message: '请选择溶剂模型' }]}
                >
                  <Select
                    placeholder="选择溶剂模型"
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
                <Form.Item
                  name="solvent_name"
                  label="溶剂"
                  rules={[
                    {
                      required: selectedSolventModel !== 'gas',
                      message: '请选择溶剂',
                    },
                  ]}
                >
                  <Select
                    placeholder="选择溶剂"
                    disabled={selectedSolventModel === 'gas'}
                  >
                    {solvents.map(s => (
                      <Select.Option key={s.value} value={s.value}>
                        {s.label}
                      </Select.Option>
                    ))}
                  </Select>
                </Form.Item>
              </Col>
            </Row>

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

