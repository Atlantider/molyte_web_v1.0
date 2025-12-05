/**
 * 热力学循环计算氧化还原电位面板
 *
 * ⚠️ 高风险警告：
 * - 结果对方法/基组/溶剂模型/构型高度敏感
 * - 计算量大，经常不收敛
 * - 数值可能存在数百 mV 的系统性偏差
 * - 仅供研究参考，不应作为定量预测
 */
import React, { useState, useEffect } from 'react';
import {
  Card,
  Button,
  Form,
  Input,
  Select,
  InputNumber,
  Table,
  Tag,
  Space,
  Modal,
  Alert,
  Spin,
  Empty,
  Popconfirm,
  message,
  Tooltip,
  Typography,
  Divider,
  Row,
  Col,
  Statistic,
  Progress,
  Checkbox,
} from 'antd';
import {
  PlusOutlined,
  DeleteOutlined,
  PlayCircleOutlined,
  WarningOutlined,
  ExperimentOutlined,
  InfoCircleOutlined,
  ReloadOutlined,
} from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import {
  RedoxJobResponse,
  RedoxJobStatus,
  SpeciesConfig,
  RedoxCalculationMode,
  ClusterTypeInfo,
  listRedoxJobs,
  createRedoxJob,
  submitRedoxJob,
  deleteRedoxJob,
  getRedoxJob,
  getAvailableClustersForRedox,
  PhysicalConstants,
} from '../api/redox';

const { Text, Title, Paragraph } = Typography;

interface RedoxPotentialPanelProps {
  mdJobId?: number;
}

const RedoxPotentialPanel: React.FC<RedoxPotentialPanelProps> = ({ mdJobId }) => {
  const [jobs, setJobs] = useState<RedoxJobResponse[]>([]);
  const [loading, setLoading] = useState(false);
  const [createModalVisible, setCreateModalVisible] = useState(false);
  const [selectedJob, setSelectedJob] = useState<RedoxJobResponse | null>(null);
  const [form] = Form.useForm();

  // 可用的 Cluster 类型
  const [availableClusters, setAvailableClusters] = useState<ClusterTypeInfo[]>([]);
  const [loadingClusters, setLoadingClusters] = useState(false);

  useEffect(() => {
    loadJobs();
    if (mdJobId) {
      loadAvailableClusters();
    }
  }, [mdJobId]);

  const loadJobs = async () => {
    setLoading(true);
    try {
      const response = await listRedoxJobs({ md_job_id: mdJobId });
      setJobs(response.jobs);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载任务列表失败');
    } finally {
      setLoading(false);
    }
  };

  const loadAvailableClusters = async () => {
    if (!mdJobId) return;
    setLoadingClusters(true);
    try {
      const response = await getAvailableClustersForRedox(mdJobId, false);
      setAvailableClusters(response.cluster_types || []);
    } catch (error: any) {
      console.error('加载可用 Cluster 失败:', error);
    } finally {
      setLoadingClusters(false);
    }
  };

  const handleCreate = async (values: any) => {
    try {
      // 从选中的 cluster 类型构建物种列表
      const selectedTypes: string[] = values.selected_clusters || [];

      if (selectedTypes.length === 0) {
        message.error('请至少选择一个 Cluster 类型');
        return;
      }

      if (selectedTypes.length > 20) {
        message.error('物种数量不能超过 20 个');
        return;
      }

      // 根据选中的类型构建 species_list，包含 smiles 用于复用匹配
      const speciesList: SpeciesConfig[] = selectedTypes.map(typeName => {
        const clusterInfo = availableClusters.find(c => c.type_name === typeName);
        return {
          name: typeName,
          smiles: clusterInfo?.example_smiles || '',  // 传递 SMILES 用于 QC 复用
          charge: clusterInfo?.charge || 0,
          multiplicity: clusterInfo?.multiplicity || 1,
          redox_type: values.redox_type || 'oxidation',
        };
      });

      await createRedoxJob({
        md_job_id: mdJobId,
        config: {
          species_list: speciesList,
          mode: values.mode,
          functional: values.functional,
          basis_set: values.basis_set,
          solvent_model: values.solvent_model,
          solvent: values.solvent,
          use_dispersion: values.use_dispersion,
          li_reference_potential: values.li_reference_potential,
          reuse_existing_qc: values.reuse_existing_qc ?? true,
        },
      });

      message.success('任务创建成功');
      setCreateModalVisible(false);
      form.resetFields();
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '创建任务失败');
    }
  };

  const handleSubmit = async (jobId: number) => {
    try {
      await submitRedoxJob(jobId);
      message.success('任务已提交');
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '提交任务失败');
    }
  };

  const handleDelete = async (jobId: number) => {
    try {
      await deleteRedoxJob(jobId);
      message.success('任务已删除');
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除任务失败');
    }
  };

  const getStatusTag = (status: RedoxJobStatus) => {
    const statusConfig: Record<RedoxJobStatus, { color: string; text: string }> = {
      CREATED: { color: 'default', text: '已创建' },
      SUBMITTED: { color: 'processing', text: '已提交' },
      RUNNING: { color: 'processing', text: '运行中' },
      COMPLETED: { color: 'success', text: '已完成' },
      FAILED: { color: 'error', text: '失败' },
    };
    const config = statusConfig[status] || { color: 'default', text: status };
    return <Tag color={config.color}>{config.text}</Tag>;
  };

  const getModeTag = (mode: RedoxCalculationMode) => {
    const modeConfig: Record<RedoxCalculationMode, { color: string; text: string }> = {
      cheap: { color: 'green', text: '快速' },
      standard: { color: 'blue', text: '标准' },
      heavy: { color: 'red', text: '精确' },
    };
    const config = modeConfig[mode] || { color: 'default', text: mode };
    return <Tag color={config.color}>{config.text}</Tag>;
  };

  const columns = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 60,
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      width: 100,
      render: (status: RedoxJobStatus) => getStatusTag(status),
    },
    {
      title: '模式',
      key: 'mode',
      width: 80,
      render: (_: any, record: RedoxJobResponse) =>
        getModeTag((record.config?.mode || 'cheap') as RedoxCalculationMode),
    },
    {
      title: '物种数',
      key: 'species_count',
      width: 80,
      render: (_: any, record: RedoxJobResponse) =>
        record.config?.species_list?.length || 0,
    },
    {
      title: '进度',
      dataIndex: 'progress',
      key: 'progress',
      width: 120,
      render: (progress: number) => (
        <Progress percent={Math.round(progress)} size="small" />
      ),
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 160,
      render: (date: string) => new Date(date).toLocaleString(),
    },
    {
      title: '操作',
      key: 'actions',
      width: 200,
      render: (_: any, record: RedoxJobResponse) => (
        <Space>
          {record.status === 'CREATED' && (
            <Button
              type="primary"
              size="small"
              icon={<PlayCircleOutlined />}
              onClick={() => handleSubmit(record.id)}
            >
              提交
            </Button>
          )}
          {record.status === 'COMPLETED' && (
            <Button
              size="small"
              onClick={() => setSelectedJob(record)}
            >
              查看结果
            </Button>
          )}
          {record.status !== 'RUNNING' && (
            <Popconfirm
              title="确定要删除此任务吗？"
              onConfirm={() => handleDelete(record.id)}
            >
              <Button size="small" danger icon={<DeleteOutlined />} />
            </Popconfirm>
          )}
        </Space>
      ),
    },
  ];

  // 结果展示组件
  const ResultView = ({ job }: { job: RedoxJobResponse }) => {
    const result = job.result;
    if (!result) return <Empty description="暂无结果" />;

    // 电位分布图
    const potentialChartOption = {
      title: { text: '氧化还原电位分布', left: 'center' },
      tooltip: { trigger: 'axis' },
      legend: { top: 30, data: ['氧化电位', '还原电位'] },
      xAxis: {
        type: 'category',
        data: result.species_results.map(s => s.name),
        axisLabel: { rotate: 45 },
      },
      yAxis: {
        type: 'value',
        name: '电位 (V vs Li+/Li)',
      },
      series: [
        {
          name: '氧化电位',
          type: 'bar',
          data: result.species_results
            .filter(s => s.redox_type === 'oxidation')
            .map(s => s.e_vs_li_v),
          itemStyle: { color: '#f5222d' },
        },
        {
          name: '还原电位',
          type: 'bar',
          data: result.species_results
            .filter(s => s.redox_type === 'reduction')
            .map(s => s.e_vs_li_v),
          itemStyle: { color: '#1890ff' },
        },
      ],
    };

    return (
      <div>
        <Alert
          type="warning"
          showIcon
          icon={<WarningOutlined />}
          message="结果仅供研究参考"
          description={result.reference_note}
          style={{ marginBottom: 16 }}
        />

        {result.electrochemical_window_v && (
          <Row gutter={16} style={{ marginBottom: 16 }}>
            <Col span={8}>
              <Statistic
                title="氧化极限 (5% 分位数)"
                value={result.oxidation_limit_v?.toFixed(2)}
                suffix="V"
              />
            </Col>
            <Col span={8}>
              <Statistic
                title="还原极限 (95% 分位数)"
                value={result.reduction_limit_v?.toFixed(2)}
                suffix="V"
              />
            </Col>
            <Col span={8}>
              <Statistic
                title="电化学窗口"
                value={result.electrochemical_window_v?.toFixed(2)}
                suffix="V"
              />
            </Col>
          </Row>
        )}

        <ReactECharts option={potentialChartOption} style={{ height: 300 }} />

        <Divider>各物种详细结果</Divider>
        <Table
          dataSource={result.species_results}
          rowKey="name"
          size="small"
          columns={[
            { title: '物种', dataIndex: 'name', width: 120 },
            { title: '类型', dataIndex: 'redox_type', width: 80,
              render: (t: string) => t === 'oxidation' ? '氧化' : '还原' },
            { title: 'E° vs Li+/Li (V)', dataIndex: 'e_vs_li_v', width: 120,
              render: (v: number) => v?.toFixed(3) },
            { title: 'ΔG(sol) (kcal/mol)', dataIndex: 'dg_sol_kcal', width: 140,
              render: (v: number) => v?.toFixed(2) },
            { title: '收敛', dataIndex: 'converged', width: 60,
              render: (v: boolean) => v ? <Tag color="green">是</Tag> : <Tag color="red">否</Tag> },
          ]}
        />

        {result.global_warnings.length > 0 && (
          <>
            <Divider>警告信息</Divider>
            {result.global_warnings.map((w, i) => (
              <Alert key={i} type="warning" message={w} style={{ marginBottom: 8 }} />
            ))}
          </>
        )}
      </div>
    );
  };

  return (
    <div>
      {/* 高风险警告 */}
      <Alert
        type="error"
        showIcon
        icon={<WarningOutlined />}
        message="⚠️ 高风险功能 - 热力学循环计算氧化还原电位"
        description={
          <ul style={{ margin: 0, paddingLeft: 20 }}>
            <li>结果对方法/基组/溶剂模型/构型<strong>高度敏感</strong></li>
            <li>计算量大，经常<strong>不收敛</strong></li>
            <li>数值可能存在<strong>数百 mV</strong> 的系统性偏差</li>
            <li>仅供<strong>研究参考</strong>，不应作为定量预测</li>
          </ul>
        }
        style={{ marginBottom: 16 }}
      />

      <Card
        title={
          <Space>
            <ExperimentOutlined />
            <span>热力学循环任务</span>
          </Space>
        }
        extra={
          <Space>
            <Button icon={<ReloadOutlined />} onClick={loadJobs}>
              刷新
            </Button>
            <Button
              type="primary"
              icon={<PlusOutlined />}
              onClick={() => setCreateModalVisible(true)}
            >
              创建任务
            </Button>
          </Space>
        }
      >
        <Table
          loading={loading}
          dataSource={jobs}
          columns={columns}
          rowKey="id"
          size="small"
          pagination={{ pageSize: 10 }}
        />
      </Card>

      {/* 创建任务弹窗 */}
      <Modal
        title="创建热力学循环任务"
        open={createModalVisible}
        onCancel={() => setCreateModalVisible(false)}
        footer={null}
        width={700}
      >
        <Alert
          type="warning"
          message="请仔细阅读以下说明后再创建任务"
          style={{ marginBottom: 16 }}
        />
        <Form
          form={form}
          layout="vertical"
          onFinish={handleCreate}
          initialValues={{
            redox_type: 'oxidation',
            mode: 'cheap',
            functional: 'B3LYP',
            basis_set: '6-31G*',
            solvent_model: 'SMD',
            solvent: 'water',
            use_dispersion: true,
            li_reference_potential: PhysicalConstants.LI_ABSOLUTE_POTENTIAL_VS_SHE,
            reuse_existing_qc: true,
          }}
        >
          {/* Cluster 选择 */}
          {availableClusters.length > 0 ? (
            <Form.Item
              name="selected_clusters"
              label={
                <Tooltip title="选择要计算氧化还原电位的 Cluster 类型">
                  选择 Cluster 类型 <InfoCircleOutlined />
                </Tooltip>
              }
              rules={[{ required: true, message: '请至少选择一个 Cluster 类型' }]}
            >
              <Select
                mode="multiple"
                placeholder="选择要计算的 Cluster 类型"
                loading={loadingClusters}
                style={{ width: '100%' }}
                optionLabelProp="label"
              >
                {availableClusters.map(cluster => (
                  <Select.Option
                    key={cluster.type_name}
                    value={cluster.type_name}
                    label={cluster.type_name}
                  >
                    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                      <span>{cluster.type_name}</span>
                      <span style={{ color: '#888' }}>
                        {cluster.count} 个样本 | 电荷: {cluster.charge}
                      </span>
                    </div>
                  </Select.Option>
                ))}
              </Select>
            </Form.Item>
          ) : (
            <Alert
              type="warning"
              message="没有可用的 Cluster"
              description={
                mdJobId
                  ? "该 MD 任务没有已完成的 QC 计算。请先运行去溶剂化能计算或 Cluster QC 计算。"
                  : "请先选择一个 MD 任务。"
              }
              style={{ marginBottom: 16 }}
            />
          )}

          <Row gutter={16}>
            <Col span={8}>
              <Form.Item name="redox_type" label="反应类型" rules={[{ required: true }]}>
                <Select>
                  <Select.Option value="oxidation">氧化 (失电子)</Select.Option>
                  <Select.Option value="reduction">还原 (得电子)</Select.Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item name="mode" label="计算模式" rules={[{ required: true }]}>
                <Select>
                  <Select.Option value="cheap">快速 (垂直近似)</Select.Option>
                  <Select.Option value="standard">标准 (gas优化)</Select.Option>
                  <Select.Option value="heavy">精确 (全优化)</Select.Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item name="functional" label="泛函">
                <Select>
                  <Select.Option value="B3LYP">B3LYP</Select.Option>
                  <Select.Option value="M062X">M06-2X</Select.Option>
                  <Select.Option value="wB97XD">ωB97X-D</Select.Option>
                </Select>
              </Form.Item>
            </Col>
          </Row>

          <Row gutter={16}>
            <Col span={8}>
              <Form.Item name="basis_set" label="基组">
                <Select>
                  <Select.Option value="6-31G*">6-31G*</Select.Option>
                  <Select.Option value="6-311+G(d,p)">6-311+G(d,p)</Select.Option>
                  <Select.Option value="def2-SVP">def2-SVP</Select.Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item name="solvent_model" label="溶剂模型">
                <Select>
                  <Select.Option value="SMD">SMD</Select.Option>
                  <Select.Option value="PCM">PCM</Select.Option>
                  <Select.Option value="CPCM">CPCM</Select.Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={8}>
              <Form.Item name="solvent" label="溶剂">
                <Input placeholder="water" />
              </Form.Item>
            </Col>
          </Row>

          <Row gutter={16}>
            <Col span={8}>
              <Form.Item name="li_reference_potential" label="Li+/Li 参考电位 (V)">
                <InputNumber style={{ width: '100%' }} step={0.01} />
              </Form.Item>
            </Col>
          </Row>

          <Form.Item name="reuse_existing_qc" valuePropName="checked">
            <Checkbox>
              复用已有 QC 结果（推荐，节省计算资源）
            </Checkbox>
          </Form.Item>

          <Form.Item>
            <Space>
              <Button type="primary" htmlType="submit">
                创建任务
              </Button>
              <Button onClick={() => setCreateModalVisible(false)}>
                取消
              </Button>
            </Space>
          </Form.Item>
        </Form>
      </Modal>

      {/* 结果查看弹窗 */}
      <Modal
        title={`任务 #${selectedJob?.id} 结果`}
        open={!!selectedJob}
        onCancel={() => setSelectedJob(null)}
        footer={null}
        width={900}
      >
        {selectedJob && <ResultView job={selectedJob} />}
      </Modal>
    </div>
  );
};

export default RedoxPotentialPanel;
