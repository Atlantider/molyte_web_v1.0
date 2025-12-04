/**
 * 重组能计算面板 (Marcus 理论)
 *
 * ⚠️ 极高风险警告：
 * - 需要 4 次几何优化 + 4 次单点，计算量极大
 * - 经常不收敛，特别是带电物种
 * - 结果对初始构型极其敏感
 * - 仅供高级研究使用
 */
import React, { useState, useEffect } from 'react';
import {
  Card,
  Button,
  Form,
  Input,
  Select,
  Table,
  Tag,
  Space,
  Modal,
  Alert,
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
} from 'antd';
import {
  PlusOutlined,
  DeleteOutlined,
  PlayCircleOutlined,
  WarningOutlined,
  ExperimentOutlined,
  InfoCircleOutlined,
  ReloadOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import {
  ReorgEnergyJobResponse,
  ReorgEnergyJobStatus,
  ReorgSpeciesConfig,
  ClusterTypeInfo,
  createReorgEnergyJob,
  submitReorgEnergyJob,
  deleteReorgEnergyJob,
  getReorgEnergyJob,
  listReorgEnergyJobs,
  getAvailableClustersForRedox,
} from '../api/redox';

const { Text, Title } = Typography;

interface ReorganizationEnergyPanelProps {
  mdJobId?: number;
}

const ReorganizationEnergyPanel: React.FC<ReorganizationEnergyPanelProps> = ({ mdJobId }) => {
  const [jobs, setJobs] = useState<ReorgEnergyJobResponse[]>([]);
  const [loading, setLoading] = useState(false);
  const [createModalVisible, setCreateModalVisible] = useState(false);
  const [selectedJob, setSelectedJob] = useState<ReorgEnergyJobResponse | null>(null);
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
      const jobList = await listReorgEnergyJobs({ md_job_id: mdJobId });
      setJobs(jobList || []);
    } catch (error: any) {
      // 如果 API 不存在，显示空列表
      setJobs([]);
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

      if (selectedTypes.length > 5) {
        message.error('重组能计算量极大，物种数量不能超过 5 个');
        return;
      }

      // 根据选中的类型构建 species_list
      const speciesList: ReorgSpeciesConfig[] = selectedTypes.map(typeName => {
        const clusterInfo = availableClusters.find(c => c.type_name === typeName);
        const charge = clusterInfo?.charge || 0;
        return {
          name: typeName,
          charge_neutral: charge,
          charge_oxidized: charge + 1,  // 氧化态 +1
          multiplicity_neutral: clusterInfo?.multiplicity || 1,
          multiplicity_oxidized: 2,  // 氧化后通常为双重态
        };
      });

      await createReorgEnergyJob({
        md_job_id: mdJobId,
        config: {
          species_list: speciesList,
          functional: values.functional,
          basis_set: values.basis_set,
          use_dispersion: values.use_dispersion,
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
      await submitReorgEnergyJob(jobId);
      message.success('任务已提交');
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '提交任务失败');
    }
  };

  const handleDelete = async (jobId: number) => {
    try {
      await deleteReorgEnergyJob(jobId);
      message.success('任务已删除');
      loadJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除任务失败');
    }
  };

  const getStatusTag = (status: ReorgEnergyJobStatus) => {
    const statusConfig: Record<ReorgEnergyJobStatus, { color: string; text: string }> = {
      CREATED: { color: 'default', text: '已创建' },
      SUBMITTED: { color: 'processing', text: '已提交' },
      RUNNING: { color: 'processing', text: '运行中' },
      COMPLETED: { color: 'success', text: '已完成' },
      FAILED: { color: 'error', text: '失败' },
    };
    const config = statusConfig[status];
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
      render: (status: ReorgEnergyJobStatus) => getStatusTag(status),
    },
    {
      title: '物种数',
      key: 'species_count',
      width: 80,
      render: (_: any, record: ReorgEnergyJobResponse) =>
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
      render: (_: any, record: ReorgEnergyJobResponse) => (
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
  const ResultView = ({ job }: { job: ReorgEnergyJobResponse }) => {
    const result = job.result;
    if (!result) return <Empty description="暂无结果" />;

    // 重组能分布图
    const lambdaChartOption = {
      title: { text: '重组能分布 (Marcus 理论)', left: 'center' },
      tooltip: { trigger: 'axis' },
      legend: { top: 30, data: ['λ_ox', 'λ_red', 'λ_total'] },
      xAxis: {
        type: 'category',
        data: result.species_results.map(s => s.name),
        axisLabel: { rotate: 45 },
      },
      yAxis: {
        type: 'value',
        name: '重组能 (eV)',
      },
      series: [
        {
          name: 'λ_ox',
          type: 'bar',
          data: result.species_results.map(s => s.lambda_ox_ev),
          itemStyle: { color: '#f5222d' },
        },
        {
          name: 'λ_red',
          type: 'bar',
          data: result.species_results.map(s => s.lambda_red_ev),
          itemStyle: { color: '#1890ff' },
        },
        {
          name: 'λ_total',
          type: 'bar',
          data: result.species_results.map(s => s.lambda_total_ev),
          itemStyle: { color: '#52c41a' },
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

        {result.lambda_total_mean_ev && (
          <Row gutter={16} style={{ marginBottom: 16 }}>
            <Col span={8}>
              <Statistic
                title="平均 λ_ox"
                value={result.lambda_ox_mean_ev?.toFixed(3)}
                suffix="eV"
              />
            </Col>
            <Col span={8}>
              <Statistic
                title="平均 λ_red"
                value={result.lambda_red_mean_ev?.toFixed(3)}
                suffix="eV"
              />
            </Col>
            <Col span={8}>
              <Statistic
                title="平均 λ_total"
                value={result.lambda_total_mean_ev?.toFixed(3)}
                suffix="eV"
              />
            </Col>
          </Row>
        )}

        <ReactECharts option={lambdaChartOption} style={{ height: 300 }} />

        <Divider>各物种详细结果</Divider>
        <Table
          dataSource={result.species_results}
          rowKey="name"
          size="small"
          columns={[
            { title: '物种', dataIndex: 'name', width: 120 },
            { title: 'λ_ox (eV)', dataIndex: 'lambda_ox_ev', width: 100,
              render: (v: number) => v?.toFixed(4) },
            { title: 'λ_red (eV)', dataIndex: 'lambda_red_ev', width: 100,
              render: (v: number) => v?.toFixed(4) },
            { title: 'λ_total (eV)', dataIndex: 'lambda_total_ev', width: 100,
              render: (v: number) => v?.toFixed(4) },
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
      {/* 极高风险警告 */}
      <Alert
        type="error"
        showIcon
        icon={<ThunderboltOutlined />}
        message="⚠️ 极高风险功能 - 重组能计算 (Marcus 理论)"
        description={
          <ul style={{ margin: 0, paddingLeft: 20 }}>
            <li>需要 <strong>4 次几何优化 + 4 次单点</strong>，计算量极大</li>
            <li>经常<strong>不收敛</strong>，特别是带电物种</li>
            <li>结果对<strong>初始构型极其敏感</strong></li>
            <li>建议<strong>限制物种数量 ≤ 5</strong></li>
            <li>仅供<strong>高级研究</strong>使用</li>
          </ul>
        }
        style={{ marginBottom: 16 }}
      />

      <Card
        title={
          <Space>
            <ExperimentOutlined />
            <span>重组能计算任务</span>
            <Tag color="red">高级研究</Tag>
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
              danger
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
        title={
          <Space>
            <WarningOutlined style={{ color: '#ff4d4f' }} />
            <span>创建重组能计算任务</span>
          </Space>
        }
        open={createModalVisible}
        onCancel={() => setCreateModalVisible(false)}
        footer={null}
        width={700}
      >
        <Alert
          type="error"
          message="⚠️ 请仔细阅读以下警告"
          description={
            <div>
              <p>重组能计算是<strong>极其耗费资源</strong>的高级功能：</p>
              <ul>
                <li>每个物种需要 4 次几何优化 + 4 次单点计算</li>
                <li>计算时间可能长达数小时甚至数天</li>
                <li>带电物种经常不收敛</li>
                <li>请确保您了解 Marcus 理论的适用范围</li>
              </ul>
            </div>
          }
          style={{ marginBottom: 16 }}
        />
        <Form
          form={form}
          layout="vertical"
          onFinish={handleCreate}
          initialValues={{
            functional: 'B3LYP',
            basis_set: '6-31G*',
            use_dispersion: true,
          }}
        >
          {/* Cluster 选择 */}
          {availableClusters.length > 0 ? (
            <Form.Item
              name="selected_clusters"
              label={
                <Tooltip title="选择要计算重组能的 Cluster 类型（最多 5 个）">
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
                maxTagCount={5}
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
            <Col span={12}>
              <Form.Item name="functional" label="泛函">
                <Select>
                  <Select.Option value="B3LYP">B3LYP</Select.Option>
                  <Select.Option value="M062X">M06-2X</Select.Option>
                  <Select.Option value="wB97XD">ωB97X-D</Select.Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item name="basis_set" label="基组">
                <Select>
                  <Select.Option value="6-31G*">6-31G*</Select.Option>
                  <Select.Option value="6-311+G(d,p)">6-311+G(d,p)</Select.Option>
                  <Select.Option value="def2-SVP">def2-SVP</Select.Option>
                </Select>
              </Form.Item>
            </Col>
          </Row>

          <Form.Item>
            <Space>
              <Button type="primary" htmlType="submit" danger>
                确认创建 (我已了解风险)
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

export default ReorganizationEnergyPanel;
