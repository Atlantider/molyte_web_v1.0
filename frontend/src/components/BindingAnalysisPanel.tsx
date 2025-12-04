/**
 * Binding 分析面板
 * 用于创建、管理和展示 Li-配体 Binding Energy 分析任务
 */
import React, { useState, useEffect, useCallback } from 'react';
import {
  Card,
  Button,
  Table,
  Space,
  Tag,
  Modal,
  Form,
  Select,
  Checkbox,
  Progress,
  Statistic,
  Row,
  Col,
  Alert,
  Spin,
  Empty,
  message,
  Tooltip,
  Popconfirm,
  Typography,
  theme,
} from 'antd';
import {
  PlusOutlined,
  PlayCircleOutlined,
  DeleteOutlined,
  ReloadOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  CloseCircleOutlined,
  ExperimentOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import {
  getBindingAnalysisJobs,
  getBindingAnalysisJob,
  createBindingAnalysisJob,
  submitBindingAnalysisJob,
  deleteBindingAnalysisJob,
  getAvailableClustersForBinding,
  BindingAnalysisJob,
  AvailableClustersInfo,
} from '../api/binding';
import { convertEnergy, getUnitOptions } from '../utils/energyUnits';

const { Text, Title } = Typography;

interface BindingAnalysisPanelProps {
  mdJobId: number;
}

const statusConfig: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
  CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
  SUBMITTED: { color: 'blue', icon: <ClockCircleOutlined />, text: '已提交' },
  RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
  COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
  FAILED: { color: 'error', icon: <CloseCircleOutlined />, text: '失败' },
};

const BindingAnalysisPanel: React.FC<BindingAnalysisPanelProps> = ({ mdJobId }) => {
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(true);
  const [jobs, setJobs] = useState<BindingAnalysisJob[]>([]);
  const [selectedJob, setSelectedJob] = useState<BindingAnalysisJob | null>(null);
  const [clustersInfo, setClustersInfo] = useState<AvailableClustersInfo | null>(null);
  const [createModalVisible, setCreateModalVisible] = useState(false);
  const [detailModalVisible, setDetailModalVisible] = useState(false);
  const [form] = Form.useForm();
  const [unit, setUnit] = useState('kcal/mol');

  const loadJobs = useCallback(async () => {
    try {
      const data = await getBindingAnalysisJobs({ md_job_id: mdJobId });
      setJobs(data.items);
    } catch (err: any) {
      message.error('加载任务列表失败');
    } finally {
      setLoading(false);
    }
  }, [mdJobId]);

  const loadClustersInfo = useCallback(async () => {
    try {
      const info = await getAvailableClustersForBinding(mdJobId);
      setClustersInfo(info);
    } catch (err: any) {
      console.error('获取 cluster 信息失败:', err);
    }
  }, [mdJobId]);

  useEffect(() => {
    loadJobs();
    loadClustersInfo();
  }, [loadJobs, loadClustersInfo]);

  // 定期刷新运行中的任务
  useEffect(() => {
    const hasRunning = jobs.some(j => j.status === 'RUNNING' || j.status === 'SUBMITTED');
    if (!hasRunning) return;
    const interval = setInterval(loadJobs, 10000);
    return () => clearInterval(interval);
  }, [jobs, loadJobs]);

  const handleCreate = async (values: any) => {
    try {
      const job = await createBindingAnalysisJob({
        md_job_id: mdJobId,
        config: {
          composition_keys: values.composition_keys,
          functional: values.functional || 'B3LYP',
          basis_set: values.basis_set || '6-31G(d)',
          solvent_model: values.solvent_model || 'gas',
          reuse_existing_qc: values.reuse_existing_qc ?? true,
        },
      });
      message.success(`任务 ${job.id} 创建成功`);
      setCreateModalVisible(false);
      form.resetFields();
      loadJobs();
    } catch (err: any) {
      message.error(err.response?.data?.detail || '创建任务失败');
    }
  };

  const handleSubmit = async (jobId: number) => {
    try {
      await submitBindingAnalysisJob(jobId);
      message.success('任务已提交');
      loadJobs();
    } catch (err: any) {
      message.error(err.response?.data?.detail || '提交任务失败');
    }
  };

  const handleDelete = async (jobId: number) => {
    try {
      await deleteBindingAnalysisJob(jobId);
      message.success('任务已删除');
      loadJobs();
    } catch (err: any) {
      message.error(err.response?.data?.detail || '删除任务失败');
    }
  };

  const handleViewDetail = async (job: BindingAnalysisJob) => {
    try {
      const detail = await getBindingAnalysisJob(job.id);
      setSelectedJob(detail);
      setDetailModalVisible(true);
    } catch (err: any) {
      message.error('获取任务详情失败');
    }
  };

  const convert = (val: number) => convertEnergy(val, 'kcal/mol', unit);
  const format = (val: number) => convert(val).toFixed(2);

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
      render: (status: string) => {
        const cfg = statusConfig[status] || { color: 'default', icon: null, text: status };
        return <Tag color={cfg.color} icon={cfg.icon}>{cfg.text}</Tag>;
      },
    },
    {
      title: '进度',
      dataIndex: 'progress',
      key: 'progress',
      width: 120,
      render: (progress: number, record: BindingAnalysisJob) => (
        record.status === 'RUNNING' ? <Progress percent={Math.round(progress)} size="small" /> : `${progress}%`
      ),
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 160,
      render: (t: string) => new Date(t).toLocaleString(),
    },
    {
      title: '操作',
      key: 'actions',
      width: 200,
      render: (_: any, record: BindingAnalysisJob) => (
        <Space size="small">
          {record.status === 'CREATED' && (
            <Button size="small" type="primary" icon={<PlayCircleOutlined />} onClick={() => handleSubmit(record.id)}>
              提交
            </Button>
          )}
          <Button size="small" onClick={() => handleViewDetail(record)}>详情</Button>
          {['CREATED', 'FAILED', 'COMPLETED'].includes(record.status) && (
            <Popconfirm title="确定删除此任务？" onConfirm={() => handleDelete(record.id)}>
              <Button size="small" danger icon={<DeleteOutlined />} />
            </Popconfirm>
          )}
        </Space>
      ),
    },
  ];

  // 结果可视化
  const renderResultChart = () => {
    if (!selectedJob?.result?.per_cluster_results?.length) return <Empty description="无结果数据" />;

    const results = selectedJob.result.per_cluster_results;
    const chartData = results
      .filter(r => r.binding_energy_kcal !== undefined)
      .map(r => ({
        name: r.composition_key,
        value: convert(r.binding_energy_kcal!),
      }));

    const option = {
      title: { text: `Binding Energy (${unit})`, left: 'center', textStyle: { fontSize: 14 } },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'category', data: chartData.map(d => d.name), axisLabel: { rotate: 45 } },
      yAxis: { type: 'value', name: unit },
      series: [{ type: 'bar', data: chartData.map(d => d.value), itemStyle: { color: token.colorPrimary } }],
      grid: { left: 60, right: 20, bottom: 80, top: 50 },
    };

    return <ReactECharts option={option} style={{ height: 300 }} />;
  };

  if (loading) {
    return <Spin tip="加载中..." style={{ display: 'block', margin: '40px auto' }} />;
  }

  return (
    <Card
      title={<Space><ExperimentOutlined /> Li-配体 Binding 分析</Space>}
      extra={
        <Space>
          <Button icon={<ReloadOutlined />} onClick={loadJobs}>刷新</Button>
          <Button type="primary" icon={<PlusOutlined />} onClick={() => setCreateModalVisible(true)}>
            创建分析
          </Button>
        </Space>
      }
    >
      <Alert
        type="info"
        showIcon
        icon={<InfoCircleOutlined />}
        message="简化版 Binding Energy 计算"
        description="E_bind = E_cluster - (E_Li+ + Σ n_j × E_ligand_j)，复用已有 QC 结果，无需完整去溶剂化流程。"
        style={{ marginBottom: 16 }}
      />

      {clustersInfo && (
        <Row gutter={16} style={{ marginBottom: 16 }}>
          <Col span={8}>
            <Statistic title="可用 Cluster 类型" value={clustersInfo.composition_keys.length} />
          </Col>
          <Col span={8}>
            <Statistic title="总溶剂化结构" value={clustersInfo.total_solvation_structures} />
          </Col>
          <Col span={8}>
            <Statistic title="已有 QC 结果类型" value={Object.keys(clustersInfo.existing_qc_by_type).length} />
          </Col>
        </Row>
      )}

      <Table columns={columns} dataSource={jobs} rowKey="id" size="small" pagination={{ pageSize: 10 }} />

      {/* 创建任务 Modal */}
      <Modal
        title="创建 Binding 分析任务"
        open={createModalVisible}
        onCancel={() => setCreateModalVisible(false)}
        onOk={() => form.submit()}
        width={600}
      >
        <Form form={form} layout="vertical" onFinish={handleCreate}>
          <Form.Item name="composition_keys" label="选择 Cluster 类型（留空则分析全部）">
            <Select mode="multiple" placeholder="选择要分析的 cluster 类型" allowClear>
              {clustersInfo?.composition_keys.map(key => (
                <Select.Option key={key} value={key}>{key}</Select.Option>
              ))}
            </Select>
          </Form.Item>
          <Row gutter={16}>
            <Col span={12}>
              <Form.Item name="functional" label="泛函" initialValue="B3LYP">
                <Select options={[{ value: 'B3LYP' }, { value: 'M06-2X' }, { value: 'wB97X-D' }]} />
              </Form.Item>
            </Col>
            <Col span={12}>
              <Form.Item name="basis_set" label="基组" initialValue="6-31G(d)">
                <Select options={[{ value: '6-31G(d)' }, { value: '6-311++G(d,p)' }, { value: 'def2-SVP' }]} />
              </Form.Item>
            </Col>
          </Row>
          <Form.Item name="reuse_existing_qc" valuePropName="checked" initialValue={true}>
            <Checkbox>复用已有 QC 结果（推荐）</Checkbox>
          </Form.Item>
        </Form>
      </Modal>

      {/* 详情 Modal */}
      <Modal
        title={`Binding 分析任务 #${selectedJob?.id}`}
        open={detailModalVisible}
        onCancel={() => setDetailModalVisible(false)}
        footer={null}
        width={800}
      >
        {selectedJob && (
          <div>
            <Row gutter={16} style={{ marginBottom: 16 }}>
              <Col span={6}><Text type="secondary">状态</Text><br /><Tag color={statusConfig[selectedJob.status]?.color}>{statusConfig[selectedJob.status]?.text}</Tag></Col>
              <Col span={6}><Text type="secondary">进度</Text><br />{selectedJob.progress}%</Col>
              <Col span={6}><Text type="secondary">创建时间</Text><br />{new Date(selectedJob.created_at).toLocaleString()}</Col>
              <Col span={6}><Text type="secondary">完成时间</Text><br />{selectedJob.finished_at ? new Date(selectedJob.finished_at).toLocaleString() : '-'}</Col>
            </Row>

            {selectedJob.error_message && (
              <Alert type="error" message="错误信息" description={selectedJob.error_message} style={{ marginBottom: 16 }} />
            )}

            {selectedJob.result && (
              <>
                <Row gutter={16} style={{ marginBottom: 16 }}>
                  <Col span={8}><Statistic title="总 Cluster" value={selectedJob.result.summary.total_clusters} /></Col>
                  <Col span={8}><Statistic title="已完成" value={selectedJob.result.summary.completed_clusters} valueStyle={{ color: token.colorSuccess }} /></Col>
                  <Col span={8}><Statistic title="失败" value={selectedJob.result.summary.failed_clusters} valueStyle={{ color: token.colorError }} /></Col>
                </Row>
                {selectedJob.result.summary.mean_total_binding_kcal !== undefined && (
                  <Row gutter={16} style={{ marginBottom: 16 }}>
                    <Col span={12}>
                      <Statistic
                        title={`平均 Binding (${unit})`}
                        value={convert(selectedJob.result.summary.mean_total_binding_kcal)}
                        precision={2}
                      />
                    </Col>
                    <Col span={12}>
                      <Space><Text>单位:</Text><Select value={unit} onChange={setUnit} style={{ width: 120 }} options={getUnitOptions()} /></Space>
                    </Col>
                  </Row>
                )}
                {renderResultChart()}
              </>
            )}
          </div>
        )}
      </Modal>
    </Card>
  );
};

export default BindingAnalysisPanel;

