/**
 * 管理员数据可见性管理页面
 */
import { useState, useEffect } from 'react';
import {
  Card,
  Table,
  Tag,
  Button,
  Modal,
  Form,
  Select,
  InputNumber,
  Input,
  message,
  Statistic,
  Row,
  Col,
  Space,
  Tooltip,
  Popconfirm,
} from 'antd';
import {
  EyeOutlined,
  EyeInvisibleOutlined,
  ClockCircleOutlined,
  LockOutlined,
  WarningOutlined,
  CheckCircleOutlined,
  UserOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import * as visibilityApi from '../../api/visibility';
import { DataVisibility, JobVisibility, VisibilityStats } from '../../api/visibility';

// 可见性标签配置
const visibilityConfig = {
  [DataVisibility.PUBLIC]: { color: 'green', icon: <EyeOutlined />, text: '公开' },
  [DataVisibility.PRIVATE]: { color: 'red', icon: <EyeInvisibleOutlined />, text: '私有' },
  [DataVisibility.DELAYED]: { color: 'orange', icon: <ClockCircleOutlined />, text: '延期公开' },
  [DataVisibility.ADMIN_ONLY]: { color: 'purple', icon: <LockOutlined />, text: '仅管理员' },
};

export default function DataVisibilityAdmin() {
  const [loading, setLoading] = useState(false);
  const [jobs, setJobs] = useState<JobVisibility[]>([]);
  const [stats, setStats] = useState<VisibilityStats | null>(null);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(20);
  const [filterVisibility, setFilterVisibility] = useState<DataVisibility | undefined>();
  const [selectedRowKeys, setSelectedRowKeys] = useState<number[]>([]);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [editingJob, setEditingJob] = useState<JobVisibility | null>(null);
  const [batchModalVisible, setBatchModalVisible] = useState(false);
  const [form] = Form.useForm();
  const [batchForm] = Form.useForm();

  // 加载数据
  const loadData = async () => {
    setLoading(true);
    try {
      const [jobsRes, statsRes] = await Promise.all([
        visibilityApi.adminGetAllJobsVisibility(filterVisibility, undefined, page, pageSize),
        visibilityApi.adminGetVisibilityStats(),
      ]);
      setJobs(jobsRes.items);
      setTotal(jobsRes.total);
      setStats(statsRes);
    } catch (error) {
      message.error('加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadData();
  }, [page, pageSize, filterVisibility]);

  // 打开编辑弹窗
  const handleEdit = (job: JobVisibility) => {
    setEditingJob(job);
    form.setFieldsValue({
      visibility: job.visibility,
      delay_days: 365,
      reason: '',
    });
    setEditModalVisible(true);
  };

  // 保存可见性设置
  const handleSave = async () => {
    if (!editingJob) return;
    try {
      const values = await form.validateFields();
      await visibilityApi.adminUpdateJobVisibility(editingJob.id, values);
      message.success('可见性设置已更新');
      setEditModalVisible(false);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新失败');
    }
  };

  // 批量更新
  const handleBatchUpdate = async () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择要更新的任务');
      return;
    }
    try {
      const values = await batchForm.validateFields();
      await visibilityApi.adminBatchUpdateVisibility(
        selectedRowKeys,
        values.visibility,
        values.delay_days,
        values.reason
      );
      message.success(`已更新 ${selectedRowKeys.length} 个任务`);
      setBatchModalVisible(false);
      setSelectedRowKeys([]);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '批量更新失败');
    }
  };

  // 表格列定义
  const columns: ColumnsType<JobVisibility> = [
    {
      title: '任务ID',
      dataIndex: 'id',
      key: 'id',
      width: 80,
    },
    {
      title: '任务名称',
      dataIndex: 'name',
      key: 'name',
      ellipsis: true,
    },
    {
      title: '用户',
      dataIndex: 'username',
      key: 'username',
      width: 120,
      render: (username: string) => (
        <Space>
          <UserOutlined />
          {username}
        </Space>
      ),
    },
    {
      title: '可见性',
      dataIndex: 'visibility',
      key: 'visibility',
      width: 120,
      render: (visibility: DataVisibility) => {
        const config = visibilityConfig[visibility];
        return (
          <Tag color={config.color} icon={config.icon}>
            {config.text}
          </Tag>
        );
      },
    },
    {
      title: '延期至',
      dataIndex: 'visibility_delay_until',
      key: 'visibility_delay_until',
      width: 120,
      render: (date: string | null, record) => {
        if (record.visibility !== DataVisibility.DELAYED || !date) return '-';
        return new Date(date).toLocaleDateString('zh-CN');
      },
    },
    {
      title: '查看/下载',
      key: 'counts',
      width: 100,
      render: (_, record) => (
        <Space>
          <Tooltip title="查看次数">
            <span><EyeOutlined /> {record.view_count}</span>
          </Tooltip>
        </Space>
      ),
    },
    {
      title: '免费核时',
      dataIndex: 'is_free_quota',
      key: 'is_free_quota',
      width: 80,
      render: (isFree: boolean) => (
        isFree ? <Tag color="blue">免费</Tag> : <Tag>付费</Tag>
      ),
    },
    {
      title: '操作',
      key: 'action',
      width: 80,
      render: (_, record) => (
        <Button type="link" size="small" onClick={() => handleEdit(record)}>
          设置
        </Button>
      ),
    },
  ];

  return (
    <div className="data-visibility-admin">
      {/* 统计卡片 */}
      {stats && (
        <Row gutter={16} style={{ marginBottom: 24 }}>
          <Col span={4}>
            <Card size="small">
              <Statistic
                title="总计"
                value={stats.total}
                prefix={<CheckCircleOutlined />}
              />
            </Card>
          </Col>
          <Col span={4}>
            <Card size="small">
              <Statistic
                title="公开"
                value={stats.public}
                valueStyle={{ color: '#52c41a' }}
                prefix={<EyeOutlined />}
              />
            </Card>
          </Col>
          <Col span={4}>
            <Card size="small">
              <Statistic
                title="延期公开"
                value={stats.delayed}
                valueStyle={{ color: '#faad14' }}
                prefix={<ClockCircleOutlined />}
              />
            </Card>
          </Col>
          <Col span={4}>
            <Card size="small">
              <Statistic
                title="私有"
                value={stats.private}
                valueStyle={{ color: '#ff4d4f' }}
                prefix={<EyeInvisibleOutlined />}
              />
            </Card>
          </Col>
          <Col span={4}>
            <Card size="small">
              <Statistic
                title="仅管理员"
                value={stats.admin_only}
                valueStyle={{ color: '#722ed1' }}
                prefix={<LockOutlined />}
              />
            </Card>
          </Col>
          <Col span={4}>
            <Card size="small">
              <Statistic
                title="即将公开(30天内)"
                value={stats.soon_public}
                valueStyle={{ color: '#eb2f96' }}
                prefix={<WarningOutlined />}
              />
            </Card>
          </Col>
        </Row>
      )}

      {/* 筛选和表格 */}
      <Card
        title="数据可见性管理"
        extra={
          <Space>
            <Select
              placeholder="筛选可见性"
              allowClear
              style={{ width: 150 }}
              value={filterVisibility}
              onChange={setFilterVisibility}
              options={[
                { value: DataVisibility.PUBLIC, label: '公开' },
                { value: DataVisibility.DELAYED, label: '延期公开' },
                { value: DataVisibility.PRIVATE, label: '私有' },
                { value: DataVisibility.ADMIN_ONLY, label: '仅管理员' },
              ]}
            />
            <Button
              type="primary"
              disabled={selectedRowKeys.length === 0}
              onClick={() => setBatchModalVisible(true)}
            >
              批量设置 ({selectedRowKeys.length})
            </Button>
          </Space>
        }
      >
        <Table
          columns={columns}
          dataSource={jobs}
          rowKey="id"
          loading={loading}
          rowSelection={{
            selectedRowKeys,
            onChange: (keys) => setSelectedRowKeys(keys as number[]),
          }}
          pagination={{
            current: page,
            pageSize,
            total,
            showSizeChanger: true,
            showTotal: (t) => `共 ${t} 条`,
            onChange: (p, ps) => {
              setPage(p);
              setPageSize(ps);
            },
          }}
        />
      </Card>

      {/* 单个编辑弹窗 */}
      <Modal
        title="设置数据可见性"
        open={editModalVisible}
        onOk={handleSave}
        onCancel={() => setEditModalVisible(false)}
        okText="保存"
        cancelText="取消"
      >
        <Form form={form} layout="vertical">
          <Form.Item
            name="visibility"
            label="可见性"
            rules={[{ required: true, message: '请选择可见性' }]}
          >
            <Select
              options={[
                { value: DataVisibility.PUBLIC, label: '🌐 公开' },
                { value: DataVisibility.DELAYED, label: '⏰ 延期公开' },
                { value: DataVisibility.PRIVATE, label: '🔒 私有' },
                { value: DataVisibility.ADMIN_ONLY, label: '👑 仅管理员' },
              ]}
            />
          </Form.Item>

          <Form.Item
            noStyle
            shouldUpdate={(prev, curr) => prev.visibility !== curr.visibility}
          >
            {({ getFieldValue }) =>
              getFieldValue('visibility') === DataVisibility.DELAYED && (
                <Form.Item
                  name="delay_days"
                  label="延期天数"
                  rules={[{ required: true, message: '请输入延期天数' }]}
                >
                  <InputNumber min={1} max={1095} addonAfter="天" style={{ width: '100%' }} />
                </Form.Item>
              )
            }
          </Form.Item>

          <Form.Item name="reason" label="修改原因">
            <Input.TextArea rows={3} placeholder="请输入修改原因（可选）" />
          </Form.Item>
        </Form>
      </Modal>

      {/* 批量编辑弹窗 */}
      <Modal
        title={`批量设置可见性 (${selectedRowKeys.length} 个任务)`}
        open={batchModalVisible}
        onOk={handleBatchUpdate}
        onCancel={() => setBatchModalVisible(false)}
        okText="确认更新"
        cancelText="取消"
      >
        <Form form={batchForm} layout="vertical">
          <Form.Item
            name="visibility"
            label="可见性"
            rules={[{ required: true, message: '请选择可见性' }]}
          >
            <Select
              options={[
                { value: DataVisibility.PUBLIC, label: '🌐 公开' },
                { value: DataVisibility.DELAYED, label: '⏰ 延期公开' },
                { value: DataVisibility.PRIVATE, label: '🔒 私有' },
                { value: DataVisibility.ADMIN_ONLY, label: '👑 仅管理员' },
              ]}
            />
          </Form.Item>

          <Form.Item
            noStyle
            shouldUpdate={(prev, curr) => prev.visibility !== curr.visibility}
          >
            {({ getFieldValue }) =>
              getFieldValue('visibility') === DataVisibility.DELAYED && (
                <Form.Item
                  name="delay_days"
                  label="延期天数"
                  rules={[{ required: true, message: '请输入延期天数' }]}
                >
                  <InputNumber min={1} max={1095} addonAfter="天" style={{ width: '100%' }} />
                </Form.Item>
              )
            }
          </Form.Item>

          <Form.Item name="reason" label="修改原因">
            <Input.TextArea rows={3} placeholder="请输入修改原因（可选）" />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  );
}

