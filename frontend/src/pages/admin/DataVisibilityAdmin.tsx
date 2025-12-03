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
  Typography,
  theme,
} from 'antd';
import {
  EyeOutlined,
  EyeInvisibleOutlined,
  ClockCircleOutlined,
  LockOutlined,
  WarningOutlined,
  CheckCircleOutlined,
  UserOutlined,
  DatabaseOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import AdminNav from '../../components/AdminNav';
import * as visibilityApi from '../../api/visibility';
import { DataVisibility, JobVisibility, VisibilityStats } from '../../api/visibility';
import { useThemeStore } from '../../stores/themeStore';

const { Title, Text } = Typography;

// 可见性标签配置
const visibilityConfig = {
  [DataVisibility.PUBLIC]: { color: 'green', icon: <EyeOutlined />, text: '公开' },
  [DataVisibility.PRIVATE]: { color: 'red', icon: <EyeInvisibleOutlined />, text: '私有' },
  [DataVisibility.DELAYED]: { color: 'orange', icon: <ClockCircleOutlined />, text: '延期公开' },
  [DataVisibility.ADMIN_ONLY]: { color: 'purple', icon: <LockOutlined />, text: '仅管理员' },
};

export default function DataVisibilityAdmin() {
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
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
      sorter: (a, b) => a.id - b.id,
      defaultSortOrder: 'descend',
    },
    {
      title: '任务名称',
      dataIndex: 'name',
      key: 'name',
      ellipsis: true,
      sorter: (a, b) => (a.name || '').localeCompare(b.name || ''),
    },
    {
      title: '用户',
      dataIndex: 'username',
      key: 'username',
      width: 120,
      sorter: (a, b) => (a.username || '').localeCompare(b.username || ''),
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
      filters: [
        { text: '公开', value: DataVisibility.PUBLIC },
        { text: '延期公开', value: DataVisibility.DELAYED },
        { text: '私有', value: DataVisibility.PRIVATE },
        { text: '仅管理员', value: DataVisibility.ADMIN_ONLY },
      ],
      onFilter: (value, record) => record.visibility === value,
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
      sorter: (a, b) => {
        if (!a.visibility_delay_until) return 1;
        if (!b.visibility_delay_until) return -1;
        return new Date(a.visibility_delay_until).getTime() - new Date(b.visibility_delay_until).getTime();
      },
      render: (date: string | null, record) => {
        if (record.visibility !== DataVisibility.DELAYED || !date) return '-';
        return new Date(date).toLocaleDateString('zh-CN');
      },
    },
    {
      title: '查看次数',
      dataIndex: 'view_count',
      key: 'view_count',
      width: 100,
      sorter: (a, b) => a.view_count - b.view_count,
      render: (count: number) => (
        <Space>
          <Tooltip title="查看次数">
            <span><EyeOutlined /> {count}</span>
          </Tooltip>
        </Space>
      ),
    },
    {
      title: '免费核时',
      dataIndex: 'is_free_quota',
      key: 'is_free_quota',
      width: 80,
      filters: [
        { text: '免费', value: true },
        { text: '付费', value: false },
      ],
      onFilter: (value, record) => record.is_free_quota === value,
      render: (isFree: boolean) => (
        isFree ? <Tag color="blue">免费</Tag> : <Tag>付费</Tag>
      ),
    },
    {
      title: '操作',
      key: 'action',
      width: 80,
      fixed: 'right' as const,
      render: (_, record) => (
        <Button type="link" size="small" onClick={() => handleEdit(record)}>
          设置
        </Button>
      ),
    },
  ];

  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <DatabaseOutlined style={{ marginRight: 12, color: token.colorPrimary }} />
          数据可见性管理
        </Title>
        <Text type="secondary">
          管理所有用户的数据可见性设置
        </Text>
      </div>

      <AdminNav />

      {/* 统计卡片 */}
      {stats && (
        <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Card
              bordered={false}
              style={{
                borderRadius: 12,
                background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                boxShadow: '0 4px 12px rgba(102, 126, 234, 0.3)',
                height: '100%',
              }}
            >
              <Statistic
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>总计</span>}
                value={stats.total}
                valueStyle={{ color: '#fff', fontSize: 28 }}
                prefix={<CheckCircleOutlined />}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Card
              bordered={false}
              style={{
                borderRadius: 12,
                background: 'linear-gradient(135deg, #11998e 0%, #38ef7d 100%)',
                boxShadow: '0 4px 12px rgba(17, 153, 142, 0.3)',
                height: '100%',
              }}
            >
              <Statistic
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>公开</span>}
                value={stats.public}
                valueStyle={{ color: '#fff', fontSize: 28 }}
                prefix={<EyeOutlined />}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Card
              bordered={false}
              style={{
                borderRadius: 12,
                background: 'linear-gradient(135deg, #fa709a 0%, #fee140 100%)',
                boxShadow: '0 4px 12px rgba(250, 112, 154, 0.3)',
                height: '100%',
              }}
            >
              <Statistic
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>延期公开</span>}
                value={stats.delayed}
                valueStyle={{ color: '#fff', fontSize: 28 }}
                prefix={<ClockCircleOutlined />}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Card
              bordered={false}
              style={{
                borderRadius: 12,
                background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
                boxShadow: '0 4px 12px rgba(240, 147, 251, 0.3)',
                height: '100%',
              }}
            >
              <Statistic
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>私有</span>}
                value={stats.private}
                valueStyle={{ color: '#fff', fontSize: 28 }}
                prefix={<EyeInvisibleOutlined />}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Card
              bordered={false}
              style={{
                borderRadius: 12,
                background: 'linear-gradient(135deg, #30cfd0 0%, #330867 100%)',
                boxShadow: '0 4px 12px rgba(48, 207, 208, 0.3)',
                height: '100%',
              }}
            >
              <Statistic
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>仅管理员</span>}
                value={stats.admin_only}
                valueStyle={{ color: '#fff', fontSize: 28 }}
                prefix={<LockOutlined />}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={8} lg={4}>
            <Card
              bordered={false}
              style={{
                borderRadius: 12,
                background: 'linear-gradient(135deg, #a8edea 0%, #fed6e3 100%)',
                boxShadow: '0 4px 12px rgba(168, 237, 234, 0.3)',
                height: '100%',
              }}
            >
              <Statistic
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>即将公开</span>}
                value={stats.soon_public}
                valueStyle={{ color: '#fff', fontSize: 28 }}
                prefix={<WarningOutlined />}
              />
            </Card>
          </Col>
        </Row>
      )}

      {/* 筛选和表格 */}
      <Card
        title={
          <Space>
            <DatabaseOutlined style={{ color: token.colorPrimary }} />
            <span>数据可见性管理</span>
          </Space>
        }
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
        style={{
          borderRadius: 12,
          border: 'none',
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
          background: token.colorBgContainer,
        }}
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

