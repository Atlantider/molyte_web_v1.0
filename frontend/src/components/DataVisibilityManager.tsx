/**
 * 数据可见性管理组件
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
  Switch,
  Input,
  message,
  Statistic,
  Row,
  Col,
  Tooltip,
  Space,
  Alert,
  Progress,
} from 'antd';
import {
  EyeOutlined,
  EyeInvisibleOutlined,
  ClockCircleOutlined,
  LockOutlined,
  GiftOutlined,
  TrophyOutlined,
  DownloadOutlined,
  SwapOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import * as visibilityApi from '../api/visibility';
import { DataVisibility, JobVisibility, VisibilityStats, ExchangeRateInfo } from '../api/visibility';

// 可见性标签配置
const visibilityConfig = {
  [DataVisibility.PUBLIC]: { color: 'green', icon: <EyeOutlined />, text: '公开' },
  [DataVisibility.PRIVATE]: { color: 'red', icon: <EyeInvisibleOutlined />, text: '私有' },
  [DataVisibility.DELAYED]: { color: 'orange', icon: <ClockCircleOutlined />, text: '延期公开' },
  [DataVisibility.ADMIN_ONLY]: { color: 'purple', icon: <LockOutlined />, text: '仅管理员' },
};

export default function DataVisibilityManager() {
  const [loading, setLoading] = useState(false);
  const [jobs, setJobs] = useState<JobVisibility[]>([]);
  const [stats, setStats] = useState<VisibilityStats | null>(null);
  const [total, setTotal] = useState(0);
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(10);
  const [filterVisibility, setFilterVisibility] = useState<DataVisibility | undefined>();
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [editingJob, setEditingJob] = useState<JobVisibility | null>(null);
  const [exchangeModalVisible, setExchangeModalVisible] = useState(false);
  const [exchangeRate, setExchangeRate] = useState<ExchangeRateInfo | null>(null);
  const [exchangePoints, setExchangePoints] = useState<number>(0);
  const [exchangeLoading, setExchangeLoading] = useState(false);
  const [form] = Form.useForm();

  // 加载数据
  const loadData = async () => {
    setLoading(true);
    try {
      const [jobsRes, statsRes] = await Promise.all([
        visibilityApi.getMyJobsVisibility(filterVisibility, page, pageSize),
        visibilityApi.getMyVisibilityStats(),
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
      anonymous_public: job.anonymous_public,
      allow_download: job.allow_download,
    });
    setEditModalVisible(true);
  };

  // 保存可见性设置
  const handleSave = async () => {
    if (!editingJob) return;
    try {
      const values = await form.validateFields();
      await visibilityApi.updateJobVisibility(editingJob.id, values);
      message.success('可见性设置已更新');
      setEditModalVisible(false);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新失败');
    }
  };

  // 领取奖励
  const handleClaimReward = async (job: JobVisibility) => {
    try {
      const result = await visibilityApi.claimPublicReward(job.id);
      message.success(result.message);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '领取失败');
    }
  };

  // 打开积分兑换弹窗
  const handleOpenExchange = async () => {
    try {
      const rate = await visibilityApi.getExchangeRate();
      setExchangeRate(rate);
      setExchangePoints(Math.min(rate.current_points, 100)); // 默认兑换100积分或全部
      setExchangeModalVisible(true);
    } catch (error: any) {
      message.error('获取兑换信息失败');
    }
  };

  // 执行积分兑换
  const handleExchange = async () => {
    if (exchangePoints <= 0) {
      message.warning('请输入要兑换的积分数量');
      return;
    }
    setExchangeLoading(true);
    try {
      const result = await visibilityApi.exchangePointsForCpuHours(exchangePoints);
      message.success(result.message);
      setExchangeModalVisible(false);
      loadData(); // 刷新数据
    } catch (error: any) {
      message.error(error.response?.data?.detail || '兑换失败');
    } finally {
      setExchangeLoading(false);
    }
  };

  // 表格列定义
  const columns: ColumnsType<JobVisibility> = [
    {
      title: '任务名称',
      dataIndex: 'name',
      key: 'name',
      ellipsis: true,
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
          <Tooltip title="下载次数">
            <span><DownloadOutlined /> {record.download_count}</span>
          </Tooltip>
        </Space>
      ),
    },
    {
      title: '奖励',
      key: 'reward',
      width: 100,
      render: (_, record) => {
        if (record.visibility !== DataVisibility.PUBLIC) return '-';
        if (record.reward_claimed) {
          return <Tag color="gold" icon={<TrophyOutlined />}>已领取</Tag>;
        }
        return (
          <Button
            type="link"
            size="small"
            icon={<GiftOutlined />}
            onClick={() => handleClaimReward(record)}
          >
            领取 +10核时
          </Button>
        );
      },
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
    <div className="data-visibility-manager">
      {/* 统计卡片 */}
      {stats && (
        <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
          <Col xs={12} sm={8} md={4}>
            <Card size="small">
              <Statistic
                title="公开数据"
                value={stats.public}
                suffix={`/ ${stats.total}`}
                valueStyle={{ color: '#52c41a' }}
                prefix={<EyeOutlined />}
              />
            </Card>
          </Col>
          <Col xs={12} sm={8} md={5}>
            <Card size="small">
              <Statistic
                title="延期公开"
                value={stats.delayed}
                valueStyle={{ color: '#faad14' }}
                prefix={<ClockCircleOutlined />}
              />
            </Card>
          </Col>
          <Col xs={12} sm={8} md={5}>
            <Card size="small">
              <Statistic
                title="私有配额"
                value={stats.private_quota_used}
                suffix={`/ ${stats.private_quota_limit}`}
                valueStyle={{ color: stats.private_quota_used! >= stats.private_quota_limit! ? '#ff4d4f' : '#1890ff' }}
                prefix={<LockOutlined />}
              />
            </Card>
          </Col>
          <Col xs={12} sm={8} md={5}>
            <Card size="small">
              <Statistic
                title="贡献积分"
                value={stats.contribution_points}
                precision={1}
                valueStyle={{ color: '#722ed1' }}
                prefix={<TrophyOutlined />}
              />
              {(stats.contribution_points || 0) >= 10 && (
                <Button
                  type="link"
                  size="small"
                  icon={<SwapOutlined />}
                  onClick={handleOpenExchange}
                  style={{ padding: 0, marginTop: 4 }}
                >
                  兑换核时
                </Button>
              )}
            </Card>
          </Col>
          <Col xs={12} sm={8} md={5}>
            <Card size="small">
              <Statistic
                title="核时余额"
                value={stats.balance_cpu_hours}
                precision={2}
                valueStyle={{ color: '#1890ff' }}
                prefix={<ThunderboltOutlined />}
                suffix="核时"
              />
            </Card>
          </Col>
        </Row>
      )}

      {/* 提示信息 */}
      <Alert
        message="数据公开奖励规则"
        description={
          <ul style={{ margin: 0, paddingLeft: 16 }}>
            <li>将数据设为公开可获得 <strong>10 核时</strong> 奖励</li>
            <li>数据被查看可获得 <strong>0.1 积分/次</strong>（每日上限 5 积分）</li>
            <li>数据被下载可获得 <strong>1 积分/次</strong></li>
            <li>积分可兑换核时：<strong>10 积分 = 1 核时</strong></li>
            <li style={{ color: '#ff4d4f' }}>⚠️ 取消公开将扣除已领取的 10 核时奖励</li>
          </ul>
        }
        type="info"
        showIcon
        style={{ marginBottom: 16 }}
      />

      {/* 筛选和表格 */}
      <Card
        title="我的数据可见性"
        extra={
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
            ]}
          />
        }
      >
        <Table
          columns={columns}
          dataSource={jobs}
          rowKey="id"
          loading={loading}
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

      {/* 编辑弹窗 */}
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
                { value: DataVisibility.PUBLIC, label: '🌐 公开 - 所有人可见，可获得奖励' },
                { value: DataVisibility.DELAYED, label: '⏰ 延期公开 - 在指定日期后自动公开' },
                { value: DataVisibility.PRIVATE, label: '🔒 私有 - 仅自己可见（需要私有配额）' },
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
                  <InputNumber
                    min={1}
                    max={stats?.max_delay_years ? stats.max_delay_years * 365 : 365}
                    addonAfter="天"
                    style={{ width: '100%' }}
                  />
                </Form.Item>
              )
            }
          </Form.Item>

          <Form.Item
            name="anonymous_public"
            label="匿名公开"
            valuePropName="checked"
            tooltip="开启后，公开时不显示您的用户名和单位信息"
          >
            <Switch checkedChildren="匿名" unCheckedChildren="实名" />
          </Form.Item>

          <Form.Item
            name="allow_download"
            label="允许下载"
            valuePropName="checked"
          >
            <Switch checkedChildren="允许" unCheckedChildren="禁止" />
          </Form.Item>

          {editingJob?.is_free_quota && (
            <Alert
              message="此任务使用免费核时"
              description="使用免费核时的任务不能设为永久私有，请选择延期公开或立即公开"
              type="warning"
              showIcon
              style={{ marginTop: 16 }}
            />
          )}

          {editingJob?.reward_claimed && form.getFieldValue('visibility') !== DataVisibility.PUBLIC && (
            <Alert
              message="取消公开将扣除奖励"
              description="您已领取 10 核时公开奖励，取消公开将从余额中扣除该奖励"
              type="error"
              showIcon
              style={{ marginTop: 16 }}
            />
          )}
        </Form>
      </Modal>

      {/* 积分兑换弹窗 */}
      <Modal
        title="积分兑换核时"
        open={exchangeModalVisible}
        onOk={handleExchange}
        onCancel={() => setExchangeModalVisible(false)}
        okText="兑换"
        cancelText="取消"
        confirmLoading={exchangeLoading}
      >
        {exchangeRate && (
          <div>
            <Row gutter={16} style={{ marginBottom: 24 }}>
              <Col span={12}>
                <Statistic
                  title="当前积分"
                  value={exchangeRate.current_points}
                  precision={1}
                  valueStyle={{ color: '#722ed1' }}
                  prefix={<TrophyOutlined />}
                />
              </Col>
              <Col span={12}>
                <Statistic
                  title="当前核时余额"
                  value={exchangeRate.current_balance}
                  precision={2}
                  valueStyle={{ color: '#1890ff' }}
                  prefix={<ThunderboltOutlined />}
                />
              </Col>
            </Row>

            <Alert
              message={`兑换比例：${exchangeRate.description}`}
              type="info"
              showIcon
              style={{ marginBottom: 16 }}
            />

            <Form layout="vertical">
              <Form.Item label="兑换积分数量">
                <InputNumber
                  value={exchangePoints}
                  onChange={(v) => setExchangePoints(v || 0)}
                  min={10}
                  max={exchangeRate.current_points}
                  step={10}
                  style={{ width: '100%' }}
                  addonAfter="积分"
                />
              </Form.Item>

              <div style={{
                padding: 16,
                background: '#f5f5f5',
                borderRadius: 8,
                textAlign: 'center'
              }}>
                <Space size="large">
                  <span style={{ fontSize: 24, color: '#722ed1' }}>
                    {exchangePoints.toFixed(1)} 积分
                  </span>
                  <SwapOutlined style={{ fontSize: 20 }} />
                  <span style={{ fontSize: 24, color: '#52c41a' }}>
                    {(exchangePoints / exchangeRate.ratio).toFixed(2)} 核时
                  </span>
                </Space>
              </div>
            </Form>
          </div>
        )}
      </Modal>
    </div>
  );
}

