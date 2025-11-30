/**
 * 充值中心页面
 */
import React, { useState, useEffect } from 'react';
import {
  Card, Row, Col, Statistic, Button, InputNumber, Radio, Space, Table, Tag, Modal,
  message, Spin, Alert, Divider, Typography, Progress, Tabs
} from 'antd';
import {
  WalletOutlined, ThunderboltOutlined, WarningOutlined, CheckCircleOutlined,
  ClockCircleOutlined, AlipayCircleOutlined, WechatOutlined, CreditCardOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import {
  getBalance, createOrder, getOrders, getTransactions, simulatePayment,
  BalanceInfo, RechargeOrder, Transaction
} from '../api/billing';

const { Title, Text } = Typography;
const { TabPane } = Tabs;

// 预设充值金额
const PRESET_AMOUNTS = [10, 50, 100, 200, 500, 1000];

const Recharge: React.FC = () => {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(true);
  const [balance, setBalance] = useState<BalanceInfo | null>(null);
  const [orders, setOrders] = useState<RechargeOrder[]>([]);
  const [transactions, setTransactions] = useState<Transaction[]>([]);
  
  // 充值表单
  const [amount, setAmount] = useState<number>(100);
  const [paymentMethod, setPaymentMethod] = useState<string>('simulated');
  const [submitting, setSubmitting] = useState(false);
  
  // 支付弹窗
  const [payModalVisible, setPayModalVisible] = useState(false);
  const [currentOrder, setCurrentOrder] = useState<RechargeOrder | null>(null);
  const [paying, setPaying] = useState(false);

  useEffect(() => {
    loadData();
  }, []);

  const loadData = async () => {
    setLoading(true);
    try {
      const [balanceData, ordersData, transData] = await Promise.all([
        getBalance(),
        getOrders(0, 10),
        getTransactions(0, 20)
      ]);
      setBalance(balanceData);
      setOrders(ordersData);
      setTransactions(transData);
    } catch (error) {
      console.error('加载数据失败:', error);
      message.error('加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  const handleCreateOrder = async () => {
    if (!amount || amount < 10) {
      message.error('最低充值金额为 ¥10');
      return;
    }
    
    setSubmitting(true);
    try {
      const order = await createOrder({ amount, payment_method: paymentMethod });
      message.success('订单创建成功');
      setCurrentOrder(order);
      setPayModalVisible(true);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '创建订单失败');
    } finally {
      setSubmitting(false);
    }
  };

  const handleSimulatePay = async () => {
    if (!currentOrder) return;
    
    setPaying(true);
    try {
      const result = await simulatePayment(currentOrder.id);
      message.success(result.message);
      setPayModalVisible(false);
      setCurrentOrder(null);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '支付失败');
    } finally {
      setPaying(false);
    }
  };

  const orderColumns = [
    { title: '订单号', dataIndex: 'order_no', key: 'order_no', width: 200 },
    { title: '金额', dataIndex: 'amount', key: 'amount', render: (v: number) => `¥${v.toFixed(2)}` },
    { title: '机时', dataIndex: 'cpu_hours', key: 'cpu_hours', render: (v: number) => `${v.toFixed(2)} 核时` },
    {
      title: '状态', dataIndex: 'payment_status', key: 'payment_status',
      render: (status: string) => {
        const config: Record<string, { color: string; text: string }> = {
          pending: { color: 'orange', text: '待支付' },
          paid: { color: 'green', text: '已支付' },
          failed: { color: 'red', text: '失败' },
          cancelled: { color: 'default', text: '已取消' },
        };
        const c = config[status] || { color: 'default', text: status };
        return <Tag color={c.color}>{c.text}</Tag>;
      }
    },
    { title: '创建时间', dataIndex: 'created_at', key: 'created_at', render: (v: string) => new Date(v).toLocaleString() },
    {
      title: '操作', key: 'action',
      render: (_: any, record: RechargeOrder) => (
        record.payment_status === 'pending' && (
          <Button type="link" size="small" onClick={() => { setCurrentOrder(record); setPayModalVisible(true); }}>
            去支付
          </Button>
        )
      )
    }
  ];

  const transColumns = [
    {
      title: '类型', dataIndex: 'type', key: 'type',
      render: (type: string) => {
        const config: Record<string, { color: string; text: string }> = {
          recharge: { color: 'green', text: '充值' },
          consume: { color: 'red', text: '消费' },
          refund: { color: 'blue', text: '退款' },
          admin_adjust: { color: 'purple', text: '管理员调整' },
          debt_repay: { color: 'orange', text: '还债' },
        };
        const c = config[type] || { color: 'default', text: type };
        return <Tag color={c.color}>{c.text}</Tag>;
      }
    },
    {
      title: '变动', dataIndex: 'amount', key: 'amount',
      render: (v: number) => (
        <Text type={v >= 0 ? 'success' : 'danger'}>
          {v >= 0 ? '+' : ''}{v.toFixed(2)} 核时
        </Text>
      )
    },
    { title: '余额', dataIndex: 'balance_after', key: 'balance_after', render: (v: number) => `${v.toFixed(2)} 核时` },
    { title: '说明', dataIndex: 'description', key: 'description', ellipsis: true },
    { title: '时间', dataIndex: 'created_at', key: 'created_at', render: (v: string) => new Date(v).toLocaleString() },
  ];

  if (loading) {
    return (
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        height: 'calc(100vh - 64px)',
        background: '#f5f7fb',
      }}>
        <Spin size="large" />
      </div>
    );
  }

  return (
    <div style={{
      padding: 24,
      background: '#f5f7fb',
      minHeight: 'calc(100vh - 64px)'
    }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <WalletOutlined style={{ marginRight: 12, color: '#1677ff' }} />
          充值中心
        </Title>
        <Text type="secondary">管理您的机时余额和充值记录</Text>
      </div>

      {/* 余额概览 */}
      <Row gutter={16} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} lg={6}>
          <Card style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none'
          }}>
            <Statistic title="可用余额" value={balance?.available || 0} suffix="核时" precision={2}
              valueStyle={{ color: (balance?.available || 0) > 10 ? '#52c41a' : '#ff4d4f' }}
              prefix={<ThunderboltOutlined />} />
          </Card>
        </Col>
        <Col xs={24} sm={12} lg={6}>
          <Card style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none'
          }}>
            <Statistic title="冻结中" value={balance?.frozen || 0} suffix="核时" precision={2}
              prefix={<ClockCircleOutlined />} />
          </Card>
        </Col>
        <Col xs={24} sm={12} lg={6}>
          <Card style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none'
          }}>
            <Statistic title="欠费" value={balance?.debt || 0} suffix="核时" precision={2}
              valueStyle={{ color: (balance?.debt || 0) > 0 ? '#ff4d4f' : '#52c41a' }}
              prefix={<WarningOutlined />} />
          </Card>
        </Col>
        <Col xs={24} sm={12} lg={6}>
          <Card style={{
            borderRadius: 12,
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            border: 'none'
          }}>
            <Statistic title="当前单价" value={balance?.price_per_hour || 0.1} prefix="¥" suffix="/核时" precision={2} />
          </Card>
        </Col>
      </Row>

      {/* 欠费警告 */}
      {(balance?.debt || 0) > 0 && (
        <Alert
          type="error"
          showIcon
          icon={<WarningOutlined />}
          message="您有欠费未还清"
          description={`当前欠费 ${balance?.debt?.toFixed(2)} 核时，请尽快充值。欠费期间无法提交新任务，已完成任务的结果将被锁定。`}
          style={{ marginBottom: 24, borderRadius: 8 }}
        />
      )}

      <Row gutter={24}>
        {/* 充值表单 */}
        <Col xs={24} lg={12}>
          <Card
            title="立即充值"
            extra={<Text type="secondary">1元 = {(1 / (balance?.price_per_hour || 0.1)).toFixed(0)} 核时</Text>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <div style={{ marginBottom: 16 }}>
              <Text strong>选择金额：</Text>
              <div style={{ marginTop: 8 }}>
                <Radio.Group value={amount} onChange={e => setAmount(e.target.value)}>
                  <Space wrap>
                    {PRESET_AMOUNTS.map(a => (
                      <Radio.Button key={a} value={a}>¥{a}</Radio.Button>
                    ))}
                  </Space>
                </Radio.Group>
              </div>
            </div>

            <div style={{ marginBottom: 16 }}>
              <Text strong>自定义金额：</Text>
              <InputNumber
                style={{ width: 200, marginLeft: 8 }}
                min={10}
                max={10000}
                value={amount}
                onChange={v => setAmount(v || 10)}
                prefix="¥"
                precision={0}
              />
            </div>

            <div style={{ marginBottom: 24 }}>
              <Text strong>支付方式：</Text>
              <div style={{ marginTop: 8 }}>
                <Radio.Group value={paymentMethod} onChange={e => setPaymentMethod(e.target.value)}>
                  <Space>
                    <Radio.Button value="simulated"><CreditCardOutlined /> 模拟支付</Radio.Button>
                    <Radio.Button value="wechat" disabled><WechatOutlined /> 微信支付</Radio.Button>
                    <Radio.Button value="alipay" disabled><AlipayCircleOutlined /> 支付宝</Radio.Button>
                  </Space>
                </Radio.Group>
              </div>
              <Text type="secondary" style={{ fontSize: 12, marginTop: 4, display: 'block' }}>
                * 微信/支付宝支付即将上线
              </Text>
            </div>

            <Divider />

            <div style={{ textAlign: 'center' }}>
              <div style={{ marginBottom: 16 }}>
                <Text>充值金额：</Text>
                <Text strong style={{ fontSize: 24, color: '#1677ff' }}>¥{amount}</Text>
                <Text style={{ marginLeft: 16 }}>可获得：</Text>
                <Text strong style={{ fontSize: 24, color: '#52c41a' }}>
                  {(amount / (balance?.price_per_hour || 0.1)).toFixed(2)} 核时
                </Text>
              </div>
              <Button
                type="primary"
                size="large"
                loading={submitting}
                onClick={handleCreateOrder}
                style={{
                  borderRadius: 8,
                  boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)',
                }}
              >
                立即充值
              </Button>
            </div>
          </Card>
        </Col>

        {/* 充值记录 */}
        <Col xs={24} lg={12}>
          <Card
            title="最近订单"
            extra={<Button type="link" size="small">查看全部</Button>}
            style={{
              borderRadius: 12,
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none'
            }}
          >
            <Table
              dataSource={orders}
              columns={orderColumns}
              rowKey="id"
              size="small"
              pagination={false}
              scroll={{ y: 300 }}
            />
          </Card>
        </Col>
      </Row>

      {/* 消费记录 */}
      <Card
        title="配额变更记录"
        style={{
          marginTop: 24,
          borderRadius: 12,
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none'
        }}
      >
        <Table
          dataSource={transactions}
          columns={transColumns}
          rowKey="id"
          size="small"
          pagination={{ pageSize: 10 }}
        />
      </Card>

      {/* 支付弹窗 */}
      <Modal
        title="确认支付"
        open={payModalVisible}
        onCancel={() => { setPayModalVisible(false); setCurrentOrder(null); }}
        footer={null}
        width={400}
      >
        {currentOrder && (
          <div style={{ textAlign: 'center', padding: 20 }}>
            <div style={{ marginBottom: 24 }}>
              <Text>订单号：{currentOrder.order_no}</Text>
            </div>
            <div style={{ marginBottom: 24 }}>
              <Text style={{ fontSize: 16 }}>支付金额</Text>
              <div style={{ fontSize: 36, fontWeight: 'bold', color: '#1890ff' }}>
                ¥{currentOrder.amount.toFixed(2)}
              </div>
              <Text type="secondary">可获得 {currentOrder.cpu_hours.toFixed(2)} 核时</Text>
            </div>

            {currentOrder.payment_method === 'simulated' && (
              <Alert
                type="info"
                message="模拟支付"
                description="点击下方按钮模拟完成支付，机时将立即到账"
                style={{ marginBottom: 24, textAlign: 'left' }}
              />
            )}

            <Button
              type="primary"
              size="large"
              block
              loading={paying}
              onClick={handleSimulatePay}
              icon={<CheckCircleOutlined />}
            >
              确认支付
            </Button>
          </div>
        )}
      </Modal>
    </div>
  );
};

export default Recharge;

