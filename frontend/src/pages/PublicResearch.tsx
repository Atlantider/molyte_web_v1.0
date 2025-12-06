import React, { useState, useEffect } from 'react';
import { useNavigate, useLocation } from 'react-router-dom';
import {
  Card,
  Form,
  Select,
  Input,
  InputNumber,
  Button,
  Table,
  Space,
  Tag,
  Typography,
  Row,
  Col,
  message,
  Tooltip,
  Empty,
  Tabs,
  Statistic,
  theme,
} from 'antd';
import {
  SearchOutlined,
  ExperimentOutlined,
  EyeOutlined,
  CheckCircleOutlined,
  LoginOutlined,
  ArrowLeftOutlined,
  DatabaseOutlined,
  ThunderboltOutlined,
  FireOutlined,
  FileSearchOutlined,
  ReloadOutlined,
  RocketOutlined,
} from '@ant-design/icons';
import { searchElectrolytes, ElectrolyteSearchResult } from '../api/research';
import { getQCJobs } from '../api/qc';
import { useAuthStore } from '../stores/authStore';
import QCDataTab from '../components/QCDataTab';
import { useThemeStore } from '../stores/themeStore';

const { Title, Text } = Typography;

export default function PublicResearch() {
  const [form] = Form.useForm();
  const navigate = useNavigate();
  const location = useLocation();
  const { isAuthenticated } = useAuthStore();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<ElectrolyteSearchResult[]>([]);
  const [total, setTotal] = useState(0);
  const [pagination, setPagination] = useState({
    current: 1,
    pageSize: 10,
  });

  // 标签页状态
  const [activeTab, setActiveTab] = useState<'md' | 'qc'>('md');

  // 统计数据
  const [mdStats, setMdStats] = useState({ total: 0, completed: 0 });
  const [qcStats, setQcStats] = useState({ total: 0, completed: 0 });

  // 加载统计数据
  useEffect(() => {
    if (isAuthenticated) {
      loadStats();
    }
  }, [isAuthenticated]);

  const loadStats = async () => {
    try {
      // 加载公开 MD 统计
      const mdResponse = await searchElectrolytes({ skip: 0, limit: 1 });
      setMdStats({
        total: mdResponse.total,
        completed: mdResponse.total,
      });

      // 加载公开 QC 统计
      const qcResponse = await getQCJobs({ skip: 0, limit: 1, visibility: 'PUBLIC', status: 'COMPLETED' });
      setQcStats({
        total: qcResponse.total,
        completed: qcResponse.total,
      });
    } catch (error) {
      console.error('Failed to load stats:', error);
    }
  };

  // 常用的阳离子、阴离子、溶剂选项
  const cationOptions = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Zn'];
  const anionOptions = ['FSI', 'TFSI', 'PF6', 'BF4', 'ClO4', 'DCA'];
  const solventOptions = [
    { label: 'EC (碳酸乙烯酯)', value: 'EC' },
    { label: 'DMC (碳酸二甲酯)', value: 'DMC' },
    { label: 'EMC (碳酸甲乙酯)', value: 'EMC' },
    { label: 'DEC (碳酸二乙酯)', value: 'DEC' },
    { label: 'PC (碳酸丙烯酯)', value: 'PC' },
    { label: 'FEC (氟代碳酸乙烯酯)', value: 'FEC' },
  ];

  // 搜索处理
  const handleSearch = async (values: any) => {
    if (!isAuthenticated) {
      message.warning('请先登录后再搜索');
      // 保存当前路径，登录后返回
      navigate('/login', { state: { from: location.pathname } });
      return;
    }

    setLoading(true);
    try {
      const params = {
        cations: values.cations,
        anions: values.anions,
        solvents: values.solvents,
        solvent_smiles: values.solvent_smiles,
        temp_min: values.temp_min,
        temp_max: values.temp_max,
        skip: (pagination.current - 1) * pagination.pageSize,
        limit: pagination.pageSize,
      };

      const response = await searchElectrolytes(params);
      setResults(response.data);
      setTotal(response.total);
      message.success(`找到 ${response.total} 个匹配的结果`);
    } catch (error: any) {
      if (error.response?.status === 401) {
        message.warning('请先登录后再搜索');
        navigate('/login', { state: { from: location.pathname } });
      } else {
        message.error(error.response?.data?.detail || '搜索失败');
      }
    } finally {
      setLoading(false);
    }
  };

  // 重置表单
  const handleReset = () => {
    form.resetFields();
    setResults([]);
    setTotal(0);
  };

  // 查看详情 - 跳转到公开结果页面
  const handleViewDetail = (record: ElectrolyteSearchResult) => {
    if (!isAuthenticated) {
      message.warning('请先登录后查看详情');
      navigate('/login', { state: { from: `/research/result/${record.job_id}` } });
      return;
    }
    navigate(`/research/result/${record.job_id}`);
  };

  // 表格列定义
  const columns = [
    {
      title: '任务 ID',
      dataIndex: 'job_id',
      key: 'job_id',
      width: 100,
      fixed: 'left' as const,
      render: (id: number) => <Text strong>#{id}</Text>,
    },
    {
      title: '配方名称',
      dataIndex: 'system_name',
      key: 'system_name',
      width: 200,
      fixed: 'left' as const,
    },
    {
      title: '阳离子',
      dataIndex: 'cations',
      key: 'cations',
      width: 150,
      render: (cations: any[]) => (
        <Space size={[0, 4]} wrap>
          {cations?.map((c, i) => (
            <Tag key={i} color="red">
              {c.name} ({c.number})
            </Tag>
          ))}
        </Space>
      ),
    },
    {
      title: '阴离子',
      dataIndex: 'anions',
      key: 'anions',
      width: 150,
      render: (anions: any[]) => (
        <Space size={[0, 4]} wrap>
          {anions?.map((a, i) => (
            <Tag key={i} color="blue">
              {a.name} ({a.number})
            </Tag>
          ))}
        </Space>
      ),
    },
    {
      title: '溶剂',
      dataIndex: 'solvents',
      key: 'solvents',
      width: 150,
      render: (solvents: any[]) => (
        <Space size={[0, 4]} wrap>
          {solvents?.map((s, i) => (
            <Tag key={i} color="green">
              {s.name} ({s.number})
            </Tag>
          ))}
        </Space>
      ),
    },
    {
      title: '温度 (K)',
      dataIndex: 'temperature',
      key: 'temperature',
      width: 100,
      render: (temp: number) => temp?.toFixed(1) || '-',
    },
    {
      title: '分析结果',
      key: 'analysis',
      width: 150,
      render: (_: any, record: ElectrolyteSearchResult) => (
        <Space size={[0, 4]} wrap>
          {record.has_rdf && <Tag color="success">RDF</Tag>}
          {record.has_msd && <Tag color="processing">MSD</Tag>}
          {record.has_solvation && <Tag color="warning">溶剂化</Tag>}
        </Space>
      ),
    },
    {
      title: '操作',
      key: 'action',
      width: 100,
      fixed: 'right' as const,
      render: (_: any, record: ElectrolyteSearchResult) => (
        <Button
          type="link"
          icon={<EyeOutlined />}
          onClick={() => handleViewDetail(record)}
        >
          查看详情
        </Button>
      ),
    },
  ];

  return (
    <div style={{ minHeight: '100vh', background: token.colorBgLayout, transition: 'background 0.3s' }}>
      {/* 顶部导航 */}
      <div style={{
        background: 'linear-gradient(135deg, #1a1f36 0%, #0d1025 100%)',
        padding: '0 24px',
        height: 64,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
        boxShadow: '0 2px 12px rgba(0,0,0,0.15)',
        position: 'sticky',
        top: 0,
        zIndex: 100,
      }}>
        <Space size={20}>
          <Button
            icon={<ArrowLeftOutlined />}
            onClick={() => navigate('/')}
            style={{
              borderRadius: 8,
              border: '1px solid rgba(255,255,255,0.2)',
              background: 'rgba(255,255,255,0.1)',
              color: '#fff',
            }}
          >
            返回首页
          </Button>
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: 12,
            padding: '8px 16px',
            background: 'rgba(255,255,255,0.08)',
            borderRadius: 10,
            border: '1px solid rgba(255,255,255,0.1)',
          }}>
            <div style={{
              width: 36,
              height: 36,
              borderRadius: 10,
              background: 'linear-gradient(135deg, #1890ff 0%, #722ed1 100%)',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              boxShadow: '0 4px 12px rgba(24, 144, 255, 0.4)',
            }}>
              <DatabaseOutlined style={{ color: '#fff', fontSize: 18 }} />
            </div>
            <span style={{
              fontSize: 18,
              fontWeight: 600,
              color: '#fff',
              letterSpacing: 1,
            }}>
              电解液研发数据库
            </span>
          </div>
        </Space>
        <Space size={12}>
          {isAuthenticated ? (
            <Button
              type="primary"
              icon={<RocketOutlined />}
              onClick={() => navigate('/workspace/dashboard')}
              style={{
                borderRadius: 8,
                background: 'linear-gradient(135deg, #1890ff 0%, #722ed1 100%)',
                border: 'none',
                boxShadow: '0 4px 12px rgba(24, 144, 255, 0.4)',
              }}
            >
              进入工作台
            </Button>
          ) : (
            <Button
              type="primary"
              icon={<LoginOutlined />}
              onClick={() => navigate('/login', { state: { from: location.pathname } })}
              style={{
                borderRadius: 8,
                background: 'linear-gradient(135deg, #1890ff 0%, #722ed1 100%)',
                border: 'none',
                boxShadow: '0 4px 12px rgba(24, 144, 255, 0.4)',
              }}
            >
              登录
            </Button>
          )}
        </Space>
      </div>

      {/* 主内容 */}
      <div style={{ padding: 24, maxWidth: 1400, margin: '0 auto' }}>
        {/* 页面介绍 */}
        <Card
          style={{
            marginBottom: 24,
            borderRadius: 12,
            border: 'none',
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
            background: 'linear-gradient(135deg, rgba(24, 144, 255, 0.05) 0%, rgba(114, 46, 209, 0.05) 100%)',
          }}
          styles={{ body: { padding: '20px 24px' } }}
        >
          <Row align="middle" gutter={24}>
            <Col>
              <div style={{
                width: 56,
                height: 56,
                borderRadius: 14,
                background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                boxShadow: '0 6px 16px rgba(102, 126, 234, 0.35)',
              }}>
                <ExperimentOutlined style={{ fontSize: 28, color: '#fff' }} />
              </div>
            </Col>
            <Col flex={1}>
              <Title level={4} style={{ margin: 0, marginBottom: 4 }}>
                公开电解液研发数据库
              </Title>
              <Text type="secondary">
                搜索和浏览已完成的分子动力学模拟和量子化学计算结果
                {!isAuthenticated && <Tag color="warning" style={{ marginLeft: 8 }}>需登录后搜索</Tag>}
              </Text>
            </Col>
          </Row>
        </Card>

        {/* 统计卡片 */}
        <Card
          style={{
            marginBottom: 24,
            borderRadius: 12,
            border: 'none',
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          }}
          styles={{ body: { padding: '24px' } }}
        >
          <Row gutter={24}>
            <Col xs={24} sm={12} md={6}>
              <Card
                bordered={false}
                style={{
                  background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                  borderRadius: 12,
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>MD数据</span>}
                  value={mdStats.total}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<DatabaseOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card
                bordered={false}
                style={{
                  background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
                  borderRadius: 12,
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>QC数据</span>}
                  value={qcStats.total}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<ExperimentOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card
                bordered={false}
                style={{
                  background: 'linear-gradient(135deg, #11998e 0%, #38ef7d 100%)',
                  borderRadius: 12,
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>已完成</span>}
                  value={mdStats.completed + qcStats.completed}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<CheckCircleOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card
                bordered={false}
                style={{
                  background: 'linear-gradient(135deg, #fa709a 0%, #fee140 100%)',
                  borderRadius: 12,
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>公开数据</span>}
                  value={mdStats.total + qcStats.total}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<FileSearchOutlined />}
                />
              </Card>
            </Col>
          </Row>
        </Card>

        {/* 标签页 */}
        <Card
          style={{
            borderRadius: 12,
            border: 'none',
            boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          }}
          styles={{ body: { padding: 0 } }}
        >
          <Tabs
            activeKey={activeTab}
            onChange={(key) => setActiveTab(key as 'md' | 'qc')}
            size="large"
            style={{ padding: '0 24px' }}
            items={[
              {
                key: 'md',
                label: (
                  <Space>
                    <DatabaseOutlined />
                    <span>MD数据库</span>
                  </Space>
                ),
                children: (
                  <div style={{ padding: '0 0 24px 0' }}>
                    {/* MD 搜索表单 */}
                    <Card
                      title={
                        <Space>
                          <SearchOutlined style={{ color: '#1677ff' }} />
                          <span>搜索条件</span>
                        </Space>
                      }
                      style={{
                        marginBottom: 24,
                        borderRadius: 12,
                        border: 'none',
                        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                      }}
                    >
          <Form form={form} onFinish={handleSearch} layout="vertical">
                        <Row gutter={16}>
                          <Col xs={24} sm={12} md={8}>
                            <Form.Item label="阳离子" name="cations" tooltip="可多选">
                              <Select
                                mode="multiple"
                                placeholder="选择阳离子"
                                options={cationOptions.map(c => ({ label: c, value: c }))}
                                allowClear
                              />
                            </Form.Item>
                          </Col>
                          <Col xs={24} sm={12} md={8}>
                            <Form.Item label="阴离子" name="anions" tooltip="可多选">
                              <Select
                                mode="multiple"
                                placeholder="选择阴离子"
                                options={anionOptions.map(a => ({ label: a, value: a }))}
                                allowClear
                              />
                            </Form.Item>
                          </Col>
                          <Col xs={24} sm={12} md={8}>
                            <Form.Item label="溶剂" name="solvents" tooltip="按名称选择，可多选">
                              <Select
                                mode="multiple"
                                placeholder="选择溶剂"
                                options={solventOptions}
                                allowClear
                              />
                            </Form.Item>
                          </Col>
                          <Col xs={24} sm={12} md={12}>
                            <Form.Item
                              label={
                                <Tooltip title="SMILES 是分子的唯一标识符，例如 EC 的 SMILES 为 C1COC(=O)O1">
                                  溶剂 SMILES
                                </Tooltip>
                              }
                              name="solvent_smiles"
                            >
                              <Input placeholder="例如: C1COC(=O)O1 (EC 的 SMILES)" allowClear />
                            </Form.Item>
                          </Col>
                          <Col xs={24} sm={12} md={6}>
                            <Form.Item label="最低温度 (K)" name="temp_min">
                              <InputNumber placeholder="273" style={{ width: '100%' }} min={0} />
                            </Form.Item>
                          </Col>
                          <Col xs={24} sm={12} md={6}>
                            <Form.Item label="最高温度 (K)" name="temp_max">
                              <InputNumber placeholder="373" style={{ width: '100%' }} min={0} />
                            </Form.Item>
                          </Col>
                          <Col xs={24} style={{ display: 'flex', alignItems: 'flex-end', marginTop: 8 }}>
                            <Form.Item style={{ marginBottom: 0, width: '100%' }}>
                              <Space size={12}>
                                <Button
                                  type="primary"
                                  htmlType="submit"
                                  icon={<SearchOutlined />}
                                  loading={loading}
                                  style={{
                                    borderRadius: 8,
                                    boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)',
                                  }}
                                >
                                  搜索
                                </Button>
                                <Button
                                  icon={<ReloadOutlined />}
                                  onClick={handleReset}
                                  style={{ borderRadius: 8 }}
                                >
                                  重置
                                </Button>
                              </Space>
                            </Form.Item>
                          </Col>
                        </Row>
                      </Form>
                    </Card>

                    {/* MD 结果表格 */}
                    <Card
                      title={
                        <Space>
                          <CheckCircleOutlined style={{ color: '#52c41a' }} />
                          <span>搜索结果</span>
                        </Space>
                      }
                      style={{
                        borderRadius: 12,
                        border: 'none',
                        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                      }}
                    >
                      {results.length > 0 ? (
                        <Table
                          columns={columns}
                          dataSource={results}
                          rowKey="job_id"
                          loading={loading}
                          pagination={{
                            current: pagination.current,
                            pageSize: pagination.pageSize,
                            total: total,
                            showSizeChanger: true,
                            showTotal: (total) => `共 ${total} 条记录`,
                            onChange: (page, pageSize) => {
                              setPagination({ current: page, pageSize: pageSize || 10 });
                            },
                          }}
                          scroll={{ x: 1200 }}
                        />
                      ) : (
                        <Empty
                          image={
                            <div style={{
                              width: 80,
                              height: 80,
                              borderRadius: '50%',
                              background: 'linear-gradient(135deg, rgba(102, 126, 234, 0.1) 0%, rgba(118, 75, 162, 0.1) 100%)',
                              display: 'flex',
                              alignItems: 'center',
                              justifyContent: 'center',
                              margin: '0 auto',
                            }}>
                              <FileSearchOutlined style={{ fontSize: 36, color: '#667eea' }} />
                            </div>
                          }
                          description={
                            <div>
                              <Text type="secondary" style={{ fontSize: 14 }}>
                                {isAuthenticated ? '请输入搜索条件查询数据' : '请登录后搜索数据'}
                              </Text>
                            </div>
                          }
                        >
                          {!isAuthenticated && (
                            <Button
                              type="primary"
                              icon={<LoginOutlined />}
                              onClick={() => navigate('/login', { state: { from: location.pathname } })}
                              style={{
                                borderRadius: 8,
                                boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)',
                              }}
                            >
                              立即登录
                            </Button>
                          )}
                        </Empty>
                      )}
                    </Card>
                  </div>
                ),
              },
              {
                key: 'qc',
                label: (
                  <Space>
                    <ExperimentOutlined />
                    <span>QC数据库</span>
                  </Space>
                ),
                children: (
                  <div style={{ padding: '0 0 24px 0' }}>
                    <QCDataTab isPublic={true} />
                  </div>
                ),
              },
            ]}
          />
        </Card>
      </div>
    </div>
  );
}

