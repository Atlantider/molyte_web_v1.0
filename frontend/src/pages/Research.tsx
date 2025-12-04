import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
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
  DatabaseOutlined,
  EyeOutlined,
  CheckCircleOutlined,
  FileSearchOutlined,
  ReloadOutlined,
  ThunderboltOutlined,
  FireOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import { searchMyElectrolytes, ElectrolyteSearchResult, getAvailableSearchOptions } from '../api/research';
import { getQCJobs } from '../api/qc';
import QCDataTab from '../components/QCDataTab';
import { useThemeStore } from '../stores/themeStore';

const { Title, Text } = Typography;

export default function Research() {
  const [form] = Form.useForm();
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
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
  const [mdStats, setMdStats] = useState({ total: 0, completed: 0, public: 0 });
  const [qcStats, setQcStats] = useState({ total: 0, completed: 0, public: 0 });

  // 动态选项数据
  const [cationOptions, setCationOptions] = useState<string[]>([]);
  const [anionOptions, setAnionOptions] = useState<string[]>([]);
  const [solventOptions, setSolventOptions] = useState<string[]>([]);
  const [optionsLoading, setOptionsLoading] = useState(false);

  // 加载统计数据和可用选项
  useEffect(() => {
    loadStats();
    loadAvailableOptions();
    // 默认加载所有 MD 结果
    loadInitialData();
  }, []);

  // 初始加载数据
  const loadInitialData = async () => {
    setLoading(true);
    try {
      const response = await searchMyElectrolytes({
        skip: 0,
        limit: pagination.pageSize,
      });
      setResults(response.data);
      setTotal(response.total);
    } catch (error) {
      console.error('Failed to load initial data:', error);
    } finally {
      setLoading(false);
    }
  };

  const loadStats = async () => {
    try {
      // 加载 MD 统计
      const mdResponse = await searchMyElectrolytes({ skip: 0, limit: 1 });
      setMdStats({
        total: mdResponse.total,
        completed: mdResponse.total, // 假设搜索结果都是已完成的
        public: 0, // 需要后端提供
      });

      // 加载 QC 统计
      const qcResponse = await getQCJobs({ skip: 0, limit: 1 });
      const qcCompleted = await getQCJobs({ skip: 0, limit: 1, status: 'COMPLETED' });
      const qcPublic = await getQCJobs({ skip: 0, limit: 1, visibility: 'PUBLIC' });
      setQcStats({
        total: qcResponse.total,
        completed: qcCompleted.total,
        public: qcPublic.total,
      });
    } catch (error) {
      console.error('Failed to load stats:', error);
    }
  };

  // 加载可用的搜索选项（从数据库中实际数据提取）
  const loadAvailableOptions = async () => {
    setOptionsLoading(true);
    try {
      const options = await getAvailableSearchOptions();
      setCationOptions(options.cations);
      setAnionOptions(options.anions);
      setSolventOptions(options.solvents);
    } catch (error) {
      console.error('Failed to load available options:', error);
      message.error('加载搜索选项失败');
    } finally {
      setOptionsLoading(false);
    }
  };

  // 搜索处理
  const handleSearch = async (values: any) => {
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

      const response = await searchMyElectrolytes(params);
      setResults(response.data);
      setTotal(response.total);
      // 不显示搜索成功消息，避免干扰用户
    } catch (error: any) {
      message.error(error.response?.data?.detail || '搜索失败');
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

  // 查看详情
  const handleViewDetail = (record: ElectrolyteSearchResult) => {
    navigate(`/workspace/jobs/${record.job_id}/detail`);
  };

  // 表格列定义
  const columns = [
    {
      title: '任务名称',
      dataIndex: 'job_name',
      key: 'job_name',
      width: 200,
      sorter: (a: ElectrolyteSearchResult, b: ElectrolyteSearchResult) =>
        (a.job_name || '').localeCompare(b.job_name || ''),
      render: (name: string, record: ElectrolyteSearchResult) => (
        <div style={{ lineHeight: 1.4 }}>
          <Text strong style={{ fontSize: 12, wordBreak: 'break-all' }}>
            {name || `#${record.job_id}`}
          </Text>
          {record.user_note && (
            <div style={{ marginTop: 2 }}>
              <Text type="secondary" style={{ fontSize: 11 }}>📝 {record.user_note}</Text>
            </div>
          )}
        </div>
      ),
    },
    {
      title: '配方',
      dataIndex: 'system_name',
      key: 'system_name',
      width: 180,
      sorter: (a: ElectrolyteSearchResult, b: ElectrolyteSearchResult) =>
        (a.system_name || '').localeCompare(b.system_name || ''),
      render: (name: string) => (
        <Text style={{ fontSize: 12, wordBreak: 'break-all', lineHeight: 1.4 }}>
          {name}
        </Text>
      ),
    },
    {
      title: '阳离子',
      dataIndex: 'cations',
      key: 'cations',
      width: 100,
      render: (cations: any[]) => (
        <Space direction="vertical" size={2}>
          {cations?.map((c, i) => (
            <Tag key={i} color="red" style={{ fontSize: 11, margin: 0 }}>
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
      width: 100,
      render: (anions: any[]) => (
        <Space direction="vertical" size={2}>
          {anions?.map((a, i) => (
            <Tag key={i} color="orange" style={{ fontSize: 11, margin: 0 }}>
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
      width: 100,
      render: (solvents: any[]) => (
        <Space direction="vertical" size={2}>
          {solvents?.map((s, i) => (
            <Tag key={i} color="green" style={{ fontSize: 11, margin: 0 }}>
              {s.name} ({s.number})
            </Tag>
          ))}
        </Space>
      ),
    },
    {
      title: '浓度',
      key: 'concentration',
      width: 80,
      render: (_: any, record: ElectrolyteSearchResult) => {
        const cationConcs = record.cations
          ?.filter(c => c.concentration !== undefined && c.concentration !== null)
          .map(c => `${c.name}: ${c.concentration}M`);
        if (cationConcs && cationConcs.length > 0) {
          return (
            <Space direction="vertical" size={0}>
              {cationConcs.map((text, i) => (
                <Text key={i} style={{ fontSize: 11, whiteSpace: 'nowrap' }}>{text}</Text>
              ))}
            </Space>
          );
        }
        return <Text type="secondary">-</Text>;
      },
    },
    {
      title: '计算',
      key: 'calc_method',
      width: 70,
      render: (_: any, record: ElectrolyteSearchResult) => (
        <Space direction="vertical" size={2}>
          {record.charge_method && (
            <Tooltip title={record.charge_method === 'resp' ? 'RESP 高精度电荷' : 'LigParGen 快速电荷'}>
              <Tag color={record.charge_method === 'resp' ? 'gold' : 'cyan'} style={{ fontSize: 10, margin: 0 }}>
                {record.charge_method === 'resp' ? 'RESP' : 'LPG'}
              </Tag>
            </Tooltip>
          )}
          {record.qc_enabled && (
            <Tooltip title="包含QC量子化学计算">
              <Tag color="purple" style={{ fontSize: 10, margin: 0 }}>QC</Tag>
            </Tooltip>
          )}
          {!record.charge_method && !record.qc_enabled && <Text type="secondary">-</Text>}
        </Space>
      ),
    },
    {
      title: '温度',
      dataIndex: 'temperature',
      key: 'temperature',
      width: 60,
      sorter: (a: ElectrolyteSearchResult, b: ElectrolyteSearchResult) =>
        (a.temperature || 0) - (b.temperature || 0),
      render: (temp: number) => <Text style={{ fontSize: 11 }}>{temp ? `${temp.toFixed(0)}K` : '-'}</Text>,
    },
    {
      title: '分析',
      key: 'analysis',
      width: 90,
      render: (_: any, record: ElectrolyteSearchResult) => (
        <Space direction="vertical" size={2}>
          <Space size={2}>
            {record.has_rdf && <Tag color="success" style={{ fontSize: 10, margin: 0 }}>RDF</Tag>}
            {record.has_msd && <Tag color="processing" style={{ fontSize: 10, margin: 0 }}>MSD</Tag>}
          </Space>
          {record.has_solvation && <Tag color="warning" style={{ fontSize: 10, margin: 0 }}>溶剂化</Tag>}
        </Space>
      ),
    },
    {
      title: '操作',
      key: 'action',
      width: 90,
      fixed: 'right' as const,
      render: (_: any, record: ElectrolyteSearchResult) => (
        <Button
          type="link"
          size="small"
          icon={<EyeOutlined />}
          onClick={() => handleViewDetail(record)}
        >
          详情
        </Button>
      ),
    },
  ];

  return (
    <div style={{
      padding: '24px',
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      transition: 'background 0.3s',
    }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <DatabaseOutlined style={{ marginRight: 12, color: token.colorPrimary }} />
          数据管理
        </Title>
        <Text type="secondary">
          搜索和浏览已完成的分子动力学模拟和量子化学计算结果
        </Text>
      </div>

      {/* 统计卡片 */}
      <Card
        style={{
          marginBottom: 24,
          borderRadius: 12,
          border: `1px solid ${token.colorBorder}`,
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
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
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>MD任务</span>}
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
                title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>QC任务</span>}
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
                value={mdStats.public + qcStats.public}
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
                  <span>MD数据</span>
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
        <Form
                      form={form}
                      onFinish={handleSearch}
                      layout="vertical"
                    >
                      <Row gutter={16}>
                        <Col xs={24} sm={12} md={8}>
                          <Form.Item
                            label="阳离子"
                            name="cations"
                            tooltip="可多选，选项来自您已完成的任务"
                          >
                            <Select
                              mode="multiple"
                              placeholder={optionsLoading ? "加载中..." : "选择阳离子"}
                              options={cationOptions.map(c => ({ label: c, value: c }))}
                              allowClear
                              loading={optionsLoading}
                              style={{ borderRadius: 8 }}
                            />
                          </Form.Item>
                        </Col>

                        <Col xs={24} sm={12} md={8}>
                          <Form.Item
                            label="阴离子"
                            name="anions"
                            tooltip="可多选，选项来自您已完成的任务"
                          >
                            <Select
                              mode="multiple"
                              placeholder={optionsLoading ? "加载中..." : "选择阴离子"}
                              options={anionOptions.map(a => ({ label: a, value: a }))}
                              allowClear
                              loading={optionsLoading}
                            />
                          </Form.Item>
                        </Col>

                        <Col xs={24} sm={12} md={8}>
                          <Form.Item
                            label="溶剂"
                            name="solvents"
                            tooltip="按名称选择，可多选，选项来自您已完成的任务"
                          >
                            <Select
                              mode="multiple"
                              placeholder={optionsLoading ? "加载中..." : "选择溶剂"}
                              options={solventOptions.map(s => ({ label: s, value: s }))}
                              allowClear
                              loading={optionsLoading}
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
                          <Form.Item
                            label="最低温度 (K)"
                            name="temp_min"
                          >
                            <InputNumber
                              placeholder="273"
                              style={{ width: '100%' }}
                              min={0}
                            />
                          </Form.Item>
                        </Col>

                        <Col xs={24} sm={12} md={6}>
                          <Form.Item
                            label="最高温度 (K)"
                            name="temp_max"
                          >
                            <InputNumber
                              placeholder="373"
                              style={{ width: '100%' }}
                              min={0}
                            />
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
                              请输入搜索条件查询数据
                            </Text>
                          </div>
                        }
                      />
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
                  <span>QC数据</span>
                </Space>
              ),
              children: (
                <div style={{ padding: '0 0 24px 0' }}>
                  <QCDataTab isPublic={false} />
                </div>
              ),
            },
          ]}
        />
      </Card>
    </div>
  );
}
