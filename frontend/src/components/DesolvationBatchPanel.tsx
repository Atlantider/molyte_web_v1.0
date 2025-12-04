/**
 * 批量去溶剂化能计算面板
 * 
 * 功能：
 * 1. 自动挑选不同配位组成的溶剂化结构
 * 2. 显示待计算的结构列表（可勾选）
 * 3. 批量提交计算任务
 * 4. 显示任务进度和结果
 */
import React, { useState, useEffect, useCallback } from 'react';
import {
  Card,
  Table,
  Button,
  Space,
  Tag,
  Progress,
  Select,
  Collapse,
  Row,
  Col,
  InputNumber,
  Typography,
  message,
  Tooltip,
  Badge,
  Spin,
  Empty,
  Divider,
  theme,
} from 'antd';
import {
  ThunderboltOutlined,
  ReloadOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  ExclamationCircleOutlined,
  BulbOutlined,
  RightOutlined,
  DownOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import { autoSelectSolvationStructures, type AutoSelectedStructure } from '../api/jobs';
import {
  batchCreateDesolvationJobs,
  getDesolvationOverview,
} from '../api/desolvation';
import type {
  DesolvationJobResponse,
  DesolvationOverviewResponse,
  SolventModel,
  SolventConfig,
} from '../types/desolvation';
import { useThemeStore } from '../stores/themeStore';
import DesolvationResultView from './DesolvationResultView';

const { Text, Title } = Typography;

interface DesolvationBatchPanelProps {
  jobId: number;  // MD Job ID
  onStructureSelect?: (structureId: number) => void;  // 选中结构时的回调
}

interface SelectedStructure extends AutoSelectedStructure {
  selected: boolean;
}

export default function DesolvationBatchPanel({ jobId, onStructureSelect }: DesolvationBatchPanelProps) {
  const { token } = theme.useToken();
  const { isDark } = useThemeStore();
  
  // 状态
  const [loading, setLoading] = useState(false);
  const [structures, setStructures] = useState<SelectedStructure[]>([]);
  const [selectedKeys, setSelectedKeys] = useState<number[]>([]);
  const [submitting, setSubmitting] = useState(false);
  const [overview, setOverview] = useState<DesolvationOverviewResponse | null>(null);
  const [expandedJobId, setExpandedJobId] = useState<number | null>(null);

  // 筛选条件
  const [cnFilter, setCnFilter] = useState<number[]>([]);  // 配位数筛选

  // 计算参数
  const [methodLevel, setMethodLevel] = useState<'fast' | 'standard' | 'accurate'>('standard');
  const [desolvationMode, setDesolvationMode] = useState<'stepwise' | 'full'>('stepwise');
  const [solventModel, setSolventModel] = useState<SolventModel>('gas');
  const [solventName, setSolventName] = useState<string>('');
  const [customParams, setCustomParams] = useState<Partial<SolventConfig>>({});

  // 获取所有可用的配位数选项
  const availableCNs = React.useMemo(() => {
    const cnSet = new Set<number>();
    structures.forEach(s => cnSet.add(s.coordination_num));
    return Array.from(cnSet).sort((a, b) => a - b);
  }, [structures]);

  // 根据筛选条件过滤后的结构
  const filteredStructures = React.useMemo(() => {
    if (cnFilter.length === 0) return structures;
    return structures.filter(s => cnFilter.includes(s.coordination_num));
  }, [structures, cnFilter]);

  // 当筛选条件变化时，更新选中的 keys
  useEffect(() => {
    if (cnFilter.length > 0) {
      const filteredIds = filteredStructures.map(s => s.id);
      setSelectedKeys(prev => prev.filter(id => filteredIds.includes(id)));
    }
  }, [cnFilter, filteredStructures]);

  // 加载自动挑选的结构
  const loadAutoSelectedStructures = useCallback(async () => {
    setLoading(true);
    try {
      const result = await autoSelectSolvationStructures(jobId);
      const selected = result.selected_structures.map(s => ({
        ...s,
        selected: true,
      }));
      setStructures(selected);
      setSelectedKeys(selected.map(s => s.id));
      message.success(`已自动挑选 ${result.unique_compositions} 种不同配位组成`);
    } catch (error) {
      message.error('加载溶剂化结构失败');
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  // 加载任务总览
  const loadOverview = useCallback(async () => {
    try {
      const data = await getDesolvationOverview(jobId);
      setOverview(data);
    } catch (error) {
      // 可能没有任务，忽略错误
    }
  }, [jobId]);

  useEffect(() => {
    loadOverview();
  }, [loadOverview]);

  // 批量提交
  const handleBatchSubmit = async () => {
    if (selectedKeys.length === 0) {
      message.warning('请选择要计算的溶剂化结构');
      return;
    }
    
    setSubmitting(true);
    try {
      const solventConfig: SolventConfig | undefined = solventModel === 'gas' ? undefined : {
        model: solventModel,
        solvent_name: solventModel !== 'custom' ? solventName : undefined,
        ...customParams,
      };
      
      const result = await batchCreateDesolvationJobs({
        md_job_id: jobId,
        structure_ids: selectedKeys,
        method_level: methodLevel,
        desolvation_mode: desolvationMode,
        solvent_config: solventConfig,
      });
      
      message.success(`已创建 ${result.created_count} 个任务，跳过 ${result.skipped_count} 个已存在任务`);
      loadOverview();
    } catch (error: any) {
      message.error(`提交失败: ${error.message || '未知错误'}`);
    } finally {
      setSubmitting(false);
    }
  };

  // 结构表格列
  const structureColumns: ColumnsType<SelectedStructure> = [
    {
      title: '配位组成',
      dataIndex: 'composition_key',
      key: 'composition_key',
      render: (key: string, record) => (
        <Space direction="vertical" size={0}>
          <Text strong style={{ fontSize: 13 }}>{key}</Text>
          <Text type="secondary" style={{ fontSize: 11 }}>
            {record.center_ion}⁺ CN={record.coordination_num}
          </Text>
        </Space>
      ),
    },
    {
      title: '分子组成',
      dataIndex: 'composition',
      key: 'composition',
      render: (composition: Record<string, number>) => (
        <Space size={4} wrap>
          {Object.entries(composition)
            .filter(([_, count]) => count > 0)
            .map(([mol, count]) => (
              <Tag key={mol} style={{ margin: 0, fontSize: 11 }}>
                {mol}: {count}
              </Tag>
            ))}
        </Space>
      ),
    },
    {
      title: '帧号',
      dataIndex: 'frame_index',
      key: 'frame_index',
      width: 80,
      render: (frame: number) => <Text type="secondary">#{frame}</Text>,
    },
  ];

  // 任务表格列
  const jobColumns: ColumnsType<DesolvationJobResponse> = [
    {
      title: '结构',
      key: 'structure',
      width: 200,
      render: (_, record) => (
        <Space direction="vertical" size={0}>
          <Text strong style={{ fontSize: 12 }}>
            {record.composition_key || `结构 #${record.solvation_structure_id}`}
          </Text>
          {record.electrolyte_name && (
            <Text type="secondary" style={{ fontSize: 11 }}>
              {record.electrolyte_name}
            </Text>
          )}
        </Space>
      ),
    },
    {
      title: '方法',
      dataIndex: 'method_level',
      key: 'method_level',
      width: 100,
      render: (level: string) => {
        const config: Record<string, { color: string; text: string }> = {
          fast: { color: 'green', text: '快速' },
          standard: { color: 'blue', text: '标准' },
          accurate: { color: 'purple', text: '精确' },
        };
        const c = config[level] || { color: 'default', text: level };
        return <Tag color={c.color}>{c.text}</Tag>;
      },
    },
    {
      title: '状态',
      key: 'status',
      width: 150,
      render: (_, record) => {
        const statusConfig: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
          CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
          SUBMITTED: { color: 'blue', icon: <ClockCircleOutlined />, text: '已提交' },
          QUEUED: { color: 'cyan', icon: <ClockCircleOutlined />, text: '排队中' },
          RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
          COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
          FAILED: { color: 'error', icon: <ExclamationCircleOutlined />, text: '失败' },
        };
        const config = statusConfig[record.status] || { color: 'default', icon: null, text: record.status };
        
        return (
          <Space direction="vertical" size={0}>
            <Tag color={config.color} icon={config.icon}>{config.text}</Tag>
            {record.qc_progress && (
              <Progress
                percent={record.qc_progress.progress_percent}
                size="small"
                style={{ width: 100 }}
                format={() => `${record.qc_progress?.completed}/${record.qc_progress?.total}`}
              />
            )}
          </Space>
        );
      },
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 140,
      render: (time: string) => new Date(time).toLocaleString('zh-CN', {
        month: '2-digit',
        day: '2-digit',
        hour: '2-digit',
        minute: '2-digit',
      }),
    },
    {
      title: '操作',
      key: 'action',
      width: 100,
      render: (_, record) => (
        <Button
          type="link"
          size="small"
          disabled={record.status !== 'COMPLETED'}
          onClick={() => setExpandedJobId(expandedJobId === record.job_id ? null : record.job_id)}
        >
          {expandedJobId === record.job_id ? '收起' : '查看结果'}
        </Button>
      ),
    },
  ];

  // 检测是否有阴离子
  const hasAnion = structures.some(s => {
    const anionPatterns = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'NO3', 'Cl', 'Br', 'I', 'OTf', 'BOB'];
    return Object.keys(s.composition).some(mol =>
      anionPatterns.some(anion => mol.toUpperCase().includes(anion.toUpperCase()))
    );
  });

  return (
    <Card
      title={
        <Space>
          <ThunderboltOutlined style={{ color: '#1890ff' }} />
          <span>去溶剂化能计算</span>
          {overview && overview.total_jobs > 0 && (
            <Badge
              count={overview.status_summary['RUNNING'] || 0}
              style={{ backgroundColor: '#1890ff' }}
              title="运行中的任务"
            />
          )}
        </Space>
      }
      extra={
        <Button
          icon={<ReloadOutlined />}
          onClick={loadOverview}
          size="small"
        >
          刷新
        </Button>
      }
      style={{
        background: isDark ? token.colorBgContainer : undefined,
        borderColor: token.colorBorder,
      }}
    >
      {/* 第一步：挑选结构 */}
      <Collapse
        defaultActiveKey={structures.length === 0 ? ['select'] : []}
        items={[{
          key: 'select',
          label: (
            <Space>
              <span>第一步：挑选溶剂化结构</span>
              {structures.length > 0 && (
                <Tag color="blue">{selectedKeys.length} 个已选</Tag>
              )}
            </Space>
          ),
          children: (
            <div>
              <div style={{ marginBottom: 16 }}>
                <Space wrap>
                  <Button
                    type="primary"
                    icon={<BulbOutlined />}
                    onClick={loadAutoSelectedStructures}
                    loading={loading}
                  >
                    自动挑选不同配位组成
                  </Button>
                  <Text type="secondary" style={{ fontSize: 12 }}>
                    系统会自动从所有溶剂化结构中挑选出不同配位组成的代表性结构
                  </Text>
                </Space>
              </div>

              {/* 配位数筛选器 */}
              {structures.length > 0 && availableCNs.length > 1 && (
                <div style={{ marginBottom: 12 }}>
                  <Space size={8} align="center">
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary }}>
                      按配位数筛选：
                    </Text>
                    <Select
                      mode="multiple"
                      placeholder="全部配位数"
                      value={cnFilter}
                      onChange={(values) => {
                        setCnFilter(values);
                        // 如果选择了筛选条件，自动选中所有符合条件的结构
                        if (values.length > 0) {
                          const filtered = structures.filter(s => values.includes(s.coordination_num));
                          setSelectedKeys(filtered.map(s => s.id));
                        } else {
                          // 清空筛选时恢复全选
                          setSelectedKeys(structures.map(s => s.id));
                        }
                      }}
                      style={{ minWidth: 200 }}
                      size="small"
                      allowClear
                      options={availableCNs.map(cn => ({
                        label: `CN = ${cn}`,
                        value: cn,
                      }))}
                    />
                    {cnFilter.length > 0 && (
                      <Tag color="orange">
                        筛选后 {filteredStructures.length} 个结构
                      </Tag>
                    )}
                    <Button
                      size="small"
                      onClick={() => setSelectedKeys(filteredStructures.map(s => s.id))}
                    >
                      全选当前
                    </Button>
                    <Button
                      size="small"
                      onClick={() => setSelectedKeys([])}
                    >
                      清空选择
                    </Button>
                  </Space>
                </div>
              )}

              {structures.length > 0 && (
                <Table
                  dataSource={filteredStructures}
                  columns={structureColumns}
                  rowKey="id"
                  size="small"
                  rowSelection={{
                    selectedRowKeys: selectedKeys,
                    onChange: (keys) => setSelectedKeys(keys as number[]),
                  }}
                  pagination={false}
                  scroll={{ y: 200 }}
                  onRow={(record) => ({
                    onClick: () => onStructureSelect?.(record.id),
                    style: { cursor: 'pointer' },
                  })}
                />
              )}
            </div>
          ),
        }]}
      />

      {/* 第二步：设置参数并提交 */}
      {structures.length > 0 && (
        <Collapse
          style={{ marginTop: 16 }}
          defaultActiveKey={['params']}
          items={[{
            key: 'params',
            label: '第二步：设置计算参数并提交',
            children: (
              <div>
                {/* 智能推荐 */}
                {hasAnion && (
                  <div style={{
                    marginBottom: 16,
                    padding: '8px 12px',
                    background: isDark ? 'rgba(250, 173, 20, 0.1)' : '#fffbe6',
                    border: `1px solid ${isDark ? 'rgba(250, 173, 20, 0.3)' : '#ffe58f'}`,
                    borderRadius: 6,
                  }}>
                    <Space size={4}>
                      <BulbOutlined style={{ color: '#faad14' }} />
                      <Text style={{ fontSize: 12 }}>
                        <strong>智能推荐：</strong>检测到阴离子，建议选择带弥散函数的基组（标准或精确）
                      </Text>
                    </Space>
                  </div>
                )}

                <Row gutter={[16, 16]}>
                  <Col span={8}>
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                      计算方法
                    </Text>
                    <Select
                      value={methodLevel}
                      onChange={setMethodLevel}
                      style={{ width: '100%' }}
                      options={[
                        { label: '快速 (B3LYP/6-31G(d))', value: 'fast' },
                        { label: '标准 (B3LYP/6-31++G(d,p))', value: 'standard' },
                        { label: '精确 (ωB97XD/6-311++G(2d,2p))', value: 'accurate' },
                      ]}
                    />
                  </Col>
                  <Col span={8}>
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                      计算模式
                    </Text>
                    <Select
                      value={desolvationMode}
                      onChange={setDesolvationMode}
                      style={{ width: '100%' }}
                      options={[
                        { label: '逐级去溶剂 (推荐)', value: 'stepwise' },
                        { label: '全部去溶剂', value: 'full' },
                      ]}
                    />
                  </Col>
                  <Col span={8}>
                    <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                      溶剂模型
                    </Text>
                    <Select
                      value={solventModel}
                      onChange={setSolventModel}
                      style={{ width: '100%' }}
                      options={[
                        { label: '气相', value: 'gas' },
                        { label: 'PCM', value: 'pcm' },
                        { label: 'SMD', value: 'smd' },
                      ]}
                    />
                  </Col>
                </Row>

                <div style={{ marginTop: 16 }}>
                  <Button
                    type="primary"
                    icon={<ThunderboltOutlined />}
                    onClick={handleBatchSubmit}
                    loading={submitting}
                    disabled={selectedKeys.length === 0}
                    size="large"
                  >
                    批量创建计算任务 ({selectedKeys.length} 个)
                  </Button>
                </div>
              </div>
            ),
          }]}
        />
      )}

      {/* 第三步：任务监控 */}
      {overview && overview.total_jobs > 0 && (
        <div style={{ marginTop: 16 }}>
          <Divider orientation="left">
            <Space>
              任务监控
              <Tag color="blue">{overview.total_jobs} 个任务</Tag>
              {overview.status_summary['COMPLETED'] > 0 && (
                <Tag color="success">{overview.status_summary['COMPLETED']} 完成</Tag>
              )}
              {overview.status_summary['RUNNING'] > 0 && (
                <Tag color="processing">{overview.status_summary['RUNNING']} 运行中</Tag>
              )}
            </Space>
          </Divider>
          
          <Table
            dataSource={overview.jobs}
            columns={jobColumns}
            rowKey="job_id"
            size="small"
            pagination={{ pageSize: 5, size: 'small' }}
            expandable={{
              expandedRowKeys: expandedJobId ? [expandedJobId] : [],
              expandedRowRender: (record) =>
                record.result ? (
                  <DesolvationResultView result={record.result} />
                ) : (
                  <Empty description="暂无结果" />
                ),
              rowExpandable: (record) => record.status === 'COMPLETED',
              expandIcon: () => null,
            }}
          />
        </div>
      )}
    </Card>
  );
}

