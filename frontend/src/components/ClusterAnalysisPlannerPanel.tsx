/**
 * Cluster 高级计算规划面板
 * 
 * 功能：
 * 1. 结构筛选和选择
 * 2. 计算类型多选（Binding/Desolvation/Redox/Reorg）
 * 3. QC 任务复用预览
 * 4. 提交和追踪任务
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
  Checkbox,
  Row,
  Col,
  Typography,
  message,
  Tooltip,
  Spin,
  Empty,
  Divider,
  Alert,
  Modal,
  Statistic,
  Badge,
  Collapse,
} from 'antd';
import {
  ThunderboltOutlined,
  ReloadOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  ExclamationCircleOutlined,
  SendOutlined,
  InfoCircleOutlined,
  RocketOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import { autoSelectSolvationStructures, type AutoSelectedStructure } from '../api/jobs';
import {
  planClusterAnalysis,
  submitClusterAnalysis,
  listClusterAnalysisJobs,
  CALC_TYPE_INFO,
  type ClusterCalcType,
  type ClusterAnalysisPlanResponse,
  type AdvancedClusterJob,
  type CalcTypeRequirements,
} from '../api/clusterAnalysis';

const { Text, Title, Paragraph } = Typography;
const { Panel } = Collapse;

interface Props {
  mdJobId: number;
}

// 计算类型选项
const CALC_TYPE_OPTIONS: { value: ClusterCalcType; label: string; risk: string }[] = [
  { value: 'BINDING_TOTAL', label: '🔗 总 Binding Energy', risk: 'low' },
  { value: 'BINDING_PAIRWISE', label: '⚛️ 分子-Li Binding', risk: 'low' },
  { value: 'DESOLVATION_STEPWISE', label: '📉 逐级去溶剂化', risk: 'medium' },
  { value: 'DESOLVATION_FULL', label: '🎯 完全去溶剂化', risk: 'low' },
  { value: 'REDOX', label: '⚡ 氧化还原电位', risk: 'high' },
  { value: 'REORGANIZATION', label: '🔄 Marcus 重组能', risk: 'high' },
];

export default function ClusterAnalysisPlannerPanel({ mdJobId }: Props) {
  // 状态
  const [loading, setLoading] = useState(false);
  const [structures, setStructures] = useState<AutoSelectedStructure[]>([]);
  const [selectedStructureIds, setSelectedStructureIds] = useState<number[]>([]);
  const [selectedCalcTypes, setSelectedCalcTypes] = useState<ClusterCalcType[]>([]);
  const [planResult, setPlanResult] = useState<ClusterAnalysisPlanResponse | null>(null);
  const [planLoading, setPlanLoading] = useState(false);
  const [submitting, setSubmitting] = useState(false);
  const [existingJobs, setExistingJobs] = useState<AdvancedClusterJob[]>([]);
  
  // QC 配置
  const [qcConfig, setQcConfig] = useState({
    functional: 'B3LYP',
    basis_set: '6-31G*',
    use_dispersion: true,
    charge_ion: 1,
  });

  // 加载溶剂化结构
  const loadStructures = useCallback(async () => {
    setLoading(true);
    try {
      const result = await autoSelectSolvationStructures(mdJobId, { max_per_composition: 100 });
      setStructures(result.selected);
    } catch (error) {
      message.error('加载溶剂化结构失败');
      console.error(error);
    } finally {
      setLoading(false);
    }
  }, [mdJobId]);

  // 加载已有任务
  const loadExistingJobs = useCallback(async () => {
    try {
      const jobs = await listClusterAnalysisJobs(mdJobId);
      setExistingJobs(jobs);
    } catch (error) {
      console.error('加载已有任务失败:', error);
    }
  }, [mdJobId]);

  useEffect(() => {
    loadStructures();
    loadExistingJobs();
  }, [loadStructures, loadExistingJobs]);

  // 规划计算
  const handlePlan = async () => {
    if (selectedStructureIds.length === 0) {
      message.warning('请先选择溶剂化结构');
      return;
    }
    if (selectedCalcTypes.length === 0) {
      message.warning('请选择至少一种计算类型');
      return;
    }

    setPlanLoading(true);
    try {
      const result = await planClusterAnalysis({
        md_job_id: mdJobId,
        solvation_structure_ids: selectedStructureIds,
        calc_types: selectedCalcTypes,
        qc_config: qcConfig,
      });
      setPlanResult(result);
      message.success('规划完成');
    } catch (error) {
      message.error('规划失败');
      console.error(error);
    } finally {
      setPlanLoading(false);
    }
  };

  // 提交任务
  const handleSubmit = async () => {
    if (!planResult) return;

    Modal.confirm({
      title: '确认提交计算任务',
      content: (
        <div>
          <p>将提交以下计算：</p>
          <ul>
            {selectedCalcTypes.map(ct => (
              <li key={ct}>{CALC_TYPE_INFO[ct].label}</li>
            ))}
          </ul>
          <p>
            <strong>新建 QC 任务：</strong> {planResult.total_new_qc_tasks} 个
          </p>
          <p>
            <strong>复用已有任务：</strong> {planResult.total_reused_qc_tasks} 个
          </p>
          <p>
            <strong>预估时间：</strong> {planResult.estimated_compute_hours.toFixed(1)} 小时
          </p>
        </div>
      ),
      okText: '提交',
      cancelText: '取消',
      onOk: async () => {
        setSubmitting(true);
        try {
          await submitClusterAnalysis({
            md_job_id: mdJobId,
            solvation_structure_ids: selectedStructureIds,
            calc_types: selectedCalcTypes,
            qc_config: qcConfig,
          });
          message.success('任务已提交');
          setPlanResult(null);
          loadExistingJobs();
        } catch (error) {
          message.error('提交失败');
          console.error(error);
        } finally {
          setSubmitting(false);
        }
      },
    });
  };

  // 结构表格列定义
  const structureColumns: ColumnsType<AutoSelectedStructure> = [
    {
      title: '结构 ID',
      dataIndex: 'structure_id',
      key: 'structure_id',
      width: 80,
    },
    {
      title: '配位数',
      dataIndex: 'coordination_number',
      key: 'coordination_number',
      width: 80,
      render: (cn: number) => <Tag color="blue">{cn}</Tag>,
    },
    {
      title: '组成',
      dataIndex: 'composition',
      key: 'composition',
      render: (comp: Record<string, number>) => (
        <Space size="small" wrap>
          {Object.entries(comp || {}).map(([mol, count]) => (
            count > 0 && <Tag key={mol}>{mol}: {count}</Tag>
          ))}
        </Space>
      ),
    },
    {
      title: '帧号',
      dataIndex: 'frame',
      key: 'frame',
      width: 80,
    },
  ];

  // 渲染计算类型复选框
  const renderCalcTypeCheckboxes = () => (
    <div style={{ marginBottom: 16 }}>
      <Text strong>选择计算类型：</Text>
      <Row gutter={[16, 8]} style={{ marginTop: 8 }}>
        {CALC_TYPE_OPTIONS.map(opt => (
          <Col key={opt.value} span={8}>
            <Checkbox
              checked={selectedCalcTypes.includes(opt.value)}
              onChange={e => {
                if (e.target.checked) {
                  setSelectedCalcTypes([...selectedCalcTypes, opt.value]);
                } else {
                  setSelectedCalcTypes(selectedCalcTypes.filter(t => t !== opt.value));
                }
                setPlanResult(null); // 清除之前的规划结果
              }}
            >
              <Space>
                {opt.label}
                {opt.risk === 'high' && (
                  <Tag color="red" style={{ marginLeft: 4 }}>⚠️ 高风险</Tag>
                )}
                {opt.risk === 'medium' && (
                  <Tag color="orange" style={{ marginLeft: 4 }}>中等</Tag>
                )}
              </Space>
            </Checkbox>
          </Col>
        ))}
      </Row>
    </div>
  );

  // 渲染规划结果
  const renderPlanResult = () => {
    if (!planResult) return null;

    return (
      <Card
        title={<><RocketOutlined /> QC 任务规划预览</>}
        style={{ marginTop: 16 }}
        extra={
          <Button
            type="primary"
            icon={<SendOutlined />}
            loading={submitting}
            onClick={handleSubmit}
          >
            提交计算
          </Button>
        }
      >
        {/* 汇总统计 */}
        <Row gutter={16} style={{ marginBottom: 16 }}>
          <Col span={6}>
            <Statistic
              title="选中结构"
              value={planResult.selected_structures_count}
              suffix="个"
            />
          </Col>
          <Col span={6}>
            <Statistic
              title="新建 QC 任务"
              value={planResult.total_new_qc_tasks}
              suffix="个"
              valueStyle={{ color: '#1890ff' }}
            />
          </Col>
          <Col span={6}>
            <Statistic
              title="复用已有任务"
              value={planResult.total_reused_qc_tasks}
              suffix="个"
              valueStyle={{ color: '#52c41a' }}
            />
          </Col>
          <Col span={6}>
            <Statistic
              title="预估时间"
              value={planResult.estimated_compute_hours.toFixed(1)}
              suffix="小时"
            />
          </Col>
        </Row>

        {/* 警告 */}
        {planResult.warnings.length > 0 && (
          <Alert
            type="warning"
            message="注意事项"
            description={
              <ul style={{ margin: 0, paddingLeft: 20 }}>
                {planResult.warnings.map((w, i) => <li key={i}>{w}</li>)}
              </ul>
            }
            style={{ marginBottom: 16 }}
          />
        )}

        {/* 各计算类型详情 */}
        <Collapse>
          {planResult.calc_requirements.map(req => (
            <Panel
              key={req.calc_type}
              header={
                <Space>
                  {CALC_TYPE_INFO[req.calc_type].icon} {CALC_TYPE_INFO[req.calc_type].label}
                  <Tag color="blue">新建 {req.new_tasks_count}</Tag>
                  <Tag color="green">复用 {req.reused_tasks_count}</Tag>
                </Space>
              }
            >
              <Paragraph type="secondary">
                公式：<code>{CALC_TYPE_INFO[req.calc_type].formula}</code>
              </Paragraph>
              <Table
                size="small"
                dataSource={req.required_qc_tasks}
                rowKey={(_, i) => `${req.calc_type}-${i}`}
                pagination={false}
                columns={[
                  { title: '类型', dataIndex: 'task_type', width: 100 },
                  { title: '描述', dataIndex: 'description' },
                  {
                    title: '状态',
                    dataIndex: 'status',
                    width: 100,
                    render: (status: string) => (
                      status === 'reused'
                        ? <Tag color="green"><CheckCircleOutlined /> 复用</Tag>
                        : <Tag color="blue"><ClockCircleOutlined /> 新建</Tag>
                    )
                  },
                ]}
              />
            </Panel>
          ))}
        </Collapse>
      </Card>
    );
  };

  // 渲染已有任务列表
  const renderExistingJobs = () => {
    if (existingJobs.length === 0) return null;

    return (
      <Card
        title="已有计算任务"
        style={{ marginTop: 16 }}
        size="small"
      >
        <Table
          size="small"
          dataSource={existingJobs}
          rowKey="id"
          pagination={false}
          columns={[
            { title: 'ID', dataIndex: 'id', width: 60 },
            {
              title: '计算类型',
              dataIndex: 'calc_types',
              render: (types: string[]) => (
                <Space size="small" wrap>
                  {types.map(t => (
                    <Tag key={t}>{CALC_TYPE_INFO[t as ClusterCalcType]?.icon} {t}</Tag>
                  ))}
                </Space>
              )
            },
            {
              title: '状态',
              dataIndex: 'status',
              width: 120,
              render: (status: string) => {
                const colors: Record<string, string> = {
                  COMPLETED: 'green',
                  RUNNING: 'blue',
                  WAITING_QC: 'orange',
                  FAILED: 'red',
                  SUBMITTED: 'cyan',
                };
                return <Tag color={colors[status] || 'default'}>{status}</Tag>;
              }
            },
            {
              title: '进度',
              dataIndex: 'progress',
              width: 100,
              render: (p: number) => <Progress percent={Math.round(p)} size="small" />
            },
            {
              title: '创建时间',
              dataIndex: 'created_at',
              width: 150,
              render: (t: string) => new Date(t).toLocaleString()
            },
          ]}
        />
      </Card>
    );
  };

  return (
    <Card
      title={
        <Space>
          <ExperimentOutlined />
          Cluster 高级计算规划
        </Space>
      }
      extra={
        <Button icon={<ReloadOutlined />} onClick={() => { loadStructures(); loadExistingJobs(); }}>
          刷新
        </Button>
      }
    >
      <Spin spinning={loading}>
        {structures.length === 0 ? (
          <Empty description="暂无溶剂化结构，请先完成 MD 计算" />
        ) : (
          <>
            {/* 步骤 1: 选择结构 */}
            <Card type="inner" title="步骤 1: 选择溶剂化结构" style={{ marginBottom: 16 }}>
              <Table
                size="small"
                dataSource={structures}
                columns={structureColumns}
                rowKey="structure_id"
                rowSelection={{
                  selectedRowKeys: selectedStructureIds,
                  onChange: keys => {
                    setSelectedStructureIds(keys as number[]);
                    setPlanResult(null);
                  },
                }}
                pagination={{ pageSize: 10 }}
              />
              <div style={{ marginTop: 8 }}>
                <Text type="secondary">
                  已选择 {selectedStructureIds.length} / {structures.length} 个结构
                </Text>
                <Button
                  type="link"
                  onClick={() => setSelectedStructureIds(structures.map(s => s.structure_id))}
                >
                  全选
                </Button>
                <Button
                  type="link"
                  onClick={() => setSelectedStructureIds([])}
                >
                  清空
                </Button>
              </div>
            </Card>

            {/* 步骤 2: 选择计算类型 */}
            <Card type="inner" title="步骤 2: 选择计算类型" style={{ marginBottom: 16 }}>
              {renderCalcTypeCheckboxes()}

              {/* 选中的计算类型说明 */}
              {selectedCalcTypes.length > 0 && (
                <Alert
                  type="info"
                  message={`已选择 ${selectedCalcTypes.length} 种计算`}
                  description={
                    <ul style={{ margin: 0, paddingLeft: 20 }}>
                      {selectedCalcTypes.map(ct => (
                        <li key={ct}>
                          <strong>{CALC_TYPE_INFO[ct].label}</strong>：{CALC_TYPE_INFO[ct].description}
                        </li>
                      ))}
                    </ul>
                  }
                />
              )}
            </Card>

            {/* 步骤 3: 规划预览 */}
            <Card type="inner" title="步骤 3: 规划与提交" style={{ marginBottom: 16 }}>
              <Space>
                <Button
                  type="primary"
                  icon={<ThunderboltOutlined />}
                  loading={planLoading}
                  onClick={handlePlan}
                  disabled={selectedStructureIds.length === 0 || selectedCalcTypes.length === 0}
                >
                  生成规划预览
                </Button>
                <Text type="secondary">
                  点击查看需要的 QC 任务和可复用的已有结果
                </Text>
              </Space>

              {renderPlanResult()}
            </Card>

            {/* 已有任务 */}
            {renderExistingJobs()}
          </>
        )}
      </Spin>
    </Card>
  );
}

