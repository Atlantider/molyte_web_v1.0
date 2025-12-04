/**
 * 去溶剂化能汇总统计面板
 * MD 任务级别的去溶剂化能统计汇总
 */
import React, { useState, useMemo } from 'react';
import { Card, Row, Col, Statistic, Table, Select, Space, Button, Tooltip, Progress, theme, Divider, Typography } from 'antd';
import { DownloadOutlined, ThunderboltOutlined, CheckCircleOutlined, ClockCircleOutlined, ExclamationCircleOutlined } from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import type { DesolvationJobResponse, TypeSummary } from '../types/desolvation';
import type { EnergyUnit } from '../utils/energyUnits';
import { getUnitOptions, convertEnergy, formatEnergy, UNIT_LABELS, UNIT_PRECISION } from '../utils/energyUnits';
import { exportBatchDesolvationToCSV } from '../utils/exportData';
import { useThemeStore } from '../stores/themeStore';

const { Text } = Typography;

interface DesolvationSummaryPanelProps {
  jobs: DesolvationJobResponse[];
  mdJobId: number;
  electrolyteName?: string;
}

export default function DesolvationSummaryPanel({ jobs, mdJobId, electrolyteName }: DesolvationSummaryPanelProps) {
  const { token } = theme.useToken();
  const { isDark } = useThemeStore();
  const [unit, setUnit] = useState<EnergyUnit>('kcal/mol');

  // 统计数据
  const stats = useMemo(() => {
    const completed = jobs.filter(j => j.status === 'completed' && j.result);
    const running = jobs.filter(j => j.status === 'running' || j.status === 'pending');
    const failed = jobs.filter(j => j.status === 'failed');
    
    // 汇总所有配体类型的去溶剂化能
    const typeAggregation: Record<string, { values: number[]; count: number }> = {};
    completed.forEach(job => {
      job.result?.per_ligand_results.forEach(lig => {
        if (!typeAggregation[lig.ligand_type]) {
          typeAggregation[lig.ligand_type] = { values: [], count: 0 };
        }
        typeAggregation[lig.ligand_type].values.push(lig.delta_e);
        typeAggregation[lig.ligand_type].count++;
      });
    });

    // 计算每种类型的统计
    const typeSummary = Object.entries(typeAggregation).map(([type, data]) => {
      const values = data.values;
      const avg = values.reduce((a, b) => a + b, 0) / values.length;
      const std = Math.sqrt(values.reduce((sum, v) => sum + Math.pow(v - avg, 2), 0) / values.length);
      const min = Math.min(...values);
      const max = Math.max(...values);
      return { ligand_type: type, avg_delta_e: avg, std_delta_e: std, min_delta_e: min, max_delta_e: max, count: data.count };
    });

    return { total: jobs.length, completed: completed.length, running: running.length, failed: failed.length, typeSummary, completedJobs: completed };
  }, [jobs]);

  // 箱线图数据
  const boxPlotOption = useMemo(() => {
    if (stats.typeSummary.length === 0) return null;

    const types = stats.typeSummary.map(t => t.ligand_type);
    const boxData = stats.typeSummary.map(t => {
      const min = convertEnergy(t.min_delta_e, unit);
      const max = convertEnergy(t.max_delta_e, unit);
      const avg = convertEnergy(t.avg_delta_e, unit);
      const q1 = avg - convertEnergy(t.std_delta_e, unit);
      const q3 = avg + convertEnergy(t.std_delta_e, unit);
      return [min, q1, avg, q3, max];
    });

    return {
      backgroundColor: 'transparent',
      tooltip: { trigger: 'item', formatter: (p: any) => `${p.name}<br/>范围: ${p.data[0].toFixed(2)} ~ ${p.data[4].toFixed(2)}<br/>平均: ${p.data[2].toFixed(2)} ${UNIT_LABELS[unit]}` },
      grid: { left: '10%', right: '10%', bottom: '15%' },
      xAxis: { type: 'category', data: types, axisLabel: { color: token.colorTextSecondary } },
      yAxis: { type: 'value', name: `ΔE (${UNIT_LABELS[unit]})`, nameTextStyle: { color: token.colorTextSecondary }, axisLabel: { color: token.colorTextSecondary }, splitLine: { lineStyle: { color: isDark ? 'rgba(255,255,255,0.1)' : 'rgba(0,0,0,0.1)' } } },
      series: [{ type: 'boxplot', data: boxData, itemStyle: { color: '#1890ff', borderColor: '#1890ff' } }]
    };
  }, [stats.typeSummary, unit, token, isDark]);

  // 导出所有数据
  const handleExportAll = () => {
    const results = stats.completedJobs.map(job => ({
      compositionKey: job.composition_key || `结构-${job.solvation_structure_id}`,
      result: job.result!
    }));
    exportBatchDesolvationToCSV(results, unit, `desolvation_md${mdJobId}_all`);
  };

  const progressPercent = stats.total > 0 ? Math.round((stats.completed / stats.total) * 100) : 0;

  return (
    <Card
      title={<Space><ThunderboltOutlined style={{ color: '#1890ff' }} /> 去溶剂化能汇总 {electrolyteName && <Text type="secondary">- {electrolyteName}</Text>}</Space>}
      extra={<Space>
        <Select value={unit} onChange={setUnit} options={getUnitOptions()} style={{ width: 120 }} size="small" />
        <Tooltip title="导出所有结果"><Button icon={<DownloadOutlined />} size="small" onClick={handleExportAll} disabled={stats.completed === 0}>导出全部</Button></Tooltip>
      </Space>}
    >
      {/* 进度概览 */}
      <Row gutter={16} style={{ marginBottom: 24 }}>
        <Col span={6}><Card size="small"><Statistic title="总任务数" value={stats.total} /></Card></Col>
        <Col span={6}><Card size="small"><Statistic title={<><CheckCircleOutlined style={{ color: '#52c41a' }} /> 已完成</>} value={stats.completed} valueStyle={{ color: '#52c41a' }} /></Card></Col>
        <Col span={6}><Card size="small"><Statistic title={<><ClockCircleOutlined style={{ color: '#1890ff' }} /> 进行中</>} value={stats.running} valueStyle={{ color: '#1890ff' }} /></Card></Col>
        <Col span={6}><Card size="small"><Statistic title={<><ExclamationCircleOutlined style={{ color: '#ff4d4f' }} /> 失败</>} value={stats.failed} valueStyle={{ color: '#ff4d4f' }} /></Card></Col>
      </Row>
      <Progress percent={progressPercent} status={stats.failed > 0 ? 'exception' : progressPercent === 100 ? 'success' : 'active'} style={{ marginBottom: 24 }} />

      {stats.typeSummary.length > 0 && (
        <>
          <Divider orientation="left">按配体类型汇总统计</Divider>
          <Row gutter={16}>
            <Col span={12}>
              {boxPlotOption && <ReactECharts option={boxPlotOption} style={{ height: 300 }} />}
            </Col>
            <Col span={12}>
              <Table
                dataSource={stats.typeSummary}
                rowKey="ligand_type"
                size="small"
                pagination={false}
                columns={[
                  { title: '配体类型', dataIndex: 'ligand_type', width: 80 },
                  { title: '样本数', dataIndex: 'count', width: 60 },
                  { title: `平均 ΔE (${UNIT_LABELS[unit]})`, dataIndex: 'avg_delta_e', render: (v: number) => <b>{formatEnergy(v, unit)}</b> },
                  { title: '标准差', dataIndex: 'std_delta_e', render: (v: number) => formatEnergy(v, unit) },
                  { title: '范围', render: (_, r) => `${formatEnergy(r.min_delta_e, unit)} ~ ${formatEnergy(r.max_delta_e, unit)}` }
                ]}
              />
            </Col>
          </Row>
        </>
      )}
    </Card>
  );
}

