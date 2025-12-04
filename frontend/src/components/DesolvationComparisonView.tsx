/**
 * 去溶剂化能结果对比组件
 * Compare desolvation energy results across different solvation structures
 */
import React, { useState, useMemo } from 'react';
import { Card, Table, Select, Space, Button, Tooltip, Empty, theme, Row, Col } from 'antd';
import { DownloadOutlined, SwapOutlined, BarChartOutlined } from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import type { DesolvationJobResponse, TypeSummary } from '../types/desolvation';
import type { EnergyUnit } from '../utils/energyUnits';
import { getUnitOptions, convertEnergy, formatEnergy, UNIT_LABELS, UNIT_PRECISION } from '../utils/energyUnits';
import { exportTypeSummaryToCSV } from '../utils/exportData';
import { useThemeStore } from '../stores/themeStore';

interface DesolvationComparisonViewProps {
  jobs: DesolvationJobResponse[];
}

export default function DesolvationComparisonView({ jobs }: DesolvationComparisonViewProps) {
  const { token } = theme.useToken();
  const { isDark } = useThemeStore();
  const [unit, setUnit] = useState<EnergyUnit>('kcal/mol');
  const [selectedTypes, setSelectedTypes] = useState<string[]>([]);

  // 过滤出有结果的任务
  const completedJobs = useMemo(() => 
    jobs.filter(j => j.status === 'completed' && j.result), 
    [jobs]
  );

  // 获取所有配体类型
  const allLigandTypes = useMemo(() => {
    const types = new Set<string>();
    completedJobs.forEach(job => {
      job.result?.per_type_summary.forEach(t => types.add(t.ligand_type));
    });
    return Array.from(types);
  }, [completedJobs]);

  // 初始化选中的类型
  React.useEffect(() => {
    if (selectedTypes.length === 0 && allLigandTypes.length > 0) {
      setSelectedTypes(allLigandTypes);
    }
  }, [allLigandTypes]);

  // 对比图表配置
  const comparisonChartOption = useMemo(() => {
    if (completedJobs.length === 0 || selectedTypes.length === 0) return null;

    const compositions = completedJobs.map(j => j.composition_key || `结构-${j.solvation_structure_id}`);
    const series = selectedTypes.map(type => {
      const data = completedJobs.map(job => {
        const typeSummary = job.result?.per_type_summary.find(t => t.ligand_type === type);
        return typeSummary ? convertEnergy(typeSummary.avg_delta_e, unit) : null;
      });
      return {
        name: type,
        type: 'bar',
        data,
        label: {
          show: completedJobs.length <= 5,
          position: 'top',
          formatter: (p: any) => p.value?.toFixed(1) || '',
          fontSize: 10
        }
      };
    });

    return {
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
        formatter: (params: any) => {
          let html = `<b>${params[0].axisValue}</b><br/>`;
          params.forEach((p: any) => {
            if (p.value !== null) {
              html += `${p.marker} ${p.seriesName}: ${p.value.toFixed(UNIT_PRECISION[unit])} ${UNIT_LABELS[unit]}<br/>`;
            }
          });
          return html;
        }
      },
      legend: {
        data: selectedTypes,
        textStyle: { color: token.colorText }
      },
      grid: { left: '3%', right: '4%', bottom: '15%', containLabel: true },
      xAxis: {
        type: 'category',
        data: compositions,
        axisLabel: { rotate: 30, color: token.colorTextSecondary, fontSize: 11 },
        axisLine: { lineStyle: { color: token.colorBorder } }
      },
      yAxis: {
        type: 'value',
        name: `平均 ΔE (${UNIT_LABELS[unit]})`,
        nameTextStyle: { color: token.colorTextSecondary },
        axisLabel: { color: token.colorTextSecondary },
        splitLine: { lineStyle: { color: isDark ? 'rgba(255,255,255,0.1)' : 'rgba(0,0,0,0.1)' } }
      },
      series
    };
  }, [completedJobs, selectedTypes, unit, token, isDark]);

  // 导出对比数据
  const handleExport = () => {
    const summaries = completedJobs.map(job => ({
      compositionKey: job.composition_key || `结构-${job.solvation_structure_id}`,
      summary: job.result?.per_type_summary || []
    }));
    exportTypeSummaryToCSV(summaries, unit, 'desolvation_comparison');
  };

  if (completedJobs.length === 0) {
    return <Empty description="暂无已完成的去溶剂化能计算结果" />;
  }

  // 对比表格数据
  const tableData = selectedTypes.map(type => {
    const row: any = { ligand_type: type };
    completedJobs.forEach(job => {
      const key = job.composition_key || `s${job.solvation_structure_id}`;
      const typeSummary = job.result?.per_type_summary.find(t => t.ligand_type === type);
      row[key] = typeSummary ? formatEnergy(typeSummary.avg_delta_e, unit) : '-';
      row[`${key}_std`] = typeSummary ? formatEnergy(typeSummary.std_delta_e, unit) : '-';
    });
    return row;
  });

  const tableColumns = [
    { title: '配体类型', dataIndex: 'ligand_type', key: 'ligand_type', fixed: 'left' as const, width: 100 },
    ...completedJobs.map(job => {
      const key = job.composition_key || `s${job.solvation_structure_id}`;
      return {
        title: job.composition_key || `结构-${job.solvation_structure_id}`,
        key,
        children: [
          { title: `平均 (${UNIT_LABELS[unit]})`, dataIndex: key, key: `${key}_avg`, width: 100 },
          { title: '标准差', dataIndex: `${key}_std`, key: `${key}_std`, width: 80 }
        ]
      };
    })
  ];

  return (
    <Card
      title={<Space><SwapOutlined /> 去溶剂化能对比</Space>}
      extra={
        <Space>
          <Select
            mode="multiple"
            value={selectedTypes}
            onChange={setSelectedTypes}
            options={allLigandTypes.map(t => ({ value: t, label: t }))}
            style={{ minWidth: 200 }}
            placeholder="选择配体类型"
            maxTagCount={2}
            size="small"
          />
          <Select value={unit} onChange={setUnit} options={getUnitOptions()} style={{ width: 120 }} size="small" />
          <Tooltip title="导出对比数据">
            <Button icon={<DownloadOutlined />} size="small" onClick={handleExport}>导出</Button>
          </Tooltip>
        </Space>
      }
    >
      {comparisonChartOption && (
        <ReactECharts option={comparisonChartOption} style={{ height: 350, marginBottom: 16 }} />
      )}
      <Table
        dataSource={tableData}
        columns={tableColumns}
        rowKey="ligand_type"
        pagination={false}
        size="small"
        scroll={{ x: 'max-content' }}
        bordered
      />
    </Card>
  );
}

