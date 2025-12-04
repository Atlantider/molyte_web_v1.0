/**
 * 去溶剂化能图表组件
 * Desolvation energy chart components using ECharts
 */
import React, { useMemo } from 'react';
import ReactECharts from 'echarts-for-react';
import { Card, Row, Col, theme } from 'antd';
import type { LigandDesolvationResult, TypeSummary } from '../types/desolvation';
import type { EnergyUnit } from '../utils/energyUnits';
import { convertEnergy, UNIT_LABELS, UNIT_PRECISION } from '../utils/energyUnits';
import { useThemeStore } from '../stores/themeStore';

interface DesolvationChartsProps {
  perLigandResults: LigandDesolvationResult[];
  perTypeSummary: TypeSummary[];
  unit: EnergyUnit;
}

export default function DesolvationCharts({
  perLigandResults,
  perTypeSummary,
  unit
}: DesolvationChartsProps) {
  const { token } = theme.useToken();
  const { isDark } = useThemeStore();

  // 柱状图配置 - 按配体展示
  const barChartOption = useMemo(() => {
    const labels = perLigandResults.map(r => r.ligand_label);
    const values = perLigandResults.map(r => convertEnergy(r.delta_e, unit));
    const colors = values.map(v => v < 0 ? '#52c41a' : '#1890ff');

    return {
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
        formatter: (params: any) => {
          const data = params[0];
          return `${data.name}<br/>ΔE: ${data.value.toFixed(UNIT_PRECISION[unit])} ${UNIT_LABELS[unit]}`;
        }
      },
      grid: {
        left: '3%',
        right: '4%',
        bottom: '15%',
        containLabel: true
      },
      xAxis: {
        type: 'category',
        data: labels,
        axisLabel: {
          rotate: 45,
          color: token.colorTextSecondary,
          fontSize: 11
        },
        axisLine: { lineStyle: { color: token.colorBorder } }
      },
      yAxis: {
        type: 'value',
        name: `ΔE (${UNIT_LABELS[unit]})`,
        nameTextStyle: { color: token.colorTextSecondary },
        axisLabel: { color: token.colorTextSecondary },
        axisLine: { lineStyle: { color: token.colorBorder } },
        splitLine: { lineStyle: { color: isDark ? 'rgba(255,255,255,0.1)' : 'rgba(0,0,0,0.1)' } }
      },
      series: [{
        type: 'bar',
        data: values.map((v, i) => ({
          value: v,
          itemStyle: { color: colors[i] }
        })),
        label: {
          show: values.length <= 10,
          position: 'top',
          formatter: (params: any) => params.value.toFixed(1),
          fontSize: 10,
          color: token.colorText
        }
      }]
    };
  }, [perLigandResults, unit, token, isDark]);

  // 按类型汇总的柱状图（带误差棒）
  const typeSummaryChartOption = useMemo(() => {
    const types = perTypeSummary.map(t => t.ligand_type);
    const avgValues = perTypeSummary.map(t => convertEnergy(t.avg_delta_e, unit));
    const stdValues = perTypeSummary.map(t => convertEnergy(t.std_delta_e, unit));
    const counts = perTypeSummary.map(t => t.count);

    return {
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        formatter: (params: any) => {
          const idx = params[0].dataIndex;
          const type = types[idx];
          const avg = avgValues[idx].toFixed(UNIT_PRECISION[unit]);
          const std = stdValues[idx].toFixed(UNIT_PRECISION[unit]);
          const count = counts[idx];
          return `<b>${type}</b><br/>平均: ${avg} ± ${std} ${UNIT_LABELS[unit]}<br/>数量: ${count}`;
        }
      },
      grid: {
        left: '3%',
        right: '4%',
        bottom: '10%',
        containLabel: true
      },
      xAxis: {
        type: 'category',
        data: types,
        axisLabel: { color: token.colorTextSecondary },
        axisLine: { lineStyle: { color: token.colorBorder } }
      },
      yAxis: {
        type: 'value',
        name: `平均 ΔE (${UNIT_LABELS[unit]})`,
        nameTextStyle: { color: token.colorTextSecondary },
        axisLabel: { color: token.colorTextSecondary },
        splitLine: { lineStyle: { color: isDark ? 'rgba(255,255,255,0.1)' : 'rgba(0,0,0,0.1)' } }
      },
      series: [
        {
          type: 'bar',
          data: avgValues,
          itemStyle: { color: '#1890ff' },
          label: {
            show: true,
            position: 'top',
            formatter: (p: any) => p.value.toFixed(1),
            color: token.colorText
          }
        },
        {
          type: 'custom',
          renderItem: (params: any, api: any) => {
            const xValue = api.value(0);
            const yValue = api.value(1);
            const highPoint = api.coord([xValue, yValue + stdValues[params.dataIndex]]);
            const lowPoint = api.coord([xValue, yValue - stdValues[params.dataIndex]]);
            const halfWidth = 6;
            return {
              type: 'group',
              children: [
                { type: 'line', shape: { x1: highPoint[0], y1: highPoint[1], x2: lowPoint[0], y2: lowPoint[1] }, style: { stroke: '#666', lineWidth: 2 } },
                { type: 'line', shape: { x1: highPoint[0] - halfWidth, y1: highPoint[1], x2: highPoint[0] + halfWidth, y2: highPoint[1] }, style: { stroke: '#666', lineWidth: 2 } },
                { type: 'line', shape: { x1: lowPoint[0] - halfWidth, y1: lowPoint[1], x2: lowPoint[0] + halfWidth, y2: lowPoint[1] }, style: { stroke: '#666', lineWidth: 2 } }
              ]
            };
          },
          data: avgValues.map((v, i) => [types[i], v]),
          z: 100
        }
      ]
    };
  }, [perTypeSummary, unit, token, isDark]);

  return (
    <Row gutter={16}>
      <Col span={12}>
        <Card size="small" title="按配体去溶剂化能" bordered={false}>
          <ReactECharts option={barChartOption} style={{ height: 300 }} />
        </Card>
      </Col>
      <Col span={12}>
        <Card size="small" title="按类型汇总 (平均值±标准差)" bordered={false}>
          <ReactECharts option={typeSummaryChartOption} style={{ height: 300 }} />
        </Card>
      </Col>
    </Row>
  );
}

