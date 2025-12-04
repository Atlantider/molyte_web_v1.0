/**
 * Binding Energy 视图组件
 * 从去溶剂化结果派生 Li-配体 Binding Energy 统计和可视化
 */
import React, { useState, useEffect, useMemo } from 'react';
import {
  Card,
  Row,
  Col,
  Statistic,
  Table,
  Tag,
  Space,
  Select,
  Spin,
  Empty,
  Alert,
  Tooltip,
  Typography,
  theme,
} from 'antd';
import {
  FireOutlined,
  LinkOutlined,
  BarChartOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import { getBindingStatisticsOverview, BindingStatisticsOverview, TypeBindingStats } from '../api/desolvation';
import { convertEnergy, formatEnergy, EnergyUnit, getUnitOptions } from '../utils/energyUnits';

const { Text, Title } = Typography;

interface BindingEnergyViewProps {
  mdJobId: number;
}

const BindingEnergyView: React.FC<BindingEnergyViewProps> = ({ mdJobId }) => {
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(true);
  const [data, setData] = useState<BindingStatisticsOverview | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [unit, setUnit] = useState<EnergyUnit>('kcal/mol');

  useEffect(() => {
    loadData();
  }, [mdJobId]);

  const loadData = async () => {
    setLoading(true);
    setError(null);
    try {
      const result = await getBindingStatisticsOverview(mdJobId);
      setData(result);
    } catch (err: any) {
      setError(err.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  // 转换能量值 (数据已经是 kcal/mol)
  const convert = (value: number) => convertEnergy(value, unit);
  const format = (value: number) => formatEnergy(value, unit);

  // 按类型的柱状图配置
  const typeBarChartOption = useMemo(() => {
    if (!data?.per_type_statistics) return {};
    
    const types = Object.keys(data.per_type_statistics);
    const means = types.map(t => convert(data.per_type_statistics[t].mean));
    const stds = types.map(t => convert(data.per_type_statistics[t].std));
    
    return {
      tooltip: {
        trigger: 'axis',
        formatter: (params: any) => {
          const idx = params[0].dataIndex;
          const stats = data.per_type_statistics[types[idx]];
          return `<b>${types[idx]}</b><br/>
            均值: ${format(stats.mean)}<br/>
            标准差: ${format(stats.std)}<br/>
            范围: ${format(stats.min)} ~ ${format(stats.max)}<br/>
            样本数: ${stats.count}`;
        }
      },
      xAxis: {
        type: 'category',
        data: types,
        axisLabel: { rotate: 45 }
      },
      yAxis: {
        type: 'value',
        name: `Binding Energy (${unit})`,
      },
      series: [
        {
          name: 'Mean',
          type: 'bar',
          data: means,
          itemStyle: {
            color: (params: any) => means[params.dataIndex] > 0 ? token.colorPrimary : token.colorSuccess
          },
          label: { show: true, position: 'top', formatter: (p: any) => p.value.toFixed(1) }
        },
        {
          name: 'Std',
          type: 'custom',
          renderItem: (params: any, api: any) => {
            const xValue = api.value(0);
            const yValue = api.value(1);
            const stdValue = api.value(2);
            const coord = api.coord([xValue, yValue]);
            const halfWidth = api.size([1, 0])[0] * 0.1;
            const highPoint = api.coord([xValue, yValue + stdValue]);
            const lowPoint = api.coord([xValue, yValue - stdValue]);
            return {
              type: 'group',
              children: [
                { type: 'line', shape: { x1: coord[0], y1: highPoint[1], x2: coord[0], y2: lowPoint[1] }, style: { stroke: '#333', lineWidth: 2 } },
                { type: 'line', shape: { x1: coord[0] - halfWidth, y1: highPoint[1], x2: coord[0] + halfWidth, y2: highPoint[1] }, style: { stroke: '#333', lineWidth: 2 } },
                { type: 'line', shape: { x1: coord[0] - halfWidth, y1: lowPoint[1], x2: coord[0] + halfWidth, y2: lowPoint[1] }, style: { stroke: '#333', lineWidth: 2 } },
              ]
            };
          },
          data: types.map((t, i) => [i, means[i], stds[i]]),
          z: 100
        }
      ]
    };
  }, [data, unit, token]);

  // 箱线图配置
  const boxPlotOption = useMemo(() => {
    if (!data?.per_type_statistics) return {};
    
    const types = Object.keys(data.per_type_statistics);
    const boxData = types.map(t => {
      const stats = data.per_type_statistics[t];
      const values = stats.values.map(v => convert(v)).sort((a, b) => a - b);
      const min = convert(stats.min);
      const max = convert(stats.max);
      const q1 = stats.percentile_25 ? convert(stats.percentile_25) : values[Math.floor(values.length * 0.25)];
      const q3 = stats.percentile_75 ? convert(stats.percentile_75) : values[Math.floor(values.length * 0.75)];
      const median = values[Math.floor(values.length * 0.5)];
      return [min, q1, median, q3, max];
    });

    return {
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: types, axisLabel: { rotate: 45 } },
      yAxis: { type: 'value', name: `Binding Energy (${unit})` },
      series: [{
        name: 'Binding Energy',
        type: 'boxplot',
        data: boxData,
        itemStyle: { color: token.colorPrimary, borderColor: token.colorPrimaryBorder }
      }]
    };
  }, [data, unit, token]);

  if (loading) {
    return <Spin tip="加载 Binding Energy 数据..." style={{ display: 'block', margin: '40px auto' }} />;
  }

  if (error) {
    return <Alert type="error" message="加载失败" description={error} />;
  }

  if (!data || data.total_completed === 0) {
    return <Empty description="暂无已完成的去溶剂化计算结果" />;
  }

  // 表格列
  const typeColumns = [
    {
      title: '配体类型',
      dataIndex: 'type',
      key: 'type',
      render: (t: string) => <Tag color="blue">{t}</Tag>
    },
    {
      title: `平均 Binding (${unit})`,
      dataIndex: 'mean',
      key: 'mean',
      render: (v: number) => <Text strong style={{ color: v > 0 ? token.colorPrimary : token.colorSuccess }}>{format(v)}</Text>,
      sorter: (a: any, b: any) => a.mean - b.mean
    },
    {
      title: `标准差 (${unit})`,
      dataIndex: 'std',
      key: 'std',
      render: (v: number) => format(v)
    },
    {
      title: `范围 (${unit})`,
      key: 'range',
      render: (_: any, record: any) => `${format(record.min)} ~ ${format(record.max)}`
    },
    { title: '样本数', dataIndex: 'count', key: 'count' }
  ];

  const typeTableData = Object.entries(data.per_type_statistics).map(([type, stats]) => ({
    key: type,
    type,
    ...stats
  }));

  return (
    <div>
      {/* 统计卡片 */}
      <Row gutter={16} style={{ marginBottom: 16 }}>
        <Col span={6}>
          <Card size="small">
            <Statistic
              title="已完成 Cluster"
              value={data.total_completed}
              prefix={<FireOutlined />}
            />
          </Card>
        </Col>
        <Col span={6}>
          <Card size="small">
            <Statistic
              title="总 Binding 数据"
              value={data.total_ligand_bindings}
              prefix={<LinkOutlined />}
            />
          </Card>
        </Col>
        <Col span={6}>
          <Card size="small">
            <Statistic
              title="配体类型数"
              value={Object.keys(data.per_type_statistics).length}
              prefix={<BarChartOutlined />}
            />
          </Card>
        </Col>
        <Col span={6}>
          <Card size="small">
            {data.last_layer_statistics && (
              <Statistic
                title={
                  <Tooltip title="最后一级去溶剂化能（最内层 binding，即最难脱除的配体）">
                    <span>最内层 Binding 均值 <InfoCircleOutlined /></span>
                  </Tooltip>
                }
                value={convert(data.last_layer_statistics.mean)}
                precision={2}
                suffix={unit}
              />
            )}
          </Card>
        </Col>
      </Row>

      {/* 单位选择 */}
      <div style={{ marginBottom: 16, textAlign: 'right' }}>
        <Space>
          <Text>能量单位:</Text>
          <Select
            value={unit}
            onChange={setUnit}
            style={{ width: 120 }}
            options={getUnitOptions()}
          />
        </Space>
      </div>

      {/* 图表 */}
      <Row gutter={16} style={{ marginBottom: 16 }}>
        <Col span={12}>
          <Card size="small" title="按配体类型的平均 Binding Energy">
            <ReactECharts option={typeBarChartOption} style={{ height: 300 }} />
          </Card>
        </Col>
        <Col span={12}>
          <Card size="small" title="Binding Energy 分布（箱线图）">
            <ReactECharts option={boxPlotOption} style={{ height: 300 }} />
          </Card>
        </Col>
      </Row>

      {/* 表格 */}
      <Card size="small" title="按类型统计">
        <Table
          columns={typeColumns}
          dataSource={typeTableData}
          pagination={false}
          size="small"
        />
      </Card>

      {/* 最后一级详情 */}
      {data.last_layer_statistics?.details && (
        <Card size="small" title="最内层 Binding 详情（各 Cluster 的最大 Binding）" style={{ marginTop: 16 }}>
          <Table
            columns={[
              { title: 'Cluster', dataIndex: 'composition_key', key: 'composition_key', render: (v: string) => <Tag>{v}</Tag> },
              { title: '配体类型', dataIndex: 'ligand_type', key: 'ligand_type' },
              { title: `Binding (${unit})`, dataIndex: 'binding_energy_kcal', key: 'binding', render: (v: number) => format(v), sorter: (a: any, b: any) => a.binding_energy_kcal - b.binding_energy_kcal }
            ]}
            dataSource={data.last_layer_statistics.details.map((d, i) => ({ ...d, key: i }))}
            pagination={{ pageSize: 10 }}
            size="small"
          />
        </Card>
      )}
    </div>
  );
};

export default BindingEnergyView;

