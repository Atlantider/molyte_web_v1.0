/**
 * Cluster QC 统计面板
 * 展示 HOMO/LUMO/gap 分布和简化的电化学窗口估计
 */
import React, { useState, useEffect, useMemo } from 'react';
import {
  Card,
  Row,
  Col,
  Statistic,
  Space,
  Spin,
  Empty,
  Alert,
  Tooltip,
  Typography,
  Switch,
  Tag,
  Divider,
  theme,
} from 'antd';
import {
  ThunderboltOutlined,
  ExperimentOutlined,
  WarningOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons';
import ReactECharts from 'echarts-for-react';
import { getClusterStatistics, ClusterStatisticsResponse, OrbitalStatistics } from '../api/qc';

const { Text, Title, Paragraph } = Typography;

interface ClusterStatisticsPanelProps {
  mdJobId: number;
}

const ClusterStatisticsPanel: React.FC<ClusterStatisticsPanelProps> = ({ mdJobId }) => {
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(true);
  const [data, setData] = useState<ClusterStatisticsResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [includeSingleMolecule, setIncludeSingleMolecule] = useState(false);

  useEffect(() => {
    loadData();
  }, [mdJobId, includeSingleMolecule]);

  const loadData = async () => {
    setLoading(true);
    setError(null);
    try {
      const result = await getClusterStatistics({
        md_job_id: mdJobId,
        include_single_molecule: includeSingleMolecule,
      });
      setData(result);
    } catch (err: any) {
      setError(err.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  // 直方图配置生成
  const createHistogramOption = (stats: OrbitalStatistics | null, title: string, color: string) => {
    if (!stats || !stats.values.length) return {};
    
    const values = stats.values;
    const binCount = Math.min(20, Math.ceil(Math.sqrt(values.length)));
    const min = Math.min(...values);
    const max = Math.max(...values);
    const binWidth = (max - min) / binCount || 1;
    
    const bins = Array(binCount).fill(0);
    const binLabels = [];
    for (let i = 0; i < binCount; i++) {
      binLabels.push((min + binWidth * (i + 0.5)).toFixed(2));
    }
    values.forEach(v => {
      const idx = Math.min(Math.floor((v - min) / binWidth), binCount - 1);
      bins[idx]++;
    });
    
    return {
      title: { text: title, left: 'center', textStyle: { fontSize: 14 } },
      tooltip: {
        trigger: 'axis',
        formatter: (params: any) => `${params[0].name} eV: ${params[0].value} 个`
      },
      xAxis: { type: 'category', data: binLabels, name: 'eV' },
      yAxis: { type: 'value', name: '频数' },
      series: [{
        type: 'bar',
        data: bins,
        itemStyle: { color },
        barWidth: '80%'
      }],
      grid: { left: 50, right: 20, bottom: 40, top: 40 }
    };
  };

  // HOMO 直方图
  const homoHistOption = useMemo(() => 
    createHistogramOption(data?.homo_statistics || null, 'HOMO 能量分布', token.colorPrimary),
    [data, token]
  );

  // LUMO 直方图
  const lumoHistOption = useMemo(() => 
    createHistogramOption(data?.lumo_statistics || null, 'LUMO 能量分布', token.colorSuccess),
    [data, token]
  );

  // Gap 直方图
  const gapHistOption = useMemo(() => 
    createHistogramOption(data?.gap_statistics || null, 'HOMO-LUMO Gap 分布', token.colorWarning),
    [data, token]
  );

  if (loading) {
    return <Spin tip="加载 Cluster 统计数据..." style={{ display: 'block', margin: '40px auto' }} />;
  }

  if (error) {
    return <Alert type="error" message="加载失败" description={error} />;
  }

  if (!data || data.total_qc_jobs === 0) {
    return <Empty description={data?.message || "暂无 Cluster QC 计算结果"} />;
  }

  const ecWindow = data.electrochemical_window_estimate;

  return (
    <div>
      {/* 开关 */}
      <div style={{ marginBottom: 16, textAlign: 'right' }}>
        <Space>
          <Text>包含单分子 QC 结果:</Text>
          <Switch checked={includeSingleMolecule} onChange={setIncludeSingleMolecule} />
        </Space>
      </div>

      {/* 统计卡片 */}
      <Row gutter={16} style={{ marginBottom: 16 }}>
        <Col span={6}>
          <Card size="small">
            <Statistic title="已完成 QC 任务" value={data.total_qc_jobs} prefix={<ExperimentOutlined />} />
          </Card>
        </Col>
        <Col span={6}>
          <Card size="small">
            <Statistic title="有轨道数据" value={data.total_with_orbital_data || 0} prefix={<ThunderboltOutlined />} />
          </Card>
        </Col>
        <Col span={6}>
          <Card size="small">
            {data.homo_statistics && (
              <Statistic
                title={<Tooltip title="最高占据分子轨道能量均值"><span>HOMO 均值 <InfoCircleOutlined /></span></Tooltip>}
                value={data.homo_statistics.mean}
                precision={2}
                suffix="eV"
                valueStyle={{ color: token.colorPrimary }}
              />
            )}
          </Card>
        </Col>
        <Col span={6}>
          <Card size="small">
            {data.gap_statistics && (
              <Statistic
                title={<Tooltip title="HOMO-LUMO 能隙均值"><span>Gap 均值 <InfoCircleOutlined /></span></Tooltip>}
                value={data.gap_statistics.mean}
                precision={2}
                suffix="eV"
                valueStyle={{ color: token.colorWarning }}
              />
            )}
          </Card>
        </Col>
      </Row>

      {/* 电化学窗口估计 */}
      {ecWindow && (
        <Alert
          type="info"
          icon={<WarningOutlined />}
          style={{ marginBottom: 16 }}
          message="简化的电化学窗口估计（基于 HOMO/LUMO）"
          description={
            <div>
              <Row gutter={16}>
                <Col span={8}>
                  <Text strong>氧化极限：</Text> {ecWindow.oxidation_limit_ev.toFixed(2)} eV
                </Col>
                <Col span={8}>
                  <Text strong>还原极限：</Text> {ecWindow.reduction_limit_ev.toFixed(2)} eV
                </Col>
                <Col span={8}>
                  <Text strong>窗口宽度：</Text> {ecWindow.window_ev.toFixed(2)} eV
                </Col>
              </Row>
              <Paragraph type="secondary" style={{ marginTop: 8, marginBottom: 0, fontSize: 12 }}>
                ⚠️ {ecWindow.note}
              </Paragraph>
            </div>
          }
        />
      )}

      {/* 直方图 */}
      <Row gutter={16} style={{ marginBottom: 16 }}>
        <Col span={8}>
          <Card size="small">
            {data.homo_statistics ? (
              <ReactECharts option={homoHistOption} style={{ height: 250 }} />
            ) : <Empty description="无 HOMO 数据" />}
          </Card>
        </Col>
        <Col span={8}>
          <Card size="small">
            {data.lumo_statistics ? (
              <ReactECharts option={lumoHistOption} style={{ height: 250 }} />
            ) : <Empty description="无 LUMO 数据" />}
          </Card>
        </Col>
        <Col span={8}>
          <Card size="small">
            {data.gap_statistics ? (
              <ReactECharts option={gapHistOption} style={{ height: 250 }} />
            ) : <Empty description="无 Gap 数据" />}
          </Card>
        </Col>
      </Row>

      {/* 分组对比：含 Li vs 不含 Li */}
      {(data.with_li_statistics || data.without_li_statistics) && (
        <Card size="small" title="含 Li vs 不含 Li 对比" style={{ marginBottom: 16 }}>
          <Row gutter={16}>
            <Col span={12}>
              <Title level={5}>含 Li 的 Cluster</Title>
              {data.with_li_statistics?.homo ? (
                <Space direction="vertical" size={4}>
                  <Text>HOMO: {data.with_li_statistics.homo.mean.toFixed(2)} ± {data.with_li_statistics.homo.std.toFixed(2)} eV ({data.with_li_statistics.homo.count} 个)</Text>
                  <Text>LUMO: {data.with_li_statistics.lumo?.mean.toFixed(2)} ± {data.with_li_statistics.lumo?.std.toFixed(2)} eV</Text>
                  <Text>Gap: {data.with_li_statistics.gap?.mean.toFixed(2)} ± {data.with_li_statistics.gap?.std.toFixed(2)} eV</Text>
                </Space>
              ) : <Text type="secondary">暂无数据</Text>}
            </Col>
            <Col span={12}>
              <Title level={5}>不含 Li 的分子</Title>
              {data.without_li_statistics?.homo ? (
                <Space direction="vertical" size={4}>
                  <Text>HOMO: {data.without_li_statistics.homo.mean.toFixed(2)} ± {data.without_li_statistics.homo.std.toFixed(2)} eV ({data.without_li_statistics.homo.count} 个)</Text>
                  <Text>LUMO: {data.without_li_statistics.lumo?.mean.toFixed(2)} ± {data.without_li_statistics.lumo?.std.toFixed(2)} eV</Text>
                  <Text>Gap: {data.without_li_statistics.gap?.mean.toFixed(2)} ± {data.without_li_statistics.gap?.std.toFixed(2)} eV</Text>
                </Space>
              ) : <Text type="secondary">暂无数据</Text>}
            </Col>
          </Row>
        </Card>
      )}

      {/* 按类型统计 */}
      {data.per_type_statistics && Object.keys(data.per_type_statistics).length > 0 && (
        <Card size="small" title="按分子类型统计">
          <Row gutter={[16, 16]}>
            {Object.entries(data.per_type_statistics).map(([type, stats]) => (
              <Col span={8} key={type}>
                <Card size="small" style={{ background: token.colorBgLayout }}>
                  <Tag color="blue" style={{ marginBottom: 8 }}>{type}</Tag>
                  <div>样本数: {stats.count}</div>
                  {stats.homo && <div>HOMO: {stats.homo.mean.toFixed(2)} eV</div>}
                  {stats.lumo && <div>LUMO: {stats.lumo.mean.toFixed(2)} eV</div>}
                  {stats.gap && <div>Gap: {stats.gap.mean.toFixed(2)} eV</div>}
                </Card>
              </Col>
            ))}
          </Row>
        </Card>
      )}
    </div>
  );
};

export default ClusterStatisticsPanel;

