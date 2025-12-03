/**
 * 溶剂化结构分析组件 - Nature 期刊风格（仪表盘布局）
 * 使用 ECharts 实现专业的科学论文级别图表
 * 支持 3D 分子结构可视化
 *
 * 布局:
 * - 上方: 统计信息卡片 + 参数设置
 * - 第二行: 三个饼图（配位数分布、阴离子配位数分布、溶剂壳组成）
 * - 第三行: 左侧溶液3D结构（带帧滑块）| 中间溶剂化列表 | 右侧选中的溶剂化结构3D
 */
import { useState, useEffect, useRef, useCallback } from 'react';
import {
  Card,
  Space,
  message,
  Spin,
  Typography,
  Tag,
  Button,
  Table,
  Tooltip,
  Slider,
  Dropdown,
  theme,
} from 'antd';
import {
  ReloadOutlined,
  DownloadOutlined,
  ExperimentOutlined,
  FileTextOutlined,
  LeftOutlined,
  RightOutlined,
  FullscreenOutlined,
  PieChartOutlined,
  TableOutlined,
  NumberOutlined,
  ApartmentOutlined,
  AppstoreOutlined,
} from '@ant-design/icons';

import ReactECharts from 'echarts-for-react';
import type { EChartsOption } from 'echarts';
import {
  getSolvationStructures,
  refreshSolvationStructures,
  getSolvationStatistics,
  getSolvationStructureContent,
  getSystemStructure,
  getFrameCount,
  exportSolvationData,
  type SolvationStructure,
  type SolvationStatistics,
  type SystemStructure,
  type SolvationStructureContent,
} from '../api/jobs';
import DesolvationAnalysisPanel from './DesolvationAnalysisPanel';
import { useThemeStore } from '../stores/themeStore';

const { Text } = Typography;

// Dashboard 样式常量
const DASHBOARD_STYLES = {
  cardBorderRadius: 12,
  cardPadding: 24,
  gutter: 24,
  titleFontSize: 16,
  titleFontWeight: 600,
  chartHeight: 280,
};

// 响应式CSS样式
const RESPONSIVE_STYLES = `
  .dashboard-card:hover {
    box-shadow: 0 8px 24px rgba(15, 23, 42, 0.12);
    transform: translateY(-2px);
  }
  .dashboard-card {
    transition: all 0.3s ease;
  }
  .stats-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 20px;
  }
  .charts-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 20px;
  }
  /* 三列布局：左侧小卡片、中间列表、右侧大卡片 */
  .structure-grid {
    display: grid;
    grid-template-columns: 340px 1fr 1fr;
    gap: 20px;
    align-items: stretch;
  }
  .structure-grid > .structure-card-left {
    min-height: 420px;
  }
  .structure-grid > .structure-card-center {
    min-height: 420px;
  }
  .structure-grid > .structure-card-right {
    min-height: 420px;
  }
  @media (max-width: 1400px) {
    .stats-grid {
      grid-template-columns: repeat(2, 1fr);
    }
    .charts-grid {
      grid-template-columns: repeat(2, 1fr);
    }
    .structure-grid {
      grid-template-columns: 300px 1fr 1fr;
    }
  }
  @media (max-width: 1200px) {
    .structure-grid {
      grid-template-columns: 1fr 1fr;
    }
    .structure-grid > .structure-card-left {
      grid-column: 1 / -1;
    }
  }
  @media (max-width: 992px) {
    .stats-grid {
      grid-template-columns: repeat(2, 1fr);
    }
    .charts-grid {
      grid-template-columns: 1fr;
    }
    .structure-grid {
      grid-template-columns: 1fr;
    }
    .structure-grid > .structure-card-left,
    .structure-grid > .structure-card-center,
    .structure-grid > .structure-card-right {
      min-height: auto;
    }
  }
  @media (max-width: 768px) {
    .stats-grid {
      grid-template-columns: 1fr;
    }
  }
  /* 表格行样式 */
  .solvation-table .ant-table-row {
    cursor: pointer;
    transition: all 0.2s ease;
  }
  .solvation-table .ant-table-row:hover {
    background: rgba(24, 144, 255, 0.06) !important;
  }
  .solvation-table .ant-table-row-selected {
    background: rgba(114, 46, 209, 0.08) !important;
    border-left: 3px solid #722ed1;
  }
  .solvation-table .ant-table-row-selected:hover {
    background: rgba(114, 46, 209, 0.12) !important;
  }
  .solvation-table .ant-table-row-selected td:first-child {
    font-weight: 600;
    color: #722ed1;
  }
  .solvation-table .ant-table-thead > tr > th {
    background: #fafafa;
    font-weight: 600;
    font-size: 12px;
    padding: 8px 12px;
    position: sticky;
    top: 0;
    z-index: 1;
  }
  .solvation-table .ant-table-tbody > tr > td {
    padding: 6px 12px;
    font-size: 12px;
  }
  .solvation-table .ant-table-body {
    max-height: 320px;
    overflow-y: auto;
  }
  /* 3D查看器容器 */
  .viewer-container {
    border: 1px solid #e8e8e8;
    border-radius: 8px;
    background: linear-gradient(135deg, #f8f9fa 0%, #f0f2f5 100%);
    overflow: hidden;
    box-shadow: inset 0 2px 4px rgba(0, 0, 0, 0.06);
    position: relative;
  }
  /* 3D结构卡片增强样式 */
  .structure-card-left,
  .structure-card-right,
  .structure-card-center {
    border: 1px solid #e8e8e8 !important;
    box-shadow: 0 4px 16px rgba(0, 0, 0, 0.08) !important;
  }
  .structure-card-left:hover,
  .structure-card-right:hover,
  .structure-card-center:hover {
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12) !important;
    border-color: #d9d9d9 !important;
  }
  /* 左右3D卡片头部样式增强 */
  .structure-card-left .ant-card-head,
  .structure-card-right .ant-card-head,
  .structure-card-center .ant-card-head {
    min-height: 52px;
    padding: 0 16px;
    border-bottom: 1px solid #e8e8e8;
    background: linear-gradient(180deg, #fafafa 0%, #f5f5f5 100%);
    flex-shrink: 0;
  }
  .structure-card-left .ant-card-head-title,
  .structure-card-right .ant-card-head-title,
  .structure-card-center .ant-card-head-title {
    padding: 14px 0;
    overflow: visible;
  }
  .structure-card-left .ant-card-head-wrapper,
  .structure-card-right .ant-card-head-wrapper {
    display: flex;
    align-items: center;
  }
  .structure-card-left .ant-card-extra,
  .structure-card-right .ant-card-extra {
    padding: 14px 0;
  }
  /* 卡片紧凑布局 */
  .structure-card-left .ant-card-body,
  .structure-card-right .ant-card-body {
    display: flex;
    flex-direction: column;
    padding: 16px !important;
  }
  .structure-card-center .ant-card-body {
    padding: 12px 16px !important;
  }
  /* 选中行高亮 - 使用 data-row-key 匹配 */
  .solvation-table .ant-table-row.row-selected {
    background: linear-gradient(90deg, rgba(114, 46, 209, 0.12) 0%, rgba(114, 46, 209, 0.06) 100%) !important;
    border-left: 3px solid #722ed1 !important;
  }
  .solvation-table .ant-table-row.row-selected:hover {
    background: linear-gradient(90deg, rgba(114, 46, 209, 0.16) 0%, rgba(114, 46, 209, 0.08) 100%) !important;
  }
  .solvation-table .ant-table-row.row-selected td:first-child {
    font-weight: 600;
    color: #722ed1;
  }
  /* 图表卡片增强样式 */
  .chart-card {
    border: 1px solid #e8e8e8 !important;
    box-shadow: 0 4px 16px rgba(0, 0, 0, 0.08) !important;
    transition: all 0.3s ease !important;
  }
  .chart-card:hover {
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12) !important;
    border-color: #d9d9d9 !important;
  }
  .chart-card .ant-card-head {
    min-height: 52px;
    padding: 0 16px;
    border-bottom: 1px solid #e8e8e8;
    background: linear-gradient(180deg, #fafafa 0%, #f5f5f5 100%);
  }
  .chart-card .ant-card-head-title {
    padding: 14px 0;
  }
  .chart-card .ant-card-body {
    padding: 16px !important;
  }
  /* 图表容器样式 */
  .chart-container {
    background: linear-gradient(135deg, #fafbfc 0%, #f0f2f5 100%);
    border: 1px solid #e8e8e8;
    border-radius: 8px;
    padding: 8px;
  }
`;

// 3Dmol.js 类型声明
declare global {
  interface Window {
    $3Dmol: any;
  }
}

interface SolvationStructureProps {
  jobId: number;
}

// Nature 期刊风格配色（低饱和度）
const NATURE_COLORS = [
  '#1f77b4', // 深蓝
  '#2ca02c', // 青绿
  '#ff7f0e', // 暗橙
  '#9467bd', // 紫色
  '#8c564b', // 棕色
  '#e377c2', // 粉色
  '#7f7f7f', // 灰色
  '#bcbd22', // 黄绿
  '#17becf', // 青色
];

export default function SolvationStructureNature({ jobId }: SolvationStructureProps) {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const [structures, setStructures] = useState<SolvationStructure[]>([]);
  const [statistics, setStatistics] = useState<SolvationStatistics | null>(null);
  const [loading, setLoading] = useState(true);
  const [refreshing, setRefreshing] = useState(false);
  const [cutoff, setCutoff] = useState(3.0);
  const [pendingCutoff, setPendingCutoff] = useState(3.0);

  // 3D 可视化状态
  const [systemStructure, setSystemStructure] = useState<SystemStructure | null>(null);
  const [currentFrame, setCurrentFrame] = useState(-1);
  const [totalFrames, setTotalFrames] = useState(0);
  const [loadingSystem, setLoadingSystem] = useState(false);

  // 当前选中的溶剂化结构ID（用于右侧面板显示和列表高亮）
  const [selectedStructureId, setSelectedStructureId] = useState<number | null>(null);
  const [sideStructureContent, setSideStructureContent] = useState<SolvationStructureContent | null>(null);
  const [loadingSideStructure, setLoadingSideStructure] = useState(false);

  const cnChartRef = useRef<any>(null);
  const pieChartRef = useRef<any>(null);
  const anionChartRef = useRef<any>(null);
  const ionPairChartRef = useRef<any>(null);
  const systemViewerRef = useRef<HTMLDivElement>(null);
  const systemViewerInstance = useRef<any>(null);
  const sideViewerRef = useRef<HTMLDivElement>(null);
  const sideViewerInstance = useRef<any>(null);

  // 加载 3Dmol.js
  useEffect(() => {
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
      script.async = true;
      document.body.appendChild(script);
    }
  }, []);

  // 加载溶剂化结构数据
  const loadData = useCallback(async () => {
    setLoading(true);
    try {
      const [structuresData, statsData] = await Promise.all([
        getSolvationStructures(jobId),
        getSolvationStatistics(jobId),
      ]);
      setStructures(structuresData);
      setStatistics(statsData);

      // 获取帧数并自动加载最后一帧
      try {
        const frameData = await getFrameCount(jobId);
        setTotalFrames(frameData.frame_count);
        if (frameData.frame_count > 0) {
          const lastFrame = frameData.frame_count - 1;
          setCurrentFrame(lastFrame);
          // 自动加载最后一帧的系统结构
          loadSystemStructure(lastFrame);
        }
      } catch (e) {
        console.warn('Failed to get frame count:', e);
      }

      // 自动加载第一个溶剂化结构用于侧边显示
      if (structuresData.length > 0) {
        setSelectedStructureId(structuresData[0].id);
        loadSideStructure(structuresData[0].id);
      }
    } catch (error: any) {
      console.error('Failed to load solvation data:', error);
      // 如果没有数据，自动开始分析
      if (error.response?.status === 404) {
        handleRefresh();
      }
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  // 刷新/重新计算溶剂化结构
  const handleRefresh = async () => {
    setRefreshing(true);
    try {
      const result = await refreshSolvationStructures(jobId, pendingCutoff);
      if (result.success) {
        message.success(result.message);
        setCutoff(pendingCutoff);
        await loadData();
      } else {
        message.warning(result.message);
      }
    } catch (error: any) {
      console.error('Failed to refresh solvation:', error);
      message.error(error.response?.data?.detail || '溶剂化结构分析失败');
    } finally {
      setRefreshing(false);
    }
  };

  // 加载体系结构
  const loadSystemStructure = async (frame: number) => {
    setLoadingSystem(true);
    try {
      const data = await getSystemStructure(jobId, frame);
      setSystemStructure(data);
      setCurrentFrame(data.frame_index);
      setTotalFrames(data.total_frames);
    } catch (error: any) {
      console.error('Failed to load system structure:', error);
      message.error('加载体系结构失败');
    } finally {
      setLoadingSystem(false);
    }
  };

  // 加载侧边显示的溶剂化结构
  const loadSideStructure = async (structureId: number) => {
    setLoadingSideStructure(true);
    try {
      const content = await getSolvationStructureContent(jobId, structureId);
      setSideStructureContent(content);
    } catch (error: any) {
      console.error('Failed to load side structure:', error);
    } finally {
      setLoadingSideStructure(false);
    }
  };

  // 获取当前选中结构的索引
  const selectedStructureIndex = structures.findIndex(s => s.id === selectedStructureId);

  // 切换到上一个/下一个溶剂化结构
  const handlePrevStructure = () => {
    if (selectedStructureIndex > 0 && structures.length > 0) {
      const prevStructure = structures[selectedStructureIndex - 1];
      setSelectedStructureId(prevStructure.id);
      loadSideStructure(prevStructure.id);
    }
  };

  const handleNextStructure = () => {
    if (selectedStructureIndex >= 0 && selectedStructureIndex < structures.length - 1) {
      const nextStructure = structures[selectedStructureIndex + 1];
      setSelectedStructureId(nextStructure.id);
      loadSideStructure(nextStructure.id);
    }
  };

  // 渲染体系 3D 结构
  useEffect(() => {
    if (!systemStructure?.xyz_content || !systemViewerRef.current || !window.$3Dmol) return;

    // 确保容器有尺寸后再创建 viewer
    const container = systemViewerRef.current;
    if (container.clientWidth === 0 || container.clientHeight === 0) {
      // 延迟重试
      const timer = setTimeout(() => {
        if (systemViewerRef.current) {
          systemViewerRef.current.dispatchEvent(new Event('resize'));
        }
      }, 100);
      return () => clearTimeout(timer);
    }

    if (systemViewerInstance.current) {
      systemViewerInstance.current.clear();
      systemViewerInstance.current = null;
    }

    const viewer = window.$3Dmol.createViewer(container, {
      backgroundColor: '#f8f9fa',
    });
    systemViewerInstance.current = viewer;

    viewer.addModel(systemStructure.xyz_content, 'xyz');
    viewer.setStyle({}, {
      sphere: { scale: 0.25, colorscheme: 'Jmol' },
    });
    // 自动居中并缩放，留一些边距（0.9倍缩放）
    viewer.zoomTo(0.9);
    // 设置略微俯视的视角（绕x轴旋转-20度，绕y轴旋转15度）
    viewer.rotate(-20, 'x');
    viewer.rotate(15, 'y');
    viewer.render();

    // 处理窗口大小变化
    const handleResize = () => {
      if (viewer) {
        viewer.resize();
        viewer.render();
      }
    };
    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, [systemStructure]);

  // 渲染侧边溶剂化结构 3D
  useEffect(() => {
    if (!sideStructureContent?.xyz_content || !sideViewerRef.current || !window.$3Dmol) return;

    const container = sideViewerRef.current;
    if (container.clientWidth === 0 || container.clientHeight === 0) {
      const timer = setTimeout(() => {
        if (sideViewerRef.current) {
          sideViewerRef.current.dispatchEvent(new Event('resize'));
        }
      }, 100);
      return () => clearTimeout(timer);
    }

    if (sideViewerInstance.current) {
      sideViewerInstance.current.clear();
      sideViewerInstance.current = null;
    }

    const viewer = window.$3Dmol.createViewer(container, {
      backgroundColor: '#f8f9fa',
    });
    sideViewerInstance.current = viewer;

    viewer.addModel(sideStructureContent.xyz_content, 'xyz');

    // 设置所有原子的默认样式（球棍模型）
    viewer.setStyle({}, {
      stick: { radius: 0.15, colorscheme: 'Jmol' },
      sphere: { scale: 0.3, colorscheme: 'Jmol' },
    });

    // 高亮中心离子（第一个原子，serial 从 0 开始）
    viewer.setStyle({ serial: 0 }, {
      sphere: { scale: 0.5, color: '#e74c3c' },
    });

    // 增强氢原子的可见性
    viewer.setStyle({ elem: 'H' }, {
      stick: { radius: 0.12, colorscheme: 'Jmol' },
      sphere: { scale: 0.25, colorscheme: 'Jmol' },
    });

    // 自动居中并缩放，留一些边距
    viewer.zoomTo(0.85);
    // 设置略微俯视的视角
    viewer.rotate(-15, 'x');
    viewer.rotate(10, 'y');
    viewer.render();

    const handleResize = () => {
      if (viewer) {
        viewer.resize();
        viewer.render();
      }
    };
    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, [sideStructureContent]);

  // 导出图片
  const exportChart = (chartRef: any, filename: string) => {
    if (chartRef.current) {
      const echartInstance = chartRef.current.getEchartsInstance();
      const url = echartInstance.getDataURL({
        type: 'png',
        pixelRatio: 3,
        backgroundColor: '#fff',
      });
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      link.click();
    }
  };

  // 导出数据
  const handleExportData = async (format: 'json' | 'csv') => {
    try {
      const data = await exportSolvationData(jobId, format);
      if (format === 'csv') {
        const blob = new Blob([data], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = `solvation_job${jobId}.csv`;
        link.click();
        URL.revokeObjectURL(url);
      } else {
        const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = `solvation_job${jobId}.json`;
        link.click();
        URL.revokeObjectURL(url);
      }
      message.success(`数据已导出为 ${format.toUpperCase()} 格式`);
    } catch (error) {
      message.error('导出数据失败');
    }
  };

  useEffect(() => {
    loadData();
  }, [loadData]);

  if (loading && structures.length === 0) {
    return (
      <div style={{ textAlign: 'center', padding: '50px' }}>
        <Spin size="large" tip="加载溶剂化结构数据..." />
      </div>
    );
  }

  // 准备配位数分布数据
  const cnDistribution = statistics?.coordination_distribution || {};
  const cnData = Object.entries(cnDistribution)
    .map(([cn, count]) => ({ cn: parseInt(cn), count: count as number }))
    .sort((a, b) => a.cn - b.cn);

  // 准备分子组成饼图数据
  const moleculeCounts = statistics?.molecule_counts || {};
  const pieData = Object.entries(moleculeCounts)
    .filter(([_, count]) => (count as number) > 0)
    .map(([name, count], index) => ({
      name,
      value: count as number,
      itemStyle: { color: NATURE_COLORS[index % NATURE_COLORS.length] },
    }));

  // 准备阴离子配位数分布数据
  const anionCnDistribution = statistics?.anion_coordination_distribution || {};
  const anionCnData = Object.entries(anionCnDistribution)
    .map(([cn, count]) => ({ cn: parseInt(cn), count: count as number }))
    .sort((a, b) => a.cn - b.cn);

  // 计算离子对分类统计 (基于 n_contact = 阴离子配位数)
  // Free Ion: Anion = 0 (实际上需要 n_shell2 = 0，但目前无法区分，暂视为 Anion = 0 中的一部分)
  // SSIP (Solvent-Separated Ion Pair): Anion = 0 且有溶剂分离的阴离子 (目前无法区分，暂计入 Anion = 0)
  // CIP (Contact Ion Pair): Anion = 1
  // AGG (Aggregates): Anion >= 2
  const ionPairStats = {
    freeIon: (anionCnDistribution['0'] as number) || 0,  // n_contact = 0 (含 Free + SSIP)
    cip: (anionCnDistribution['1'] as number) || 0,       // n_contact = 1
    agg: Object.entries(anionCnDistribution)
      .filter(([cn]) => parseInt(cn) >= 2)
      .reduce((sum, [_, count]) => sum + (count as number), 0), // n_contact >= 2
  };
  const totalIonPairs = ionPairStats.freeIon + ionPairStats.cip + ionPairStats.agg;

  // Nature 期刊统一样式常量
  const CHART_FONT_FAMILY = 'Arial, Helvetica, sans-serif';
  const CHART_TITLE_SIZE = 13;
  const CHART_LABEL_SIZE = 11;
  const CHART_LEGEND_SIZE = 11;

  // 配位数分布饼图数据
  const cnPieData = cnData.map((d, index) => ({
    name: `CN = ${d.cn}`,
    value: d.count,
    itemStyle: { color: NATURE_COLORS[index % NATURE_COLORS.length] },
  }));

  // 配位数分布饼图配置 (Nature 风格)
  const cnPieOption: EChartsOption = {
    tooltip: {
      trigger: 'item',
      backgroundColor: 'rgba(255, 255, 255, 0.96)',
      borderColor: '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
      formatter: '{b}: {c} ({d}%)',
    },
    legend: {
      orient: 'horizontal',
      bottom: 8,
      itemWidth: 12,
      itemHeight: 12,
      itemGap: 16,
      textStyle: {
        fontSize: CHART_LEGEND_SIZE,
        fontFamily: CHART_FONT_FAMILY,
        color: '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['28%', '58%'],
        center: ['50%', '42%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: '#fff',
          borderWidth: 2,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: '#333',
          formatter: '{d}%',
        },
        labelLine: {
          show: true,
          length: 8,
          length2: 6,
        },
        emphasis: {
          label: {
            show: true,
            fontSize: CHART_LABEL_SIZE + 1,
            fontWeight: 600,
          },
        },
        data: cnPieData,
      },
    ],
  };

  // 溶剂壳组成饼图配置 (Nature 风格)
  const pieOption: EChartsOption = {
    tooltip: {
      trigger: 'item',
      backgroundColor: 'rgba(255, 255, 255, 0.96)',
      borderColor: '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
      formatter: '{b}: {c} ({d}%)',
    },
    legend: {
      orient: 'horizontal',
      bottom: 8,
      itemWidth: 12,
      itemHeight: 12,
      itemGap: 16,
      textStyle: {
        fontSize: CHART_LEGEND_SIZE,
        fontFamily: CHART_FONT_FAMILY,
        color: '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['28%', '58%'],
        center: ['50%', '42%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: '#fff',
          borderWidth: 2,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: '#333',
          formatter: '{d}%',
        },
        labelLine: {
          show: true,
          length: 8,
          length2: 6,
        },
        emphasis: {
          label: {
            show: true,
            fontSize: CHART_LABEL_SIZE + 1,
            fontWeight: 600,
          },
        },
        data: pieData,
      },
    ],
  };

  // 阴离子配位数分布饼图数据（使用不同色系）- 使用英文标签 Anion = n
  const ANION_COLORS = ['#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#1f77b4', '#ff7f0e'];
  const anionCnPieData = anionCnData.map((d, index) => ({
    name: `Anion = ${d.cn}`,
    value: d.count,
    itemStyle: { color: ANION_COLORS[index % ANION_COLORS.length] },
  }));

  // 阴离子配位数分布饼图配置 (Nature 风格)
  const anionCnPieOption: EChartsOption = {
    tooltip: {
      trigger: 'item',
      backgroundColor: 'rgba(255, 255, 255, 0.96)',
      borderColor: '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
      formatter: '{b}: {c} ({d}%)',
    },
    legend: {
      orient: 'horizontal',
      bottom: 8,
      itemWidth: 12,
      itemHeight: 12,
      itemGap: 16,
      textStyle: {
        fontSize: CHART_LEGEND_SIZE,
        fontFamily: CHART_FONT_FAMILY,
        color: '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['28%', '58%'],
        center: ['50%', '42%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: '#fff',
          borderWidth: 2,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: '#333',
          formatter: '{d}%',
        },
        labelLine: {
          show: true,
          length: 8,
          length2: 6,
        },
        emphasis: {
          label: {
            show: true,
            fontSize: CHART_LABEL_SIZE + 1,
            fontWeight: 600,
          },
        },
        data: anionCnPieData,
      },
    ],
  };

  // 离子对分类饼图数据 (CIP/AGG/Free Ion)
  const ION_PAIR_COLORS = {
    'Free Ion / SSIP': '#52c41a',  // 绿色 - 自由离子/溶剂分离离子对
    'CIP': '#1890ff',               // 蓝色 - 接触离子对
    'AGG': '#ff4d4f',               // 红色 - 聚集体
  };
  const ionPairPieData = [
    { name: 'Free Ion / SSIP', value: ionPairStats.freeIon, itemStyle: { color: ION_PAIR_COLORS['Free Ion / SSIP'] } },
    { name: 'CIP', value: ionPairStats.cip, itemStyle: { color: ION_PAIR_COLORS['CIP'] } },
    { name: 'AGG', value: ionPairStats.agg, itemStyle: { color: ION_PAIR_COLORS['AGG'] } },
  ].filter(d => d.value > 0);

  // 离子对分类饼图配置
  const ionPairPieOption: EChartsOption = {
    tooltip: {
      trigger: 'item',
      backgroundColor: 'rgba(255, 255, 255, 0.96)',
      borderColor: '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
      formatter: (params: any) => {
        const percent = totalIonPairs > 0 ? ((params.value / totalIonPairs) * 100).toFixed(1) : '0';
        return `${params.name}: ${params.value} (${percent}%)`;
      },
    },
    legend: {
      orient: 'horizontal',
      bottom: 8,
      itemWidth: 12,
      itemHeight: 12,
      itemGap: 16,
      textStyle: {
        fontSize: CHART_LEGEND_SIZE,
        fontFamily: CHART_FONT_FAMILY,
        color: '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['28%', '58%'],
        center: ['50%', '42%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: '#fff',
          borderWidth: 2,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: '#333',
          formatter: '{d}%',
        },
        labelLine: {
          show: true,
          length: 8,
          length2: 6,
        },
        emphasis: {
          label: {
            show: true,
            fontSize: CHART_LABEL_SIZE + 1,
            fontWeight: 600,
          },
        },
        data: ionPairPieData,
      },
    ],
  };

  // 表格列定义 - 中文（英文）格式
  const columns = [
    {
      title: '#',
      dataIndex: 'id',
      key: 'id',
      width: 50,
      render: (_: any, __: any, index: number) => index + 1,
    },
    {
      title: '中心离子',
      dataIndex: 'center_ion',
      key: 'center_ion',
      width: 80,
      render: (ion: string) => <Tag color="blue">{ion}⁺</Tag>,
    },
    {
      title: 'CN',
      dataIndex: 'coordination_num',
      key: 'coordination_num',
      width: 60,
      sorter: (a: SolvationStructure, b: SolvationStructure) =>
        a.coordination_num - b.coordination_num,
      render: (cn: number) => (
        <Tag color={cn > 4 ? 'green' : cn > 2 ? 'orange' : 'red'}>{cn}</Tag>
      ),
    },
    {
      title: '壳层组成 (Shell Composition)',
      dataIndex: 'composition',
      key: 'composition',
      render: (comp: Record<string, number>) => (
        <Space size="small" wrap>
          {Object.entries(comp || {})
            .filter(([_, count]) => count > 0)
            .map(([mol, count], idx) => (
              <Tag key={mol} color={NATURE_COLORS[idx % NATURE_COLORS.length]}>
                {mol}: {count}
              </Tag>
            ))}
        </Space>
      ),
    },
    {
      title: '',
      key: 'action',
      width: 40,
      align: 'right' as const,
      render: (_: any, record: SolvationStructure) =>
        record.file_path ? (
          <Tooltip title="下载 XYZ">
            <Button
              type="text"
              size="small"
              icon={<DownloadOutlined style={{ color: '#52c41a' }} />}
              href={`/api/v1/jobs/${jobId}/solvation/structure/${record.id}`}
              target="_blank"
              onClick={(e) => e.stopPropagation()}
            />
          </Tooltip>
        ) : null,
    },
  ];

  // Dashboard 卡片样式
  const dashboardCardStyle: React.CSSProperties = {
    background: token.colorBgContainer,
    borderRadius: DASHBOARD_STYLES.cardBorderRadius,
    boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(15, 23, 42, 0.08)',
    border: `1px solid ${token.colorBorder}`,
    transition: 'all 0.3s ease',
  };

  return (
    <div style={{
      background: token.colorBgLayout,
      padding: DASHBOARD_STYLES.gutter,
      minHeight: '100%',
      borderRadius: 8,
    }}>
      <style>{RESPONSIVE_STYLES}</style>

      {/* 第一部分：统计信息卡片 */}
      <div className="stats-grid" style={{ marginBottom: DASHBOARD_STYLES.gutter }}>
        {/* 总结构数 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #1890ff 0%, #096dd9 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <NumberOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>总结构数 (Total)</div>
                <div style={{ fontSize: 28, fontWeight: 700, color: '#1890ff' }}>
                  {statistics?.total_count || 0}
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* 平均配位数 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #52c41a 0%, #389e0d 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <ApartmentOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>平均配位数 (Avg CN)</div>
                <div style={{ fontSize: 28, fontWeight: 700, color: '#52c41a' }}>
                  {(statistics?.average_coordination_number || 0).toFixed(2)}
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* 结构类型 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
              <div style={{
                width: 48, height: 48, borderRadius: 12,
                background: 'linear-gradient(135deg, #722ed1 0%, #531dab 100%)',
                display: 'flex', alignItems: 'center', justifyContent: 'center',
              }}>
                <AppstoreOutlined style={{ fontSize: 24, color: '#fff' }} />
              </div>
              <div>
                <div style={{ fontSize: 12, color: '#6b7280' }}>结构类型 (Types)</div>
                <div style={{ fontSize: 28, fontWeight: 700, color: '#722ed1' }}>
                  {Object.keys(statistics?.composition_distribution || {}).length}
                </div>
              </div>
            </div>
          </div>
        </div>

        {/* 截断距离 */}
        <div className="dashboard-card" style={dashboardCardStyle}>
          <div style={{ padding: DASHBOARD_STYLES.cardPadding }}>
            <div style={{ fontSize: 12, color: '#6b7280', marginBottom: 8 }}>截断距离 (Cutoff)</div>
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <Slider
                style={{ flex: 1 }}
                min={2.0}
                max={6.0}
                step={0.1}
                value={pendingCutoff}
                onChange={(v) => setPendingCutoff(v)}
                tooltip={{ formatter: (v) => `${v?.toFixed(1)} Å` }}
              />
              <Tag color="blue" style={{ margin: 0 }}>{pendingCutoff.toFixed(1)} Å</Tag>
            </div>
            <div style={{ marginTop: 12, display: 'flex', gap: 8 }}>
              <Button
                type="primary"
                size="small"
                icon={<ReloadOutlined />}
                onClick={handleRefresh}
                loading={refreshing}
                disabled={pendingCutoff === cutoff}
                style={{ borderRadius: 6 }}
              >
                应用 (Apply)
              </Button>
              <Dropdown menu={{ items: [
                { key: 'json', icon: <FileTextOutlined />, label: '导出 JSON', onClick: () => handleExportData('json') },
                { key: 'csv', icon: <FileTextOutlined />, label: '导出 CSV', onClick: () => handleExportData('csv') },
              ]}}>
                <Button size="small" icon={<DownloadOutlined />} style={{ borderRadius: 6 }}>导出 (Export)</Button>
              </Dropdown>
            </div>
          </div>
        </div>
      </div>

      {/* 第二部分：三个饼图 */}
      <div className="charts-grid" style={{ marginBottom: DASHBOARD_STYLES.gutter }}>
        {/* 配位数分布饼图 */}
        <Card
          className="dashboard-card chart-card"
          style={dashboardCardStyle}
          title={
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <PieChartOutlined style={{ color: '#1f77b4', fontSize: 16 }} />
              <span style={{ fontSize: DASHBOARD_STYLES.titleFontSize, fontWeight: DASHBOARD_STYLES.titleFontWeight, color: token.colorText }}>
                配位数分布 (CN Distribution)
              </span>
            </div>
          }
          extra={
            <Button
              icon={<DownloadOutlined />}
              size="small"
              type="text"
              style={{ borderRadius: 6 }}
              onClick={() => exportChart(cnChartRef, `solvation_cn_job${jobId}.png`)}
            />
          }
        >
          <div className="chart-container">
            <ReactECharts
              ref={cnChartRef}
              option={cnPieOption}
              style={{ height: DASHBOARD_STYLES.chartHeight }}
              notMerge={true}
              theme={isDark ? 'dark' : undefined}
            />
          </div>
        </Card>

        {/* 阴离子配位分布 */}
        <Card
          className="dashboard-card chart-card"
          style={dashboardCardStyle}
          title={
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <PieChartOutlined style={{ color: '#9467bd', fontSize: 16 }} />
              <span style={{ fontSize: DASHBOARD_STYLES.titleFontSize, fontWeight: DASHBOARD_STYLES.titleFontWeight, color: token.colorText }}>
                阴离子配位分布 (Anion CN)
              </span>
            </div>
          }
          extra={
            <Button
              icon={<DownloadOutlined />}
              size="small"
              type="text"
              style={{ borderRadius: 6 }}
              onClick={() => exportChart(anionChartRef, `solvation_anion_cn_job${jobId}.png`)}
            />
          }
        >
          <div className="chart-container">
            <ReactECharts
              ref={anionChartRef}
              option={anionCnPieOption}
              style={{ height: DASHBOARD_STYLES.chartHeight }}
              notMerge={true}
              theme={isDark ? 'dark' : undefined}
            />
          </div>
        </Card>

        {/* 溶剂化壳层组成 */}
        <Card
          className="dashboard-card chart-card"
          style={dashboardCardStyle}
          title={
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <PieChartOutlined style={{ color: '#2ca02c', fontSize: 16 }} />
              <span style={{ fontSize: DASHBOARD_STYLES.titleFontSize, fontWeight: DASHBOARD_STYLES.titleFontWeight, color: token.colorText }}>
                溶剂化壳层组成 (Shell Composition)
              </span>
            </div>
          }
          extra={
            <Button
              icon={<DownloadOutlined />}
              size="small"
              type="text"
              style={{ borderRadius: 6 }}
              onClick={() => exportChart(pieChartRef, `solvation_composition_job${jobId}.png`)}
            />
          }
        >
          <div className="chart-container">
            <ReactECharts
              ref={pieChartRef}
              option={pieOption}
              style={{ height: DASHBOARD_STYLES.chartHeight }}
              notMerge={true}
              theme={isDark ? 'dark' : undefined}
            />
          </div>
        </Card>

        {/* 离子对分类 (CIP/AGG/Free Ion) */}
        <Card
          className="dashboard-card chart-card"
          style={dashboardCardStyle}
          title={
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <PieChartOutlined style={{ color: '#fa8c16', fontSize: 16 }} />
              <span style={{ fontSize: DASHBOARD_STYLES.titleFontSize, fontWeight: DASHBOARD_STYLES.titleFontWeight, color: token.colorText }}>
                离子对分类 (Ion Pair Classification)
              </span>
            </div>
          }
          extra={
            <Button
              icon={<DownloadOutlined />}
              size="small"
              type="text"
              style={{ borderRadius: 6 }}
              onClick={() => exportChart(ionPairChartRef, `solvation_ion_pair_job${jobId}.png`)}
            />
          }
        >
          <div className="chart-container">
            <ReactECharts
              ref={ionPairChartRef}
              option={ionPairPieOption}
              style={{ height: DASHBOARD_STYLES.chartHeight }}
              notMerge={true}
              theme={isDark ? 'dark' : undefined}
            />
          </div>
        </Card>
      </div>

      {/* 第三部分：三列布局 - 溶液结构 | 列表 | 溶剂化结构 */}
      <div className="structure-grid" style={{ marginTop: 16 }}>
        {/* 左侧：整体溶液结构（紧凑方卡片） */}
        <Card
          className="dashboard-card structure-card-left"
          style={{ ...dashboardCardStyle, maxWidth: 360, overflow: 'hidden' }}
          title={
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <FullscreenOutlined style={{ color: token.colorPrimary, fontSize: 16 }} />
              <span
                style={{
                  fontSize: 14,
                  fontWeight: 600,
                  color: token.colorText,
                }}
              >
                整体溶液结构 (System)
              </span>
            </div>
          }
        >
          {/* 3D 视窗区域 */}
          {loadingSystem ? (
            <div
              style={{
                height: 280,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                background: isDark ? 'rgba(255,255,255,0.04)' : '#f5f5f5',
                borderRadius: 8,
              }}
            >
              <Spin size="small" tip="加载中..." />
            </div>
          ) : systemStructure ? (
            <>
              <div
                ref={systemViewerRef}
                className="viewer-container"
                style={{
                  width: '100%',
                  height: 280,
                  minHeight: 280,
                  borderRadius: 8,
                  overflow: 'hidden',
                  background: 'linear-gradient(135deg, #f5f7fa 0%, #e8ecf1 100%)',
                  border: '1px solid #e0e0e0',
                }}
              />
              <div style={{ marginTop: 6, textAlign: 'center', lineHeight: 1.4 }}>
                <Text type="secondary" style={{ fontSize: 10 }}>
                  盒子 (Box): {systemStructure.box.map((b) => b.toFixed(1)).join(' × ')} Å |{' '}
                  {systemStructure.atom_count} 原子 | 拖动旋转
                </Text>
              </div>
            </>
          ) : (
            <div
              style={{
                height: 280,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                background: '#f5f5f5',
                borderRadius: 8,
                border: '1px dashed #d9d9d9',
              }}
            >
              <Spin size="small" tip="加载体系..." />
            </div>
          )}

          {/* 帧数滑块 */}
          {totalFrames > 1 && (
            <div
              style={{
                marginTop: 8,
                display: 'flex',
                alignItems: 'center',
                gap: 4,
              }}
            >
              <Button
                size="small"
                type="text"
                icon={<LeftOutlined style={{ fontSize: 9 }} />}
                disabled={currentFrame <= 0}
                onClick={() => loadSystemStructure(currentFrame - 1)}
                style={{ padding: '0 2px', minWidth: 20, height: 20 }}
              />
              <Slider
                style={{ flex: 1, margin: '0 4px' }}
                min={0}
                max={totalFrames - 1}
                value={currentFrame}
                onChange={(v) => loadSystemStructure(v)}
                tooltip={{ formatter: (v) => `帧 (Frame) ${(v || 0) + 1}` }}
              />
              <Button
                size="small"
                type="text"
                icon={<RightOutlined style={{ fontSize: 9 }} />}
                disabled={currentFrame >= totalFrames - 1}
                onClick={() => loadSystemStructure(currentFrame + 1)}
                style={{ padding: '0 2px', minWidth: 20, height: 20 }}
              />
            </div>
          )}
        </Card>

        {/* 中间：溶剂化结构列表 */}
        <Card
          className="dashboard-card structure-card-center"
          style={{ ...dashboardCardStyle, overflow: 'hidden' }}
          title={
            <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <TableOutlined style={{ color: '#fa8c16', fontSize: 16 }} />
              <span
                style={{
                  fontSize: 14,
                  fontWeight: 600,
                  color: token.colorText,
                }}
              >
                溶剂化结构列表 (Structures)
              </span>
              <Tag color="blue" style={{ marginLeft: 4, fontSize: 11 }}>
                {structures.length}
              </Tag>
            </div>
          }
        >
          <Table
            className="solvation-table"
            dataSource={structures}
            columns={columns}
            rowKey="id"
            size="small"
            pagination={{
              pageSize: 8,
              showSizeChanger: false,
              showTotal: (total) => (
                <span style={{ fontSize: 11 }}>共 {total} 个</span>
              ),
              size: 'small',
            }}
            scroll={{ y: 300 }}
            rowClassName={(record) =>
              record.id === selectedStructureId ? 'row-selected' : ''
            }
            onRow={(record) => ({
              onClick: () => {
                setSelectedStructureId(record.id);
                loadSideStructure(record.id);
              },
            })}
          />
        </Card>

        {/* 右侧：当前选中溶剂化结构 */}
        <Card
          className="dashboard-card structure-card-right"
          style={{ ...dashboardCardStyle, overflow: 'hidden' }}
          title={
            <div
              style={{
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'space-between',
                width: '100%',
              }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                <ExperimentOutlined style={{ color: '#722ed1', fontSize: 16 }} />
                <span
                  style={{
                    fontSize: 14,
                    fontWeight: 600,
                    color: token.colorText,
                  }}
                >
                  典型溶剂化结构 (Solvation)
                </span>
              </div>
              {structures.length > 0 && selectedStructureIndex >= 0 && (
                <Tag color="purple" style={{ margin: 0, fontSize: 11 }}>
                  {selectedStructureIndex + 1}/{structures.length}
                </Tag>
              )}
            </div>
          }
          extra={
            structures.length > 0 && (
              <Space size={4}>
                <Button
                  size="small"
                  type="text"
                  icon={<LeftOutlined style={{ fontSize: 10 }} />}
                  disabled={selectedStructureIndex <= 0}
                  onClick={handlePrevStructure}
                  style={{ padding: '2px 4px', minWidth: 24, height: 24 }}
                />
                <Button
                  size="small"
                  type="text"
                  icon={<RightOutlined style={{ fontSize: 10 }} />}
                  disabled={selectedStructureIndex < 0 || selectedStructureIndex >= structures.length - 1}
                  onClick={handleNextStructure}
                  style={{ padding: '2px 4px', minWidth: 24, height: 24 }}
                />
              </Space>
            )
          }
        >
          {loadingSideStructure ? (
            <div
              style={{
                height: 340,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                background: '#f5f5f5',
                borderRadius: 8,
              }}
            >
              <Spin size="small" tip="加载中..." />
            </div>
          ) : sideStructureContent ? (
            <>
              {/* 中心离子信息 - 紧凑显示 */}
              <div style={{ marginBottom: 8, textAlign: 'center' }}>
                <Tag color="purple" style={{ fontSize: 12, padding: '2px 8px' }}>
                  {sideStructureContent.center_ion}⁺ 配位数 (CN) ={' '}
                  {sideStructureContent.coordination_num}
                </Tag>
              </div>
              {/* 3D 结构查看器 */}
              <div
                ref={sideViewerRef}
                className="viewer-container"
                style={{
                  width: '100%',
                  height: 280,
                  minHeight: 280,
                  borderRadius: 8,
                  overflow: 'hidden',
                  background: 'linear-gradient(135deg, #f5f7fa 0%, #e8ecf1 100%)',
                  border: '1px solid #e0e0e0',
                }}
              />
              {/* 溶剂壳组成标签 */}
              <div style={{ marginTop: 10, textAlign: 'center' }}>
                <Space wrap size={4} style={{ justifyContent: 'center' }}>
                  {Object.entries(sideStructureContent.composition || {})
                    .filter(([_, count]) => count > 0)
                    .map(([mol, count], idx) => (
                      <Tag
                        key={mol}
                        color={
                          NATURE_COLORS[idx % NATURE_COLORS.length]
                        }
                        style={{ margin: 0, fontSize: 11, padding: '1px 6px' }}
                      >
                        {mol}: {count}
                      </Tag>
                    ))}
                </Space>
              </div>
              {/* 操作提示 */}
              <div style={{ marginTop: 6, textAlign: 'center' }}>
                <Text type="secondary" style={{ fontSize: 10 }}>
                  红色球体为中心离子 (Red sphere: center ion) | 拖动旋转
                </Text>
              </div>
            </>
          ) : (
            <div
              style={{
                height: 340,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                background: '#f5f5f5',
                borderRadius: 8,
                border: '1px dashed #d9d9d9',
              }}
            >
              <Text type="secondary" style={{ fontSize: 12 }}>
                点击列表查看结构
              </Text>
            </div>
          )}
        </Card>
      </div>

      {/* 去溶剂化能分析面板 */}
      <div style={{ marginTop: DASHBOARD_STYLES.gutter }}>
        <DesolvationAnalysisPanel
          jobId={jobId}
          clusterId={selectedStructureId || undefined}
        />
      </div>

    </div>
  );
}

