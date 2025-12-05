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
  Modal,
  Divider,
  theme,
  Select,
  InputNumber,
  Row,
  Col,
  Collapse,
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
  ThunderboltOutlined,
  BulbOutlined,
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
  autoSelectSolvationStructures,
  type SolvationStructure,
  type SolvationStatistics,
  type SystemStructure,
  type SolvationStructureContent,
  type AutoSelectResponse,
} from '../api/jobs';
import { useThemeStore } from '../stores/themeStore';
import {
  createDesolvationJob,
  getDesolvationJob,
  listClusterDesolvationJobs,
} from '../api/desolvation';
import type { DesolvationJobResponse, SolventModel, SolventConfig } from '../types/desolvation';
import DesolvationResultView from './DesolvationResultView';

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
    gap: 12px;
  }
  .charts-grid {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 12px;
  }
  /* 三列布局：左侧小卡片、中间列表、右侧大卡片 */
  .structure-grid {
    display: grid;
    grid-template-columns: 340px 1fr 1fr;
    gap: 12px;
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
  .solvation-table .ant-table-thead > tr > th {
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
    border-radius: 8px;
    overflow: hidden;
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
    overflow: hidden !important;
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
    padding: 8px !important;
    overflow: hidden !important;
  }
  /* 饼图容器样式 - 使用特定类名避免影响其他组件 */
  .solvation-pie-chart-container {
    background: transparent;
    border: none;
    border-radius: 0;
    padding: 0;
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
  onGoToDesolvation?: () => void;  // 跳转到去溶剂化能计算 Tab（已废弃，保留兼容）
  onGoToPostProcess?: () => void;  // 跳转到后处理分析页面
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

export default function SolvationStructureNature({ jobId, onGoToDesolvation, onGoToPostProcess }: SolvationStructureProps) {
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

  // 去溶剂化能计算状态
  const [desolvationModalVisible, setDesolvationModalVisible] = useState(false);
  const [desolvationJobs, setDesolvationJobs] = useState<DesolvationJobResponse[]>([]);
  const [selectedDesolvationJob, setSelectedDesolvationJob] = useState<DesolvationJobResponse | null>(null);
  const [loadingDesolvation, setLoadingDesolvation] = useState(false);
  const [creatingDesolvation, setCreatingDesolvation] = useState(false);
  const [selectedMethodLevel, setSelectedMethodLevel] = useState<'fast' | 'standard' | 'accurate'>('standard');
  const [selectedDesolvationMode, setSelectedDesolvationMode] = useState<'stepwise' | 'full'>('stepwise');
  const [selectedSolventModel, setSelectedSolventModel] = useState<SolventModel>('gas');
  const [selectedSolventName, setSelectedSolventName] = useState<string>('');
  const [customSolventParams, setCustomSolventParams] = useState<Partial<SolventConfig>>({});

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
      backgroundColor: isDark ? '#1a1a1a' : '#f8f9fa',
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
  }, [systemStructure, isDark]);

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
      backgroundColor: isDark ? '#1a1a1a' : '#f8f9fa',
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
  }, [sideStructureContent, isDark]);

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

  // 自动挑选不同配位组成的溶剂化结构
  const handleAutoSelect = async () => {
    try {
      const result: AutoSelectResponse = await autoSelectSolvationStructures(jobId);

      if (result.selected_structures.length === 0) {
        message.warning('没有找到溶剂化结构');
        return;
      }

      // 只显示自动挑选的结构
      const selectedIds = new Set(result.selected_structures.map(s => s.id));
      const filteredStructures = structures.filter(s => selectedIds.has(s.id));
      setStructures(filteredStructures);

      message.success(
        `已自动挑选 ${result.unique_compositions} 种不同配位组成的结构（共 ${result.total_structures} 个结构）`
      );

      // 自动选择第一个结构
      if (filteredStructures.length > 0) {
        setSelectedStructureId(filteredStructures[0].id);
        loadSideStructure(filteredStructures[0].id);
      }
    } catch (error) {
      message.error('自动挑选失败');
    }
  };

  // 加载去溶剂化能任务列表
  const loadDesolvationJobs = async (clusterId: number) => {
    setLoadingDesolvation(true);
    try {
      const jobs = await listClusterDesolvationJobs(clusterId);
      setDesolvationJobs(jobs);
    } catch (error: any) {
      message.error(`加载去溶剂化能任务失败: ${error.message || '未知错误'}`);
    } finally {
      setLoadingDesolvation(false);
    }
  };

  // 创建去溶剂化能任务
  const handleCreateDesolvationJob = async (clusterId: number) => {
    setCreatingDesolvation(true);
    try {
      // 构建溶剂配置
      let solventConfig: SolventConfig | undefined;
      if (selectedSolventModel !== 'gas') {
        solventConfig = {
          model: selectedSolventModel,
          solvent_name: selectedSolventName || undefined,
          ...(selectedSolventModel === 'custom' ? customSolventParams : {}),
        };
      }

      await createDesolvationJob({
        md_job_id: jobId,
        solvation_structure_id: clusterId,
        method_level: selectedMethodLevel,
        desolvation_mode: selectedDesolvationMode,
        solvent_config: solventConfig,
      });
      message.success('去溶剂化能任务已创建');
      await loadDesolvationJobs(clusterId);
    } catch (error: any) {
      message.error(`创建任务失败: ${error.response?.data?.detail || error.message || '未知错误'}`);
    } finally {
      setCreatingDesolvation(false);
    }
  };

  // 查看去溶剂化能结果
  const handleViewDesolvationResult = async (jobId: number) => {
    try {
      const job = await getDesolvationJob(jobId);
      setSelectedDesolvationJob(job);
    } catch (error: any) {
      message.error(`获取任务详情失败: ${error.message || '未知错误'}`);
    }
  };

  // 打开去溶剂化能 Modal
  const handleOpenDesolvationModal = (clusterId: number) => {
    setSelectedStructureId(clusterId);
    setDesolvationModalVisible(true);
    loadDesolvationJobs(clusterId);
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
      backgroundColor: isDark ? 'rgba(31, 31, 31, 0.96)' : 'rgba(255, 255, 255, 0.96)',
      borderColor: isDark ? '#404040' : '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: isDark ? '#E8E8E8' : '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
      formatter: '{b}: {c} ({d}%)',
    },
    legend: {
      orient: 'horizontal',
      bottom: 4,
      left: 'center',
      itemWidth: 12,
      itemHeight: 12,
      itemGap: 12,
      textStyle: {
        fontSize: CHART_LEGEND_SIZE,
        fontFamily: CHART_FONT_FAMILY,
        color: isDark ? '#E8E8E8' : '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['35%', '70%'],
        center: ['50%', '45%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: isDark ? '#1F1F1F' : '#fff',
          borderWidth: 0.5,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: isDark ? '#E8E8E8' : '#333',
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
      backgroundColor: isDark ? 'rgba(31, 31, 31, 0.96)' : 'rgba(255, 255, 255, 0.96)',
      borderColor: isDark ? '#404040' : '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: isDark ? '#E8E8E8' : '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
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
        color: isDark ? '#E8E8E8' : '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['35%', '70%'],
        center: ['50%', '45%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: isDark ? '#1F1F1F' : '#fff',
          borderWidth: 0.5,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: isDark ? '#E8E8E8' : '#333',
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
      backgroundColor: isDark ? 'rgba(31, 31, 31, 0.96)' : 'rgba(255, 255, 255, 0.96)',
      borderColor: isDark ? '#404040' : '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: isDark ? '#E8E8E8' : '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
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
        color: isDark ? '#E8E8E8' : '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['35%', '70%'],
        center: ['50%', '45%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: isDark ? '#1F1F1F' : '#fff',
          borderWidth: 0.5,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: isDark ? '#E8E8E8' : '#333',
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
      backgroundColor: isDark ? 'rgba(31, 31, 31, 0.96)' : 'rgba(255, 255, 255, 0.96)',
      borderColor: isDark ? '#404040' : '#e0e0e0',
      borderWidth: 1,
      textStyle: { color: isDark ? '#E8E8E8' : '#333', fontSize: CHART_LABEL_SIZE, fontFamily: CHART_FONT_FAMILY },
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
        color: isDark ? '#E8E8E8' : '#333',
      },
    },
    series: [
      {
        type: 'pie',
        radius: ['35%', '70%'],
        center: ['50%', '45%'],
        avoidLabelOverlap: true,
        itemStyle: {
          borderRadius: 3,
          borderColor: isDark ? '#1F1F1F' : '#fff',
          borderWidth: 0.5,
        },
        label: {
          show: true,
          fontSize: CHART_LABEL_SIZE,
          fontFamily: CHART_FONT_FAMILY,
          color: isDark ? '#E8E8E8' : '#333',
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
      title: '操作',
      key: 'action',
      width: 120,
      align: 'center' as const,
      render: (_: any, record: SolvationStructure) => (
        <Space size="small">
          <Tooltip title="计算去溶剂化能">
            <Button
              type="text"
              size="small"
              icon={<ThunderboltOutlined style={{ color: '#1890ff' }} />}
              onClick={(e) => {
                e.stopPropagation();
                handleOpenDesolvationModal(record.id);
              }}
            />
          </Tooltip>
          {record.file_path && (
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
          )}
        </Space>
      ),
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
              <Button
                size="small"
                icon={<AppstoreOutlined />}
                onClick={handleAutoSelect}
                style={{ borderRadius: 6 }}
              >
                自动挑选
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
          <div className="solvation-pie-chart-container">
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
          <div className="solvation-pie-chart-container">
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
          <div className="solvation-pie-chart-container">
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
          <div className="solvation-pie-chart-container">
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
                  border: `1px solid ${token.colorBorder}`,
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
                background: isDark ? '#1F1F1F' : '#f5f5f5',
                borderRadius: 8,
                border: `1px dashed ${token.colorBorder}`,
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
            <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', width: '100%' }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                <TableOutlined style={{ color: '#fa8c16', fontSize: 16 }} />
                <span
                  style={{
                    fontSize: 14,
                    fontWeight: 600,
                    color: token.colorText,
                  }}
                >
                  溶剂化结构列表
                </span>
                <Tag color="blue" style={{ marginLeft: 4, fontSize: 11 }}>
                  {structures.length}
                </Tag>
              </div>
              {structures.length > 0 && (onGoToPostProcess || onGoToDesolvation) && (
                <Button
                  type="primary"
                  size="small"
                  icon={<ThunderboltOutlined />}
                  onClick={onGoToPostProcess || onGoToDesolvation}
                  style={{ fontSize: 12 }}
                >
                  前往后处理分析
                </Button>
              )}
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
                background: isDark ? 'rgba(255,255,255,0.04)' : '#f5f5f5',
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
                  border: `1px solid ${token.colorBorder}`,
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
                background: isDark ? '#1F1F1F' : '#f5f5f5',
                borderRadius: 8,
                border: `1px dashed ${token.colorBorder}`,
              }}
            >
              <Text type="secondary" style={{ fontSize: 12 }}>
                点击列表查看结构
              </Text>
            </div>
          )}
        </Card>
      </div>

      {/* 去溶剂化能计算 Modal（单个结构详细操作） */}
      <Modal
        title={
          <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
            <ThunderboltOutlined style={{ color: '#1890ff' }} />
            <span>去溶剂化能计算 (Desolvation Energy)</span>
          </div>
        }
        open={desolvationModalVisible}
        onCancel={() => {
          setDesolvationModalVisible(false);
          setSelectedDesolvationJob(null);
        }}
        footer={null}
        width={1000}
        destroyOnClose
      >
        <div style={{ marginTop: 16 }}>
          {/* 溶剂化结构信息 */}
          {selectedStructureId && sideStructureContent && (
            <Card
              size="small"
              style={{ marginBottom: 16, background: isDark ? 'rgba(255,255,255,0.04)' : '#f5f5f5' }}
            >
              <div style={{ display: 'flex', alignItems: 'center', gap: 16 }}>
                <div>
                  <Text strong style={{ fontSize: 14, color: token.colorText }}>
                    溶剂化结构 #{selectedStructureId}
                  </Text>
                </div>
                <Tag color="purple">
                  {sideStructureContent.center_ion}⁺ CN = {sideStructureContent.coordination_num}
                </Tag>
                <div style={{ flex: 1 }}>
                  <Space size={4} wrap>
                    {Object.entries(sideStructureContent.composition || {})
                      .filter(([_, count]) => count > 0)
                      .map(([mol, count]) => (
                        <Tag key={mol} style={{ margin: 0, fontSize: 11 }}>
                          {mol}: {count}
                        </Tag>
                      ))}
                  </Space>
                </div>
              </div>
            </Card>
          )}

          {/* 计算说明 */}
          <Card
            title={
              <Space>
                <ExperimentOutlined style={{ color: '#1890ff' }} />
                <span>计算说明</span>
              </Space>
            }
            size="small"
            style={{ marginBottom: 16 }}
          >
            <div style={{ fontSize: 13, lineHeight: 1.8, color: token.colorText }}>
              <p style={{ margin: '0 0 8px 0' }}>
                <strong>去溶剂化能定义：</strong> ΔE<sub>i</sub> = E<sub>cluster</sub> - (E<sub>cluster-i</sub> + E<sub>i</sub>)
              </p>
              <p style={{ margin: '0 0 8px 0' }}>
                <strong>计算流程：</strong>
              </p>
              <ol style={{ margin: '0 0 8px 0', paddingLeft: 20 }}>
                <li>创建完整溶剂化团簇的 QC 任务 (E<sub>cluster</sub>)</li>
                <li>为每个配体分子创建 2 个 QC 任务：
                  <ul style={{ marginTop: 4 }}>
                    <li>单独配体分子 (E<sub>i</sub>)</li>
                    <li>移除该配体后的团簇 (E<sub>cluster-i</sub>)</li>
                  </ul>
                </li>
                <li>等待所有 QC 任务完成</li>
                <li>计算每个配体的去溶剂化能</li>
              </ol>
              <p style={{ margin: 0 }}>
                <strong>预计创建 QC 任务数：</strong>
                {sideStructureContent && (
                  <Tag color="blue" style={{ marginLeft: 8 }}>
                    {1 + Object.values(sideStructureContent.composition || {}).reduce((sum, count) => sum + count, 0) * 2} 个
                  </Tag>
                )}
              </p>
            </div>
          </Card>

          {/* 方法选择和创建任务 */}
          <Card
            title={
              <Space>
                <ThunderboltOutlined style={{ color: '#52c41a' }} />
                <span>创建计算任务</span>
              </Space>
            }
            size="small"
            style={{ marginBottom: 16 }}
          >
            {/* 智能推荐提示 */}
            {(() => {
              // 检测阴离子：常见阴离子包括 PF6, TFSI, FSI, BF4, ClO4, NO3 等
              const ANION_PATTERNS = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'NO3', 'Cl', 'Br', 'I', 'OTf', 'BOB', 'DFOB'];
              const composition = sideStructureContent?.composition || {};
              const hasAnion = Object.keys(composition).some(mol =>
                ANION_PATTERNS.some(anion => mol.toUpperCase().includes(anion.toUpperCase()))
              );

              if (hasAnion) {
                return (
                  <div style={{
                    marginBottom: 12,
                    padding: '8px 12px',
                    background: isDark ? 'rgba(250, 173, 20, 0.1)' : '#fffbe6',
                    border: `1px solid ${isDark ? 'rgba(250, 173, 20, 0.3)' : '#ffe58f'}`,
                    borderRadius: 6,
                  }}>
                    <Space size={4}>
                      <BulbOutlined style={{ color: '#faad14' }} />
                      <Text style={{ fontSize: 12, color: token.colorText }}>
                        <strong>智能推荐：</strong>检测到阴离子，建议选择带弥散函数的基组（标准或精确），以提高阴离子电子云描述精度
                      </Text>
                    </Space>
                  </div>
                );
              }
              return null;
            })()}

            <div style={{ marginBottom: 12 }}>
              <Text strong style={{ fontSize: 13, color: token.colorText }}>选择计算方法：</Text>
            </div>
            <Space direction="vertical" style={{ width: '100%' }} size={8}>
              <div
                onClick={() => setSelectedMethodLevel('fast')}
                style={{
                  padding: 12,
                  border: `2px solid ${selectedMethodLevel === 'fast' ? '#1890ff' : token.colorBorder}`,
                  borderRadius: 8,
                  cursor: 'pointer',
                  background: selectedMethodLevel === 'fast' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                  transition: 'all 0.3s',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <div>
                    <Tag color="green">快速 (Fast)</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>B3LYP/6-31G(d)</Text>
                  </div>
                  <Text type="secondary" style={{ fontSize: 12 }}>~1-2 小时</Text>
                </div>
                <Text type="secondary" style={{ fontSize: 11, marginTop: 4, display: 'block' }}>
                  适用于快速预览，无弥散函数
                </Text>
              </div>
              <div
                onClick={() => setSelectedMethodLevel('standard')}
                style={{
                  padding: 12,
                  border: `2px solid ${selectedMethodLevel === 'standard' ? '#1890ff' : token.colorBorder}`,
                  borderRadius: 8,
                  cursor: 'pointer',
                  background: selectedMethodLevel === 'standard' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                  transition: 'all 0.3s',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <div>
                    <Tag color="blue">标准 (Standard)</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>B3LYP/6-31++G(d,p)</Text>
                    <Tag color="orange" style={{ marginLeft: 8, fontSize: 10 }}>含弥散</Tag>
                  </div>
                  <Text type="secondary" style={{ fontSize: 12 }}>~2-4 小时</Text>
                </div>
                <Text type="secondary" style={{ fontSize: 11, marginTop: 4, display: 'block' }}>
                  推荐用于含阴离子体系，++G 提供重/轻原子弥散函数
                </Text>
              </div>
              <div
                onClick={() => setSelectedMethodLevel('accurate')}
                style={{
                  padding: 12,
                  border: `2px solid ${selectedMethodLevel === 'accurate' ? '#1890ff' : token.colorBorder}`,
                  borderRadius: 8,
                  cursor: 'pointer',
                  background: selectedMethodLevel === 'accurate' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                  transition: 'all 0.3s',
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <div>
                    <Tag color="purple">精确 (Accurate)</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>ωB97XD/6-311++G(2d,2p)</Text>
                    <Tag color="orange" style={{ marginLeft: 8, fontSize: 10 }}>含弥散+色散</Tag>
                  </div>
                  <Text type="secondary" style={{ fontSize: 12 }}>~4-8 小时</Text>
                </div>
                <Text type="secondary" style={{ fontSize: 11, marginTop: 4, display: 'block' }}>
                  高精度计算，ωB97XD 含长程校正和色散修正，适合弱相互作用
                </Text>
              </div>
            </Space>

            {/* 去溶剂化模式选择 */}
            <Divider orientation="left" style={{ margin: '16px 0 12px' }}>
              <Text strong style={{ color: token.colorText }}>去溶剂化模式</Text>
            </Divider>
            <Space direction="horizontal" size={12}>
              <div
                onClick={() => setSelectedDesolvationMode('stepwise')}
                style={{
                  padding: 12,
                  border: `2px solid ${selectedDesolvationMode === 'stepwise' ? '#1890ff' : token.colorBorder}`,
                  borderRadius: 8,
                  cursor: 'pointer',
                  background: selectedDesolvationMode === 'stepwise' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                  transition: 'all 0.3s',
                  flex: 1,
                }}
              >
                <div>
                  <Tag color="cyan">逐级去溶剂</Tag>
                  <Text style={{ fontSize: 13, color: token.colorText }}>Stepwise</Text>
                </div>
                <Text type="secondary" style={{ fontSize: 12, display: 'block', marginTop: 4 }}>
                  逐个移除配体，计算每个配体的去溶剂化能
                </Text>
              </div>
              <div
                onClick={() => setSelectedDesolvationMode('full')}
                style={{
                  padding: 12,
                  border: `2px solid ${selectedDesolvationMode === 'full' ? '#1890ff' : token.colorBorder}`,
                  borderRadius: 8,
                  cursor: 'pointer',
                  background: selectedDesolvationMode === 'full' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                  transition: 'all 0.3s',
                  flex: 1,
                }}
              >
                <div>
                  <Tag color="orange">全部去溶剂</Tag>
                  <Text style={{ fontSize: 13, color: token.colorText }}>Full</Text>
                </div>
                <Text type="secondary" style={{ fontSize: 12, display: 'block', marginTop: 4 }}>
                  一次性移除所有配体，计算总去溶剂化能
                </Text>
              </div>
            </Space>

            {/* 溶剂模型选择 */}
            <Divider orientation="left" style={{ margin: '16px 0 12px' }}>
              <Text strong style={{ color: token.colorText }}>溶剂模型</Text>
            </Divider>
            <Row gutter={[12, 12]}>
              <Col span={24}>
                <Space direction="horizontal" size={12} wrap>
                  <div
                    onClick={() => setSelectedSolventModel('gas')}
                    style={{
                      padding: 12,
                      border: `2px solid ${selectedSolventModel === 'gas' ? '#1890ff' : token.colorBorder}`,
                      borderRadius: 8,
                      cursor: 'pointer',
                      background: selectedSolventModel === 'gas' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                      transition: 'all 0.3s',
                      minWidth: 120,
                    }}
                  >
                    <Tag color="default">气相</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>Gas</Text>
                    <Text type="secondary" style={{ fontSize: 11, display: 'block', marginTop: 4 }}>
                      无溶剂效应
                    </Text>
                  </div>
                  <div
                    onClick={() => setSelectedSolventModel('pcm')}
                    style={{
                      padding: 12,
                      border: `2px solid ${selectedSolventModel === 'pcm' ? '#1890ff' : token.colorBorder}`,
                      borderRadius: 8,
                      cursor: 'pointer',
                      background: selectedSolventModel === 'pcm' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                      transition: 'all 0.3s',
                      minWidth: 120,
                    }}
                  >
                    <Tag color="blue">PCM</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>极化连续介质</Text>
                    <Text type="secondary" style={{ fontSize: 11, display: 'block', marginTop: 4 }}>
                      IEFPCM 模型
                    </Text>
                  </div>
                  <div
                    onClick={() => setSelectedSolventModel('smd')}
                    style={{
                      padding: 12,
                      border: `2px solid ${selectedSolventModel === 'smd' ? '#1890ff' : token.colorBorder}`,
                      borderRadius: 8,
                      cursor: 'pointer',
                      background: selectedSolventModel === 'smd' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                      transition: 'all 0.3s',
                      minWidth: 120,
                    }}
                  >
                    <Tag color="green">SMD</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>溶剂密度模型</Text>
                    <Text type="secondary" style={{ fontSize: 11, display: 'block', marginTop: 4 }}>
                      更精确的溶剂化
                    </Text>
                  </div>
                  <div
                    onClick={() => setSelectedSolventModel('custom')}
                    style={{
                      padding: 12,
                      border: `2px solid ${selectedSolventModel === 'custom' ? '#1890ff' : token.colorBorder}`,
                      borderRadius: 8,
                      cursor: 'pointer',
                      background: selectedSolventModel === 'custom' ? (isDark ? 'rgba(24,144,255,0.1)' : '#e6f7ff') : 'transparent',
                      transition: 'all 0.3s',
                      minWidth: 120,
                    }}
                  >
                    <Tag color="purple">自定义</Tag>
                    <Text style={{ fontSize: 13, color: token.colorText }}>Custom</Text>
                    <Text type="secondary" style={{ fontSize: 11, display: 'block', marginTop: 4 }}>
                      自定义溶剂参数
                    </Text>
                  </div>
                </Space>
              </Col>
            </Row>

            {/* 溶剂名称选择（PCM/SMD 模式） */}
            {(selectedSolventModel === 'pcm' || selectedSolventModel === 'smd') && (
              <Row gutter={[12, 12]} style={{ marginTop: 12 }}>
                <Col span={12}>
                  <Text style={{ fontSize: 12, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                    溶剂名称
                  </Text>
                  <Select
                    value={selectedSolventName || undefined}
                    onChange={(value) => setSelectedSolventName(value)}
                    placeholder="选择溶剂"
                    style={{ width: '100%' }}
                    options={[
                      { label: 'Water (水)', value: 'water' },
                      { label: 'Acetonitrile (乙腈)', value: 'acetonitrile' },
                      { label: 'Methanol (甲醇)', value: 'methanol' },
                      { label: 'Ethanol (乙醇)', value: 'ethanol' },
                      { label: 'DMSO (二甲基亚砜)', value: 'dmso' },
                      { label: 'DMF (二甲基甲酰胺)', value: 'dmf' },
                      { label: 'THF (四氢呋喃)', value: 'thf' },
                      { label: 'Dichloromethane (二氯甲烷)', value: 'dichloromethane' },
                      { label: 'Chloroform (氯仿)', value: 'chloroform' },
                      { label: 'Acetone (丙酮)', value: 'acetone' },
                      { label: 'Toluene (甲苯)', value: 'toluene' },
                      { label: 'Hexane (正己烷)', value: 'hexane' },
                      { label: 'Diethyl Ether (乙醚)', value: 'diethylether' },
                      { label: 'Propylene Carbonate (碳酸丙烯酯)', value: 'propylenecarbonate' },
                      { label: 'Ethylene Carbonate (碳酸乙烯酯)', value: 'ethylenecarbonate' },
                    ]}
                  />
                </Col>
              </Row>
            )}

            {/* 自定义溶剂参数（Custom 模式） */}
            {selectedSolventModel === 'custom' && (
              <Collapse
                size="small"
                style={{ marginTop: 12 }}
                items={[{
                  key: 'custom-params',
                  label: <Text style={{ fontSize: 12 }}>自定义溶剂参数</Text>,
                  children: (
                    <Row gutter={[12, 12]}>
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          介电常数 ε
                        </Text>
                        <InputNumber
                          value={customSolventParams.eps}
                          onChange={(value) => setCustomSolventParams(prev => ({ ...prev, eps: value ?? undefined }))}
                          placeholder="如: 78.4"
                          style={{ width: '100%' }}
                          min={1}
                          step={0.1}
                        />
                      </Col>
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          光学介电常数 n²
                        </Text>
                        <InputNumber
                          value={customSolventParams.eps_inf}
                          onChange={(value) => setCustomSolventParams(prev => ({ ...prev, eps_inf: value ?? undefined }))}
                          placeholder="如: 1.78"
                          style={{ width: '100%' }}
                          min={1}
                          step={0.01}
                        />
                      </Col>
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          表面张力 γ
                        </Text>
                        <InputNumber
                          value={customSolventParams.surface_tension}
                          onChange={(value) => setCustomSolventParams(prev => ({ ...prev, surface_tension: value ?? undefined }))}
                          placeholder="如: 72.0"
                          style={{ width: '100%' }}
                          min={0}
                          step={0.1}
                        />
                      </Col>
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          氢键酸度 α
                        </Text>
                        <InputNumber
                          value={customSolventParams.hbond_acidity}
                          onChange={(value) => setCustomSolventParams(prev => ({ ...prev, hbond_acidity: value ?? undefined }))}
                          placeholder="如: 0.82"
                          style={{ width: '100%' }}
                          min={0}
                          max={1}
                          step={0.01}
                        />
                      </Col>
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          氢键碱度 β
                        </Text>
                        <InputNumber
                          value={customSolventParams.hbond_basicity}
                          onChange={(value) => setCustomSolventParams(prev => ({ ...prev, hbond_basicity: value ?? undefined }))}
                          placeholder="如: 0.35"
                          style={{ width: '100%' }}
                          min={0}
                          max={1}
                          step={0.01}
                        />
                      </Col>
                      <Col span={8}>
                        <Text style={{ fontSize: 11, color: token.colorTextSecondary, display: 'block', marginBottom: 4 }}>
                          芳香碳比例 φ
                        </Text>
                        <InputNumber
                          value={customSolventParams.carbon_aromaticity}
                          onChange={(value) => setCustomSolventParams(prev => ({ ...prev, carbon_aromaticity: value ?? undefined }))}
                          placeholder="如: 0.0"
                          style={{ width: '100%' }}
                          min={0}
                          max={1}
                          step={0.01}
                        />
                      </Col>
                    </Row>
                  ),
                }]}
              />
            )}

            <div style={{ marginTop: 16 }}>
              <Space>
                <Button
                  type="primary"
                  icon={<ThunderboltOutlined />}
                  onClick={() => selectedStructureId && handleCreateDesolvationJob(selectedStructureId)}
                  loading={creatingDesolvation}
                  disabled={!selectedStructureId}
                  size="large"
                >
                  创建计算任务
                </Button>
                <Button
                  icon={<ReloadOutlined />}
                  onClick={() => selectedStructureId && loadDesolvationJobs(selectedStructureId)}
                  loading={loadingDesolvation}
                  disabled={!selectedStructureId}
                >
                  刷新任务列表
                </Button>
              </Space>
            </div>
          </Card>

          {/* 任务列表 */}
          <Card title="任务列表" size="small" style={{ marginBottom: 16 }}>
            <Table
              dataSource={desolvationJobs}
              rowKey="job_id"
              loading={loadingDesolvation}
              columns={[
                {
                  title: 'ID',
                  dataIndex: 'job_id',
                  key: 'job_id',
                  width: 60,
                },
                {
                  title: '方法',
                  dataIndex: 'method_level',
                  key: 'method_level',
                  width: 100,
                },
                {
                  title: '模式',
                  dataIndex: 'desolvation_mode',
                  key: 'desolvation_mode',
                  width: 100,
                  render: (mode: string) => (
                    <Tag color={mode === 'stepwise' ? 'cyan' : 'orange'}>
                      {mode === 'stepwise' ? '逐级' : '全部'}
                    </Tag>
                  ),
                },
                {
                  title: '状态',
                  dataIndex: 'status',
                  key: 'status',
                  width: 100,
                  render: (status: string) => {
                    const statusConfig: Record<string, { color: string; text: string }> = {
                      CREATED: { color: 'default', text: '已创建' },
                      SUBMITTED: { color: 'blue', text: '已提交' },
                      QUEUED: { color: 'cyan', text: '排队中' },
                      RUNNING: { color: 'processing', text: '运行中' },
                      POSTPROCESSING: { color: 'purple', text: '后处理' },
                      COMPLETED: { color: 'success', text: '已完成' },
                      FAILED: { color: 'error', text: '失败' },
                      CANCELLED: { color: 'default', text: '已取消' },
                    };
                    const config = statusConfig[status] || { color: 'default', text: status };
                    return <Tag color={config.color}>{config.text}</Tag>;
                  },
                },
                {
                  title: '创建时间',
                  dataIndex: 'created_at',
                  key: 'created_at',
                  width: 160,
                  render: (time: string) => new Date(time).toLocaleString('zh-CN'),
                },
                {
                  title: '操作',
                  key: 'action',
                  width: 100,
                  render: (_: any, record: DesolvationJobResponse) => (
                    <Button
                      type="link"
                      size="small"
                      disabled={record.status !== 'COMPLETED'}
                      onClick={() => handleViewDesolvationResult(record.job_id)}
                    >
                      查看结果
                    </Button>
                  ),
                },
              ]}
              pagination={false}
              size="small"
            />
          </Card>

          {/* 结果展示 */}
          {selectedDesolvationJob && selectedDesolvationJob.result && (
            <DesolvationResultView result={selectedDesolvationJob.result} />
          )}
        </div>
      </Modal>

    </div>
  );
}

