/**
 * 后处理详情页面
 * - 创建模式：选择 MD Job → 选择结构 → 选择计算类型 → 提交
 * - 查看模式：显示任务状态、QC 进度、计算结果
 */
import { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import { useParams, useNavigate, useSearchParams } from 'react-router-dom';
import {
  Card,
  Steps,
  Button,
  Space,
  Select,
  Table,
  Checkbox,
  Tag,
  Typography,
  Row,
  Col,
  Statistic,
  Alert,
  Spin,
  Empty,
  message,
  Modal,
  Descriptions,
  Progress,
  Tooltip,
  Divider,
  Collapse,
  theme,
  Switch,
  InputNumber,
} from 'antd';
import {
  ArrowLeftOutlined,
  RocketOutlined,
  ExperimentOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  CloseCircleOutlined,
  SendOutlined,
  InfoCircleOutlined,
  BulbOutlined,
  AppstoreOutlined,
  UnorderedListOutlined,
  CalculatorOutlined,
  SettingOutlined,
  EyeOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import {
  listClusterAnalysisJobs,
  getClusterAnalysisJob,
  planClusterAnalysis,
  submitClusterAnalysis,
  getClusterAnalysisResults,
  CALC_TYPE_INFO,
  type AdvancedClusterJob,
  type ClusterCalcType,
  type ClusterAnalysisPlanResponse,
  type CalcTypeRequirements,
  type PlannedQCTask,
} from '../api/clusterAnalysis';
import { getMDJobs, getMDJob, getSolvationStructures, autoSelectSolvationStructures, type SolvationStructure, type AutoSelectResponse } from '../api/jobs';
import { previewDesolvationStructures, type DesolvationPreviewResponse } from '../api/desolvation';
import type { MDJob } from '../types';
import { JobStatus } from '../types';
import ClusterAnalysisResultsPanel from '../components/ClusterAnalysisResultsPanel';
import { useThemeStore } from '../stores/themeStore';
import dayjs from 'dayjs';

// 3Dmol.js 类型声明
declare global {
  interface Window {
    $3Dmol: any;
  }
}

const { Title, Text, Paragraph } = Typography;

// 状态配置
const STATUS_CONFIG: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
  CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
  SUBMITTED: { color: 'processing', icon: <SyncOutlined spin />, text: '已提交' },
  RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
  WAITING_QC: { color: 'warning', icon: <ThunderboltOutlined />, text: '等待 QC' },
  CALCULATING: { color: 'processing', icon: <SyncOutlined spin />, text: '计算中' },
  COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
  FAILED: { color: 'error', icon: <CloseCircleOutlined />, text: '失败' },
  CANCELLED: { color: 'default', icon: <CloseCircleOutlined />, text: '已取消' },
};

// 计算类型详细信息（补充物理意义、复用逻辑、路径图）
const CALC_TYPE_EXTRA: Record<string, {
  meaning: string;
  reuse: string;
  diagram: string;
  diagramTitle: string;
}> = {
  'BINDING_TOTAL': {
    meaning: '评估整个溶剂化壳层的稳定性，值越负表示离子与溶剂结合越强',
    reuse: '分子能量(Li⁺, EC, DMC等)跨结构/跨类型共享',
    diagramTitle: '能量分解',
    diagram: `┌─────────────────────────────────────┐
│  E(Li⁺·EC₂·DMC₂)  ← 簇能量(1次计算)  │
└──────────────────┬──────────────────┘
                   │
     ┌─────────────┼─────────────┐
     ↓             ↓             ↓
 ┌───────┐    ┌────────┐    ┌────────┐
 │ E(Li⁺)│    │ E(EC)  │    │ E(DMC) │
 │ (共享)│    │ ×2(共享)│    │ ×2(共享)│
 └───────┘    └────────┘    └────────┘
     └─────────────┴─────────────┘
                   ↓
         ΔE = E_簇 - E_Li - Σ(n×E_配体)`,
  },
  'BINDING_PAIRWISE': {
    meaning: '比较不同配体与离子的亲和力强弱，指导电解液配方优化',
    reuse: '二聚体能量按离子-配体对复用，单分子能量全局共享',
    diagramTitle: '配对能量计算',
    diagram: `┌──────────────────────────────────────────┐
│            二聚体能量计算                 │
├────────────────────┬─────────────────────┤
│    Li-EC 二聚体    │    Li-DMC 二聚体    │
│   E(Li⁺·EC)       │   E(Li⁺·DMC)        │
└─────────┬──────────┴──────────┬──────────┘
          ↓                     ↓
    ┌─────┴─────┐         ┌─────┴─────┐
    ↓     ↓     ↓         ↓     ↓     ↓
 E(Li⁺) E(EC) E(Li-EC)  E(Li⁺) E(DMC) E(Li-DMC)
    │     │     │         │     │       │
    └──┬──┘     │         └──┬──┘       │
       ↓        ↓            ↓          ↓
  ΔE(Li-EC)=E(Li-EC)-E(Li)-E(EC)
  ΔE(Li-DMC)=E(Li-DMC)-E(Li)-E(DMC)`,
  },
  'DESOLVATION_STEPWISE': {
    meaning: '分析离子迁移时脱溶剂化能垒，影响离子电导率',
    reuse: '所有中间态组成的能量复用，单分子能量共享',
    diagramTitle: '脱溶剂化路径树',
    diagram: `            Li⁺·EC₂·DMC₂ (完整簇)
                ↙            ↘
       Li⁺·EC₁·DMC₂      Li⁺·EC₂·DMC₁
         ↙      ↘          ↙      ↘
   Li⁺·DMC₂  Li⁺·EC₁·DMC₁  Li⁺·EC₂
         ↘       ↓    ↓       ↙
          Li⁺·DMC₁   Li⁺·EC₁
                ↘    ↙
                 Li⁺ (裸离子)

中间态数 = (n₁+1)×(n₂+1)×... - 1
例: (2+1)×(2+1)-1 = 8 种中间态`,
  },
  'DESOLVATION_FULL': {
    meaning: '与溶剂化能本质相同，计算完全脱溶剂化的总能量',
    reuse: '与 BINDING_TOTAL 共享计算',
    diagramTitle: '完全脱溶剂化',
    diagram: `Li⁺·EC₂·DMC₂ ──────→ Li⁺ + 2×EC + 2×DMC

      │                         │
      ↓                         ↓
   E_cluster              E_ion + Σ E_ligand

ΔE_desolvation = E_cluster - E_ion - Σ E_ligand`,
  },
  'REDOX': {
    meaning: '预测电解液的电化学稳定窗口（氧化/还原电位）',
    reuse: '每个唯一组成需独立计算氧化态和还原态',
    diagramTitle: '热力学循环',
    diagram: `                氧化过程
    M(gas) ─────────────────→ M⁺(gas) + e⁻
       │                          │
  ΔG_solv(M)                 ΔG_solv(M⁺)
       ↓                          ↓
    M(sol) ─────────────────→ M⁺(sol) + e⁻
                ΔG_ox(sol)

E°_ox = -ΔG_ox(sol) / nF

需要计算:
├─ 中性态气相优化 E(M,gas)
├─ 氧化态气相优化 E(M⁺,gas)
├─ 中性态溶剂化 E(M,sol)
└─ 氧化态溶剂化 E(M⁺,sol)`,
  },
  'REORGANIZATION': {
    meaning: 'Marcus理论电子转移速率常数，评估电极/电解液界面反应动力学',
    reuse: '每个唯一组成需4个计算(2优化+2单点)',
    diagramTitle: 'Marcus 4点方案',
    diagram: `        势能面示意图
    E↑
     │    ╱╲         ╱╲
     │   ╱  ╲       ╱  ╲
     │  ╱    ╲     ╱    ╲
     │ ╱  R₁  ╲   ╱  R₂  ╲
     │╱        ╲ ╱        ╲
     └──────────────────────→ Q

四点计算方案:
┌────────────────────────────────────┐
│ 1. 优化态1几何 R₁ → E(R₁,Q₁)      │
│ 2. 优化态2几何 R₂ → E(R₂,Q₂)      │
│ 3. 态1几何+态2波函 → E(R₁,Q₂)     │
│ 4. 态2几何+态1波函 → E(R₂,Q₁)     │
└────────────────────────────────────┘

λ = ½[E(R₁,Q₂) + E(R₂,Q₁)] - ½[E(R₁,Q₁) + E(R₂,Q₂)]`,
  },
};

// 计算类型选项
const CALC_TYPE_OPTIONS: { value: ClusterCalcType; label: string; description: string; riskLevel: string }[] = [
  { value: 'BINDING_TOTAL', label: '溶剂化能', description: '整个溶剂化簇的形成/脱溶剂化能', riskLevel: 'low' },
  { value: 'BINDING_PAIRWISE', label: '分子配位能', description: '单分子与离子的结合能对比', riskLevel: 'low' },
  { value: 'DESOLVATION_STEPWISE', label: '逐级脱溶剂化', description: '逐个移除配体的能量路径', riskLevel: 'medium' },
  { value: 'REDOX', label: '氧化还原电位', description: '热力学循环法计算电化学稳定性', riskLevel: 'high' },
  { value: 'REORGANIZATION', label: 'Marcus 重组能', description: 'Marcus 理论计算电子转移', riskLevel: 'high' },
];

export default function PostProcessDetail() {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const { token } = theme.useToken();

  const isCreateMode = !id || id === 'create';

  // 通用状态
  const [loading, setLoading] = useState(false);
  const [job, setJob] = useState<AdvancedClusterJob | null>(null);

  // 创建模式状态
  const [currentStep, setCurrentStep] = useState(0);
  const [mdJobs, setMdJobs] = useState<MDJob[]>([]);
  const [selectedMdJobId, setSelectedMdJobId] = useState<number | undefined>(
    searchParams.get('md_job_id') ? Number(searchParams.get('md_job_id')) : undefined
  );
  const [selectedMdJob, setSelectedMdJob] = useState<MDJob | null>(null);
  const [structures, setStructures] = useState<SolvationStructure[]>([]);
  const [selectedStructureIds, setSelectedStructureIds] = useState<number[]>([]);
  const [selectedCalcTypes, setSelectedCalcTypes] = useState<ClusterCalcType[]>([]);
  const [planResult, setPlanResult] = useState<ClusterAnalysisPlanResponse | null>(null);
  const [planLoading, setPlanLoading] = useState(false);
  const [submitLoading, setSubmitLoading] = useState(false);

  // REDOX 和 REORGANIZATION 子选项
  const [redoxOptions, setRedoxOptions] = useState({ include_molecule: true, include_dimer: true, include_cluster: false });
  const [reorganizationOptions, setReorganizationOptions] = useState({ include_molecule: true, include_cluster: true });

  // QC 配置
  const [qcConfig, setQcConfig] = useState({
    functional: 'B3LYP',
    basis_set: '6-31G*',        // 使用 6-31G* 格式以匹配已有任务
    use_dispersion: true,
    charge_ion: 1,
    solvent_model: 'pcm',       // 默认 PCM 溶剂模型（电解液计算需要溶液环境）
    solvent: 'Water',
    custom_eps: 78.4,           // 自定义介电常数
    custom_eps_inf: 1.78,       // 自定义光学介电常数
  });

  // 智能推荐状态
  const [useSmartRecommend, setUseSmartRecommend] = useState(true);
  const [recommendReason, setRecommendReason] = useState('');

  // 筛选状态
  const [filterCoordNums, setFilterCoordNums] = useState<number[]>([]);
  const [filterAnions, setFilterAnions] = useState<string[]>([]);
  const [filterSolvents, setFilterSolvents] = useState<string[]>([]);

  // 视图模式: list = 列表视图, grouped = 分组视图
  const [viewMode, setViewMode] = useState<'list' | 'grouped'>('list');
  // 智能选择加载状态
  const [autoSelectLoading, setAutoSelectLoading] = useState(false);
  // 展开的分组（结构列表）
  const [expandedGroups, setExpandedGroups] = useState<string[]>([]);
  // Step 2 Collapse 展开状态
  const [planCollapseKeys, setPlanCollapseKeys] = useState<string[]>([]);
  // 计算类型详情弹窗
  const [calcTypeDetailVisible, setCalcTypeDetailVisible] = useState(false);
  const [selectedCalcTypeForDetail, setSelectedCalcTypeForDetail] = useState<ClusterCalcType | null>(null);

  // 结构预览状态
  const [previewVisible, setPreviewVisible] = useState(false);
  const [previewLoading, setPreviewLoading] = useState(false);
  const [previewData, setPreviewData] = useState<DesolvationPreviewResponse | null>(null);
  const [selectedPreviewTab, setSelectedPreviewTab] = useState<string>('cluster');
  const [previewCalcType, setPreviewCalcType] = useState<string>('');  // 当前预览的计算类型
  const previewViewerRef = useRef<HTMLDivElement>(null);
  const previewViewerInstance = useRef<any>(null);
  const { isDark } = useThemeStore();

  // 数字转下标
  const toSubscript = (num: number): string => {
    const subscripts = '₀₁₂₃₄₅₆₇₈₉';
    return String(num).split('').map(d => subscripts[parseInt(d)]).join('');
  };

  // 从 composition 生成簇名称
  // 旧格式: Na1MPN4PF6 → 新格式: Na⁺·MPN·(PF6)₄
  const generateClusterName = (centerIon: string, composition: Record<string, number>): string => {
    if (!composition || Object.keys(composition).length === 0) return centerIon ? `${centerIon}⁺` : '-';

    // 阴离子列表（需要括号包裹）
    const anionList = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'DFOB', 'BOB', 'NO3'];

    const parts = Object.entries(composition)
      .filter(([_, count]) => count > 0)
      .sort(([a], [b]) => a.localeCompare(b))
      .map(([mol, count]) => {
        const isAnion = anionList.some(a => mol.toUpperCase().includes(a.toUpperCase()));
        if (count === 1) {
          return mol;
        } else if (isAnion) {
          // 阴离子用括号: (PF6)₄
          return `(${mol})${toSubscript(count)}`;
        } else {
          // 溶剂直接加下标: MPN₂
          return `${mol}${toSubscript(count)}`;
        }
      });

    return `${centerIon || ''}⁺·${parts.join('·')}`;
  };

  // 提取所有可用的筛选选项
  const filterOptions = useMemo(() => {
    const coordNums = new Set<number>();
    const anions = new Set<string>();
    const solvents = new Set<string>();

    // 常见阴离子列表
    const anionList = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'DFOB', 'BOB', 'Cl', 'Br', 'I', 'NO3'];

    structures.forEach(s => {
      if (s.coordination_num) coordNums.add(s.coordination_num);
      if (s.composition) {
        Object.keys(s.composition).forEach(mol => {
          if (anionList.some(a => mol.toUpperCase().includes(a.toUpperCase()))) {
            anions.add(mol);
          } else {
            solvents.add(mol);
          }
        });
      }
    });

    return {
      coordNums: Array.from(coordNums).sort((a, b) => a - b),
      anions: Array.from(anions).sort(),
      solvents: Array.from(solvents).sort(),
    };
  }, [structures]);

  // 筛选后的结构
  const filteredStructures = useMemo(() => {
    const filtered = structures.filter(s => {
      // 配位数筛选
      if (filterCoordNums.length > 0 && !filterCoordNums.includes(s.coordination_num)) {
        return false;
      }
      // 阴离子筛选
      if (filterAnions.length > 0) {
        const hasAnion = filterAnions.some(anion =>
          s.composition && s.composition[anion] && s.composition[anion] > 0
        );
        if (!hasAnion) return false;
      }
      // 溶剂筛选
      if (filterSolvents.length > 0) {
        const hasSolvent = filterSolvents.some(solvent =>
          s.composition && s.composition[solvent] && s.composition[solvent] > 0
        );
        if (!hasSolvent) return false;
      }
      return true;
    });
    // 选中的结构排到前面
    return filtered.sort((a, b) => {
      const aSelected = selectedStructureIds.includes(a.id) ? 0 : 1;
      const bSelected = selectedStructureIds.includes(b.id) ? 0 : 1;
      return aSelected - bSelected;
    });
  }, [structures, filterCoordNums, filterAnions, filterSolvents, selectedStructureIds]);

  // 重置筛选
  const resetFilters = () => {
    setFilterCoordNums([]);
    setFilterAnions([]);
    setFilterSolvents([]);
  };

  // 按组成分组的结构
  const groupedStructures = useMemo(() => {
    const groups: Record<string, { structures: SolvationStructure[]; count: number; percentage: number }> = {};
    const total = structures.length;

    structures.forEach(s => {
      const key = generateClusterName(s.center_ion, s.composition);
      if (!groups[key]) {
        groups[key] = { structures: [], count: 0, percentage: 0 };
      }
      groups[key].structures.push(s);
      groups[key].count++;
    });

    // 计算百分比
    Object.keys(groups).forEach(key => {
      groups[key].percentage = total > 0 ? (groups[key].count / total) * 100 : 0;
    });

    return groups;
  }, [structures]);

  // 分组列表（按数量排序）
  const sortedGroupKeys = useMemo(() => {
    return Object.keys(groupedStructures).sort(
      (a, b) => groupedStructures[b].count - groupedStructures[a].count
    );
  }, [groupedStructures]);

  // 获取结构的组成键
  const getStructureGroupKey = useCallback((s: SolvationStructure) => {
    return generateClusterName(s.center_ion, s.composition);
  }, []);

  // 智能选择：占比前3
  const handleSelectTop3 = useCallback(async () => {
    if (!selectedMdJobId) return;
    setAutoSelectLoading(true);
    try {
      const result = await autoSelectSolvationStructures(selectedMdJobId, 'top3');
      const ids = result.selected_structures.map(s => s.id);
      setSelectedStructureIds(ids);
      const percentages = result.selected_structures.map(
        (s: { percentage?: number; composition_key?: string }) =>
          `${s.composition_key?.split('-').slice(1).join('-') || '?'} (${s.percentage || 0}%)`
      ).join(', ');
      message.success(`占比前3: ${percentages}`);
    } catch (err) {
      console.error('Auto select failed:', err);
      message.error('选择失败');
    } finally {
      setAutoSelectLoading(false);
    }
  }, [selectedMdJobId]);

  // 智能选择：每种1个（全覆盖）
  const handleSelectAll = useCallback(async () => {
    if (!selectedMdJobId) return;
    setAutoSelectLoading(true);
    try {
      const result = await autoSelectSolvationStructures(selectedMdJobId, 'all');
      const ids = result.selected_structures.map(s => s.id);
      setSelectedStructureIds(ids);
      message.success(`每种1个: 已选择 ${ids.length} 种组成`);
    } catch (err) {
      console.error('Auto select failed:', err);
      message.error('选择失败');
    } finally {
      setAutoSelectLoading(false);
    }
  }, [selectedMdJobId]);

  // 每种组成选N个
  const selectNPerGroup = useCallback((n: number) => {
    const ids: number[] = [];
    Object.values(groupedStructures).forEach(group => {
      const selected = group.structures.slice(0, n).map(s => s.id);
      ids.push(...selected);
    });
    setSelectedStructureIds(ids);
    message.success(`已选择 ${ids.length} 个结构（每种组成最多 ${n} 个）`);
  }, [groupedStructures]);

  // 选择某个分组的全部结构
  const selectGroup = useCallback((groupKey: string, select: boolean) => {
    const groupStructures = groupedStructures[groupKey]?.structures || [];
    const groupIds = groupStructures.map(s => s.id);
    if (select) {
      setSelectedStructureIds(prev => [...new Set([...prev, ...groupIds])]);
    } else {
      setSelectedStructureIds(prev => prev.filter(id => !groupIds.includes(id)));
    }
  }, [groupedStructures]);

  // 检查分组是否全选
  const isGroupAllSelected = useCallback((groupKey: string) => {
    const groupStructures = groupedStructures[groupKey]?.structures || [];
    return groupStructures.length > 0 && groupStructures.every(s => selectedStructureIds.includes(s.id));
  }, [groupedStructures, selectedStructureIds]);

  // 检查分组是否部分选中
  const isGroupPartiallySelected = useCallback((groupKey: string) => {
    const groupStructures = groupedStructures[groupKey]?.structures || [];
    const selectedCount = groupStructures.filter(s => selectedStructureIds.includes(s.id)).length;
    return selectedCount > 0 && selectedCount < groupStructures.length;
  }, [groupedStructures, selectedStructureIds]);

  // 计算预估 QC 任务数（跨类型复用分子能量）
  const estimatedQCTasks = useMemo(() => {
    if (selectedStructureIds.length === 0 || selectedCalcTypes.length === 0) {
      return { total: 0, details: {} as Record<string, number>, breakdown: '', baseMonomerCount: 0 };
    }

    // 1. 统计唯一组成（簇）及其配体组成
    const uniqueCompositions = new Map<string, Record<string, number>>(); // key -> composition
    // 2. 统计唯一分子类型（所有结构共享）
    const uniqueMolecules = new Set<string>();
    // 3. 统计唯一中心离子
    const uniqueCenterIons = new Set<string>();
    // 4. 统计唯一离子-配体对（用于 pairwise）
    const uniquePairs = new Set<string>();

    selectedStructureIds.forEach(id => {
      const structure = structures.find(s => s.id === id);
      if (structure) {
        const compositionKey = generateClusterName(structure.center_ion, structure.composition);
        if (!uniqueCompositions.has(compositionKey)) {
          uniqueCompositions.set(compositionKey, structure.composition);
        }
        uniqueCenterIons.add(structure.center_ion);
        Object.entries(structure.composition).forEach(([mol, count]) => {
          if (count > 0) {
            uniqueMolecules.add(mol);
            uniquePairs.add(`${structure.center_ion}-${mol}`);
          }
        });
      }
    });

    const numClusters = uniqueCompositions.size;
    const numMolecules = uniqueMolecules.size;
    const numCenterIons = uniqueCenterIons.size;
    const numPairs = uniquePairs.size;

    /**
     * 计算所有中间态数量（方案1：计算所有可能的中间态组成）
     *
     * 对于组成 Li + n₁×A + n₂×B + n₃×C：
     * 唯一中间态数 = (n₁+1) × (n₂+1) × (n₃+1) × ... - 1
     *
     * 例如 Li·EC₂·DMC₂：
     * = (2+1) × (2+1) - 1 = 8 种中间态（包括完整簇，不含裸离子）
     */
    let totalIntermediates = 0;
    uniqueCompositions.forEach((composition) => {
      const counts = Object.values(composition).filter(c => c > 0);
      if (counts.length > 0) {
        // 中间态数 = ∏(n_i + 1) - 1（排除裸离子）
        const intermediateCount = counts.reduce((acc, n) => acc * (n + 1), 1) - 1;
        totalIntermediates += intermediateCount;
      }
    });

    /**
     * 跨类型复用策略：
     *
     * 共享的基础分子能量（只算一次）:
     *   - 离子能量: numCenterIons（如 Li⁺）
     *   - 配体能量: numMolecules（如 EC, DMC, PF6）
     *
     * 各类型独有计算:
     *   - BINDING_TOTAL: 完整簇能量 (numClusters)
     *   - BINDING_PAIRWISE: 二聚体能量 (numPairs)
     *   - DESOLVATION_STEPWISE: 所有中间态能量 (totalIntermediates)
     *   - REDOX: 每簇 4 个计算
     *   - REORGANIZATION: 每簇 4 个计算
     */

    // 判断哪些类型需要基础分子能量
    const needsMonomerEnergy = selectedCalcTypes.some(t =>
      ['BINDING_TOTAL', 'BINDING_PAIRWISE', 'DESOLVATION_STEPWISE'].includes(t)
    );

    // 计算基础分子能量（跨类型共享，只算一次）
    const baseMonomerCount = needsMonomerEnergy ? (numCenterIons + numMolecules) : 0;

    // 计算各类型独有的任务
    const details: Record<string, number> = {};
    let typeSpecificTotal = 0;

    selectedCalcTypes.forEach(calcType => {
      let count = 0;
      switch (calcType) {
        case 'BINDING_TOTAL':
          // 只需要完整簇能量
          count = numClusters;
          break;
        case 'BINDING_PAIRWISE':
          // 只需要二聚体能量
          count = numPairs;
          break;
        case 'DESOLVATION_STEPWISE':
          // 所有中间态能量（包括完整簇，不含裸离子）
          count = totalIntermediates;
          break;
        case 'REDOX':
          // 每个簇: 还原态优化 + 氧化态优化 + 溶剂化校正
          count = numClusters * 4;
          break;
        case 'REORGANIZATION':
          // Marcus 4点方案
          count = numClusters * 4;
          break;
      }
      details[calcType] = count;
      typeSpecificTotal += count;
    });

    // 总任务 = 基础分子能量(共享) + 各类型独有任务
    const total = baseMonomerCount + typeSpecificTotal;

    return { total, details, breakdown: '', baseMonomerCount };
  }, [selectedStructureIds, selectedCalcTypes, structures]);

  // 加载已完成的 MD Jobs
  const loadMdJobs = useCallback(async () => {
    try {
      const data = await getMDJobs();
      setMdJobs(data.filter(j => j.status === JobStatus.COMPLETED));
    } catch (err) {
      console.error('Failed to load MD jobs:', err);
    }
  }, []);

  // 加载单个 MD Job 详情
  const loadMdJobDetail = useCallback(async (mdJobId: number) => {
    try {
      const data = await getMDJob(mdJobId);
      setSelectedMdJob(data);
    } catch (err) {
      console.error('Failed to load MD job detail:', err);
    }
  }, []);

  // 加载溶剂化结构
  const loadStructures = useCallback(async (mdJobId: number) => {
    setLoading(true);
    try {
      const data = await getSolvationStructures(mdJobId);
      setStructures(data);
    } catch (err) {
      console.error('Failed to load structures:', err);
      message.error('加载溶剂化结构失败');
    } finally {
      setLoading(false);
    }
  }, []);

  // 加载现有任务详情
  const loadJob = useCallback(async () => {
    if (!id || isCreateMode) return;
    setLoading(true);
    try {
      const data = await getClusterAnalysisJob(Number(id));
      setJob(data);
    } catch (err) {
      console.error('Failed to load job:', err);
      message.error('加载任务详情失败');
    } finally {
      setLoading(false);
    }
  }, [id, isCreateMode]);

  useEffect(() => {
    if (isCreateMode) {
      loadMdJobs();
    } else {
      loadJob();
    }
  }, [isCreateMode, loadMdJobs, loadJob]);

  useEffect(() => {
    if (selectedMdJobId && isCreateMode) {
      loadMdJobDetail(selectedMdJobId);
      loadStructures(selectedMdJobId);
      setCurrentStep(1);
    }
  }, [selectedMdJobId, isCreateMode, loadMdJobDetail, loadStructures]);

  // 加载 3Dmol.js
  useEffect(() => {
    if (!window.$3Dmol) {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
      script.async = true;
      document.body.appendChild(script);
    }
  }, []);

  // 智能推荐：当选择的结构或计算类型变化时自动更新推荐参数
  useEffect(() => {
    if (!useSmartRecommend || selectedStructureIds.length === 0) return;

    // 检测是否包含阴离子（需要弥散基组）
    const selectedStructs = structures.filter(s => selectedStructureIds.includes(s.id));
    const hasAnion = selectedStructs.some(s => {
      const comp = s.composition || {};
      // 常见阴离子名称
      const anionPatterns = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'NO3', 'Cl', 'Br', 'I'];
      return Object.keys(comp).some(mol =>
        anionPatterns.some(pattern => mol.toUpperCase().includes(pattern))
      );
    });

    const hasRedox = selectedCalcTypes.includes('REDOX');
    const hasReorg = selectedCalcTypes.includes('REORGANIZATION');

    // 根据分子特性选择基组
    const recommendedBasis = hasAnion ? '6-31+G(d,p)' : '6-31G*';
    const basisReason = hasAnion ? '(阴离子需弥散函数)' : '';

    if (hasRedox || hasReorg) {
      setQcConfig(prev => ({
        ...prev,
        functional: 'B3LYP',
        basis_set: recommendedBasis,
        use_dispersion: true,
        solvent_model: 'pcm',
        solvent: prev.solvent || 'Water',
      }));
      setRecommendReason(`Redox/重组能推荐 B3LYP/${recommendedBasis} + PCM溶剂 ${basisReason}`);
    } else if (selectedCalcTypes.length > 0) {
      setQcConfig(prev => ({
        ...prev,
        functional: 'B3LYP',
        basis_set: recommendedBasis,
        use_dispersion: true,
        solvent_model: 'pcm',
        solvent: prev.solvent || 'Water',
      }));
      setRecommendReason(`Binding 推荐 B3LYP-D3BJ/${recommendedBasis} + PCM溶剂 ${basisReason}`);
    }
  }, [useSmartRecommend, selectedStructureIds, selectedCalcTypes, structures]);

  // 打开结构预览
  const handlePreviewStructure = async (structureId: number, calcType: string = '') => {
    setPreviewLoading(true);
    setPreviewVisible(true);
    setPreviewCalcType(calcType);
    // 根据计算类型设置默认选中的 Tab
    if (calcType === 'BINDING_PAIRWISE') {
      setSelectedPreviewTab('dimer_0');  // Pairwise 默认显示第一个 dimer
    } else {
      setSelectedPreviewTab('cluster');  // 其他默认显示 cluster
    }
    try {
      const data = await previewDesolvationStructures(structureId);
      setPreviewData(data);
    } catch (error: any) {
      message.error(`加载预览失败: ${error.message || '未知错误'}`);
      setPreviewVisible(false);
    } finally {
      setPreviewLoading(false);
    }
  };

  // 渲染 3D 预览
  const renderMolecule = useCallback((xyzContent: string, highlightCenterIon: boolean = true) => {
    if (!previewViewerRef.current || !window.$3Dmol || !xyzContent) return;

    // 清除旧的 viewer
    if (previewViewerInstance.current) {
      previewViewerInstance.current.clear();
      previewViewerInstance.current = null;
    }

    const container = previewViewerRef.current;
    const viewer = window.$3Dmol.createViewer(container, {
      backgroundColor: isDark ? '#1a1a1a' : '#f8f9fa',
    });
    previewViewerInstance.current = viewer;

    viewer.addModel(xyzContent, 'xyz');
    viewer.setStyle({}, {
      stick: { radius: 0.15, colorscheme: 'Jmol' },
      sphere: { scale: 0.3, colorscheme: 'Jmol' },
    });
    // 只对 cluster 高亮中心离子（第一个原子）
    if (highlightCenterIon) {
      viewer.setStyle({ serial: 0 }, {
        sphere: { scale: 0.5, color: '#e74c3c' },
      });
    }
    viewer.zoomTo(0.85);
    viewer.rotate(-15, 'x');
    viewer.rotate(10, 'y');
    viewer.render();
  }, [isDark]);

  // 当预览 Tab 变化或数据变化时重新渲染
  useEffect(() => {
    if (!previewData || !previewVisible) return;

    let xyzContent = '';
    let highlightCenterIon = true;

    if (selectedPreviewTab === 'cluster') {
      xyzContent = previewData.cluster.xyz_content;
      highlightCenterIon = true;
    } else if (selectedPreviewTab === 'center_ion') {
      xyzContent = previewData.center_ion_structure.xyz_content;
      highlightCenterIon = true;
    } else if (selectedPreviewTab.startsWith('dimer_')) {
      // dimer 结构（Li + 配体）
      const idx = parseInt(selectedPreviewTab.replace('dimer_', ''), 10);
      const dimer = previewData.dimer_structures?.[idx];
      if (dimer) {
        xyzContent = dimer.xyz_content;
      }
      highlightCenterIon = true;  // 高亮中心离子
    } else if (selectedPreviewTab.startsWith('cluster_minus_')) {
      // cluster-minus 结构
      const idx = parseInt(selectedPreviewTab.replace('cluster_minus_', ''), 10);
      const cm = previewData.cluster_minus_structures?.[idx];
      if (cm) {
        xyzContent = cm.xyz_content;
      }
      highlightCenterIon = true;
    } else if (selectedPreviewTab.startsWith('ligand_')) {
      const idx = parseInt(selectedPreviewTab.replace('ligand_', ''), 10);
      const ligand = previewData.ligands[idx];
      if (ligand) {
        xyzContent = ligand.xyz_content;
      }
      highlightCenterIon = false;
    }

    // 延迟渲染，确保 Modal 已完全展开
    const timer = setTimeout(() => {
      renderMolecule(xyzContent, highlightCenterIon);
    }, 100);

    return () => clearTimeout(timer);
  }, [previewData, previewVisible, selectedPreviewTab, renderMolecule]);

  // 生成规划预览
  const handlePlan = async () => {
    if (!selectedMdJobId || selectedStructureIds.length === 0 || selectedCalcTypes.length === 0) {
      message.warning('请选择结构和计算类型');
      return;
    }
    setPlanLoading(true);
    try {
      const result = await planClusterAnalysis({
        md_job_id: selectedMdJobId,
        solvation_structure_ids: selectedStructureIds,
        calc_types: selectedCalcTypes,
        redox_options: selectedCalcTypes.includes('REDOX') ? redoxOptions : undefined,
        reorganization_options: selectedCalcTypes.includes('REORGANIZATION') ? reorganizationOptions : undefined,
        qc_config: qcConfig,
      });
      setPlanResult(result);
      // 初始化 Collapse 展开状态：默认全部展开
      setPlanCollapseKeys(result.calc_requirements.map(r => r.calc_type));
      setCurrentStep(2);  // Step 2: 确认提交页面
    } catch (err) {
      console.error('Failed to plan:', err);
      message.error('生成规划失败');
    } finally {
      setPlanLoading(false);
    }
  };

  // 提交任务
  const handleSubmit = async () => {
    if (!selectedMdJobId || selectedStructureIds.length === 0 || selectedCalcTypes.length === 0) {
      return;
    }
    setSubmitLoading(true);
    try {
      const job = await submitClusterAnalysis({
        md_job_id: selectedMdJobId,
        solvation_structure_ids: selectedStructureIds,
        calc_types: selectedCalcTypes,
        qc_config: qcConfig,
      });
      message.success('分析任务已提交');
      navigate(`/workspace/liquid-electrolyte/analysis/${job.id}`);
    } catch (err) {
      console.error('Failed to submit:', err);
      message.error('提交任务失败');
    } finally {
      setSubmitLoading(false);
    }
  };

  // 结构选择表格列
  const structureColumns: ColumnsType<SolvationStructure> = [
    {
      title: '',
      width: 40,
      render: (_, record) => (
        <Checkbox
          checked={selectedStructureIds.includes(record.id)}
          onChange={(e) => {
            if (e.target.checked) {
              setSelectedStructureIds([...selectedStructureIds, record.id]);
            } else {
              setSelectedStructureIds(selectedStructureIds.filter(id => id !== record.id));
            }
          }}
        />
      ),
    },
    {
      title: '#',
      dataIndex: 'id',
      width: 60,
      render: (id: number) => <Text type="secondary">#{id}</Text>,
    },
    {
      title: '簇名称',
      key: 'cluster_name',
      width: 150,
      render: (_, record) => (
        <Text strong style={{ fontFamily: 'monospace' }}>
          {generateClusterName(record.center_ion, record.composition)}
        </Text>
      ),
    },
    {
      title: '组成',
      dataIndex: 'composition',
      width: 200,
      render: (comp: Record<string, number>) => {
        if (!comp || Object.keys(comp).length === 0) return '-';
        return (
          <Space size={4} wrap>
            {Object.entries(comp)
              .filter(([_, count]) => count > 0)
              .map(([mol, count]) => (
                <Tag key={mol} color="cyan" style={{ margin: 0 }}>{mol}×{count}</Tag>
              ))}
          </Space>
        );
      },
    },
    {
      title: 'CN',
      dataIndex: 'coordination_num',
      width: 55,
      align: 'center' as const,
      render: (cn: number) => <Tag color="blue">{cn ?? '-'}</Tag>,
    },
    {
      title: '帧号',
      dataIndex: 'snapshot_frame',
      width: 70,
      align: 'center' as const,
      render: (frame: number) => <Text type="secondary">{frame ?? '-'}</Text>,
    },
    {
      title: '占比',
      key: 'percentage',
      width: 80,
      align: 'center' as const,
      render: (_, record) => {
        const groupKey = getStructureGroupKey(record);
        const group = groupedStructures[groupKey];
        const pct = group?.percentage || 0;
        return (
          <Tooltip title={`该组成共 ${group?.count || 0} 个`}>
            <Text type="secondary">{pct.toFixed(1)}%</Text>
          </Tooltip>
        );
      },
    },
    {
      title: '组内',
      key: 'group_count',
      width: 60,
      align: 'center' as const,
      render: (_, record) => {
        const groupKey = getStructureGroupKey(record);
        const group = groupedStructures[groupKey];
        return <Tag>{group?.count || 0}</Tag>;
      },
    },
  ];

  // 计算当前步骤的有效性
  const canProceedToStep2 = selectedMdJobId && selectedStructureIds.length > 0 && selectedCalcTypes.length > 0;
  const canSubmit = canProceedToStep2 && planResult;

  // 计算类型的详细说明和示意图
  const CALC_TYPE_DETAILS: Record<ClusterCalcType, { diagram: string; explanation: string }> = {
    BINDING_TOTAL: {
      diagram: 'Li⁺ + A + B + C → [Li·ABC]  (ΔE = 溶剂化能)',
      explanation: '计算离子与所有配体结合释放的能量，等于完全脱溶剂化能的负值',
    },
    BINDING_PAIRWISE: {
      diagram: 'Li⁺ + A → [Li-A]，Li⁺ + B → [Li-B] ...',
      explanation: '分别计算每个配体的结合能，对比不同分子（阴离子 vs 溶剂）的配位能力',
    },
    DESOLVATION_STEPWISE: {
      diagram: '[Li·ABC] → [Li·AB] + C → [Li·A] + B → Li⁺ + A',
      explanation: '模拟逐步脱溶剂化过程，找出能量最优的脱离顺序',
    },
    DESOLVATION_FULL: {
      diagram: '[Li·ABC] → Li⁺ + A + B + C',
      explanation: '与溶剂化能相同（符号相反）',
    },
    REDOX: {
      diagram: 'M → M⁺ + e⁻  (氧化电位)  |  M + e⁻ → M⁻  (还原电位)',
      explanation: '预测电解液的电化学稳定窗口',
    },
    REORGANIZATION: {
      diagram: 'λ = [E(R₁,G₂) - E(R₁,G₁)] + [E(R₂,G₁) - E(R₂,G₂)]',
      explanation: 'Marcus 4点法计算电子转移活化能',
    },
  };

  // 渲染创建模式 - 单页布局，无 Step 0
  const renderCreateMode = () => {
    const coveredCompositions = Object.keys(groupedStructures).filter(k =>
      groupedStructures[k].structures.some(s => selectedStructureIds.includes(s.id))
    ).length;

    // Step 2: 确认提交 - 详细任务规划页面
    if (currentStep === 2) {
      // QC 任务表格列定义
      const taskColumns: ColumnsType<PlannedQCTask> = [
        {
          title: 'Name',
          dataIndex: 'task_type',
          key: 'task_type',
          width: 200,
          render: (taskType: string) => {
            // 从 task_type 提取英文名称
            if (taskType === 'cluster') return <Tag color="purple">Cluster</Tag>;
            if (taskType === 'ion') return <Tag color="gold">Li⁺ Ion</Tag>;
            if (taskType?.startsWith('ligand_')) {
              const name = taskType.replace('ligand_', '');
              return <Tag color="blue">Ligand: {name}</Tag>;
            }
            if (taskType?.startsWith('dimer_')) {
              const name = taskType.replace('dimer_', '');
              return <Tag color="green">Dimer: Li-{name}</Tag>;
            }
            if (taskType?.startsWith('intermediate_')) {
              const parts = taskType.replace('intermediate_', '').split('_');
              return <Tag color="orange">Intermediate #{parts[0]}</Tag>;
            }
            if (taskType?.startsWith('redox_mol_')) {
              const parts = taskType.replace('redox_mol_', '').split('_');
              const mol = parts[0];
              const state = parts.slice(1).join('_');
              const stateMap: Record<string, string> = {
                'neutral_gas': 'N/Gas', 'charged_gas': 'Ox/Gas',
                'neutral_sol': 'N/Sol', 'charged_sol': 'Ox/Sol',
              };
              return <Tag color="magenta">Redox-Mol: {mol} ({stateMap[state] || state})</Tag>;
            }
            if (taskType?.startsWith('redox_dimer_')) {
              const parts = taskType.replace('redox_dimer_', '').split('_');
              const mol = parts[0];
              const state = parts.slice(1).join('_');
              const stateMap: Record<string, string> = {
                'neutral_gas': 'N/Gas', 'charged_gas': 'Ox/Gas',
                'neutral_sol': 'N/Sol', 'charged_sol': 'Ox/Sol',
              };
              return <Tag color="geekblue">Redox-Dimer: Li-{mol} ({stateMap[state] || state})</Tag>;
            }
            if (taskType?.startsWith('reorg_mol_')) {
              const parts = taskType.replace('reorg_mol_', '').split('_');
              const mol = parts[0];
              const mode = parts.slice(1).join('_');
              const modeMap: Record<string, string> = {
                'opt_neutral': 'Opt-N', 'opt_charged': 'Opt-Ox',
                'sp_charged_at_neutral': 'SP-Ox@N', 'sp_neutral_at_charged': 'SP-N@Ox',
              };
              return <Tag color="volcano">Reorg-Mol: {mol} ({modeMap[mode] || mode})</Tag>;
            }
            if (taskType?.startsWith('reorg_cluster_')) {
              const parts = taskType.replace('reorg_cluster_', '').split('_');
              const id = parts[0];
              const mode = parts.slice(1).join('_');
              const modeMap: Record<string, string> = {
                'opt_neutral': 'Opt-N', 'opt_charged': 'Opt-Ox',
                'sp_charged_at_neutral': 'SP-Ox@N', 'sp_neutral_at_charged': 'SP-N@Ox',
              };
              return <Tag color="purple">Reorg-Cluster #{id} ({modeMap[mode] || mode})</Tag>;
            }
            return <Tag>{taskType}</Tag>;
          },
        },
        {
          title: '描述',
          dataIndex: 'description',
          key: 'description',
          ellipsis: true,
          render: (desc: string) => (
            <Tooltip title={desc}>
              <span style={{ fontSize: 12 }}>{desc}</span>
            </Tooltip>
          ),
        },
        {
          title: 'SMILES',
          dataIndex: 'smiles',
          key: 'smiles',
          width: 140,
          render: (smiles: string) => smiles ? (
            <Tooltip title={smiles}>
              <Text code style={{ fontSize: 10 }}>{smiles.length > 18 ? smiles.slice(0, 18) + '...' : smiles}</Text>
            </Tooltip>
          ) : '-',
        },
        {
          title: '电荷',
          dataIndex: 'charge',
          key: 'charge',
          width: 50,
          align: 'center',
          render: (c: number) => c > 0 ? `+${c}` : c,
        },
        {
          title: '多重度',
          dataIndex: 'multiplicity',
          key: 'multiplicity',
          width: 55,
          align: 'center',
        },
        {
          title: '状态',
          dataIndex: 'status',
          key: 'status',
          width: 95,
          align: 'center',
          render: (status: string) => {
            if (status === 'reused') {
              return <Tag color="success" icon={<CheckCircleOutlined />}>全局复用</Tag>;
            }
            if (status === 'local_reused') {
              return <Tag color="cyan" icon={<CheckCircleOutlined />}>局部复用</Tag>;
            }
            return <Tag color="warning" icon={<ThunderboltOutlined />}>新建</Tag>;
          },
        },
      ];

      // 获取第一个有 structure_id 的任务，用于预览
      const getFirstStructureId = (tasks: PlannedQCTask[]): number | null => {
        for (const task of tasks) {
          if (task.structure_id) return task.structure_id;
        }
        return null;
      };

      // 过滤任务：根据子选项筛选 REDOX 和 REORGANIZATION 的任务（仅影响显示，不影响实际提交）
      const filterTasks = (tasks: PlannedQCTask[], calcType: string): PlannedQCTask[] => {
        if (calcType === 'REDOX') {
          return tasks.filter(t => {
            if (t.task_type?.startsWith('redox_mol_')) return redoxOptions.include_molecule;
            if (t.task_type?.startsWith('redox_dimer_')) return redoxOptions.include_dimer;
            if (t.task_type?.startsWith('redox_cluster_')) return redoxOptions.include_cluster;
            return true;
          });
        }
        if (calcType === 'REORGANIZATION') {
          return tasks.filter(t => {
            if (t.task_type?.startsWith('reorg_mol_')) return reorganizationOptions.include_molecule;
            if (t.task_type?.startsWith('reorg_cluster_')) return reorganizationOptions.include_cluster;
            return true;
          });
        }
        return tasks;
      };

      // 构建 Collapse 项
      const collapseItems = planResult?.calc_requirements.map((req: CalcTypeRequirements) => {
        const info = CALC_TYPE_INFO[req.calc_type as ClusterCalcType];
        const firstStructureId = getFirstStructureId(req.required_qc_tasks);
        const filteredTasks = filterTasks(req.required_qc_tasks, req.calc_type);
        const newCount = filteredTasks.filter(t => t.status === 'new').length;
        const reusedCount = filteredTasks.filter(t => t.status === 'reused' || t.status === 'local_reused').length;

        return {
          key: req.calc_type,
          label: (
            <div style={{
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'space-between',
              width: '100%',
              flexWrap: 'wrap',
              gap: 8
            }}>
              <span style={{ flex: '1 1 auto', minWidth: 200 }}>
                {info?.icon} <Text strong>{info?.label || req.calc_type}</Text>
                <br />
                <Text type="secondary" style={{ fontSize: 12 }}>{req.description}</Text>
              </span>
              <span style={{ display: 'flex', alignItems: 'center', gap: 8, flexShrink: 0 }}>
                {/* REDOX 子选项 - 仅用于前端过滤显示 */}
                {req.calc_type === 'REDOX' && (
                  <>
                    <Checkbox
                      checked={redoxOptions.include_molecule}
                      onChange={(e) => { e.stopPropagation(); setRedoxOptions({ ...redoxOptions, include_molecule: e.target.checked }); }}
                      onClick={(e) => e.stopPropagation()}
                    >
                      <Text style={{ fontSize: 11 }}>Molecule</Text>
                    </Checkbox>
                    <Checkbox
                      checked={redoxOptions.include_dimer}
                      onChange={(e) => { e.stopPropagation(); setRedoxOptions({ ...redoxOptions, include_dimer: e.target.checked }); }}
                      onClick={(e) => e.stopPropagation()}
                    >
                      <Text style={{ fontSize: 11 }}>Li-Dimer</Text>
                    </Checkbox>
                    <Checkbox
                      checked={redoxOptions.include_cluster}
                      onChange={(e) => { e.stopPropagation(); setRedoxOptions({ ...redoxOptions, include_cluster: e.target.checked }); }}
                      onClick={(e) => e.stopPropagation()}
                    >
                      <Text style={{ fontSize: 11 }}>Cluster</Text>
                    </Checkbox>
                    <span style={{ width: 1, height: 16, background: token.colorBorder, margin: '0 4px' }} />
                  </>
                )}
                {/* REORGANIZATION 子选项 - 仅用于前端过滤显示 */}
                {req.calc_type === 'REORGANIZATION' && (
                  <>
                    <Checkbox
                      checked={reorganizationOptions.include_molecule}
                      onChange={(e) => { e.stopPropagation(); setReorganizationOptions({ ...reorganizationOptions, include_molecule: e.target.checked }); }}
                      onClick={(e) => e.stopPropagation()}
                    >
                      <Text style={{ fontSize: 11 }}>Molecule</Text>
                    </Checkbox>
                    <Checkbox
                      checked={reorganizationOptions.include_cluster}
                      onChange={(e) => { e.stopPropagation(); setReorganizationOptions({ ...reorganizationOptions, include_cluster: e.target.checked }); }}
                      onClick={(e) => e.stopPropagation()}
                    >
                      <Text style={{ fontSize: 11 }}>Cluster</Text>
                    </Checkbox>
                    <span style={{ width: 1, height: 16, background: token.colorBorder, margin: '0 4px' }} />
                  </>
                )}
                {firstStructureId && (
                  <Tooltip title="Preview structures">
                    <Button
                      size="small"
                      icon={<EyeOutlined />}
                      onClick={(e) => {
                        e.stopPropagation();
                        handlePreviewStructure(firstStructureId, req.calc_type);
                      }}
                    >
                      Preview
                    </Button>
                  </Tooltip>
                )}
                <Tag color="warning">New {newCount}</Tag>
                <Tag color="success">Reuse {reusedCount}</Tag>
              </span>
            </div>
          ),
          children: (
            <Table
              dataSource={filteredTasks}
              columns={taskColumns}
              rowKey={(record, index) => `${record.task_type}-${record.smiles || record.structure_id}-${index}`}
              size="small"
              pagination={false}
              scroll={{ y: 300 }}
              rowClassName={(record) =>
                record.status === 'reused' ? 'row-reused' :
                record.status === 'local_reused' ? 'row-local-reused' : 'row-new'
              }
            />
          ),
        };
      }) || [];

      return (
        <div style={{ padding: '24px 48px', maxHeight: 'calc(100vh - 120px)', overflowY: 'auto' }}>
          {/* 返回按钮 - 更醒目的样式 */}
          <Button
            type="primary"
            ghost
            icon={<ArrowLeftOutlined />}
            onClick={() => setCurrentStep(1)}
            style={{ marginBottom: 16 }}
          >
            返回修改配置
          </Button>

          {planResult ? (
            <>
              {/* 统计概览 */}
              <Card title="📋 任务规划概览" style={{ marginBottom: 16 }}>
                <Row gutter={24}>
                  <Col span={6}>
                    <Statistic
                      title="选中结构"
                      value={planResult.selected_structures_count}
                      suffix="个"
                      prefix={<AppstoreOutlined />}
                    />
                  </Col>
                  <Col span={6}>
                    <Statistic
                      title="新建 QC 任务"
                      value={planResult.total_new_qc_tasks}
                      valueStyle={{ color: token.colorWarning }}
                      prefix={<ThunderboltOutlined />}
                    />
                  </Col>
                  <Col span={6}>
                    <Statistic
                      title="复用 QC 任务"
                      value={planResult.total_reused_qc_tasks}
                      valueStyle={{ color: token.colorSuccess }}
                      prefix={<CheckCircleOutlined />}
                    />
                  </Col>
                  <Col span={6}>
                    <Statistic
                      title="预估计算时间"
                      value={(planResult.estimated_compute_hours ?? 0).toFixed(1)}
                      suffix="核时"
                      prefix={<ClockCircleOutlined />}
                    />
                  </Col>
                </Row>
              </Card>

              {/* 计算参数配置 */}
              <Card
                title={<><SettingOutlined /> 计算参数配置</>}
                size="small"
                style={{ marginBottom: 16 }}
              >
                <Descriptions column={4} size="small">
                  <Descriptions.Item label="泛函">
                    <Tag color="blue">{qcConfig.functional}{qcConfig.use_dispersion ? '-D3BJ' : ''}</Tag>
                  </Descriptions.Item>
                  <Descriptions.Item label="基组">
                    <Tag color="green">{qcConfig.basis_set}</Tag>
                  </Descriptions.Item>
                  <Descriptions.Item label="溶剂模型">
                    <Tag color="orange">
                      {qcConfig.solvent_model === 'gas' ? '气相' : `${qcConfig.solvent_model.toUpperCase()}/${qcConfig.solvent}`}
                    </Tag>
                  </Descriptions.Item>
                  <Descriptions.Item label="色散校正">
                    <Tag color={qcConfig.use_dispersion ? 'success' : 'default'}>
                      {qcConfig.use_dispersion ? '已启用' : '未启用'}
                    </Tag>
                  </Descriptions.Item>
                </Descriptions>
              </Card>

              {/* 按计算类型分组的任务详情 */}
              <Card
                title={<><UnorderedListOutlined /> QC 任务详情 (按计算类型分组)</>}
                style={{ marginBottom: 16 }}
              >
                <Collapse
                  items={collapseItems}
                  activeKey={planCollapseKeys}
                  onChange={(keys) => setPlanCollapseKeys(keys as string[])}
                  style={{ background: 'transparent' }}
                  className="no-animation-collapse"
                />
              </Card>

              {/* 警告信息 */}
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

              {/* 提交按钮 */}
              <div style={{ textAlign: 'center', padding: '16px 0' }}>
                <Button
                  type="primary"
                  size="large"
                  icon={<SendOutlined />}
                  loading={submitLoading}
                  onClick={handleSubmit}
                  style={{ paddingLeft: 48, paddingRight: 48, height: 48 }}
                >
                  提交分析任务
                </Button>
              </div>
            </>
          ) : (
            <Card>
              <div style={{ textAlign: 'center', padding: 48 }}>
                <Spin size="large" />
                <div style={{ marginTop: 16 }}><Text type="secondary">正在生成规划...</Text></div>
              </div>
            </Card>
          )}

          {/* 自定义样式：区分新建和复用任务行，支持 dark mode */}
          <style>{`
            .row-reused {
              background-color: ${token.colorSuccessBg} !important;
            }
            .row-reused td {
              color: ${token.colorText} !important;
            }
            .row-local-reused {
              background-color: ${token.colorInfoBg} !important;
            }
            .row-local-reused td {
              color: ${token.colorText} !important;
            }
            .row-new {
              background-color: ${token.colorWarningBg} !important;
            }
            .row-new td {
              color: ${token.colorText} !important;
            }
            .ant-collapse-header {
              padding: 12px 16px !important;
            }
            .ant-table-cell {
              word-break: break-word !important;
            }
            /* 禁用 Collapse 内容变化时的动画 */
            .no-animation-collapse .ant-collapse-content {
              transition: none !important;
            }
            .no-animation-collapse .ant-collapse-content-box {
              transition: none !important;
            }
          `}</style>
        </div>
      );
    }

    // Step 1: 配置分析（三栏布局）
    return (
      <div style={{ display: 'flex', flexDirection: 'column', height: 'calc(100vh - 160px)' }}>
        {/* 顶部工具栏：数据源选择 + 进度 */}
        <div style={{
          background: token.colorBgContainer,
          borderRadius: 8,
          padding: '12px 16px',
          marginBottom: 12,
          display: 'flex',
          justifyContent: 'space-between',
          alignItems: 'center',
        }}>
          <Space>
            <Text strong>数据源:</Text>
            <Select
              style={{ width: 350 }}
              placeholder="选择 MD 模拟任务"
              value={selectedMdJobId}
              onChange={(v) => {
                setSelectedMdJobId(v);
                setSelectedStructureIds([]);
                setSelectedCalcTypes([]);
                setPlanResult(null);
              }}
              options={mdJobs.map(j => ({
                value: j.id,
                label: `${j.config?.job_name || 'MD Job'} (#${j.id})`,
              }))}
              showSearch
              filterOption={(input, option) =>
                (option?.label ?? '').toLowerCase().includes(input.toLowerCase())
              }
            />
            {selectedMdJobId && (
              <Tag color="green">{structures.length} 个溶剂化结构</Tag>
            )}
          </Space>
          <Space>
            <Text type="secondary">已选 {selectedStructureIds.length} 结构 · {selectedCalcTypes.length} 种计算</Text>
            <Button
              type="primary"
              disabled={!canProceedToStep2}
              loading={planLoading}
              onClick={handlePlan}
            >
              生成规划 →
            </Button>
          </Space>
        </div>

          {/* 主内容区 - 三栏布局 */}
          <Row gutter={16} style={{ flex: 1, minHeight: 0 }}>
            {/* 左栏：结构选择 */}
            <Col span={10} style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
              <Card
                size="small"
                title={<Space><ExperimentOutlined /> 溶剂化结构</Space>}
                extra={
                  <Space.Compact size="small">
                    <Button
                      icon={<UnorderedListOutlined />}
                      type={viewMode === 'list' ? 'primary' : 'default'}
                      onClick={() => setViewMode('list')}
                    />
                    <Button
                      icon={<AppstoreOutlined />}
                      type={viewMode === 'grouped' ? 'primary' : 'default'}
                      onClick={() => setViewMode('grouped')}
                    />
                  </Space.Compact>
                }
                bodyStyle={{ padding: 12, display: 'flex', flexDirection: 'column', flex: 1, minHeight: 0, overflow: 'hidden' }}
                style={{ height: '100%', display: 'flex', flexDirection: 'column' }}
              >
                <Spin spinning={loading} style={{ flex: 1, display: 'flex', flexDirection: 'column' }}>
                  {/* 快速选择 */}
                  <Space wrap style={{ marginBottom: 12 }}>
                    <Tooltip title="选择出现频率最高的前3种配位组成">
                      <Button size="small" icon={<BulbOutlined />} loading={autoSelectLoading} onClick={handleSelectTop3}>
                        占比前3
                      </Button>
                    </Tooltip>
                    <Tooltip title="每种配位组成选1个代表性结构">
                      <Button size="small" loading={autoSelectLoading} onClick={handleSelectAll}>
                        每种1个
                      </Button>
                    </Tooltip>
                    <Button size="small" onClick={() => selectNPerGroup(3)}>每种3个</Button>
                    <Button size="small" type={filteredStructures.length !== structures.length ? 'primary' : 'default'}
                      onClick={() => setSelectedStructureIds(filteredStructures.map(s => s.id))}>
                      全选 ({filteredStructures.length})
                    </Button>
                    <Button size="small" danger onClick={() => setSelectedStructureIds([])}>清空</Button>
                  </Space>

                  {/* 筛选器 */}
                  <Space wrap style={{ marginBottom: 12 }}>
                    <Select size="small" mode="multiple" style={{ minWidth: 80 }} placeholder="CN"
                      value={filterCoordNums} onChange={setFilterCoordNums} allowClear maxTagCount={1}
                      options={filterOptions.coordNums.map(n => ({ value: n, label: `CN=${n}` }))} />
                    <Select size="small" mode="multiple" style={{ minWidth: 100 }} placeholder="阴离子"
                      value={filterAnions} onChange={setFilterAnions} allowClear maxTagCount={1}
                      options={filterOptions.anions.map(a => ({ value: a, label: a }))} />
                    <Select size="small" mode="multiple" style={{ minWidth: 100 }} placeholder="溶剂"
                      value={filterSolvents} onChange={setFilterSolvents} allowClear maxTagCount={1}
                      options={filterOptions.solvents.map(s => ({ value: s, label: s }))} />
                    {(filterCoordNums.length > 0 || filterAnions.length > 0 || filterSolvents.length > 0) && (
                      <Button size="small" type="link" onClick={resetFilters}>重置</Button>
                    )}
                  </Space>

                  {/* 结构列表 - 带滚动 */}
                  <div style={{ flex: 1, overflow: 'auto', minHeight: 200 }}>
                    {viewMode === 'grouped' ? (
                      <div style={{ paddingRight: 4 }}>
                        {sortedGroupKeys.filter(key => groupedStructures[key].structures.some(s => filteredStructures.includes(s))).map(groupKey => {
                          const group = groupedStructures[groupKey];
                          const isAllSelected = isGroupAllSelected(groupKey);
                          const isPartial = isGroupPartiallySelected(groupKey);
                          const selectedInGroup = group.structures.filter(s => selectedStructureIds.includes(s.id)).length;
                          return (
                            <div key={groupKey} style={{
                              padding: '8px 12px', marginBottom: 4, borderRadius: 8,
                              background: isAllSelected ? token.colorPrimaryBg : (isPartial ? token.colorWarningBg : token.colorBgLayout),
                              border: `1px solid ${isAllSelected ? token.colorPrimary : token.colorBorder}`,
                              cursor: 'pointer',
                            }}
                              onClick={() => selectGroup(groupKey, !isAllSelected)}
                            >
                              <Row justify="space-between" align="middle">
                                <Col>
                                  <Space>
                                    <Checkbox checked={isAllSelected} indeterminate={isPartial} />
                                    <Text strong style={{ fontFamily: 'monospace', fontSize: 13 }}>{groupKey}</Text>
                                  </Space>
                                </Col>
                                <Col>
                                  <Space size={4}>
                                    <Tag color="blue" style={{ margin: 0 }}>{group.count}</Tag>
                                    <Text type="secondary" style={{ fontSize: 12 }}>{group.percentage.toFixed(1)}%</Text>
                                    {selectedInGroup > 0 && <Tag color="green" style={{ margin: 0 }}>{selectedInGroup}✓</Tag>}
                                  </Space>
                                </Col>
                              </Row>
                            </div>
                          );
                        })}
                      </div>
                    ) : (
                      <Table
                        columns={structureColumns}
                        dataSource={filteredStructures}
                        rowKey="id"
                        size="small"
                        pagination={{ pageSize: 8, size: 'small', showTotal: (t) => `共${t}条`, showSizeChanger: true, pageSizeOptions: ['8', '15', '30', '50'] }}
                        scroll={{ x: 600, y: 'calc(100vh - 500px)' }}
                      />
                    )}
                  </div>
                </Spin>
              </Card>
            </Col>

            {/* 中栏：计算类型选择 */}
            <Col span={8} style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
              <Card
                size="small"
                title={<Space><ThunderboltOutlined /> 计算类型</Space>}
                bodyStyle={{ padding: 12, overflow: 'auto' }}
                style={{ height: '100%' }}
              >
                <Space direction="vertical" style={{ width: '100%' }} size={4}>
                  {CALC_TYPE_OPTIONS.map(opt => {
                    const isSelected = selectedCalcTypes.includes(opt.value);
                    const info = CALC_TYPE_INFO[opt.value];
                    return (
                      <div
                        key={opt.value}
                        style={{
                          padding: '8px 12px',
                          borderRadius: 6,
                          border: `1px solid ${isSelected ? token.colorPrimary : token.colorBorder}`,
                          background: isSelected ? token.colorPrimaryBg : 'transparent',
                          cursor: 'pointer',
                          transition: 'all 0.2s',
                        }}
                        onClick={() => {
                          if (isSelected) {
                            setSelectedCalcTypes(selectedCalcTypes.filter(t => t !== opt.value));
                          } else {
                            setSelectedCalcTypes([...selectedCalcTypes, opt.value]);
                          }
                        }}
                      >
                        <Row justify="space-between" align="middle">
                          <Col>
                            <Space>
                              <Checkbox checked={isSelected} />
                              <Text strong style={{ fontSize: 13 }}>{info.icon} {opt.label}</Text>
                            </Space>
                          </Col>
                          <Col>
                            <Space size={4}>
                              <Tag color={opt.riskLevel === 'high' ? 'red' : opt.riskLevel === 'medium' ? 'orange' : 'green'} style={{ margin: 0 }}>
                                {opt.riskLevel === 'high' ? '高' : opt.riskLevel === 'medium' ? '中' : '低'}
                              </Tag>
                              <Button
                                type="link"
                                size="small"
                                style={{ padding: 0, height: 'auto', fontSize: 11 }}
                                onClick={(e) => {
                                  e.stopPropagation();
                                  setSelectedCalcTypeForDetail(opt.value);
                                  setCalcTypeDetailVisible(true);
                                }}
                              >
                                详情
                              </Button>
                            </Space>
                          </Col>
                        </Row>
                        {/* 关键公式 - 使用 CALC_TYPE_INFO 中已有的 formula */}
                        <div style={{
                          marginTop: 4,
                          marginLeft: 24,
                          padding: '2px 6px',
                          background: token.colorBgLayout,
                          borderRadius: 4,
                          fontFamily: 'monospace',
                          fontSize: 11,
                          color: token.colorTextSecondary,
                          display: 'inline-block',
                        }}>
                          {info.formula}
                        </div>
                        {/* REDOX 子选项 */}
                        {opt.value === 'REDOX' && isSelected && (
                          <div style={{ marginTop: 8, marginLeft: 24, display: 'flex', gap: 12, flexWrap: 'wrap' }}
                               onClick={(e) => e.stopPropagation()}>
                            <Checkbox
                              checked={redoxOptions.include_molecule}
                              onChange={(e) => setRedoxOptions({ ...redoxOptions, include_molecule: e.target.checked })}
                            >
                              <Text style={{ fontSize: 12 }}>分子 Redox</Text>
                            </Checkbox>
                            <Checkbox
                              checked={redoxOptions.include_dimer}
                              onChange={(e) => setRedoxOptions({ ...redoxOptions, include_dimer: e.target.checked })}
                            >
                              <Text style={{ fontSize: 12 }}>Li-Dimer Redox</Text>
                            </Checkbox>
                            <Checkbox
                              checked={redoxOptions.include_cluster}
                              onChange={(e) => setRedoxOptions({ ...redoxOptions, include_cluster: e.target.checked })}
                            >
                              <Text style={{ fontSize: 12 }}>Cluster Redox</Text>
                            </Checkbox>
                          </div>
                        )}
                        {/* REORGANIZATION 子选项 */}
                        {opt.value === 'REORGANIZATION' && isSelected && (
                          <div style={{ marginTop: 8, marginLeft: 24, display: 'flex', gap: 12, flexWrap: 'wrap' }}
                               onClick={(e) => e.stopPropagation()}>
                            <Checkbox
                              checked={reorganizationOptions.include_molecule}
                              onChange={(e) => setReorganizationOptions({ ...reorganizationOptions, include_molecule: e.target.checked })}
                            >
                              <Text style={{ fontSize: 12 }}>分子重组能</Text>
                            </Checkbox>
                            <Checkbox
                              checked={reorganizationOptions.include_cluster}
                              onChange={(e) => setReorganizationOptions({ ...reorganizationOptions, include_cluster: e.target.checked })}
                            >
                              <Text style={{ fontSize: 12 }}>Cluster 重组能</Text>
                            </Checkbox>
                          </div>
                        )}
                      </div>
                    );
                  })}
                </Space>

                {selectedCalcTypes.some(t => ['REDOX', 'REORGANIZATION'].includes(t)) && (
                  <Alert
                    type="warning"
                    message="高风险计算，结果不确定性较大"
                    style={{ marginTop: 12 }}
                    showIcon
                  />
                )}
              </Card>
            </Col>

            {/* 右栏：计算参数 + 分析预览 */}
            <Col span={6} style={{ display: 'flex', flexDirection: 'column', height: '100%', gap: 12 }}>
              {/* 计算参数配置 */}
              <Card
                size="small"
                title={
                  <Row justify="space-between" align="middle" style={{ width: '100%' }}>
                    <Space><SettingOutlined /> 计算参数</Space>
                    <Tooltip title="根据计算类型和分子特性自动推荐最佳参数">
                      <Switch
                        size="small"
                        checked={useSmartRecommend}
                        onChange={(checked) => {
                          setUseSmartRecommend(checked);
                          if (checked) {
                            // 检测是否包含阴离子（需要弥散基组）
                            const selectedStructs = structures.filter(s => selectedStructureIds.includes(s.id));
                            const hasAnion = selectedStructs.some(s => {
                              const comp = s.composition || {};
                              // 常见阴离子名称
                              const anionPatterns = ['PF6', 'TFSI', 'FSI', 'BF4', 'ClO4', 'NO3', 'Cl', 'Br', 'I'];
                              return Object.keys(comp).some(mol =>
                                anionPatterns.some(pattern => mol.toUpperCase().includes(pattern))
                              );
                            });

                            const hasRedox = selectedCalcTypes.includes('REDOX');
                            const hasReorg = selectedCalcTypes.includes('REORGANIZATION');

                            // 根据分子特性选择基组
                            const recommendedBasis = hasAnion ? '6-31+G(d,p)' : '6-31G*';
                            const basisReason = hasAnion ? '(阴离子需弥散函数)' : '';

                            if (hasRedox || hasReorg) {
                              setQcConfig(prev => ({
                                ...prev,
                                functional: 'B3LYP',
                                basis_set: recommendedBasis,
                                use_dispersion: true,
                                solvent_model: 'pcm',
                                solvent: prev.solvent || 'Water',
                              }));
                              setRecommendReason(`Redox/重组能推荐 B3LYP/${recommendedBasis} + PCM溶剂 ${basisReason}`);
                            } else {
                              setQcConfig(prev => ({
                                ...prev,
                                functional: 'B3LYP',
                                basis_set: recommendedBasis,
                                use_dispersion: true,
                                solvent_model: 'pcm',
                                solvent: prev.solvent || 'Water',
                              }));
                              setRecommendReason(`Binding 推荐 B3LYP-D3BJ/${recommendedBasis} + PCM溶剂 ${basisReason}`);
                            }
                          } else {
                            setRecommendReason('');
                          }
                        }}
                        checkedChildren="智能"
                        unCheckedChildren="手动"
                      />
                    </Tooltip>
                  </Row>
                }
                bodyStyle={{ padding: 12 }}
              >
                {recommendReason && (
                  <Alert
                    type="info"
                    message={recommendReason}
                    style={{ marginBottom: 8, padding: '4px 8px', fontSize: 11 }}
                    showIcon
                  />
                )}
                <Row gutter={[8, 6]}>
                  <Col span={12}>
                    <Text type="secondary" style={{ fontSize: 11 }}>泛函</Text>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={qcConfig.functional}
                      onChange={(v) => setQcConfig(prev => ({ ...prev, functional: v }))}
                      showSearch
                    >
                      <Select.OptGroup label="常用">
                        <Select.Option value="B3LYP">B3LYP (推荐)</Select.Option>
                        <Select.Option value="PBE0">PBE0</Select.Option>
                        <Select.Option value="M06-2X">M06-2X (非共价)</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="长程校正">
                        <Select.Option value="wB97X-D">ωB97X-D</Select.Option>
                        <Select.Option value="CAM-B3LYP">CAM-B3LYP</Select.Option>
                        <Select.Option value="LC-wPBE">LC-ωPBE</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="Meta-GGA">
                        <Select.Option value="M06">M06</Select.Option>
                        <Select.Option value="M06-L">M06-L</Select.Option>
                        <Select.Option value="TPSS">TPSS</Select.Option>
                      </Select.OptGroup>
                    </Select>
                  </Col>
                  <Col span={12}>
                    <Text type="secondary" style={{ fontSize: 11 }}>基组</Text>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={qcConfig.basis_set}
                      onChange={(v) => setQcConfig(prev => ({ ...prev, basis_set: v }))}
                      showSearch
                    >
                      <Select.OptGroup label="Pople (推荐)">
                        <Select.Option value="6-31G*">6-31G* (标准)</Select.Option>
                        <Select.Option value="6-31+G(d,p)">6-31+G(d,p) (阴离子)</Select.Option>
                        <Select.Option value="6-311++G(d,p)">6-311++G(d,p) (高精度)</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="Karlsruhe">
                        <Select.Option value="def2-SVP">def2-SVP (快速)</Select.Option>
                        <Select.Option value="def2-TZVP">def2-TZVP (精确)</Select.Option>
                        <Select.Option value="def2-QZVP">def2-QZVP (极高精度)</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="Dunning">
                        <Select.Option value="cc-pVDZ">cc-pVDZ</Select.Option>
                        <Select.Option value="cc-pVTZ">cc-pVTZ</Select.Option>
                        <Select.Option value="aug-cc-pVDZ">aug-cc-pVDZ</Select.Option>
                      </Select.OptGroup>
                    </Select>
                  </Col>
                  <Col span={12}>
                    <Text type="secondary" style={{ fontSize: 11 }}>溶剂模型</Text>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={qcConfig.solvent_model}
                      onChange={(v) => setQcConfig(prev => ({ ...prev, solvent_model: v }))}
                    >
                      <Select.Option value="gas">气相 (推荐)</Select.Option>
                      <Select.Option value="pcm">PCM</Select.Option>
                      <Select.Option value="smd">SMD (精确)</Select.Option>
                    </Select>
                  </Col>
                  <Col span={12}>
                    <Text type="secondary" style={{ fontSize: 11 }}>溶剂</Text>
                    <Select
                      size="small"
                      style={{ width: '100%' }}
                      value={qcConfig.solvent}
                      onChange={(v) => setQcConfig(prev => ({ ...prev, solvent: v }))}
                      disabled={qcConfig.solvent_model === 'gas'}
                      showSearch
                      optionFilterProp="children"
                    >
                      <Select.OptGroup label="📌 水系电解液 (ε>50)">
                        <Select.Option value="Water">水 Water ε=78.4</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="📌 高介电 (ε=40-90)">
                        <Select.Option value="DiMethylSulfoxide">DMSO ε=46.8</Select.Option>
                        <Select.Option value="1,2-EthaneDiol">乙二醇 ε=40.2</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="📌 中介电 (ε=15-40)">
                        <Select.Option value="Acetonitrile">乙腈 ACN ε=35.7</Select.Option>
                        <Select.Option value="Methanol">甲醇 ε=32.6</Select.Option>
                        <Select.Option value="Ethanol">乙醇 ε=24.9</Select.Option>
                        <Select.Option value="Acetone">丙酮 ε=20.5</Select.Option>
                        <Select.Option value="1-Propanol">正丙醇 ε=20.5</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="📌 低介电 (ε<15) - 碳酸酯">
                        <Select.Option value="DiChloroEthane">二氯乙烷 ε=10.1</Select.Option>
                        <Select.Option value="Dichloromethane">二氯甲烷 ε=8.9</Select.Option>
                        <Select.Option value="TetraHydroFuran">THF ε=7.4</Select.Option>
                        <Select.Option value="Chloroform">氯仿 ε=4.7 (DMC参考)</Select.Option>
                        <Select.Option value="DiethylEther">乙醚 ε=4.2</Select.Option>
                        <Select.Option value="Toluene">甲苯 ε=2.4</Select.Option>
                        <Select.Option value="Benzene">苯 ε=2.3</Select.Option>
                      </Select.OptGroup>
                      <Select.OptGroup label="🔧 自定义">
                        <Select.Option value="custom">自定义溶剂参数...</Select.Option>
                      </Select.OptGroup>
                    </Select>
                  </Col>
                  {/* 自定义溶剂参数 */}
                  {qcConfig.solvent_model !== 'gas' && qcConfig.solvent === 'custom' && (
                    <>
                      <Col span={12}>
                        <Text type="secondary" style={{ fontSize: 11 }}>介电常数 ε</Text>
                        <InputNumber
                          size="small"
                          style={{ width: '100%' }}
                          value={qcConfig.custom_eps}
                          onChange={(v) => setQcConfig(prev => ({ ...prev, custom_eps: v || 78.4 }))}
                          min={1}
                          max={200}
                          step={0.1}
                          placeholder="如 78.4"
                        />
                      </Col>
                      <Col span={12}>
                        <Text type="secondary" style={{ fontSize: 11 }}>光学介电 ε∞</Text>
                        <InputNumber
                          size="small"
                          style={{ width: '100%' }}
                          value={qcConfig.custom_eps_inf}
                          onChange={(v) => setQcConfig(prev => ({ ...prev, custom_eps_inf: v || 1.78 }))}
                          min={1}
                          max={10}
                          step={0.01}
                          placeholder="如 1.78"
                        />
                      </Col>
                    </>
                  )}
                  <Col span={24}>
                    <Checkbox
                      checked={qcConfig.use_dispersion}
                      onChange={(e) => setQcConfig(prev => ({ ...prev, use_dispersion: e.target.checked }))}
                    >
                      <Text style={{ fontSize: 12 }}>色散校正 D3BJ (推荐)</Text>
                    </Checkbox>
                  </Col>
                </Row>
              </Card>

              {/* 分析预览 */}
              <Card
                size="small"
                title={<Space><CalculatorOutlined /> 分析预览</Space>}
                bodyStyle={{ padding: 12 }}
                style={{ flex: 1 }}
              >
                {/* 选择统计 - 更紧凑 */}
                <Row gutter={8} style={{ marginBottom: 12 }}>
                  <Col span={12}>
                    <div style={{ background: token.colorPrimaryBg, borderRadius: 6, padding: '8px 12px', textAlign: 'center' }}>
                      <Text type="secondary" style={{ fontSize: 11 }}>已选结构</Text>
                      <div><Text strong style={{ fontSize: 18 }}>{selectedStructureIds.length}</Text><Text type="secondary"> / {structures.length}</Text></div>
                    </div>
                  </Col>
                  <Col span={12}>
                    <div style={{ background: token.colorBgLayout, borderRadius: 6, padding: '8px 12px', textAlign: 'center' }}>
                      <Text type="secondary" style={{ fontSize: 11 }}>覆盖组成</Text>
                      <div><Text strong style={{ fontSize: 18 }}>{coveredCompositions}</Text><Text type="secondary"> / {sortedGroupKeys.length}</Text></div>
                    </div>
                  </Col>
                </Row>

                <Divider style={{ margin: '8px 0' }}>预估 QC 任务</Divider>

                {selectedCalcTypes.length > 0 && selectedStructureIds.length > 0 ? (
                  <>
                    {/* 共享的基础分子能量 */}
                    {estimatedQCTasks.baseMonomerCount > 0 && (
                      <Row justify="space-between" style={{ marginBottom: 4, background: token.colorSuccessBg, padding: '2px 6px', borderRadius: 4 }}>
                        <Col><Text style={{ fontSize: 11 }}>🔗 共享分子能量</Text></Col>
                        <Col><Text style={{ fontSize: 11 }}>{estimatedQCTasks.baseMonomerCount}</Text></Col>
                      </Row>
                    )}
                    {/* 各类型独有任务 */}
                    {Object.entries(estimatedQCTasks.details).map(([calcType, count]) => (
                      <Row key={calcType} justify="space-between" style={{ marginBottom: 4 }}>
                        <Col><Text style={{ fontSize: 12 }}>{CALC_TYPE_INFO[calcType as ClusterCalcType]?.icon} {CALC_TYPE_INFO[calcType as ClusterCalcType]?.label}</Text></Col>
                        <Col><Text strong>{count}</Text></Col>
                      </Row>
                    ))}
                    <Divider style={{ margin: '6px 0' }} />
                    <Row justify="space-between">
                      <Col><Text strong>总计</Text></Col>
                      <Col><Text strong style={{ color: token.colorPrimary, fontSize: 16 }}>{estimatedQCTasks.total}</Text></Col>
                    </Row>
                  </>
                ) : (
                  <Empty description="选择结构和计算类型后显示" image={Empty.PRESENTED_IMAGE_SIMPLE} />
                )}
              </Card>
            </Col>
          </Row>

        </div>
      );
  };

  // 渲染查看模式
  const renderViewMode = () => {
    if (!job) return <Spin spinning={loading}><Empty description="加载中..." /></Spin>;

    const statusCfg = STATUS_CONFIG[job.status] || STATUS_CONFIG.CREATED;

    return (
      <>
        {/* 任务状态卡片 */}
        <Card style={{ marginBottom: 24 }}>
          <Row gutter={24} align="middle">
            <Col flex="auto">
              <Space direction="vertical">
                <Space>
                  <Title level={4} style={{ margin: 0 }}>分析任务 #{job.id}</Title>
                  <Tag icon={statusCfg.icon} color={statusCfg.color}>{statusCfg.text}</Tag>
                </Space>
                <Space>
                  <Text type="secondary">MD Job:</Text>
                  <Button
                    type="link"
                    size="small"
                    icon={<RocketOutlined />}
                    onClick={() => navigate(`/workspace/liquid-electrolyte/md/${job.md_job_id}`)}
                  >
                    #{job.md_job_id}
                  </Button>
                  <Divider type="vertical" />
                  <Text type="secondary">创建时间: {dayjs(job.created_at).format('YYYY-MM-DD HH:mm')}</Text>
                </Space>
              </Space>
            </Col>
            <Col>
              <Row gutter={16}>
                <Col>
                  <Statistic title="选中结构" value={job.selected_structures?.count || 0} />
                </Col>
                <Col>
                  <Statistic
                    title="QC 进度"
                    value={job.qc_task_plan?.completed_qc_tasks || 0}
                    suffix={`/ ${job.qc_task_plan?.total_qc_tasks || 0}`}
                  />
                </Col>
                <Col>
                  <Statistic
                    title="总进度"
                    value={job.progress}
                    suffix="%"
                    valueStyle={{ color: job.progress === 100 ? '#52c41a' : undefined }}
                  />
                </Col>
              </Row>
            </Col>
          </Row>

          {job.error_message && (
            <Alert
              type="error"
              message="错误信息"
              description={job.error_message}
              style={{ marginTop: 16 }}
            />
          )}
        </Card>

        {/* 计算类型 */}
        <Card title="计算类型" size="small" style={{ marginBottom: 24 }}>
          <Space wrap>
            {job.calc_types.map(ct => {
              const info = CALC_TYPE_INFO[ct as ClusterCalcType];
              return (
                <Tag key={ct} color={info?.riskLevel === 'high' ? 'red' : info?.riskLevel === 'medium' ? 'orange' : 'blue'}>
                  {info?.icon} {info?.label || ct}
                </Tag>
              );
            })}
          </Space>
        </Card>

        {/* 结果展示 */}
        {job.status === 'COMPLETED' && (
          <ClusterAnalysisResultsPanel jobId={job.id} />
        )}

        {/* 等待中提示 */}
        {['SUBMITTED', 'RUNNING', 'WAITING_QC', 'CALCULATING'].includes(job.status) && (
          <Card>
            <div style={{ textAlign: 'center', padding: 40 }}>
              <Spin size="large" />
              <Title level={4} style={{ marginTop: 16 }}>任务运行中...</Title>
              <Text type="secondary">QC 任务进度: {job.qc_task_plan?.completed_qc_tasks || 0} / {job.qc_task_plan?.total_qc_tasks || 0}</Text>
              <div style={{ marginTop: 16 }}>
                <Progress
                  percent={job.progress}
                  status="active"
                  strokeColor={{ '0%': '#108ee9', '100%': '#87d068' }}
                />
              </div>
            </div>
          </Card>
        )}
      </>
    );
  };

  return (
    <div style={{ padding: 24 }}>
      {/* 返回按钮 */}
      <div style={{ marginBottom: 16 }}>
        <Button
          icon={<ArrowLeftOutlined />}
          onClick={() => navigate('/workspace/liquid-electrolyte/analysis')}
        >
          返回列表
        </Button>
      </div>

      {/* 页面标题 */}
      <Title level={4} style={{ marginBottom: 24 }}>
        {isCreateMode ? (
          <>
            <ExperimentOutlined style={{ marginRight: 8 }} />
            新建后处理分析
          </>
        ) : (
          <>
            <ExperimentOutlined style={{ marginRight: 8 }} />
            分析任务详情
          </>
        )}
      </Title>

      {isCreateMode ? renderCreateMode() : renderViewMode()}

      {/* 计算类型详情弹窗 */}
      <Modal
        title={selectedCalcTypeForDetail && (
          <Space>
            {CALC_TYPE_INFO[selectedCalcTypeForDetail]?.icon}
            {CALC_TYPE_OPTIONS.find(o => o.value === selectedCalcTypeForDetail)?.label}
            <Tag color={
              CALC_TYPE_OPTIONS.find(o => o.value === selectedCalcTypeForDetail)?.riskLevel === 'high' ? 'red' :
              CALC_TYPE_OPTIONS.find(o => o.value === selectedCalcTypeForDetail)?.riskLevel === 'medium' ? 'orange' : 'green'
            }>
              {CALC_TYPE_OPTIONS.find(o => o.value === selectedCalcTypeForDetail)?.riskLevel === 'high' ? '高风险' :
               CALC_TYPE_OPTIONS.find(o => o.value === selectedCalcTypeForDetail)?.riskLevel === 'medium' ? '中风险' : '低风险'}
            </Tag>
          </Space>
        )}
        open={calcTypeDetailVisible}
        onCancel={() => setCalcTypeDetailVisible(false)}
        footer={null}
        width={600}
      >
        {selectedCalcTypeForDetail && (() => {
          const info = CALC_TYPE_INFO[selectedCalcTypeForDetail];
          const extra = CALC_TYPE_EXTRA[selectedCalcTypeForDetail];
          return (
          <div>
            {/* 核心公式 */}
            <div style={{
              background: token.colorPrimaryBg,
              padding: 16,
              borderRadius: 8,
              marginBottom: 16,
              textAlign: 'center',
            }}>
              <Text style={{ fontFamily: 'monospace', fontSize: 16 }}>
                {info.formula}
              </Text>
            </div>

            {/* 物理意义 */}
            <div style={{ marginBottom: 16 }}>
              <Text strong>📖 物理意义</Text>
              <Paragraph style={{ marginTop: 8, marginBottom: 0 }}>
                {extra?.meaning || info.description}
              </Paragraph>
            </div>

            {/* 复用逻辑 */}
            <div style={{ marginBottom: 16 }}>
              <Text strong>🔗 复用逻辑</Text>
              <Paragraph style={{ marginTop: 8, marginBottom: 0 }}>
                {extra?.reuse || '按唯一组成复用计算结果'}
              </Paragraph>
            </div>

            {/* 计算路径/流程图 */}
            {extra?.diagram && (
              <div style={{ marginBottom: 16 }}>
                <Text strong>🌳 {extra.diagramTitle || '计算流程'}</Text>
                <pre style={{
                  marginTop: 8,
                  padding: 12,
                  background: token.colorBgLayout,
                  borderRadius: 6,
                  fontSize: 11,
                  fontFamily: 'monospace',
                  whiteSpace: 'pre',
                  overflow: 'auto',
                  lineHeight: 1.4,
                }}>
                  {extra.diagram}
                </pre>
              </div>
            )}

            {/* QC 任务估算说明 */}
            <div style={{
              background: token.colorInfoBg,
              padding: 12,
              borderRadius: 6,
              border: `1px solid ${token.colorInfoBorder}`,
            }}>
              <Text strong style={{ display: 'block', marginBottom: 8 }}>💡 QC 任务估算</Text>
              {selectedCalcTypeForDetail === 'BINDING_TOTAL' && (
                <Text style={{ fontSize: 12 }}>
                  任务数 = 唯一组成数 + 共享分子能量（离子 + 配体种类）<br/>
                  例：2种组成 + Li⁺ + EC + DMC = 2 + 3 = 5 个任务
                </Text>
              )}
              {selectedCalcTypeForDetail === 'BINDING_PAIRWISE' && (
                <Text style={{ fontSize: 12 }}>
                  任务数 = 离子-配体对数 + 共享分子能量<br/>
                  例：Li-EC + Li-DMC + Li⁺ + EC + DMC = 2 + 3 = 5 个任务
                </Text>
              )}
              {selectedCalcTypeForDetail === 'DESOLVATION_STEPWISE' && (
                <Text style={{ fontSize: 12 }}>
                  任务数 = Σ(每个组成的中间态数) + 共享分子能量<br/>
                  例：Li·EC₂·DMC₂(8) + Li·EC₁·DMC₃(7) + 3 = 18 个任务
                </Text>
              )}
              {selectedCalcTypeForDetail === 'REDOX' && (
                <Text style={{ fontSize: 12 }}>
                  任务数 = 唯一组成数 × 4（氧化态优化 + 还原态优化 + 2×溶剂化）<br/>
                  例：2种组成 × 4 = 8 个任务
                </Text>
              )}
              {selectedCalcTypeForDetail === 'REORGANIZATION' && (
                <Text style={{ fontSize: 12 }}>
                  任务数 = 唯一组成数 × 4（2个几何优化 + 2个交叉单点）<br/>
                  例：2种组成 × 4 = 8 个任务
                </Text>
              )}
            </div>
          </div>
          );
        })()}
      </Modal>

      {/* 结构预览 Modal */}
      <Modal
        title={
          <Space>
            <EyeOutlined />
            <span>结构预览 - {previewData?.cluster_name || ''}</span>
          </Space>
        }
        open={previewVisible}
        onCancel={() => {
          setPreviewVisible(false);
          setPreviewData(null);
        }}
        width={800}
        footer={null}
        destroyOnClose
      >
        {previewLoading ? (
          <div style={{ textAlign: 'center', padding: 60 }}>
            <Spin tip="正在加载结构..." />
          </div>
        ) : previewData ? (
          <div>
            {/* 结构信息摘要 */}
            <Alert
              type="info"
              style={{ marginBottom: 16 }}
              message={
                <Space split={<Divider type="vertical" />}>
                  <span>中心离子: <strong>{previewData.center_ion}</strong></span>
                  <span>配位数: <strong>{previewData.ligands.length}</strong></span>
                  <span>总电荷: <strong>{previewData.total_charge}</strong></span>
                  <span>
                    配位组成: {Object.entries(previewData.composition)
                      .filter(([, v]) => v > 0)
                      .map(([k, v]) => `${k}×${v}`)
                      .join(', ')}
                  </span>
                </Space>
              }
            />

            <Row gutter={16}>
              {/* 左侧：结构列表 - 根据计算类型显示不同结构 */}
              <Col span={8}>
                <Card
                  size="small"
                  title={
                    <span>
                      可查看的结构
                      {previewCalcType && (
                        <Text type="secondary" style={{ fontSize: 12, marginLeft: 8 }}>
                          ({CALC_TYPE_INFO[previewCalcType as ClusterCalcType]?.label || previewCalcType})
                        </Text>
                      )}
                    </span>
                  }
                  style={{ height: 420, overflow: 'auto' }}
                >
                  <Space direction="vertical" style={{ width: '100%' }}>
                    {/* BINDING_TOTAL, DESOLVATION, REDOX, REORGANIZATION 需要完整 Cluster */}
                    {(!previewCalcType || ['BINDING_TOTAL', 'DESOLVATION_STEPWISE', 'REDOX', 'REORGANIZATION'].includes(previewCalcType)) && (
                      <Button
                        block
                        type={selectedPreviewTab === 'cluster' ? 'primary' : 'default'}
                        onClick={() => setSelectedPreviewTab('cluster')}
                      >
                        Full Cluster ({previewData.cluster.atom_count} atoms)
                      </Button>
                    )}

                    {/* 所有 Binding 类型都需要中心离子 */}
                    {(!previewCalcType || ['BINDING_TOTAL', 'BINDING_PAIRWISE', 'DESOLVATION_STEPWISE'].includes(previewCalcType)) && (
                      <Button
                        block
                        type={selectedPreviewTab === 'center_ion' ? 'primary' : 'default'}
                        onClick={() => setSelectedPreviewTab('center_ion')}
                      >
                        Center Ion ({previewData.center_ion})
                      </Button>
                    )}

                    {/* BINDING_PAIRWISE 和 REDOX 需要 Dimer 结构（Li + 配体）*/}
                    {['BINDING_PAIRWISE', 'REDOX'].includes(previewCalcType) && previewData.dimer_structures && previewData.dimer_structures.length > 0 && (
                      <>
                        <Divider style={{ margin: '8px 0' }}>Li-Ligand Dimer</Divider>
                        {previewData.dimer_structures.map((dimer: any, idx: number) => (
                          <Button
                            key={`dimer_${idx}`}
                            block
                            type={selectedPreviewTab === `dimer_${idx}` ? 'primary' : 'default'}
                            onClick={() => setSelectedPreviewTab(`dimer_${idx}`)}
                            style={{ background: selectedPreviewTab === `dimer_${idx}` ? undefined : '#e6f7e6' }}
                          >
                            {dimer.name} ({dimer.atom_count} atoms)
                          </Button>
                        ))}
                      </>
                    )}

                    {/* DESOLVATION_STEPWISE 需要 Cluster-minus 结构 */}
                    {previewCalcType === 'DESOLVATION_STEPWISE' && previewData.cluster_minus_structures && previewData.cluster_minus_structures.length > 0 && (
                      <>
                        <Divider style={{ margin: '8px 0' }}>Cluster-minus</Divider>
                        {previewData.cluster_minus_structures.map((cm: any, idx: number) => (
                          <Button
                            key={`cluster_minus_${idx}`}
                            block
                            type={selectedPreviewTab === `cluster_minus_${idx}` ? 'primary' : 'default'}
                            onClick={() => setSelectedPreviewTab(`cluster_minus_${idx}`)}
                            style={{ background: selectedPreviewTab === `cluster_minus_${idx}` ? undefined : '#fff7e6' }}
                          >
                            -{cm.removed_ligand} ({cm.atom_count} atoms)
                          </Button>
                        ))}
                      </>
                    )}

                    {/* 单独配体 - Binding、REDOX、REORGANIZATION 类型都需要 */}
                    {(!previewCalcType || ['BINDING_TOTAL', 'BINDING_PAIRWISE', 'DESOLVATION_STEPWISE', 'REDOX', 'REORGANIZATION'].includes(previewCalcType)) && (
                      <>
                        <Divider style={{ margin: '8px 0' }}>Ligands</Divider>
                        {previewData.ligands.map((ligand: any, idx: number) => (
                          <Button
                            key={`ligand_${idx}`}
                            block
                            type={selectedPreviewTab === `ligand_${idx}` ? 'primary' : 'default'}
                            onClick={() => setSelectedPreviewTab(`ligand_${idx}`)}
                          >
                            {ligand.ligand_label} ({ligand.atom_count} atoms)
                          </Button>
                        ))}
                      </>
                    )}
                  </Space>
                </Card>
              </Col>

              {/* 右侧：3D 视图 */}
              <Col span={16}>
                <Card
                  size="small"
                  title={
                    selectedPreviewTab === 'cluster'
                      ? 'Full Cluster'
                      : selectedPreviewTab === 'center_ion'
                      ? 'Center Ion'
                      : selectedPreviewTab.startsWith('dimer_')
                      ? `Dimer: ${previewData.dimer_structures?.[parseInt(selectedPreviewTab.replace('dimer_', ''), 10)]?.name}`
                      : selectedPreviewTab.startsWith('cluster_minus_')
                      ? `Cluster-minus: -${previewData.cluster_minus_structures?.[parseInt(selectedPreviewTab.replace('cluster_minus_', ''), 10)]?.removed_ligand}`
                      : selectedPreviewTab.startsWith('ligand_')
                      ? `Ligand: ${previewData.ligands[parseInt(selectedPreviewTab.replace('ligand_', ''), 10)]?.ligand_label}`
                      : '3D Structure'
                  }
                  style={{ height: 420 }}
                >
                  <div
                    ref={previewViewerRef}
                    style={{
                      width: '100%',
                      height: 340,
                      background: isDark ? '#1a1a1a' : '#f8f9fa',
                      borderRadius: 8,
                    }}
                  />
                  <div style={{
                    marginTop: 8,
                    textAlign: 'center',
                    color: token.colorTextSecondary,
                    fontSize: 12,
                  }}>
                    拖动旋转 | 滚轮缩放 | 右键平移 | 红色球体为中心离子
                  </div>
                </Card>
              </Col>
            </Row>
          </div>
        ) : null}
      </Modal>
    </div>
  );
}
