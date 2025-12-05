/**
 * QC 数据标签页组件
 * 用于数据管理和公开搜索页面
 */
import { useState, useEffect } from 'react';
import {
  Card,
  Form,
  Input,
  Select,
  Button,
  Table,
  Space,
  Tag,
  message,
  Modal,
  Tooltip,
  Row,
  Col,
  Radio,
  Alert,
  Descriptions,
  Statistic,
  Divider,
  Image,
  Spin,
  Tabs,
  Typography,
  theme,
} from 'antd';

const { Text } = Typography;
import {
  SearchOutlined,
  ReloadOutlined,
  EyeOutlined,
  ThunderboltOutlined,
  DownloadOutlined,
  EditOutlined,
  ExperimentOutlined,
  BarChartOutlined,
  PlusOutlined,
  DeleteOutlined,
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import type { QCJob } from '../types/qc';
import { getQCJobs, getESPImage, getHOMOImage, getLUMOImage } from '../api/qc';
import { renderSolventModel } from '../utils/qc';
import QCRecalculateModal from './QCRecalculateModal';
import { useThemeStore } from '../stores/themeStore';

interface QCDataTabProps {
  isPublic?: boolean; // 是否为公开搜索
}

export default function QCDataTab({ isPublic = false }: QCDataTabProps) {
  const [form] = Form.useForm();
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(false);
  const [data, setData] = useState<QCJob[]>([]);
  const [total, setTotal] = useState(0);
  const [currentPage, setCurrentPage] = useState(1);
  const [pageSize, setPageSize] = useState(10);
  const [recalculateModalVisible, setRecalculateModalVisible] = useState(false);
  const [selectedJob, setSelectedJob] = useState<QCJob | null>(null);

  // 详情模态框状态（用于公开数据查看）
  const [detailModalVisible, setDetailModalVisible] = useState(false);
  const [detailJob, setDetailJob] = useState<QCJob | null>(null);
  const [detailImages, setDetailImages] = useState<{esp?: string; homo?: string; lumo?: string}>({});
  const [detailImagesLoading, setDetailImagesLoading] = useState(false);

  // 对比分析相关状态
  const [comparisonModalVisible, setComparisonModalVisible] = useState(false);
  const [comparisonType, setComparisonType] = useState<'molecule' | 'functional' | 'basis_set' | 'solvent'>('molecule');
  const [comparisonStep, setComparisonStep] = useState(1); // 当前步骤

  // 对比参数
  const [selectedMolecules, setSelectedMolecules] = useState<string[]>([]); // 选择的分子（SMILES）
  const [selectedFunctionals, setSelectedFunctionals] = useState<string[]>([]);
  const [selectedBasisSets, setSelectedBasisSets] = useState<string[]>([]);
  const [selectedSolvents, setSelectedSolvents] = useState<Array<{model: string, solvent?: string}>>([]);

  // 可用选项（从已有任务中提取）
  const [availableMolecules, setAvailableMolecules] = useState<Array<{smiles: string, name: string}>>([]);
  const [availableFunctionals, setAvailableFunctionals] = useState<string[]>([]);
  const [availableBasisSets, setAvailableBasisSets] = useState<string[]>([]);
  const [availableSolvents, setAvailableSolvents] = useState<Array<{model: string, solvent?: string, label: string}>>([]);

  // 对比结果
  const [comparisonResults, setComparisonResults] = useState<QCJob[]>([]);

  // 加载数据
  const loadData = async (page = 1, size = 10, searchParams = {}) => {
    setLoading(true);
    try {
      const params: any = {
        skip: (page - 1) * size,
        limit: size,
        ...searchParams,
      };

      // 如果是公开搜索，只显示公开数据
      if (isPublic) {
        params.visibility = 'PUBLIC';
        params.status = 'COMPLETED'; // 只显示已完成的
      }

      const response = await getQCJobs(params);
      setData(response.jobs);
      setTotal(response.total);
      setCurrentPage(page);
      setPageSize(size);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  // 初始加载
  useEffect(() => {
    loadData();
  }, [isPublic]);

  // 从已有任务中提取可用选项（仅已完成的任务）
  useEffect(() => {
    if (data.length === 0) return;

    // 只考虑已完成的任务
    const completedJobs = data.filter(job => job.status === 'COMPLETED');

    // 提取分子
    const moleculesMap = new Map<string, string>();
    completedJobs.forEach(job => {
      if (job.smiles && job.molecule_name) {
        moleculesMap.set(job.smiles, job.molecule_name);
      }
    });
    setAvailableMolecules(Array.from(moleculesMap.entries()).map(([smiles, name]) => ({ smiles, name })));
  }, [data]);

  // 根据已选分子动态更新可用的泛函、基组、溶剂选项
  useEffect(() => {
    if (selectedMolecules.length === 0) {
      setAvailableFunctionals([]);
      setAvailableBasisSets([]);
      setAvailableSolvents([]);
      return;
    }

    const completedJobs = data.filter(job =>
      job.status === 'COMPLETED' && selectedMolecules.includes(job.smiles)
    );

    // 提取该分子的所有泛函
    const functionalsSet = new Set<string>();
    completedJobs.forEach(job => {
      if (job.functional) functionalsSet.add(job.functional);
    });
    setAvailableFunctionals(Array.from(functionalsSet));

    // 提取该分子的所有基组
    const basisSetsSet = new Set<string>();
    completedJobs.forEach(job => {
      if (job.basis_set) basisSetsSet.add(job.basis_set);
    });
    setAvailableBasisSets(Array.from(basisSetsSet));

    // 提取该分子的所有溶剂
    const solventsMap = new Map<string, {model: string, solvent?: string, label: string}>();
    completedJobs.forEach(job => {
      const solventConfig = job.solvent_config || job.config?.solvent_config;
      const model = solventConfig?.model || 'gas';
      const solvent = solventConfig?.solvent_name;
      const key = `${model}_${solvent || ''}`;
      const label = solvent ? `${model.toUpperCase()} / ${solvent}` : model.toUpperCase();
      if (!solventsMap.has(key)) {
        solventsMap.set(key, { model, solvent, label });
      }
    });
    setAvailableSolvents(Array.from(solventsMap.values()));
  }, [selectedMolecules, data]);

  // 搜索
  const handleSearch = () => {
    const values = form.getFieldsValue();
    const searchParams: any = {};

    if (values.molecule_name) searchParams.molecule_name = values.molecule_name;
    if (values.smiles) searchParams.smiles = values.smiles;
    if (values.functional) searchParams.functional = values.functional;
    if (values.basis_set) searchParams.basis_set = values.basis_set;
    if (values.status) searchParams.status = values.status;
    if (values.visibility) searchParams.visibility = values.visibility;

    loadData(1, pageSize, searchParams);
  };

  // 开始对比分析
  const handleStartComparison = () => {
    // 重置状态
    setComparisonStep(1);
    setSelectedMolecules([]);
    setSelectedFunctionals([]);
    setSelectedBasisSets([]);
    setSelectedSolvents([]);
    setComparisonResults([]);
    setComparisonModalVisible(true);
  };

  // 切换对比类型
  const handleComparisonTypeChange = (type: 'molecule' | 'functional' | 'basis_set' | 'solvent') => {
    setComparisonType(type);
    setComparisonStep(1);
    setSelectedMolecules([]);
    setSelectedFunctionals([]);
    setSelectedBasisSets([]);
    setSelectedSolvents([]);
    setComparisonResults([]);
  };

  // 执行对比查询
  const handleExecuteComparison = () => {
    let filteredJobs: QCJob[] = [];

    switch (comparisonType) {
      case 'molecule':
        // 分子对比：固定泛函、基组、溶剂，对比不同分子
        if (selectedFunctionals.length === 0 || selectedBasisSets.length === 0 || selectedSolvents.length === 0) {
          message.error('请选择固定参数');
          return;
        }
        if (selectedMolecules.length < 2) {
          message.error('请至少选择 2 个分子进行对比');
          return;
        }

        filteredJobs = data.filter(job => {
          if (job.status !== 'COMPLETED') return false;
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const jobSolventKey = `${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          const selectedSolventKey = `${selectedSolvents[0].model}_${selectedSolvents[0].solvent || ''}`;

          return selectedMolecules.includes(job.smiles) &&
                 job.functional === selectedFunctionals[0] &&
                 job.basis_set === selectedBasisSets[0] &&
                 jobSolventKey === selectedSolventKey;
        });
        break;

      case 'functional':
        // 泛函对比：固定分子、基组、溶剂，对比不同泛函
        if (selectedMolecules.length === 0 || selectedBasisSets.length === 0 || selectedSolvents.length === 0) {
          message.error('请选择固定参数');
          return;
        }
        if (selectedFunctionals.length < 2) {
          message.error('请至少选择 2 个泛函进行对比');
          return;
        }

        filteredJobs = data.filter(job => {
          if (job.status !== 'COMPLETED') return false;
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const jobSolventKey = `${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          const selectedSolventKey = `${selectedSolvents[0].model}_${selectedSolvents[0].solvent || ''}`;

          return job.smiles === selectedMolecules[0] &&
                 selectedFunctionals.includes(job.functional) &&
                 job.basis_set === selectedBasisSets[0] &&
                 jobSolventKey === selectedSolventKey;
        });
        break;

      case 'basis_set':
        // 基组对比：固定分子、泛函、溶剂，对比不同基组
        if (selectedMolecules.length === 0 || selectedFunctionals.length === 0 || selectedSolvents.length === 0) {
          message.error('请选择固定参数');
          return;
        }
        if (selectedBasisSets.length < 2) {
          message.error('请至少选择 2 个基组进行对比');
          return;
        }

        filteredJobs = data.filter(job => {
          if (job.status !== 'COMPLETED') return false;
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const jobSolventKey = `${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          const selectedSolventKey = `${selectedSolvents[0].model}_${selectedSolvents[0].solvent || ''}`;

          return job.smiles === selectedMolecules[0] &&
                 job.functional === selectedFunctionals[0] &&
                 selectedBasisSets.includes(job.basis_set) &&
                 jobSolventKey === selectedSolventKey;
        });
        break;

      case 'solvent':
        // 溶剂对比：固定分子、泛函、基组，对比不同溶剂
        if (selectedMolecules.length === 0 || selectedFunctionals.length === 0 || selectedBasisSets.length === 0) {
          message.error('请选择固定参数');
          return;
        }
        if (selectedSolvents.length < 2) {
          message.error('请至少选择 2 个溶剂进行对比');
          return;
        }

        filteredJobs = data.filter(job => {
          if (job.status !== 'COMPLETED') return false;
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const jobSolventKey = `${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;

          return job.smiles === selectedMolecules[0] &&
                 job.functional === selectedFunctionals[0] &&
                 job.basis_set === selectedBasisSets[0] &&
                 selectedSolvents.some(s => `${s.model}_${s.solvent || ''}` === jobSolventKey);
        });
        break;
    }

    if (filteredJobs.length === 0) {
      message.warning('没有找到符合条件的已完成任务');
      return;
    }

    console.log('对比任务数据:', filteredJobs);
    console.log('第一个任务的 results:', filteredJobs[0]?.results);
    console.log('第一个任务的 result[0]:', filteredJobs[0]?.results?.[0]);

    setComparisonResults(filteredJobs);
    setComparisonStep(comparisonType === 'molecule' ? 3 : 4);
    message.success(`找到 ${filteredJobs.length} 个任务进行对比`);
  };

  // 检查是否有可对比的任务
  const getComparisonAvailability = () => {
    const completedJobs = data.filter(job => job.status === 'COMPLETED');

    if (completedJobs.length === 0) {
      return {
        available: false,
        message: '暂无已完成的 QC 计算任务',
        suggestion: '请先创建并运行 QC 计算任务'
      };
    }

    switch (comparisonType) {
      case 'molecule':
        // 检查是否有相同参数的不同分子
        const paramGroups = new Map<string, Set<string>>();
        completedJobs.forEach(job => {
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const key = `${job.functional}_${job.basis_set}_${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          if (!paramGroups.has(key)) {
            paramGroups.set(key, new Set());
          }
          paramGroups.get(key)!.add(job.smiles);
        });
        const hasComparable = Array.from(paramGroups.values()).some(molecules => molecules.size >= 2);
        if (!hasComparable) {
          return {
            available: false,
            message: '没有找到可对比的分子',
            suggestion: '需要至少 2 个分子使用相同的泛函、基组和溶剂参数'
          };
        }
        break;

      case 'functional':
        // 检查是否有相同分子的不同泛函
        const molGroups = new Map<string, Set<string>>();
        completedJobs.forEach(job => {
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const key = `${job.smiles}_${job.basis_set}_${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          if (!molGroups.has(key)) {
            molGroups.set(key, new Set());
          }
          molGroups.get(key)!.add(job.functional);
        });
        const hasFunctionals = Array.from(molGroups.values()).some(functionals => functionals.size >= 2);
        if (!hasFunctionals) {
          return {
            available: false,
            message: '没有找到可对比的泛函',
            suggestion: '需要至少 2 个泛函计算相同的分子（基组和溶剂也要相同）'
          };
        }
        break;

      case 'basis_set':
        // 检查是否有相同分子的不同基组
        const basisGroups = new Map<string, Set<string>>();
        completedJobs.forEach(job => {
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const key = `${job.smiles}_${job.functional}_${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          if (!basisGroups.has(key)) {
            basisGroups.set(key, new Set());
          }
          basisGroups.get(key)!.add(job.basis_set);
        });
        const hasBasisSets = Array.from(basisGroups.values()).some(basisSets => basisSets.size >= 2);
        if (!hasBasisSets) {
          return {
            available: false,
            message: '没有找到可对比的基组',
            suggestion: '需要至少 2 个基组计算相同的分子（泛函和溶剂也要相同）'
          };
        }
        break;

      case 'solvent':
        // 检查是否有相同分子的不同溶剂
        const solventGroups = new Map<string, Set<string>>();
        completedJobs.forEach(job => {
          const solventConfig = job.solvent_config || job.config?.solvent_config;
          const key = `${job.smiles}_${job.functional}_${job.basis_set}`;
          if (!solventGroups.has(key)) {
            solventGroups.set(key, new Set());
          }
          const solventKey = `${solventConfig?.model || 'gas'}_${solventConfig?.solvent_name || ''}`;
          solventGroups.get(key)!.add(solventKey);
        });
        const hasSolvents = Array.from(solventGroups.values()).some(solvents => solvents.size >= 2);
        if (!hasSolvents) {
          return {
            available: false,
            message: '没有找到可对比的溶剂',
            suggestion: '需要至少 2 个溶剂计算相同的分子（泛函和基组也要相同）'
          };
        }
        break;
    }

    return { available: true, message: '', suggestion: '' };
  };

  // 渲染参数选择界面
  const renderParameterSelection = () => {
    const availability = getComparisonAvailability();

    if (!availability.available) {
      return (
        <Space direction="vertical" size="large" style={{ width: '100%', textAlign: 'center', padding: '40px 0' }}>
          <Alert
            type="warning"
            showIcon
            message={availability.message}
            description={
              <Space direction="vertical" size="small">
                <div>{availability.suggestion}</div>
                <div style={{ marginTop: 16 }}>
                  <Button type="primary" onClick={() => window.location.href = '/workspace/job-create'}>
                    创建批量 QC 计算任务
                  </Button>
                </div>
              </Space>
            }
          />
          <Button onClick={() => setComparisonStep(1)}>返回</Button>
        </Space>
      );
    }

    // 统一的参数选择流程
    return (
      <Space direction="vertical" size="large" style={{ width: '100%' }}>
        {/* 步骤指示 */}
        <Alert
          type="info"
          showIcon
          message={
            comparisonType === 'molecule' ? '分子对比：选择多个分子，使用相同的计算参数进行对比' :
            comparisonType === 'functional' ? '泛函对比：选择分子，对比不同泛函的计算结果' :
            comparisonType === 'basis_set' ? '基组对比：选择分子，对比不同基组的计算结果' :
            '溶剂对比：选择分子，对比不同溶剂环境的计算结果'
          }
        />

        {/* 步骤 1: 选择分子 */}
        <Card
          size="small"
          title={
            <Space>
              <span style={{
                display: 'inline-block',
                width: 24,
                height: 24,
                lineHeight: '24px',
                textAlign: 'center',
                borderRadius: '50%',
                background: '#1890ff',
                color: 'white',
                fontSize: 12,
                fontWeight: 'bold'
              }}>1</span>
              <span>选择分子</span>
            </Space>
          }
          style={{ borderColor: selectedMolecules.length > 0 ? '#52c41a' : '#d9d9d9' }}
        >
          <Select
            mode={comparisonType === 'molecule' ? 'multiple' : undefined}
            style={{ width: '100%' }}
            placeholder={comparisonType === 'molecule' ? '选择要对比的分子（至少 2 个）' : '选择分子'}
            value={comparisonType === 'molecule' ? selectedMolecules : selectedMolecules[0]}
            onChange={(value) => {
              if (comparisonType === 'molecule') {
                setSelectedMolecules(value as string[]);
              } else {
                setSelectedMolecules([value as string]);
              }
              // 清空后续选择
              setSelectedFunctionals([]);
              setSelectedBasisSets([]);
              setSelectedSolvents([]);
            }}
            options={availableMolecules.map(m => ({
              label: `${m.name} (${m.smiles})`,
              value: m.smiles
            }))}
            showSearch
            filterOption={(input, option) =>
              (option?.label ?? '').toLowerCase().includes(input.toLowerCase())
            }
          />
          {selectedMolecules.length > 0 && (
            <div style={{ marginTop: 8, color: '#52c41a', fontSize: 12 }}>
              ✓ 已选择 {selectedMolecules.length} 个分子
            </div>
          )}
        </Card>

        {/* 步骤 2: 根据对比类型选择参数 */}
        {selectedMolecules.length > 0 && renderComparisonTypeSpecificParams()}

        {/* 操作按钮 */}
        <div style={{ textAlign: 'right', marginTop: 16 }}>
          <Space>
            <Button onClick={() => setComparisonStep(1)}>上一步</Button>
            <Button
              type="primary"
              onClick={handleExecuteComparison}
              disabled={!canExecuteComparison()}
            >
              开始对比
            </Button>
          </Space>
        </div>
      </Space>
    );
  };

  // 检查是否可以执行对比
  const canExecuteComparison = () => {
    if (selectedMolecules.length === 0) return false;

    switch (comparisonType) {
      case 'molecule':
        return selectedMolecules.length >= 2 &&
               selectedFunctionals.length > 0 &&
               selectedBasisSets.length > 0 &&
               selectedSolvents.length > 0;
      case 'functional':
        return selectedMolecules.length > 0 &&
               selectedFunctionals.length >= 2 &&
               selectedBasisSets.length > 0 &&
               selectedSolvents.length > 0;
      case 'basis_set':
        return selectedMolecules.length > 0 &&
               selectedFunctionals.length > 0 &&
               selectedBasisSets.length >= 2 &&
               selectedSolvents.length > 0;
      case 'solvent':
        return selectedMolecules.length > 0 &&
               selectedFunctionals.length > 0 &&
               selectedBasisSets.length > 0 &&
               selectedSolvents.length >= 2;
      default:
        return false;
    }
  };

  // 渲染对比类型特定的参数选择
  const renderComparisonTypeSpecificParams = () => {
    switch (comparisonType) {
      case 'molecule':
        // 分子对比：选择固定的泛函、基组、溶剂
        return (
          <>
            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>2</span>
                  <span>选择固定的计算参数</span>
                </Space>
              }
              style={{ borderColor: selectedFunctionals.length > 0 && selectedBasisSets.length > 0 && selectedSolvents.length > 0 ? '#52c41a' : '#d9d9d9' }}
            >
              <Row gutter={[16, 16]}>
                <Col span={8}>
                  <div style={{ marginBottom: 8 }}><strong>泛函</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择泛函"
                    value={selectedFunctionals[0]}
                    onChange={(value) => setSelectedFunctionals([value])}
                    options={availableFunctionals.map(f => ({ label: f, value: f }))}
                    disabled={availableFunctionals.length === 0}
                  />
                  {availableFunctionals.length === 0 && (
                    <div style={{ marginTop: 4, color: '#ff4d4f', fontSize: 12 }}>
                      所选分子没有可用的泛函
                    </div>
                  )}
                </Col>
                <Col span={8}>
                  <div style={{ marginBottom: 8 }}><strong>基组</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择基组"
                    value={selectedBasisSets[0]}
                    onChange={(value) => setSelectedBasisSets([value])}
                    options={availableBasisSets.map(b => ({ label: b, value: b }))}
                    disabled={availableBasisSets.length === 0}
                  />
                  {availableBasisSets.length === 0 && (
                    <div style={{ marginTop: 4, color: '#ff4d4f', fontSize: 12 }}>
                      所选分子没有可用的基组
                    </div>
                  )}
                </Col>
                <Col span={8}>
                  <div style={{ marginBottom: 8 }}><strong>溶剂</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择溶剂"
                    value={selectedSolvents[0] ? `${selectedSolvents[0].model}_${selectedSolvents[0].solvent || ''}` : undefined}
                    onChange={(value) => {
                      const solvent = availableSolvents.find(s => `${s.model}_${s.solvent || ''}` === value);
                      if (solvent) setSelectedSolvents([solvent]);
                    }}
                    options={availableSolvents.map(s => ({
                      label: s.label,
                      value: `${s.model}_${s.solvent || ''}`
                    }))}
                    disabled={availableSolvents.length === 0}
                  />
                  {availableSolvents.length === 0 && (
                    <div style={{ marginTop: 4, color: '#ff4d4f', fontSize: 12 }}>
                      所选分子没有可用的溶剂
                    </div>
                  )}
                </Col>
              </Row>
            </Card>
          </>
        );

      case 'functional':
        // 泛函对比：选择固定的基组、溶剂，选择多个泛函
        return (
          <>
            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>2</span>
                  <span>选择固定的计算参数</span>
                </Space>
              }
              style={{ borderColor: selectedBasisSets.length > 0 && selectedSolvents.length > 0 ? '#52c41a' : '#d9d9d9' }}
            >
              <Row gutter={[16, 16]}>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}><strong>基组</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择基组"
                    value={selectedBasisSets[0]}
                    onChange={(value) => setSelectedBasisSets([value])}
                    options={availableBasisSets.map(b => ({ label: b, value: b }))}
                    disabled={availableBasisSets.length === 0}
                  />
                </Col>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}><strong>溶剂</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择溶剂"
                    value={selectedSolvents[0] ? `${selectedSolvents[0].model}_${selectedSolvents[0].solvent || ''}` : undefined}
                    onChange={(value) => {
                      const solvent = availableSolvents.find(s => `${s.model}_${s.solvent || ''}` === value);
                      if (solvent) setSelectedSolvents([solvent]);
                    }}
                    options={availableSolvents.map(s => ({
                      label: s.label,
                      value: `${s.model}_${s.solvent || ''}`
                    }))}
                    disabled={availableSolvents.length === 0}
                  />
                </Col>
              </Row>
            </Card>

            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>3</span>
                  <span>选择要对比的泛函（至少 2 个）</span>
                </Space>
              }
              style={{ borderColor: selectedFunctionals.length >= 2 ? '#52c41a' : '#d9d9d9' }}
            >
              <Select
                mode="multiple"
                style={{ width: '100%' }}
                placeholder="选择泛函"
                value={selectedFunctionals}
                onChange={setSelectedFunctionals}
                options={availableFunctionals.map(f => ({ label: f, value: f }))}
                disabled={availableFunctionals.length === 0}
              />
              {selectedFunctionals.length > 0 && (
                <div style={{ marginTop: 8, color: selectedFunctionals.length >= 2 ? '#52c41a' : '#faad14', fontSize: 12 }}>
                  {selectedFunctionals.length >= 2 ? '✓' : '⚠'} 已选择 {selectedFunctionals.length} 个泛函
                  {selectedFunctionals.length < 2 && '（至少需要 2 个）'}
                </div>
              )}
            </Card>
          </>
        );

      case 'basis_set':
        // 基组对比：选择固定的泛函、溶剂，选择多个基组
        return (
          <>
            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>2</span>
                  <span>选择固定的计算参数</span>
                </Space>
              }
              style={{ borderColor: selectedFunctionals.length > 0 && selectedSolvents.length > 0 ? '#52c41a' : '#d9d9d9' }}
            >
              <Row gutter={[16, 16]}>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}><strong>泛函</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择泛函"
                    value={selectedFunctionals[0]}
                    onChange={(value) => setSelectedFunctionals([value])}
                    options={availableFunctionals.map(f => ({ label: f, value: f }))}
                    disabled={availableFunctionals.length === 0}
                  />
                </Col>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}><strong>溶剂</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择溶剂"
                    value={selectedSolvents[0] ? `${selectedSolvents[0].model}_${selectedSolvents[0].solvent || ''}` : undefined}
                    onChange={(value) => {
                      const solvent = availableSolvents.find(s => `${s.model}_${s.solvent || ''}` === value);
                      if (solvent) setSelectedSolvents([solvent]);
                    }}
                    options={availableSolvents.map(s => ({
                      label: s.label,
                      value: `${s.model}_${s.solvent || ''}`
                    }))}
                    disabled={availableSolvents.length === 0}
                  />
                </Col>
              </Row>
            </Card>

            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>3</span>
                  <span>选择要对比的基组（至少 2 个）</span>
                </Space>
              }
              style={{ borderColor: selectedBasisSets.length >= 2 ? '#52c41a' : '#d9d9d9' }}
            >
              <Select
                mode="multiple"
                style={{ width: '100%' }}
                placeholder="选择基组"
                value={selectedBasisSets}
                onChange={setSelectedBasisSets}
                options={availableBasisSets.map(b => ({ label: b, value: b }))}
                disabled={availableBasisSets.length === 0}
              />
              {selectedBasisSets.length > 0 && (
                <div style={{ marginTop: 8, color: selectedBasisSets.length >= 2 ? '#52c41a' : '#faad14', fontSize: 12 }}>
                  {selectedBasisSets.length >= 2 ? '✓' : '⚠'} 已选择 {selectedBasisSets.length} 个基组
                  {selectedBasisSets.length < 2 && '（至少需要 2 个）'}
                </div>
              )}
            </Card>
          </>
        );

      case 'solvent':
        // 溶剂对比：选择固定的泛函、基组，选择多个溶剂
        return (
          <>
            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>2</span>
                  <span>选择固定的计算参数</span>
                </Space>
              }
              style={{ borderColor: selectedFunctionals.length > 0 && selectedBasisSets.length > 0 ? '#52c41a' : '#d9d9d9' }}
            >
              <Row gutter={[16, 16]}>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}><strong>泛函</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择泛函"
                    value={selectedFunctionals[0]}
                    onChange={(value) => setSelectedFunctionals([value])}
                    options={availableFunctionals.map(f => ({ label: f, value: f }))}
                    disabled={availableFunctionals.length === 0}
                  />
                </Col>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}><strong>基组</strong></div>
                  <Select
                    style={{ width: '100%' }}
                    placeholder="选择基组"
                    value={selectedBasisSets[0]}
                    onChange={(value) => setSelectedBasisSets([value])}
                    options={availableBasisSets.map(b => ({ label: b, value: b }))}
                    disabled={availableBasisSets.length === 0}
                  />
                </Col>
              </Row>
            </Card>

            <Card
              size="small"
              title={
                <Space>
                  <span style={{
                    display: 'inline-block',
                    width: 24,
                    height: 24,
                    lineHeight: '24px',
                    textAlign: 'center',
                    borderRadius: '50%',
                    background: '#1890ff',
                    color: 'white',
                    fontSize: 12,
                    fontWeight: 'bold'
                  }}>3</span>
                  <span>选择要对比的溶剂（至少 2 个）</span>
                </Space>
              }
              style={{ borderColor: selectedSolvents.length >= 2 ? '#52c41a' : '#d9d9d9' }}
            >
              <Select
                mode="multiple"
                style={{ width: '100%' }}
                placeholder="选择溶剂"
                value={selectedSolvents.map(s => `${s.model}_${s.solvent || ''}`)}
                onChange={(values) => {
                  const solvents = values.map(v => availableSolvents.find(s => `${s.model}_${s.solvent || ''}` === v)).filter(Boolean) as typeof availableSolvents;
                  setSelectedSolvents(solvents);
                }}
                options={availableSolvents.map(s => ({
                  label: s.label,
                  value: `${s.model}_${s.solvent || ''}`
                }))}
                disabled={availableSolvents.length === 0}
              />
              {selectedSolvents.length > 0 && (
                <div style={{ marginTop: 8, color: selectedSolvents.length >= 2 ? '#52c41a' : '#faad14', fontSize: 12 }}>
                  {selectedSolvents.length >= 2 ? '✓' : '⚠'} 已选择 {selectedSolvents.length} 个溶剂
                  {selectedSolvents.length < 2 && '（至少需要 2 个）'}
                </div>
              )}
            </Card>
          </>
        );

      default:
        return null;
    }
  };

  // 渲染单个对比项的图片
  const ComparisonItemImages = ({ job }: { job: QCJob }) => {
    const [espImageUrl, setEspImageUrl] = useState<string | null>(null);
    const [homoImageUrl, setHomoImageUrl] = useState<string | null>(null);
    const [lumoImageUrl, setLumoImageUrl] = useState<string | null>(null);
    const [loading, setLoading] = useState(false);

    const result = job.results && job.results.length > 0 ? job.results[0] : null;

    useEffect(() => {
      if (!result?.id) return;

      const loadImages = async () => {
        setLoading(true);
        try {
          // 加载 ESP 图片
          if (result.esp_image_path) {
            try {
              const blob = await getESPImage(result.id);
              setEspImageUrl(URL.createObjectURL(blob));
            } catch (err) {
              console.error('加载ESP图片失败:', err);
            }
          }

          // 加载 HOMO 图片
          if (result.homo_image_path) {
            try {
              const blob = await getHOMOImage(result.id);
              setHomoImageUrl(URL.createObjectURL(blob));
            } catch (err) {
              console.error('加载HOMO图片失败:', err);
            }
          }

          // 加载 LUMO 图片
          if (result.lumo_image_path) {
            try {
              const blob = await getLUMOImage(result.id);
              setLumoImageUrl(URL.createObjectURL(blob));
            } catch (err) {
              console.error('加载LUMO图片失败:', err);
            }
          }
        } finally {
          setLoading(false);
        }
      };

      loadImages();

      // 清理 URL
      return () => {
        if (espImageUrl) URL.revokeObjectURL(espImageUrl);
        if (homoImageUrl) URL.revokeObjectURL(homoImageUrl);
        if (lumoImageUrl) URL.revokeObjectURL(lumoImageUrl);
      };
    }, [result?.id]);

    if (loading) {
      return (
        <div style={{ textAlign: 'center', padding: 20 }}>
          <Spin tip="加载图片中..." />
        </div>
      );
    }

    return (
      <Tabs
        defaultActiveKey="structure"
        size="small"
        items={[
          {
            key: 'structure',
            label: '分子结构',
            children: (
              <div style={{ textAlign: 'center', padding: 12, background: token.colorBgContainer, borderRadius: 8 }}>
                {job.smiles ? (
                  <img
                    src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(job.smiles)}/PNG?image_size=200x200`}
                    alt={job.molecule_name}
                    style={{ maxWidth: '100%', maxHeight: 200, objectFit: 'contain' }}
                    onError={(e) => {
                      const target = e.target as HTMLImageElement;
                      target.style.display = 'none';
                      const parent = target.parentElement;
                      if (parent) {
                        parent.innerHTML = `<div style="color: #999; padding: 40px;">${job.molecule_name}</div>`;
                      }
                    }}
                  />
                ) : (
                  <div style={{ color: '#999', padding: 40 }}>{job.molecule_name}</div>
                )}
              </div>
            ),
          },
          {
            key: 'homo',
            label: 'HOMO',
            children: homoImageUrl ? (
              <div style={{ textAlign: 'center', padding: 12, background: token.colorBgContainer, borderRadius: 8 }}>
                <Image
                  src={homoImageUrl}
                  alt="HOMO Orbital"
                  style={{ maxWidth: '100%', maxHeight: 200 }}
                  placeholder
                />
              </div>
            ) : (
              <div style={{ textAlign: 'center', padding: 40, color: token.colorTextSecondary, background: token.colorBgContainer, borderRadius: 8 }}>
                HOMO 图片暂未生成
              </div>
            ),
          },
          {
            key: 'lumo',
            label: 'LUMO',
            children: lumoImageUrl ? (
              <div style={{ textAlign: 'center', padding: 12, background: token.colorBgContainer, borderRadius: 8 }}>
                <Image
                  src={lumoImageUrl}
                  alt="LUMO Orbital"
                  style={{ maxWidth: '100%', maxHeight: 200 }}
                  placeholder
                />
              </div>
            ) : (
              <div style={{ textAlign: 'center', padding: 40, color: token.colorTextSecondary, background: token.colorBgContainer, borderRadius: 8 }}>
                LUMO 图片暂未生成
              </div>
            ),
          },
          {
            key: 'esp',
            label: 'ESP',
            children: espImageUrl ? (
              <div style={{ textAlign: 'center', padding: 12, background: token.colorBgContainer, borderRadius: 8 }}>
                <Image
                  src={espImageUrl}
                  alt="ESP Surface"
                  style={{ maxWidth: '100%', maxHeight: 200 }}
                  placeholder
                />
              </div>
            ) : (
              <div style={{ textAlign: 'center', padding: 40, color: token.colorTextSecondary, background: token.colorBgContainer, borderRadius: 8 }}>
                ESP 图片暂未生成
              </div>
            ),
          },
        ]}
      />
    );
  };

  // 渲染对比结果
  const renderComparisonResults = () => {
    if (comparisonResults.length === 0) return null;

    // 统计信息
    const getResultValue = (job: QCJob, field: 'homo_ev' | 'lumo_ev') => {
      const result = job.results && job.results.length > 0 ? job.results[0] : null;
      return result?.[field];
    };

    const homoValues = comparisonResults.map(j => getResultValue(j, 'homo_ev')).filter(v => v !== undefined) as number[];
    const lumoValues = comparisonResults.map(j => getResultValue(j, 'lumo_ev')).filter(v => v !== undefined) as number[];
    const gapValues = comparisonResults.map(j => {
      const homo = getResultValue(j, 'homo_ev');
      const lumo = getResultValue(j, 'lumo_ev');
      if (homo !== undefined && lumo !== undefined) {
        return lumo - homo;
      }
      return undefined;
    }).filter(v => v !== undefined) as number[];

    const stats = {
      homo: {
        min: homoValues.length > 0 ? Math.min(...homoValues) : 0,
        max: homoValues.length > 0 ? Math.max(...homoValues) : 0,
        avg: homoValues.length > 0 ? homoValues.reduce((a, b) => a + b, 0) / homoValues.length : 0,
      },
      lumo: {
        min: lumoValues.length > 0 ? Math.min(...lumoValues) : 0,
        max: lumoValues.length > 0 ? Math.max(...lumoValues) : 0,
        avg: lumoValues.length > 0 ? lumoValues.reduce((a, b) => a + b, 0) / lumoValues.length : 0,
      },
      gap: {
        min: gapValues.length > 0 ? Math.min(...gapValues) : 0,
        max: gapValues.length > 0 ? Math.max(...gapValues) : 0,
        avg: gapValues.length > 0 ? gapValues.reduce((a, b) => a + b, 0) / gapValues.length : 0,
      },
    };

    // 确定对比维度的列
    const variableColumn = comparisonType === 'molecule' ? {
      title: '分子名称',
      dataIndex: 'molecule_name',
      key: 'molecule_name',
      width: 150,
      fixed: 'left' as const,
      render: (text: string) => <Tag color="blue" style={{ fontSize: 14 }}>{text}</Tag>,
    } : comparisonType === 'functional' ? {
      title: '泛函',
      dataIndex: 'functional',
      key: 'functional',
      width: 120,
      fixed: 'left' as const,
      render: (text: string) => <Tag color="blue" style={{ fontSize: 14 }}>{text}</Tag>,
    } : comparisonType === 'basis_set' ? {
      title: '基组',
      dataIndex: 'basis_set',
      key: 'basis_set',
      width: 150,
      fixed: 'left' as const,
      render: (text: string) => <Tag color="blue" style={{ fontSize: 14 }}>{text}</Tag>,
    } : comparisonType === 'solvent' ? {
      title: '溶剂',
      key: 'solvent',
      width: 200,
      fixed: 'left' as const,
      render: (_: any, record: QCJob) => {
        const solventConfig = record.solvent_config || record.config?.solvent_config;
        return (
          <Space>
            <Tag color="purple">{solventConfig?.model || 'gas'}</Tag>
            {solventConfig?.solvent_name && <Tag color="green">{solventConfig.solvent_name}</Tag>}
          </Space>
        );
      },
    } : {
      title: '任务',
      dataIndex: 'id',
      key: 'id',
      width: 80,
      fixed: 'left' as const,
    };

    return (
      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
        {/* 统计信息 */}
        <Row gutter={16}>
          <Col span={8}>
            <Card size="small">
              <Statistic
                title="HOMO 能量范围 (eV)"
                value={`${stats.homo.min.toFixed(4)} ~ ${stats.homo.max.toFixed(4)}`}
                valueStyle={{ fontSize: 14 }}
                suffix={
                  <Tooltip title="平均值">
                    <span style={{ fontSize: 12, color: '#666' }}>
                      (avg: {stats.homo.avg.toFixed(4)})
                    </span>
                  </Tooltip>
                }
              />
            </Card>
          </Col>
          <Col span={8}>
            <Card size="small">
              <Statistic
                title="LUMO 能量范围 (eV)"
                value={`${stats.lumo.min.toFixed(4)} ~ ${stats.lumo.max.toFixed(4)}`}
                valueStyle={{ fontSize: 14 }}
                suffix={
                  <Tooltip title="平均值">
                    <span style={{ fontSize: 12, color: '#666' }}>
                      (avg: {stats.lumo.avg.toFixed(4)})
                    </span>
                  </Tooltip>
                }
              />
            </Card>
          </Col>
          <Col span={8}>
            <Card size="small">
              <Statistic
                title="能隙范围 (eV)"
                value={`${stats.gap.min.toFixed(4)} ~ ${stats.gap.max.toFixed(4)}`}
                valueStyle={{ fontSize: 14 }}
                suffix={
                  <Tooltip title="平均值">
                    <span style={{ fontSize: 12, color: '#666' }}>
                      (avg: {stats.gap.avg.toFixed(4)})
                    </span>
                  </Tooltip>
                }
              />
            </Card>
          </Col>
        </Row>

        <Divider style={{ margin: '12px 0' }} />

        {/* 对比卡片网格 */}
        <Row gutter={[16, 16]}>
          {comparisonResults.map((job) => {
            const result = job.results && job.results.length > 0 ? job.results[0] : null;
            const homoValue = result?.homo_ev;
            const lumoValue = result?.lumo_ev;
            const gapValue = homoValue !== undefined && lumoValue !== undefined ? lumoValue - homoValue : undefined;

            // 获取对比维度的显示值
            let comparisonLabel = '';
            let comparisonValue: any = '';

            switch (comparisonType) {
              case 'molecule':
                comparisonLabel = '分子';
                comparisonValue = <Tag color="blue" style={{ fontSize: 14 }}>{job.molecule_name}</Tag>;
                break;
              case 'functional':
                comparisonLabel = '泛函';
                comparisonValue = <Tag color="blue" style={{ fontSize: 14 }}>{job.functional}</Tag>;
                break;
              case 'basis_set':
                comparisonLabel = '基组';
                comparisonValue = <Tag color="blue" style={{ fontSize: 14 }}>{job.basis_set}</Tag>;
                break;
              case 'solvent':
                comparisonLabel = '溶剂';
                const solventConfig = job.solvent_config || job.config?.solvent_config;
                const model = solventConfig?.model || 'gas';
                const solvent = solventConfig?.solvent_name;
                comparisonValue = (
                  <Space size="small">
                    <Tag color="purple">{model.toUpperCase()}</Tag>
                    {solvent && <Tag color="cyan">{solvent}</Tag>}
                  </Space>
                );
                break;
            }

            return (
              <Col xs={24} sm={12} lg={8} xl={6} key={job.id}>
                <Card
                  size="small"
                  title={comparisonValue}
                  styles={{
                    header: { background: '#f0f5ff', borderBottom: '1px solid #d9d9d9' },
                    body: { padding: 12 }
                  }}
                >
                  {/* 图片展示 */}
                  <ComparisonItemImages job={job} />

                  <Divider style={{ margin: '12px 0' }} />

                  {/* 计算参数 */}
                  <Space direction="vertical" size="small" style={{ width: '100%', fontSize: 12 }}>
                    {comparisonType !== 'molecule' && (
                      <div>
                        <Text type="secondary">分子：</Text>
                        <Text strong>{job.molecule_name}</Text>
                      </div>
                    )}
                    {comparisonType !== 'functional' && (
                      <div>
                        <Text type="secondary">泛函：</Text>
                        <Tag color="blue" style={{ fontSize: 12 }}>{job.functional}</Tag>
                      </div>
                    )}
                    {comparisonType !== 'basis_set' && (
                      <div>
                        <Text type="secondary">基组：</Text>
                        <Tag color="green" style={{ fontSize: 12 }}>{job.basis_set}</Tag>
                      </div>
                    )}
                    {comparisonType !== 'solvent' && (() => {
                      const solventConfig = job.solvent_config || job.config?.solvent_config;
                      const model = solventConfig?.model || 'gas';
                      const solvent = solventConfig?.solvent_name;
                      return (
                        <div>
                          <Text type="secondary">溶剂：</Text>
                          <Space size="small">
                            <Tag color="purple" style={{ fontSize: 12 }}>{model.toUpperCase()}</Tag>
                            {solvent && <Tag color="cyan" style={{ fontSize: 12 }}>{solvent}</Tag>}
                          </Space>
                        </div>
                      );
                    })()}
                  </Space>

                  <Divider style={{ margin: '12px 0' }} />

                  {/* 计算结果 */}
                  <Space direction="vertical" size="small" style={{ width: '100%' }}>
                    <Row gutter={8}>
                      <Col span={12}>
                        <Statistic
                          title="HOMO"
                          value={homoValue?.toFixed(4) || '-'}
                          suffix="eV"
                          valueStyle={{
                            fontSize: 14,
                            fontWeight: homoValue === stats.homo.min || homoValue === stats.homo.max ? 'bold' : 'normal',
                            color: homoValue === stats.homo.min ? '#52c41a' : homoValue === stats.homo.max ? '#ff4d4f' : 'inherit'
                          }}
                        />
                      </Col>
                      <Col span={12}>
                        <Statistic
                          title="LUMO"
                          value={lumoValue?.toFixed(4) || '-'}
                          suffix="eV"
                          valueStyle={{
                            fontSize: 14,
                            fontWeight: lumoValue === stats.lumo.min || lumoValue === stats.lumo.max ? 'bold' : 'normal',
                            color: lumoValue === stats.lumo.min ? '#ff4d4f' : lumoValue === stats.lumo.max ? '#52c41a' : 'inherit'
                          }}
                        />
                      </Col>
                    </Row>
                    <Row gutter={8}>
                      <Col span={12}>
                        <Statistic
                          title="Gap"
                          value={gapValue?.toFixed(4) || '-'}
                          suffix="eV"
                          valueStyle={{
                            fontSize: 14,
                            fontWeight: gapValue && (Math.abs(gapValue - stats.gap.min) < 0.0001 || Math.abs(gapValue - stats.gap.max) < 0.0001) ? 'bold' : 'normal',
                            color: gapValue && Math.abs(gapValue - stats.gap.min) < 0.0001 ? '#ff4d4f' : gapValue && Math.abs(gapValue - stats.gap.max) < 0.0001 ? '#52c41a' : 'inherit'
                          }}
                        />
                      </Col>
                      <Col span={12}>
                        <Statistic
                          title="偶极矩"
                          value={result?.dipole_moment?.toFixed(3) || '-'}
                          suffix="D"
                          valueStyle={{ fontSize: 14 }}
                        />
                      </Col>
                    </Row>
                  </Space>
                </Card>
              </Col>
            );
          })}
        </Row>

        {/* 对比说明 */}
        <Alert
          type="success"
          showIcon
          message="对比说明"
          description={
            <ul style={{ margin: 0, paddingLeft: 20 }}>
              <li>绿色 ↑ 表示该参数的最大值（HOMO 越高越好，表示更容易失去电子）</li>
              <li>红色 ↓ 表示该参数的最小值（LUMO 越低越好，表示更容易接受电子）</li>
              <li>能隙 (Gap) 越大表示化学稳定性越好</li>
              <li>点击列标题可以按该列排序</li>
            </ul>
          }
        />
      </Space>
    );
  };

  // 导出对比数据
  const handleExportComparison = () => {
    const csv = [
      ['分子名称', 'SMILES', '泛函', '基组', '溶剂模型', '溶剂', 'HOMO (eV)', 'LUMO (eV)', 'Gap (eV)', '总能量 (Hartree)', '偶极矩 (Debye)'],
      ...comparisonResults.map(job => {
        const solventConfig = job.solvent_config || job.config?.solvent_config;
        const result = job.results && job.results.length > 0 ? job.results[0] : null;
        return [
          job.molecule_name,
          job.smiles,
          job.functional,
          job.basis_set,
          solventConfig?.model || 'gas',
          solventConfig?.solvent_name || '',
          result?.homo_ev?.toFixed(4) || '',
          result?.lumo_ev?.toFixed(4) || '',
          (result?.homo_ev !== undefined && result?.lumo_ev !== undefined) ? (result.lumo_ev - result.homo_ev).toFixed(4) : '',
          result?.energy_au?.toFixed(6) || '',
          result?.dipole_moment?.toFixed(4) || '',
        ];
      })
    ].map(row => row.join(',')).join('\n');

    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = `qc_comparison_${new Date().getTime()}.csv`;
    link.click();
    message.success('数据已导出');
  };

  // 重置
  const handleReset = () => {
    form.resetFields();
    loadData(1, pageSize);
  };

  // 查看详情
  const handleViewDetail = async (record: QCJob) => {
    if (isPublic) {
      // 公开数据：显示模态框
      setDetailJob(record);
      setDetailModalVisible(true);
      setDetailImagesLoading(true);
      setDetailImages({});

      // 加载图片
      try {
        const resultId = record.results && record.results.length > 0 ? record.results[0].id : null;
        if (resultId) {
          const [espImg, homoImg, lumoImg] = await Promise.allSettled([
            getESPImage(resultId),
            getHOMOImage(resultId),
            getLUMOImage(resultId),
          ]);

          setDetailImages({
            esp: espImg.status === 'fulfilled' && espImg.value ? URL.createObjectURL(espImg.value) : undefined,
            homo: homoImg.status === 'fulfilled' && homoImg.value ? URL.createObjectURL(homoImg.value) : undefined,
            lumo: lumoImg.status === 'fulfilled' && lumoImg.value ? URL.createObjectURL(lumoImg.value) : undefined,
          });
        }
      } catch (e) {
        console.warn('Failed to load some images:', e);
      } finally {
        setDetailImagesLoading(false);
      }
    } else {
      // 私有数据：跳转到详情页
      navigate(`/workspace/liquid-electrolyte/qc/${record.id}`);
    }
  };

  // 重新计算
  const handleRecalculate = (record: QCJob) => {
    setSelectedJob(record);
    setRecalculateModalVisible(true);
  };

  // 重新计算成功
  const handleRecalculateSuccess = (newJob: QCJob) => {
    setRecalculateModalVisible(false);
    setSelectedJob(null);
    message.success('重新计算任务已创建，正在跳转...');
    setTimeout(() => {
      navigate(`/workspace/liquid-electrolyte/qc/${newJob.id}`);
    }, 1000);
  };

  // 表格列定义
  const columns: any[] = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 80,
      fixed: 'left',
    },
    {
      title: '分子名称',
      dataIndex: 'molecule_name',
      key: 'molecule_name',
      width: 150,
      ellipsis: true,
    },
    {
      title: 'SMILES',
      dataIndex: 'smiles',
      key: 'smiles',
      width: 200,
      ellipsis: true,
      render: (text: string) => (
        <Tooltip title={text}>
          <code style={{ fontSize: 12 }}>{text}</code>
        </Tooltip>
      ),
    },
    {
      title: '泛函',
      dataIndex: 'functional',
      key: 'functional',
      width: 100,
    },
    {
      title: '基组',
      dataIndex: 'basis_set',
      key: 'basis_set',
      width: 120,
    },
    {
      title: '溶剂模型',
      dataIndex: 'config',
      key: 'solvent_model',
      width: 150,
      render: (config: any) => renderSolventModel(config),
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      width: 100,
      filters: [
        { text: '已创建', value: 'CREATED' },
        { text: '已提交', value: 'SUBMITTED' },
        { text: '运行中', value: 'RUNNING' },
        { text: '已完成', value: 'COMPLETED' },
        { text: '失败', value: 'FAILED' },
        { text: '已取消', value: 'CANCELLED' },
      ],
      onFilter: (value: any, record: any) => record.status === value,
      render: (status: string) => {
        const statusConfig: any = {
          CREATED: { color: 'default', text: '已创建' },
          SUBMITTED: { color: 'processing', text: '已提交' },
          RUNNING: { color: 'processing', text: '运行中' },
          COMPLETED: { color: 'success', text: '已完成' },
          FAILED: { color: 'error', text: '失败' },
          CANCELLED: { color: 'default', text: '已取消' },
        };
        const config = statusConfig[status] || { color: 'default', text: status };
        return <Tag color={config.color}>{config.text}</Tag>;
      },
    },
  ];

  // 如果不是公开搜索，添加可见性列
  if (!isPublic) {
    columns.push({
      title: '可见性',
      dataIndex: 'visibility',
      key: 'visibility',
      width: 100,
      render: (visibility: string) => {
        const visibilityConfig: any = {
          PUBLIC: { color: 'green', text: '公开' },
          DELAYED: { color: 'orange', text: '延迟公开' },
          PRIVATE: { color: 'red', text: '私有' },
        };
        const config = visibilityConfig[visibility] || { color: 'default', text: visibility };
        return <Tag color={config.color}>{config.text}</Tag>;
      },
    });
  }

  // 操作列
  columns.push({
    title: '操作',
    key: 'action',
    width: isPublic ? 150 : 250,
    fixed: 'right',
    render: (_: any, record: QCJob) => (
      <Space size="small">
        <Tooltip title="查看详情">
          <Button
            type="link"
            size="small"
            icon={<EyeOutlined />}
            onClick={() => handleViewDetail(record)}
          >
            详情
          </Button>
        </Tooltip>
        {!isPublic && record.status === 'COMPLETED' && (
          <Tooltip title="基于此任务重新计算">
            <Button
              type="link"
              size="small"
              icon={<ThunderboltOutlined />}
              onClick={() => handleRecalculate(record)}
              style={{ color: '#722ed1' }}
            >
              重算
            </Button>
          </Tooltip>
        )}
      </Space>
    ),
  });

  return (
    <Space direction="vertical" size="middle" style={{ width: '100%' }}>
      {/* 搜索表单 */}
      <Card size="small">
        <Form form={form} layout="vertical">
          <Row gutter={16}>
            <Col span={6}>
              <Form.Item name="molecule_name" label="分子名称">
                <Input placeholder="输入分子名称" allowClear />
              </Form.Item>
            </Col>
            <Col span={6}>
              <Form.Item name="smiles" label="SMILES">
                <Input placeholder="输入SMILES" allowClear />
              </Form.Item>
            </Col>
            <Col span={6}>
              <Form.Item name="functional" label="泛函">
                <Select placeholder="选择泛函" allowClear>
                  <Select.Option value="HF">HF</Select.Option>
                  <Select.Option value="B3LYP">B3LYP</Select.Option>
                  <Select.Option value="M062X">M06-2X</Select.Option>
                  <Select.Option value="wB97XD">ωB97X-D</Select.Option>
                  <Select.Option value="PBE0">PBE0</Select.Option>
                </Select>
              </Form.Item>
            </Col>
            <Col span={6}>
              <Form.Item name="basis_set" label="基组">
                <Select placeholder="选择基组" allowClear>
                  <Select.Option value="STO-3G">STO-3G</Select.Option>
                  <Select.Option value="6-31G(d)">6-31G(d)</Select.Option>
                  <Select.Option value="6-31G(d,p)">6-31G(d,p)</Select.Option>
                  <Select.Option value="6-31++G(d,p)">6-31++G(d,p)</Select.Option>
                  <Select.Option value="6-311G(d,p)">6-311G(d,p)</Select.Option>
                  <Select.Option value="6-311++G(d,p)">6-311++G(d,p)</Select.Option>
                  <Select.Option value="Def2-TZVP">Def2-TZVP</Select.Option>
                </Select>
              </Form.Item>
            </Col>
          </Row>

          {!isPublic && (
            <Row gutter={16}>
              <Col span={6}>
                <Form.Item name="status" label="状态">
                  <Select placeholder="选择状态" allowClear>
                    <Select.Option value="CREATED">已创建</Select.Option>
                    <Select.Option value="SUBMITTED">已提交</Select.Option>
                    <Select.Option value="RUNNING">运行中</Select.Option>
                    <Select.Option value="COMPLETED">已完成</Select.Option>
                    <Select.Option value="FAILED">失败</Select.Option>
                    <Select.Option value="CANCELLED">已取消</Select.Option>
                  </Select>
                </Form.Item>
              </Col>
              <Col span={6}>
                <Form.Item name="visibility" label="可见性">
                  <Select placeholder="选择可见性" allowClear>
                    <Select.Option value="PUBLIC">公开</Select.Option>
                    <Select.Option value="DELAYED">延迟公开</Select.Option>
                    <Select.Option value="PRIVATE">私有</Select.Option>
                  </Select>
                </Form.Item>
              </Col>
            </Row>
          )}

          <Row>
            <Col span={24} style={{ textAlign: 'right' }}>
              <Space>
                <Button onClick={handleReset}>
                  <ReloadOutlined /> 重置
                </Button>
                <Button type="primary" onClick={handleSearch} icon={<SearchOutlined />}>
                  搜索
                </Button>
                <Button
                  type="primary"
                  icon={<BarChartOutlined />}
                  onClick={handleStartComparison}
                >
                  对比分析
                </Button>
              </Space>
            </Col>
          </Row>
        </Form>
      </Card>

      {/* 对比分析 Modal */}
      <Modal
        title={
          <Space>
            <BarChartOutlined style={{ color: '#1890ff' }} />
            <span>QC 计算对比分析</span>
          </Space>
        }
        open={comparisonModalVisible}
        onCancel={() => {
          setComparisonModalVisible(false);
          setComparisonStep(1);
          setSelectedMolecules([]);
          setSelectedFunctionals([]);
          setSelectedBasisSets([]);
          setSelectedSolvents([]);
          setComparisonResults([]);
        }}
        footer={null}
        width={1200}
        centered
        styles={{
          body: { maxHeight: '80vh', overflowY: 'auto', padding: '24px' }
        }}
      >
        <Space direction="vertical" size="large" style={{ width: '100%' }}>
          {/* 步骤 1: 选择对比类型 */}
          {comparisonStep === 1 && (
            <Space direction="vertical" size="large" style={{ width: '100%' }}>
              <div style={{ textAlign: 'center', padding: '20px 0' }}>
                <h3 style={{ fontSize: 18, marginBottom: 8 }}>选择对比类型</h3>
                <Text type="secondary">根据您的研究需求，选择合适的对比维度</Text>
              </div>

              <Row gutter={[16, 16]}>
                {[
                  {
                    key: 'molecule',
                    icon: '🧪',
                    title: '分子对比',
                    description: '在相同计算条件下，对比不同分子的电子结构性质',
                    example: '例如：对比 Li+、Na+、K+ 在 B3LYP/6-31++G(d,p) 下的 HOMO/LUMO',
                  },
                  {
                    key: 'functional',
                    icon: '📊',
                    title: '泛函对比',
                    description: '评估不同泛函对同一体系计算结果的影响',
                    example: '例如：对比 Li+ 在 B3LYP、M06-2X、ωB97X-D 下的能隙差异',
                  },
                  {
                    key: 'basis_set',
                    icon: '📐',
                    title: '基组对比',
                    description: '评估基组大小对计算精度的影响',
                    example: '例如：对比 EC 在 6-31G(d,p)、6-31++G(d,p)、Def2-TZVP 下的收敛性',
                  },
                  {
                    key: 'solvent',
                    icon: '💧',
                    title: '溶剂对比',
                    description: '研究溶剂效应对分子性质的影响',
                    example: '例如：对比 EC 在气相、PCM/Water、SMD/Acetonitrile 下的偶极矩',
                  },
                ].map((type) => (
                  <Col xs={24} sm={12} key={type.key}>
                    <Card
                      hoverable
                      style={{
                        borderColor: comparisonType === type.key ? '#1890ff' : '#d9d9d9',
                        borderWidth: comparisonType === type.key ? 2 : 1,
                        background: comparisonType === type.key ? '#f0f5ff' : 'white',
                        height: '100%',
                        cursor: 'pointer',
                        transition: 'all 0.3s',
                      }}
                      onClick={() => handleComparisonTypeChange(type.key as any)}
                      styles={{
                        body: { padding: 20 }
                      }}
                    >
                      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
                        <div style={{ display: 'flex', alignItems: 'center', gap: 12 }}>
                          <span style={{ fontSize: 32 }}>{type.icon}</span>
                          <div>
                            <div style={{ fontSize: 16, fontWeight: 'bold', marginBottom: 4 }}>
                              {type.title}
                              {comparisonType === type.key && (
                                <Tag color="blue" style={{ marginLeft: 8 }}>已选择</Tag>
                              )}
                            </div>
                            <Text type="secondary" style={{ fontSize: 13 }}>
                              {type.description}
                            </Text>
                          </div>
                        </div>
                        <div style={{
                          padding: 12,
                          background: token.colorBgContainer,
                          borderRadius: 4,
                          borderLeft: `3px solid ${token.colorPrimary}`
                        }}>
                          <Text style={{ fontSize: 12, color: token.colorTextSecondary }}>
                            {type.example}
                          </Text>
                        </div>
                      </Space>
                    </Card>
                  </Col>
                ))}
              </Row>

              <div style={{ textAlign: 'center', marginTop: 16 }}>
                <Button
                  type="primary"
                  size="large"
                  onClick={() => setComparisonStep(2)}
                  style={{ minWidth: 120 }}
                >
                  下一步
                </Button>
              </div>
            </Space>
          )}

          {/* 步骤 2: 选择参数 */}
          {comparisonStep === 2 && renderParameterSelection()}

          {/* 步骤 3/4: 显示对比结果 */}
          {(comparisonStep === 3 || comparisonStep === 4) && comparisonResults.length > 0 && (
            <Card
              size="small"
              title={
                <Space>
                  <BarChartOutlined style={{ color: '#52c41a' }} />
                  <span>对比结果 ({comparisonResults.length} 个任务)</span>
                </Space>
              }
              extra={
                <Space>
                  <Button size="small" onClick={() => setComparisonStep(2)}>
                    重新选择
                  </Button>
                  <Button
                    size="small"
                    type="primary"
                    icon={<DownloadOutlined />}
                    onClick={handleExportComparison}
                  >
                    导出数据
                  </Button>
                </Space>
              }
            >
              {renderComparisonResults()}
            </Card>
          )}
        </Space>
      </Modal>

      {/* 数据表格 */}
      <Card size="small">
        <Table
          columns={columns}
          dataSource={data}
          rowKey="id"
          loading={loading}
          scroll={{ x: 1200 }}
          pagination={{
            current: currentPage,
            pageSize: pageSize,
            total: total,
            showSizeChanger: true,
            showQuickJumper: true,
            showTotal: (total) => `共 ${total} 条`,
            onChange: (page, size) => {
              const values = form.getFieldsValue();
              loadData(page, size, values);
            },
          }}
        />
      </Card>



      {/* 重新计算对话框 */}
      <QCRecalculateModal
        visible={recalculateModalVisible}
        job={selectedJob}
        onClose={() => {
          setRecalculateModalVisible(false);
          setSelectedJob(null);
        }}
        onSuccess={handleRecalculateSuccess}
      />

      {/* QC详情模态框（公开数据） */}
      <Modal
        title={
          <Space>
            <ExperimentOutlined style={{ color: '#1890ff' }} />
            <span>QC 计算结果详情</span>
          </Space>
        }
        open={detailModalVisible}
        onCancel={() => {
          setDetailModalVisible(false);
          setDetailJob(null);
          setDetailImages({});
        }}
        footer={null}
        width={900}
        destroyOnClose
      >
        {detailJob && (
          <Space direction="vertical" style={{ width: '100%' }} size="middle">
            {/* 基本信息 */}
            <Descriptions bordered size="small" column={2}>
              <Descriptions.Item label="分子名称">{detailJob.molecule_name || '-'}</Descriptions.Item>
              <Descriptions.Item label="SMILES">
                <Text copyable style={{ maxWidth: 300 }}>{detailJob.smiles || '-'}</Text>
              </Descriptions.Item>
              <Descriptions.Item label="泛函">{detailJob.functional || '-'}</Descriptions.Item>
              <Descriptions.Item label="基组">{detailJob.basis_set || '-'}</Descriptions.Item>
              <Descriptions.Item label="溶剂模型">{renderSolventModel(detailJob)}</Descriptions.Item>
              <Descriptions.Item label="任务类型">OPT+FREQ</Descriptions.Item>
            </Descriptions>

            {/* 计算结果 */}
            {detailJob.results && detailJob.results.length > 0 && (
              <Card title="计算结果" size="small">
                <Row gutter={16}>
                  <Col span={8}>
                    <Statistic
                      title="总能量"
                      value={detailJob.results[0].energy_au?.toFixed(6) || '-'}
                      suffix="Hartree"
                      valueStyle={{ fontSize: 16 }}
                    />
                  </Col>
                  <Col span={8}>
                    <Statistic
                      title="HOMO 能量"
                      value={detailJob.results[0].homo?.toFixed(4) || '-'}
                      suffix="eV"
                      valueStyle={{ fontSize: 16, color: '#1890ff' }}
                    />
                  </Col>
                  <Col span={8}>
                    <Statistic
                      title="LUMO 能量"
                      value={detailJob.results[0].lumo?.toFixed(4) || '-'}
                      suffix="eV"
                      valueStyle={{ fontSize: 16, color: '#52c41a' }}
                    />
                  </Col>
                </Row>
                {(detailJob.results[0].homo && detailJob.results[0].lumo) && (
                  <Row style={{ marginTop: 16 }}>
                    <Col span={8}>
                      <Statistic
                        title="HOMO-LUMO Gap"
                        value={(detailJob.results[0].lumo - detailJob.results[0].homo).toFixed(4)}
                        suffix="eV"
                        valueStyle={{ fontSize: 16, color: '#722ed1' }}
                      />
                    </Col>
                  </Row>
                )}
              </Card>
            )}

            {/* 轨道图片 */}
            <Card title="分子轨道可视化" size="small">
              {detailImagesLoading ? (
                <div style={{ textAlign: 'center', padding: 40 }}>
                  <Spin tip="加载图片中..." />
                </div>
              ) : (
                <Row gutter={16}>
                  <Col span={8}>
                    <Card size="small" title="ESP 静电势">
                      {detailImages.esp ? (
                        <Image src={detailImages.esp} alt="ESP" style={{ width: '100%' }} />
                      ) : (
                        <div style={{ textAlign: 'center', padding: 20, color: '#999' }}>暂无图片</div>
                      )}
                    </Card>
                  </Col>
                  <Col span={8}>
                    <Card size="small" title="HOMO 轨道">
                      {detailImages.homo ? (
                        <Image src={detailImages.homo} alt="HOMO" style={{ width: '100%' }} />
                      ) : (
                        <div style={{ textAlign: 'center', padding: 20, color: '#999' }}>暂无图片</div>
                      )}
                    </Card>
                  </Col>
                  <Col span={8}>
                    <Card size="small" title="LUMO 轨道">
                      {detailImages.lumo ? (
                        <Image src={detailImages.lumo} alt="LUMO" style={{ width: '100%' }} />
                      ) : (
                        <div style={{ textAlign: 'center', padding: 20, color: '#999' }}>暂无图片</div>
                      )}
                    </Card>
                  </Col>
                </Row>
              )}
            </Card>
          </Space>
        )}
      </Modal>
    </Space>
  );
}

