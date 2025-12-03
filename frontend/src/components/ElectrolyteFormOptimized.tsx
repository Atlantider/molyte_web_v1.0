/**
 * 优化版电解质配方表单组件
 * 采用双栏布局 + 分类折叠 + 精简UI + 用户自定义常用组合
 */
import { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Form,
  Select,
  Input,
  InputNumber,
  Button,
  Space,
  Row,
  Col,
  Card,
  message,
  Spin,
  Alert,
  Collapse,
  Tag,
  Typography,
  Statistic,
  Modal,
  Divider,
  Popconfirm,
} from 'antd';
import {
  PlusOutlined,
  DeleteOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ExperimentOutlined,
  FireOutlined,
  StarOutlined,
  StarFilled,
  EditOutlined,
  FolderAddOutlined,
} from '@ant-design/icons';
import type { IonInfo, Project } from '../types';
import { getAvailableIons } from '../api/electrolytes';
import {
  getUserSolventCombinations,
  createSolventCombination,
  deleteSolventCombination,
  type CustomSolventCombinationResponse,
} from '../api/userPreferences';

const { Text, Title } = Typography;
const { Panel } = Collapse;

// 溶剂分类（精简版，按使用频率分类）
const SOLVENT_CATEGORIES = [
  {
    category: '碳酸酯类',
    icon: '🔋',
    color: '#1890ff',
    description: '锂离子电池最常用',
    solvents: [
      { name: 'EC', label: 'EC', fullName: '碳酸乙烯酯', smiles: 'C1COC(=O)O1' },
      { name: 'PC', label: 'PC', fullName: '碳酸丙烯酯', smiles: 'CC1COC(=O)O1' },
      { name: 'DMC', label: 'DMC', fullName: '碳酸二甲酯', smiles: 'COC(=O)OC' },
      { name: 'DEC', label: 'DEC', fullName: '碳酸二乙酯', smiles: 'CCOC(=O)OCC' },
      { name: 'EMC', label: 'EMC', fullName: '碳酸甲乙酯', smiles: 'CCOC(=O)OC' },
      { name: 'VC', label: 'VC', fullName: '碳酸亚乙烯酯', smiles: 'C1=COC(=O)O1' },
    ],
  },
  {
    category: '氟代碳酸酯',
    icon: '⚡',
    color: '#722ed1',
    description: '高电压/锂金属电池',
    solvents: [
      { name: 'FEC', label: 'FEC', fullName: '氟代碳酸乙烯酯', smiles: 'C1C(OC(=O)O1)F' },
      { name: 'DFEC', label: 'DFEC', fullName: '二氟碳酸乙烯酯', smiles: 'FC1OC(=O)OC1F' },
      { name: 'TFPC', label: 'TFPC', fullName: '三氟碳酸丙烯酯', smiles: 'CC(F)(F)C1COC(=O)O1' },
    ],
  },
  {
    category: '醚类',
    icon: '💧',
    color: '#13c2c2',
    description: '锂硫/锂空电池',
    solvents: [
      { name: 'DME', label: 'DME', fullName: '乙二醇二甲醚', smiles: 'COCCOC' },
      { name: 'DOL', label: 'DOL', fullName: '1,3-二氧戊环', smiles: 'C1COCO1' },
      { name: 'DIOX', label: 'DIOX', fullName: '1,4-二氧六环', smiles: 'C1COCCO1' },
      { name: 'THF', label: 'THF', fullName: '四氢呋喃', smiles: 'C1CCOC1' },
      { name: 'DEGDME', label: 'DEGDME', fullName: '二乙二醇二甲醚', smiles: 'COCCOCCOC' },
      { name: 'TEGDME', label: 'TEGDME', fullName: '四乙二醇二甲醚', smiles: 'COCCOCCOCCOCCOC' },
    ],
  },
  {
    category: '氟代醚类',
    icon: '🧪',
    color: '#fa8c16',
    description: '局部高浓度稀释剂',
    solvents: [
      { name: 'TTE', label: 'TTE', fullName: '四氟乙基-四氟丙基醚', smiles: 'FC(F)C(F)(F)OCC(F)(F)C(F)(F)F' },
      { name: 'BTFE', label: 'BTFE', fullName: '双三氟乙基醚', smiles: 'FC(F)(F)COCC(F)(F)F' },
      { name: 'TFEE', label: 'TFEE', fullName: '四氟乙基-三氟乙基醚', smiles: 'FC(F)C(F)(F)OCC(F)(F)F' },
      { name: 'HFE', label: 'HFE', fullName: '八氟戊基-四氟乙基醚', smiles: 'FC(F)C(F)(F)OCC(F)(F)C(F)(F)C(F)(F)CF' },
    ],
  },
  {
    category: '其他常用',
    icon: '🔬',
    color: '#52c41a',
    description: '腈类、砜类等',
    solvents: [
      { name: 'ACN', label: 'ACN', fullName: '乙腈', smiles: 'CC#N' },
      { name: 'SN', label: 'SN', fullName: '琥珀腈', smiles: 'N#CCCC#N' },
      { name: 'Water', label: 'Water', fullName: '水', smiles: 'O' },
      { name: 'Sulfolane', label: 'Sulfolane', fullName: '环丁砜', smiles: 'C1CCS(=O)(=O)C1' },
      { name: 'GBL', label: 'GBL', fullName: 'γ-丁内酯', smiles: 'C1CC(=O)OC1' },
      { name: 'NMP', label: 'NMP', fullName: 'N-甲基吡咯烷酮', smiles: 'CN1CCCC1=O' },
    ],
  },
];

// 常用组合（精简到4个最常用的）
const COMMON_COMBINATIONS = [
  { name: 'EC/DMC (1:1)', solvents: [{ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 }, { name: 'DMC', smiles: 'COC(=O)OC', molar_ratio: 1.0 }] },
  { name: 'EC/DEC (1:1)', solvents: [{ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 }, { name: 'DEC', smiles: 'CCOC(=O)OCC', molar_ratio: 1.0 }] },
  { name: 'EC/EMC (1:1)', solvents: [{ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 }, { name: 'EMC', smiles: 'CCOC(=O)OC', molar_ratio: 1.0 }] },
  { name: 'DME/DOL (1:1)', solvents: [{ name: 'DME', smiles: 'COCCOC', molar_ratio: 1.0 }, { name: 'DOL', smiles: 'C1COCO1', molar_ratio: 1.0 }] },
];

interface ElectrolyteFormOptimizedProps {
  form: any;
  projects: Project[];
  onValuesChange?: (changedValues: any, allValues: any) => void;
  onIonsChange?: (cations: SelectedIon[], anions: SelectedIon[]) => void;
  initialCations?: SelectedIon[];
  initialAnions?: SelectedIon[];
  onCreateProject?: () => void; // 新增：创建项目的回调
}

interface SelectedIon {
  name: string;
  charge: number;
  concentration: number;
}

export default function ElectrolyteFormOptimized({
  form,
  projects,
  onValuesChange,
  onIonsChange,
  initialCations = [],
  initialAnions = [],
  onCreateProject
}: ElectrolyteFormOptimizedProps) {
  const navigate = useNavigate();
  const [availableCations, setAvailableCations] = useState<IonInfo[]>([]);
  const [availableAnions, setAvailableAnions] = useState<IonInfo[]>([]);
  const [selectedCations, setSelectedCations] = useState<SelectedIon[]>(initialCations);
  const [selectedAnions, setSelectedAnions] = useState<SelectedIon[]>(initialAnions);
  const [loading, setLoading] = useState(false);
  const [isElectricallyNeutral, setIsElectricallyNeutral] = useState<boolean | null>(null);
  const [solvents, setSolvents] = useState<any[]>([]);
  const [userModifiedName, setUserModifiedName] = useState(false);
  const [lastGeneratedName, setLastGeneratedName] = useState('');

  // 用户自定义常用溶剂组合
  const [userSolventCombinations, setUserSolventCombinations] = useState<CustomSolventCombinationResponse[]>([]);
  const [saveModalVisible, setSaveModalVisible] = useState(false);
  const [saveName, setSaveName] = useState('');
  const [saveDescription, setSaveDescription] = useState('');

  useEffect(() => {
    loadAvailableIons();
    loadUserPreferences();
    const initialSolvents = form.getFieldValue('solvents') || [];
    setSolvents(initialSolvents);
  }, []);

  useEffect(() => {
    if (initialCations.length > 0) setSelectedCations(initialCations);
  }, [initialCations]);

  useEffect(() => {
    if (initialAnions.length > 0) setSelectedAnions(initialAnions);
  }, [initialAnions]);

  const loadAvailableIons = async () => {
    setLoading(true);
    try {
      const data = await getAvailableIons();
      setAvailableCations(data.cations);
      setAvailableAnions(data.anions);
    } catch (error: any) {
      message.error('加载可用离子列表失败');
    } finally {
      setLoading(false);
    }
  };

  // 加载用户自定义偏好
  const loadUserPreferences = async () => {
    try {
      const combinations = await getUserSolventCombinations();
      setUserSolventCombinations(combinations);
    } catch (error: any) {
      console.error('加载用户偏好失败:', error);
      // 不显示错误消息，因为这不是关键功能
    }
  };

  // 保存当前溶剂配置为自定义组合
  const handleSaveAsCombination = async () => {
    if (!saveName.trim()) {
      message.error('请输入组合名称');
      return;
    }

    const currentSolvents = form.getFieldValue('solvents') || [];
    if (currentSolvents.length === 0) {
      message.error('请先添加溶剂');
      return;
    }

    try {
      await createSolventCombination({
        name: saveName.trim(),
        description: saveDescription.trim() || undefined,
        solvents: currentSolvents.map((s: any) => ({
          name: s.name,
          smiles: s.smiles,
          molar_ratio: s.molar_ratio || 1.0
        }))
      });
      message.success('保存成功！');
      setSaveModalVisible(false);
      setSaveName('');
      setSaveDescription('');
      loadUserPreferences(); // 重新加载
    } catch (error: any) {
      message.error(error.response?.data?.detail || '保存失败');
    }
  };

  // 删除自定义组合
  const handleDeleteCombination = async (id: number) => {
    try {
      await deleteSolventCombination(id);
      message.success('删除成功');
      loadUserPreferences(); // 重新加载
    } catch (error: any) {
      message.error('删除失败');
    }
  };

  // 应用自定义组合
  const handleApplyUserCombination = (combination: CustomSolventCombinationResponse) => {
    const currentSolvents = form.getFieldValue('solvents') || [];
    const newSolvents = [...currentSolvents, ...combination.solvents];
    form.setFieldsValue({ solvents: newSolvents });
    setSolvents(newSolvents);
    message.success(`已添加 "${combination.name}"`);
  };

  const checkElectricalNeutrality = () => {
    if (selectedCations.length === 0 || selectedAnions.length === 0) {
      setIsElectricallyNeutral(null);
      return;
    }
    const totalPositiveCharge = selectedCations.reduce((sum, ion) => sum + ion.concentration * ion.charge, 0);
    const totalNegativeCharge = selectedAnions.reduce((sum, ion) => sum + ion.concentration * Math.abs(ion.charge), 0);
    const isNeutral = Math.abs(totalPositiveCharge - totalNegativeCharge) < 0.01;
    setIsElectricallyNeutral(isNeutral);
  };

  useEffect(() => {
    checkElectricalNeutrality();
    if (onIonsChange) {
      onIonsChange(selectedCations, selectedAnions);
    }
  }, [selectedCations, selectedAnions]);

  const addCation = (ionName: string) => {
    const ion = availableCations.find(i => i.name === ionName);
    if (ion && !selectedCations.find(c => c.name === ionName)) {
      setSelectedCations([...selectedCations, { ...ion, concentration: 1.0 }]);
    }
  };

  const addAnion = (ionName: string) => {
    const ion = availableAnions.find(i => i.name === ionName);
    if (ion && !selectedAnions.find(a => a.name === ionName)) {
      setSelectedAnions([...selectedAnions, { ...ion, concentration: 1.0 }]);
    }
  };

  const removeCation = (ionName: string) => {
    setSelectedCations(selectedCations.filter(c => c.name !== ionName));
  };

  const removeAnion = (ionName: string) => {
    setSelectedAnions(selectedAnions.filter(a => a.name !== ionName));
  };

  const updateCationConcentration = (name: string, concentration: number) => {
    setSelectedCations(selectedCations.map(c => c && c.name === name ? { ...c, concentration } : c));
  };

  const updateAnionConcentration = (name: string, concentration: number) => {
    setSelectedAnions(selectedAnions.map(a => a && a.name === name ? { ...a, concentration } : a));
  };

  const generateDefaultName = () => {
    if (selectedCations.length === 0 && selectedAnions.length === 0 && solvents.length === 0) return '';
    const parts: string[] = [];
    const validCations = selectedCations.filter(c => c && c.name);
    const sortedCations = [...validCations].sort((a, b) => (b.concentration || 0) - (a.concentration || 0));
    // 清理名称中的特殊字符，将/替换为-
    const topCations = sortedCations.slice(0, 2).map(c => c.name?.replace(/\//g, '-')).filter(Boolean);
    if (topCations.length > 0) parts.push(topCations.join('-'));
    const validAnions = selectedAnions.filter(a => a && a.name);
    const sortedAnions = [...validAnions].sort((a, b) => (b.concentration || 0) - (a.concentration || 0));
    const topAnions = sortedAnions.slice(0, 2).map(a => a.name?.replace(/\//g, '-')).filter(Boolean);
    if (topAnions.length > 0) parts.push(topAnions.join('-'));
    const validSolvents = solvents.filter((s: any) => s && s.name);
    const sortedSolvents = [...validSolvents].sort((a: any, b: any) => (b.molar_ratio || 0) - (a.molar_ratio || 0));
    // 清理溶剂名称中的特殊字符，将/替换为-
    const topSolvents = sortedSolvents.slice(0, 3).map((s: any) => s.name?.replace(/\//g, '-')).filter(Boolean);
    if (topSolvents.length > 0) parts.push(topSolvents.join('-'));
    return parts.join('-');
  };

  useEffect(() => {
    if (!form) return;
    const currentName = form.getFieldValue('name');
    if (userModifiedName && currentName !== lastGeneratedName) return;
    try {
      const newName = generateDefaultName();
      if (newName && newName !== currentName) {
        form.setFieldsValue({ name: newName });
        setLastGeneratedName(newName);
        setUserModifiedName(false);
      }
    } catch (error) {
      console.error('生成默认名称时出错:', error);
    }
  }, [selectedCations, selectedAnions, solvents]);

  const handleFormChange = (changedValues: any, allValues: any) => {
    if (changedValues.solvents) {
      setSolvents(allValues.solvents || []);
    }
    if (changedValues.name !== undefined) {
      setUserModifiedName(true);
    }
    if (onValuesChange) {
      onValuesChange(changedValues, allValues);
    }
  };

  if (loading) {
    return <Spin tip="加载离子列表..." />;
  }

  return (
    <>
    <Row gutter={24}>
      {/* 左侧：表单配置 */}
      <Col span={16}>
        <Form form={form} layout="vertical" onValuesChange={handleFormChange}>
          {/* 基本信息 - 三栏布局 */}
          <Card size="small" title="📋 基本信息" style={{ marginBottom: 16 }}>
            <Row gutter={16}>
              <Col span={10}>
                <Form.Item name="project_id" label="所属项目" rules={[{ required: true, message: '请选择所属项目' }]}>
                  <Select
                    placeholder="选择项目"
                    notFoundContent={
                      <div style={{ textAlign: 'center', padding: '16px 0' }}>
                        <div style={{ color: '#999', marginBottom: 8 }}>暂无项目</div>
                        <Button
                          type="primary"
                          icon={<FolderAddOutlined />}
                          onClick={() => {
                            if (onCreateProject) {
                              onCreateProject();
                            } else {
                              navigate('/workspace/projects?action=create');
                            }
                          }}
                        >
                          新建项目
                        </Button>
                      </div>
                    }
                    dropdownRender={(menu) => {
                      // 只在有项目时显示下拉菜单底部的"新建项目"按钮
                      const hasProjects = projects && projects.length > 0;
                      return (
                        <>
                          {menu}
                          {hasProjects && (
                            <>
                              <Divider style={{ margin: '8px 0' }} />
                              <div style={{ padding: '4px 8px' }}>
                                <Button
                                  type="link"
                                  icon={<FolderAddOutlined />}
                                  onClick={() => {
                                    if (onCreateProject) {
                                      onCreateProject();
                                    } else {
                                      navigate('/workspace/projects?action=create');
                                    }
                                  }}
                                  style={{ width: '100%', textAlign: 'left' }}
                                >
                                  新建项目
                                </Button>
                              </div>
                            </>
                          )}
                        </>
                      );
                    }}
                  >
                    {projects?.filter(p => p && p.id && p.name).map(p => <Select.Option key={p.id} value={p.id}>{p.name}</Select.Option>)}
                  </Select>
                </Form.Item>
              </Col>
              <Col span={8}>
                <Form.Item
                  name="name"
                  label="配方名称"
                  rules={[{ max: 100, message: '名称不能超过100个字符' }]}
                  tooltip="留空将自动生成"
                >
                  <Input placeholder="留空自动生成" allowClear />
                </Form.Item>
              </Col>
              <Col span={6}>
                <Form.Item name="temperature" label="温度 (K)" rules={[{ required: true }]} initialValue={298.15}>
                  <InputNumber min={0} max={1000} step={0.01} style={{ width: '100%' }} />
                </Form.Item>
              </Col>
            </Row>
          </Card>

          {/* 离子配置 - 双栏布局 */}
          <Card size="small" title={<><ExperimentOutlined /> 离子配置</>} style={{ marginBottom: 16 }}>
            <Row gutter={16}>
              {/* 阳离子 */}
              <Col span={12}>
                <div>
                  <Text strong style={{ color: '#1890ff' }}>阳离子</Text>
                  <Select
                    placeholder="选择阳离子添加"
                    style={{ width: '100%', marginTop: 8 }}
                    onChange={addCation}
                    value={undefined}
                    showSearch
                    size="small"
                    filterOption={(input, option) =>
                      (option?.label as string)?.toLowerCase().includes(input.toLowerCase())
                    }
                  >
                    {availableCations
                      .filter(ion => !selectedCations.find(c => c.name === ion.name))
                      .map(ion => (
                        <Select.Option key={ion.name} value={ion.name}>
                          {ion.name} (+{ion.charge})
                        </Select.Option>
                      ))}
                  </Select>
                  <div style={{ marginTop: 12 }}>
                    {selectedCations.filter(ion => ion && ion.name).map((ion, index) => (
                      <Card key={ion.name} size="small" style={{ marginBottom: 8, borderColor: index === 0 ? '#1890ff' : undefined, borderWidth: index === 0 ? 2 : 1 }}>
                        <div style={{ marginBottom: 8 }}>
                          <Tag color="blue">{ion.name}</Tag>
                          <Text type="secondary" style={{ fontSize: 12 }}>+{ion.charge}</Text>
                          {index === 0 && <Tag color="gold" style={{ marginLeft: 4 }}>第一种</Tag>}
                        </div>
                        <div style={{ marginBottom: 4 }}>
                          <Text strong style={{ color: '#1890ff', fontSize: 13 }}>浓度 (mol/L):</Text>
                        </div>
                        <Space.Compact style={{ width: '100%' }}>
                          <InputNumber
                            size="middle"
                            min={0.001}
                            max={10}
                            step={0.1}
                            value={ion.concentration}
                            onChange={(value) => updateCationConcentration(ion.name, value || 0)}
                            style={{
                              width: 'calc(100% - 32px)',
                              fontWeight: 'bold',
                              fontSize: 14
                            }}
                            addonAfter={<span style={{ fontWeight: 'bold' }}>M</span>}
                          />
                          <Button
                            type="text"
                            danger
                            size="small"
                            icon={<DeleteOutlined />}
                            onClick={() => removeCation(ion.name)}
                          />
                        </Space.Compact>
                      </Card>
                    ))}
                    {selectedCations.length === 0 && (
                      <Alert message="请选择阳离子" type="warning" showIcon style={{ padding: '4px 8px' }} />
                    )}
                  </div>
                </div>
              </Col>

              {/* 阴离子 */}
              <Col span={12}>
                <div>
                  <Text strong style={{ color: '#f5222d' }}>阴离子</Text>
                  <Select
                    placeholder="选择阴离子添加"
                    style={{ width: '100%', marginTop: 8 }}
                    onChange={addAnion}
                    value={undefined}
                    showSearch
                    size="small"
                  >
                    {availableAnions
                      .filter(ion => !selectedAnions.find(a => a.name === ion.name))
                      .map(ion => (
                        <Select.Option key={ion.name} value={ion.name}>
                          {ion.name} ({ion.charge})
                        </Select.Option>
                      ))}
                  </Select>
                  <div style={{ marginTop: 12 }}>
                    {selectedAnions.filter(ion => ion && ion.name).map(ion => (
                      <Card key={ion.name} size="small" style={{ marginBottom: 8 }}>
                        <div style={{ marginBottom: 8 }}>
                          <Tag color="red">{ion.name}</Tag>
                          <Text type="secondary" style={{ fontSize: 12 }}>{ion.charge}</Text>
                        </div>
                        <div style={{ marginBottom: 4 }}>
                          <Text strong style={{ color: '#f5222d', fontSize: 13 }}>浓度 (mol/L):</Text>
                        </div>
                        <Space.Compact style={{ width: '100%' }}>
                          <InputNumber
                            size="middle"
                            min={0.001}
                            max={10}
                            step={0.1}
                            value={ion.concentration}
                            onChange={(value) => updateAnionConcentration(ion.name, value || 0)}
                            style={{
                              width: 'calc(100% - 32px)',
                              fontWeight: 'bold',
                              fontSize: 14
                            }}
                            addonAfter={<span style={{ fontWeight: 'bold' }}>M</span>}
                          />
                          <Button
                            type="text"
                            danger
                            size="small"
                            icon={<DeleteOutlined />}
                            onClick={() => removeAnion(ion.name)}
                          />
                        </Space.Compact>
                      </Card>
                    ))}
                    {selectedAnions.length === 0 && (
                      <Alert message="请选择阴离子" type="warning" showIcon style={{ padding: '4px 8px' }} />
                    )}
                  </div>
                </div>
              </Col>
            </Row>
          </Card>

          {/* 溶剂配置 - 使用折叠面板 */}
          <Card size="small" title={<><FireOutlined /> 溶剂配置</>} style={{ marginBottom: 16 }}>
            {/* 醒目提示：摩尔比说明 */}
            <Alert
              message={
                <span>
                  <Text strong style={{ color: '#fa8c16' }}>⚠️ 重要提示：</Text>
                  <Text style={{ marginLeft: 8 }}>
                    溶剂的<Text strong style={{ color: '#fa8c16' }}>摩尔比</Text>是相对于
                    <Text strong style={{ color: '#1890ff' }}>第一种阳离子</Text>的比例
                  </Text>
                </span>
              }
              description={
                <div style={{ fontSize: 12 }}>
                  例如：第一种阳离子浓度为 1.0 M，溶剂摩尔比为 5.0，则该溶剂的实际浓度为 5.0 M
                </div>
              }
              type="warning"
              showIcon
              style={{ marginBottom: 16 }}
            />

            <Form.List name="solvents" initialValue={[]}>
              {(fields, { add, remove }) => (
                <>
                  {/* 用户自定义组合 */}
                  {userSolventCombinations.length > 0 && (
                    <div style={{ marginBottom: 12 }}>
                      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
                        <Text type="secondary" style={{ fontSize: 12 }}>
                          <StarFilled style={{ color: '#faad14', marginRight: 4 }} />
                          我的常用组合：
                        </Text>
                      </div>
                      <div style={{ marginTop: 8 }}>
                        <Space wrap size="small">
                          {userSolventCombinations.map((combo) => (
                            <Popconfirm
                              key={combo.id}
                              title={
                                <div>
                                  <div><strong>{combo.name}</strong></div>
                                  {combo.description && <div style={{ fontSize: 12, color: '#666' }}>{combo.description}</div>}
                                  <div style={{ marginTop: 8, fontSize: 12 }}>
                                    溶剂：{combo.solvents?.filter(s => s && s.name).map(s => `${s.name}(${s.molar_ratio})`).join(' + ')}
                                  </div>
                                </div>
                              }
                              description={
                                <Space>
                                  <Button
                                    size="small"
                                    type="primary"
                                    onClick={() => {
                                      handleApplyUserCombination(combo);
                                    }}
                                  >
                                    应用
                                  </Button>
                                  <Button
                                    size="small"
                                    danger
                                    onClick={() => handleDeleteCombination(combo.id)}
                                  >
                                    删除
                                  </Button>
                                </Space>
                              }
                              icon={null}
                            >
                              <Button
                                size="small"
                                type="primary"
                                ghost
                                icon={<StarFilled />}
                                title={combo.description}
                              >
                                {combo.name}
                              </Button>
                            </Popconfirm>
                          ))}
                        </Space>
                      </div>
                    </div>
                  )}

                  {/* 常用组合快速添加 */}
                  <div style={{ marginBottom: 12 }}>
                    <Text type="secondary" style={{ fontSize: 12 }}>常用组合：</Text>
                    <div style={{ marginTop: 8 }}>
                      <Space wrap size="small">
                        {COMMON_COMBINATIONS.map((combo, idx) => (
                          <Button
                            key={idx}
                            size="small"
                            onClick={() => {
                              combo.solvents.forEach(s => add(s));
                            }}
                          >
                            {combo.name}
                          </Button>
                        ))}
                      </Space>
                    </div>
                  </div>

                  {/* 分类选择器 - 使用折叠面板 */}
                  <Collapse
                    size="small"
                    style={{ marginBottom: 12 }}
                    items={SOLVENT_CATEGORIES.map((cat, catIdx) => ({
                      key: catIdx,
                      label: (
                        <span>
                          <span style={{ marginRight: 8 }}>{cat.icon}</span>
                          <Text strong style={{ color: cat.color }}>{cat.category}</Text>
                          <Text type="secondary" style={{ fontSize: 12, marginLeft: 8 }}>
                            {cat.description}
                          </Text>
                        </span>
                      ),
                      children: (
                        <Space wrap size={[4, 4]}>
                          {cat.solvents.map((solvent, idx) => (
                            <Button
                              key={idx}
                              size="small"
                              type="default"
                              onClick={() => add({ name: solvent.name, smiles: solvent.smiles, molar_ratio: 1.0 })}
                              title={solvent.fullName}
                            >
                              {solvent.label}
                            </Button>
                          ))}
                        </Space>
                      ),
                    }))}
                  />

                  {/* 已添加的溶剂列表 - 双栏布局 */}
                  {fields.length > 0 && (
                    <div style={{ marginBottom: 12 }}>
                      <Text type="secondary" style={{ fontSize: 12 }}>已添加的溶剂：</Text>
                      <Row gutter={[8, 8]} style={{ marginTop: 8 }}>
                        {fields.map((field, index) => (
                          <Col span={12} key={field.key}>
                            <Card
                              size="small"
                              title={
                                <Space size="small">
                                  <Text strong style={{ fontSize: 12 }}>溶剂 {index + 1}</Text>
                                </Space>
                              }
                              extra={
                                <Button
                                  type="text"
                                  danger
                                  size="small"
                                  icon={<DeleteOutlined />}
                                  onClick={() => remove(index)}
                                />
                              }
                              style={{ height: '100%' }}
                            >
                              <Form.Item
                                {...field}
                                name={[field.name, 'name']}
                                label="名称"
                                rules={[{ required: true, message: '请输入名称' }]}
                                style={{ marginBottom: 8 }}
                              >
                                <Input placeholder="例如: EC" size="small" />
                              </Form.Item>
                              <Form.Item
                                {...field}
                                name={[field.name, 'smiles']}
                                label="SMILES"
                                rules={[{ required: true, message: '请输入 SMILES' }]}
                                style={{ marginBottom: 8 }}
                              >
                                <Input placeholder="例如: C1COC(=O)O1" size="small" />
                              </Form.Item>
                              <Form.Item
                                {...field}
                                name={[field.name, 'molar_ratio']}
                                label={
                                  <span>
                                    <Text strong style={{ color: '#fa8c16', fontSize: 13 }}>摩尔比</Text>
                                    <Text type="secondary" style={{ fontSize: 11, marginLeft: 4 }}>
                                      (相对第一种阳离子)
                                    </Text>
                                  </span>
                                }
                                rules={[{ required: true, message: '请输入摩尔比' }]}
                                initialValue={1.0}
                                style={{ marginBottom: 0 }}
                                tooltip={{
                                  title: '溶剂的摩尔数量相对于第一种阳离子的比例。例如：第一种阳离子为 1M，摩尔比为 5，则该溶剂为 5M',
                                  color: '#fa8c16'
                                }}
                              >
                                <InputNumber
                                  min={0.1}
                                  max={100}
                                  step={0.1}
                                  style={{ width: '100%', fontWeight: 'bold' }}
                                  size="middle"
                                />
                              </Form.Item>
                            </Card>
                          </Col>
                        ))}
                      </Row>
                    </div>
                  )}

                  {/* 自定义添加按钮 */}
                  <Space direction="vertical" style={{ width: '100%' }}>
                    <Button type="dashed" onClick={() => add()} block icon={<PlusOutlined />}>
                      自定义添加溶剂
                    </Button>

                    {/* 保存为常用组合按钮 */}
                    {fields.length > 0 && (
                      <Button
                        type="primary"
                        ghost
                        block
                        icon={<StarOutlined />}
                        onClick={() => setSaveModalVisible(true)}
                      >
                        保存为我的常用组合
                      </Button>
                    )}
                  </Space>
                </>
              )}
            </Form.List>
          </Card>

          {/* 盒子尺寸配置 - 双栏布局 */}
          <Card size="small" title="📦 模拟盒子" style={{ marginBottom: 16 }}>
            <Row gutter={16}>
              <Col span={12}>
                <Form.Item name="box_type" label="盒子类型" initialValue="cubic">
                  <Select>
                    <Select.Option value="cubic">立方体</Select.Option>
                    <Select.Option value="rectangular">长方体</Select.Option>
                  </Select>
                </Form.Item>
              </Col>
              <Col span={12}>
                <Form.Item noStyle shouldUpdate={(prev, curr) => prev.box_type !== curr.box_type}>
                  {({ getFieldValue }) => {
                    const boxType = getFieldValue('box_type');
                    if (boxType === 'cubic' || boxType === undefined) {
                      return (
                        <Form.Item
                          name="box_size"
                          label="边长 (Å)"
                          rules={[{ required: true, message: '请输入边长' }]}
                          initialValue={40}
                          extra="建议: 30-50 Å"
                        >
                          <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                        </Form.Item>
                      );
                    } else {
                      return null;
                    }
                  }}
                </Form.Item>
              </Col>
            </Row>

            <Form.Item noStyle shouldUpdate={(prev, curr) => prev.box_type !== curr.box_type}>
              {({ getFieldValue }) => {
                const boxType = getFieldValue('box_type');
                if (boxType === 'rectangular') {
                  return (
                    <Row gutter={8}>
                      <Col span={8}>
                        <Form.Item
                          name={['box_dimensions', 0]}
                          label="长 (Å)"
                          rules={[{ required: true }]}
                          initialValue={40}
                        >
                          <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                        </Form.Item>
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          name={['box_dimensions', 1]}
                          label="宽 (Å)"
                          rules={[{ required: true }]}
                          initialValue={40}
                        >
                          <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                        </Form.Item>
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          name={['box_dimensions', 2]}
                          label="高 (Å)"
                          rules={[{ required: true }]}
                          initialValue={40}
                        >
                          <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                        </Form.Item>
                      </Col>
                    </Row>
                  );
                }
                return null;
              }}
            </Form.Item>
          </Card>

          <Form.Item name="description" label="备注">
            <Input.TextArea rows={2} placeholder="可选的配方说明" />
          </Form.Item>
        </Form>
      </Col>

      {/* 右侧：配置预览 */}
      <Col span={8}>
        <Card
          title="配置预览"
          size="small"
          style={{
            position: 'sticky',
            top: 24,
            background: '#fafafa',
          }}
        >
          <Space direction="vertical" style={{ width: '100%' }} size="middle">
            {/* 电中性检查 */}
            {isElectricallyNeutral !== null && (
              <Alert
                message={isElectricallyNeutral ? '电荷平衡 ✓' : '电荷不平衡 ✗'}
                description={
                  isElectricallyNeutral
                    ? '体系满足电中性条件'
                    : '请调整离子浓度以满足电中性'
                }
                type={isElectricallyNeutral ? 'success' : 'error'}
                showIcon
                icon={isElectricallyNeutral ? <CheckCircleOutlined /> : <CloseCircleOutlined />}
              />
            )}

            {/* 统计信息 */}
            <Row gutter={8}>
              <Col span={12}>
                <Statistic
                  title="阳离子"
                  value={selectedCations.length}
                  suffix="种"
                  valueStyle={{ fontSize: 20, color: '#1890ff' }}
                />
              </Col>
              <Col span={12}>
                <Statistic
                  title="阴离子"
                  value={selectedAnions.length}
                  suffix="种"
                  valueStyle={{ fontSize: 20, color: '#f5222d' }}
                />
              </Col>
            </Row>

            <Statistic
              title="溶剂"
              value={solvents.length}
              suffix="种"
              valueStyle={{ fontSize: 20, color: '#52c41a' }}
            />

            {/* 已选离子列表 */}
            {selectedCations.length > 0 && (
              <div>
                <Text strong style={{ fontSize: 12 }}>阳离子：</Text>
                <div style={{ marginTop: 4 }}>
                  {selectedCations.filter(ion => ion && ion.name).map(ion => (
                    <Tag key={ion.name} color="blue" style={{ marginBottom: 4 }}>
                      {ion.name}: {ion.concentration} M
                    </Tag>
                  ))}
                </div>
              </div>
            )}

            {selectedAnions.length > 0 && (
              <div>
                <Text strong style={{ fontSize: 12 }}>阴离子：</Text>
                <div style={{ marginTop: 4 }}>
                  {selectedAnions.filter(ion => ion && ion.name).map(ion => (
                    <Tag key={ion.name} color="red" style={{ marginBottom: 4 }}>
                      {ion.name}: {ion.concentration} M
                    </Tag>
                  ))}
                </div>
              </div>
            )}

            {/* 已选溶剂列表 */}
            {solvents.length > 0 && (
              <div>
                <Text strong style={{ fontSize: 12 }}>溶剂：</Text>
                <div style={{ marginTop: 4 }}>
                  {solvents.filter((s: any) => s && s.name).map((s: any, idx: number) => (
                    <Tag key={idx} color="green" style={{ marginBottom: 4 }}>
                      {s.name}: {s.molar_ratio}
                    </Tag>
                  ))}
                </div>
              </div>
            )}
          </Space>
        </Card>
      </Col>
    </Row>

    {/* 保存为常用组合的对话框 */}
    <Modal
      title={<><StarOutlined /> 保存为我的常用组合</>}
      open={saveModalVisible}
      onOk={handleSaveAsCombination}
      onCancel={() => {
        setSaveModalVisible(false);
        setSaveName('');
        setSaveDescription('');
      }}
      okText="保存"
      cancelText="取消"
    >
      <Form layout="vertical">
        <Form.Item label="组合名称" required>
          <Input
            value={saveName}
            onChange={(e) => setSaveName(e.target.value)}
            placeholder="例如：我的EC/DMC 3:7"
            maxLength={100}
          />
        </Form.Item>
        <Form.Item label="描述（可选）">
          <Input.TextArea
            value={saveDescription}
            onChange={(e) => setSaveDescription(e.target.value)}
            placeholder="例如：常用于锂离子电池，低温性能好"
            rows={3}
            maxLength={500}
          />
        </Form.Item>
        <Alert
          message="提示"
          description={
            <div>
              <div>当前溶剂配置：</div>
              <div style={{ marginTop: 8 }}>
                {(form.getFieldValue('solvents') || []).filter((s: any) => s && s.name).map((s: any, idx: number) => (
                  <Tag key={idx} color="green" style={{ marginBottom: 4 }}>
                    {s.name}: {s.molar_ratio}
                  </Tag>
                ))}
              </div>
            </div>
          }
          type="info"
          showIcon
        />
      </Form>
    </Modal>
  </>
  );
}
