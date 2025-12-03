/**
 * 新版电解质配方表单组件
 * 使用勾选框选择离子，输入浓度，自动验证电中性
 */
import { useState, useEffect } from 'react';
import {
  Form,
  Select,
  Input,
  InputNumber,
  Checkbox,
  Button,
  Space,
  Divider,
  Row,
  Col,
  Card,
  message,
  Spin,
  Alert,
  AutoComplete,
  theme,
} from 'antd';
import { PlusOutlined, MinusCircleOutlined, CheckCircleOutlined, CloseCircleOutlined } from '@ant-design/icons';
import type { IonInfo, Project } from '../types';
import { getAvailableIons } from '../api/electrolytes';
import { useThemeStore } from '../stores/themeStore';

// 按分类的常用溶剂列表
const SOLVENT_CATEGORIES = [
  {
    category: '碳酸酯类',
    description: '锂离子电池最常用溶剂',
    solvents: [
      { name: 'EC', label: '碳酸乙烯酯 (EC)', smiles: 'C1COC(=O)O1' },
      { name: 'PC', label: '碳酸丙烯酯 (PC)', smiles: 'CC1COC(=O)O1' },
      { name: 'DMC', label: '碳酸二甲酯 (DMC)', smiles: 'COC(=O)OC' },
      { name: 'DEC', label: '碳酸二乙酯 (DEC)', smiles: 'CCOC(=O)OCC' },
      { name: 'EMC', label: '碳酸甲乙酯 (EMC)', smiles: 'CCOC(=O)OC' },
      { name: 'VC', label: '碳酸亚乙烯酯 (VC)', smiles: 'C1=COC(=O)O1' },
    ],
  },
  {
    category: '氟代碳酸酯',
    description: '高电压/锂金属电池用',
    solvents: [
      { name: 'FEC', label: '氟代碳酸乙烯酯 (FEC)', smiles: 'C1C(OC(=O)O1)F' },
      { name: 'DFEC', label: '二氟碳酸乙烯酯 (DFEC)', smiles: 'FC1OC(=O)OC1F' },
      { name: 'TFPC', label: '三氟碳酸丙烯酯 (TFPC)', smiles: 'CC(F)(F)C1COC(=O)O1' },
    ],
  },
  {
    category: '醚类',
    description: '锂硫/锂空电池常用',
    solvents: [
      { name: 'DME', label: '乙二醇二甲醚 (DME)', smiles: 'COCCOC' },
      { name: 'DOL', label: '1,3-二氧戊环 (DOL)', smiles: 'C1COCO1' },
      { name: 'DIOX', label: '1,4-二氧六环 (DIOX)', smiles: 'C1COCCO1' },
      { name: 'THF', label: '四氢呋喃 (THF)', smiles: 'C1CCOC1' },
      { name: 'DEGDME', label: '二乙二醇二甲醚 (DEGDME)', smiles: 'COCCOCCOC' },
      { name: 'TEGDME', label: '四乙二醇二甲醚 (TEGDME)', smiles: 'COCCOCCOCCOCCOC' },
      { name: '2-MeTHF', label: '2-甲基四氢呋喃 (2-MeTHF)', smiles: 'CC1CCCO1' },
    ],
  },
  {
    category: '氟代醚类(稀释剂)',
    description: '局部高浓度电解液稀释剂',
    solvents: [
      { name: 'TTE', label: '1,1,2,2-四氟乙基-2,2,3,3-四氟丙基醚 (TTE)', smiles: 'FC(F)C(F)(F)OCC(F)(F)C(F)(F)F' },
      { name: 'BTFE', label: '双(2,2,2-三氟乙基)醚 (BTFE)', smiles: 'FC(F)(F)COCC(F)(F)F' },
      { name: 'TFEE', label: '1,1,2,2-四氟乙基-2,2,2-三氟乙基醚 (TFEE)', smiles: 'FC(F)C(F)(F)OCC(F)(F)F' },
      { name: 'HFE', label: '1H,1H,5H-八氟戊基-1,1,2,2-四氟乙基醚 (HFE)', smiles: 'FC(F)C(F)(F)OCC(F)(F)C(F)(F)C(F)(F)CF' },
      { name: 'FDEE', label: '1,2-二氟乙基乙基醚 (FDEE)', smiles: 'CCOCC(F)F' },
    ],
  },
  {
    category: '腈类',
    description: '高电压电解液添加剂',
    solvents: [
      { name: 'ACN', label: '乙腈 (ACN)', smiles: 'CC#N' },
      { name: 'SN', label: '琥珀腈 (SN)', smiles: 'N#CCCC#N' },
      { name: 'AN', label: '己二腈 (AN)', smiles: 'N#CCCCCC#N' },
      { name: 'GN', label: '戊二腈 (GN)', smiles: 'N#CCCCC#N' },
    ],
  },
  {
    category: '砜类',
    description: '高电压/高安全性溶剂',
    solvents: [
      { name: 'Sulfolane', label: '环丁砜 (Sulfolane)', smiles: 'C1CCS(=O)(=O)C1' },
      { name: 'DMS', label: '二甲基砜 (DMS)', smiles: 'CS(=O)(=O)C' },
      { name: 'EMS', label: '乙基甲基砜 (EMS)', smiles: 'CCS(=O)(=O)C' },
      { name: 'DMSO', label: '二甲亚砜 (DMSO)', smiles: 'CS(=O)C' },
    ],
  },
  {
    category: '酯类',
    description: '低温电解液用',
    solvents: [
      { name: 'EA', label: '乙酸乙酯 (EA)', smiles: 'CCOC(=O)C' },
      { name: 'MA', label: '乙酸甲酯 (MA)', smiles: 'COC(=O)C' },
      { name: 'EP', label: '丙酸乙酯 (EP)', smiles: 'CCOC(=O)CC' },
      { name: 'MP', label: '丙酸甲酯 (MP)', smiles: 'COC(=O)CC' },
      { name: 'MB', label: '丁酸甲酯 (MB)', smiles: 'CCCC(=O)OC' },
      { name: 'GBL', label: 'γ-丁内酯 (GBL)', smiles: 'C1CC(=O)OC1' },
    ],
  },
  {
    category: '酰胺类',
    description: '高溶解性溶剂',
    solvents: [
      { name: 'DMF', label: 'N,N-二甲基甲酰胺 (DMF)', smiles: 'CN(C)C=O' },
      { name: 'DMAc', label: 'N,N-二甲基乙酰胺 (DMAc)', smiles: 'CC(=O)N(C)C' },
      { name: 'NMP', label: 'N-甲基吡咯烷酮 (NMP)', smiles: 'CN1CCCC1=O' },
      { name: 'TMU', label: '四甲基脲 (TMU)', smiles: 'CN(C)C(=O)N(C)C' },
    ],
  },
  {
    category: '磷酸酯类',
    description: '阻燃溶剂',
    solvents: [
      { name: 'TMP', label: '磷酸三甲酯 (TMP)', smiles: 'COP(=O)(OC)OC' },
      { name: 'TEP', label: '磷酸三乙酯 (TEP)', smiles: 'CCOP(=O)(OCC)OCC' },
      { name: 'DMMP', label: '二甲基甲基膦酸酯 (DMMP)', smiles: 'COP(=O)(C)OC' },
      { name: 'TFP', label: '三(2,2,2-三氟乙基)磷酸酯 (TFP)', smiles: 'FC(F)(F)COP(=O)(OCC(F)(F)F)OCC(F)(F)F' },
    ],
  },
  {
    category: '其他常用溶剂',
    description: '水系及其他',
    solvents: [
      { name: 'Water', label: '水 (Water)', smiles: 'O' },
      { name: 'Methanol', label: '甲醇 (Methanol)', smiles: 'CO' },
      { name: 'Ethanol', label: '乙醇 (Ethanol)', smiles: 'CCO' },
      { name: 'Acetone', label: '丙酮 (Acetone)', smiles: 'CC(=O)C' },
      { name: 'IPA', label: '异丙醇 (IPA)', smiles: 'CC(C)O' },
    ],
  },
];

// 扁平化的溶剂列表（用于AutoComplete）
const ALL_SOLVENTS = SOLVENT_CATEGORIES.flatMap(cat =>
  cat.solvents.map(s => ({ ...s, category: cat.category }))
);

interface ElectrolyteFormNewProps {
  form: any;
  projects: Project[];
  onValuesChange?: (changedValues: any, allValues: any) => void;
  onIonsChange?: (cations: SelectedIon[], anions: SelectedIon[]) => void;
  initialCations?: SelectedIon[];
  initialAnions?: SelectedIon[];
}

interface SelectedIon {
  name: string;
  charge: number;
  concentration: number; // mol/L
}

interface Solvent {
  name: string;
  smiles: string;
  molar_ratio: number;
  validated?: boolean;
  molecule_info?: any;
}

export default function ElectrolyteFormNew({
  form,
  projects,
  onValuesChange,
  onIonsChange,
  initialCations = [],
  initialAnions = []
}: ElectrolyteFormNewProps) {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [availableCations, setAvailableCations] = useState<IonInfo[]>([]);
  const [availableAnions, setAvailableAnions] = useState<IonInfo[]>([]);
  const [selectedCations, setSelectedCations] = useState<SelectedIon[]>(initialCations);
  const [selectedAnions, setSelectedAnions] = useState<SelectedIon[]>(initialAnions);
  const [loading, setLoading] = useState(false);
  const [isElectricallyNeutral, setIsElectricallyNeutral] = useState<boolean | null>(null);
  const [solvents, setSolvents] = useState<any[]>([]);
  const [userModifiedName, setUserModifiedName] = useState(false); // 跟踪用户是否手动修改过名称
  const [lastGeneratedName, setLastGeneratedName] = useState(''); // 记录上次自动生成的名称

  // 加载可用离子列表
  useEffect(() => {
    loadAvailableIons();
    // 初始化溶剂状态
    const initialSolvents = form.getFieldValue('solvents') || [];
    setSolvents(initialSolvents);
  }, []);

  // 当初始离子改变时更新选中状态
  useEffect(() => {
    if (initialCations.length > 0) {
      setSelectedCations(initialCations);
    }
  }, [initialCations]);

  useEffect(() => {
    if (initialAnions.length > 0) {
      setSelectedAnions(initialAnions);
    }
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

  // 检查电中性
  const checkElectricalNeutrality = () => {
    if (selectedCations.length === 0 || selectedAnions.length === 0) {
      setIsElectricallyNeutral(null);
      return;
    }

    // 计算总正电荷 = Σ(浓度 × 电荷)
    const totalPositiveCharge = selectedCations.reduce(
      (sum, ion) => sum + ion.concentration * ion.charge,
      0
    );

    // 计算总负电荷 = Σ(浓度 × |电荷|)
    const totalNegativeCharge = selectedAnions.reduce(
      (sum, ion) => sum + ion.concentration * Math.abs(ion.charge),
      0
    );

    // 允许 1% 的误差
    const isNeutral = Math.abs(totalPositiveCharge - totalNegativeCharge) < 0.01;
    setIsElectricallyNeutral(isNeutral);
  };

  // 当离子选择或浓度改变时，检查电中性
  useEffect(() => {
    checkElectricalNeutrality();
    // 通知父组件离子变化
    if (onIonsChange) {
      onIonsChange(selectedCations, selectedAnions);
    }
  }, [selectedCations, selectedAnions]);

  // 处理阳离子勾选
  const handleCationCheck = (checked: boolean, ion: IonInfo) => {
    if (checked) {
      setSelectedCations([...selectedCations, { ...ion, concentration: 1.0 }]);
    } else {
      setSelectedCations(selectedCations.filter((c) => c.name !== ion.name));
    }
  };

  // 处理阴离子勾选
  const handleAnionCheck = (checked: boolean, ion: IonInfo) => {
    if (checked) {
      setSelectedAnions([...selectedAnions, { ...ion, concentration: 1.0 }]);
    } else {
      setSelectedAnions(selectedAnions.filter((a) => a.name !== ion.name));
    }
  };

  // 更新阳离子浓度
  const updateCationConcentration = (name: string, concentration: number) => {
    setSelectedCations(
      selectedCations.map((c) => (c.name === name ? { ...c, concentration } : c))
    );
  };

  // 更新阴离子浓度
  const updateAnionConcentration = (name: string, concentration: number) => {
    setSelectedAnions(
      selectedAnions.map((a) => (a.name === name ? { ...a, concentration } : a))
    );
  };

  // 生成默认配方名称
  const generateDefaultName = () => {
    // 如果没有选择离子，返回空
    if (selectedCations.length === 0 && selectedAnions.length === 0 && solvents.length === 0) {
      return '';
    }

    const parts: string[] = [];

    // 添加阳离子（按浓度降序排列，取前2个）
    // 清理名称中的特殊字符，将/替换为-
    const sortedCations = [...selectedCations].sort((a, b) => b.concentration - a.concentration);
    const topCations = sortedCations.slice(0, 2).map(c => c.name?.replace(/\//g, '-')).filter(Boolean);
    if (topCations.length > 0) {
      parts.push(topCations.join('-'));
    }

    // 添加阴离子（按浓度降序排列，取前2个）
    const sortedAnions = [...selectedAnions].sort((a, b) => b.concentration - a.concentration);
    const topAnions = sortedAnions.slice(0, 2).map(a => a.name?.replace(/\//g, '-')).filter(Boolean);
    if (topAnions.length > 0) {
      parts.push(topAnions.join('-'));
    }

    // 添加溶剂（按摩尔比降序排列，取前3个）
    // 安全检查：确保溶剂对象存在且有 name 属性
    // 清理溶剂名称中的特殊字符，将/替换为-
    const validSolvents = solvents.filter((s: any) => s && s.name);
    const sortedSolvents = [...validSolvents].sort((a: any, b: any) => (b.molar_ratio || 0) - (a.molar_ratio || 0));
    const topSolvents = sortedSolvents.slice(0, 3).map((s: any) => s.name?.replace(/\//g, '-')).filter(Boolean);
    if (topSolvents.length > 0) {
      parts.push(topSolvents.join('-'));
    }

    return parts.join('-');
  };

  // 当离子或溶剂变化时，自动更新配方名称
  useEffect(() => {
    // 安全检查：确保所有数据都已初始化
    if (!form) return;

    const currentName = form.getFieldValue('name');

    // 如果用户手动修改过名称，且当前名称不是上次自动生成的名称，则不再自动更新
    if (userModifiedName && currentName !== lastGeneratedName) {
      return;
    }

    // 生成新名称
    try {
      const newName = generateDefaultName();
      if (newName && newName !== currentName) {
        form.setFieldsValue({ name: newName });
        setLastGeneratedName(newName);
        setUserModifiedName(false); // 重置标记，因为这是自动生成的
      }
    } catch (error) {
      console.error('生成默认名称时出错:', error);
    }
  }, [selectedCations, selectedAnions, solvents]);

  if (loading) {
    return <Spin tip="加载离子列表..." />;
  }

  // 处理表单值变化，用于监听溶剂变化
  const handleFormChange = (changedValues: any, allValues: any) => {
    // 如果溶剂发生变化，更新 solvents 状态
    if (changedValues.solvents) {
      setSolvents(allValues.solvents || []);
    }

    // 如果名称发生变化，标记为用户手动修改
    if (changedValues.name !== undefined) {
      setUserModifiedName(true);
    }

    // 调用父组件的 onValuesChange
    if (onValuesChange) {
      onValuesChange(changedValues, allValues);
    }
  };

  return (
    <div>
      <Form form={form} layout="vertical" onValuesChange={handleFormChange}>
        {/* 基本信息 */}
        <Form.Item
          name="project_id"
          label="所属项目"
          rules={[{ required: true, message: '请选择所属项目' }]}
        >
          <Select placeholder="选择项目">
            {projects.map((p) => (
              <Select.Option key={p.id} value={p.id}>
                {p.name}
              </Select.Option>
            ))}
          </Select>
        </Form.Item>

        <Form.Item
          name="name"
          label="配方名称（可选）"
          rules={[
            { max: 100, message: '名称不能超过100个字符' },
          ]}
          tooltip="留空将根据离子和溶剂自动生成，例如：Li-PF6-EC/DMC"
          extra="留空自动生成，或手动输入自定义名称"
        >
          <Input placeholder="留空自动生成（推荐），或手动输入" allowClear />
        </Form.Item>

        <Form.Item
          name="temperature"
          label="温度 (K)"
          rules={[{ required: true, message: '请输入温度' }]}
          initialValue={298.15}
        >
          <InputNumber min={0} max={1000} step={0.01} style={{ width: '100%' }} />
        </Form.Item>

        {/* 阳离子选择 */}
        <Divider>阳离子配置</Divider>
        <Card size="small" style={{ marginBottom: 16 }}>
          <Alert
            message="请设置各离子的浓度（mol/L），浓度决定了模拟体系中离子的数量"
            type="info"
            showIcon
            style={{ marginBottom: 12 }}
          />
          <Row gutter={[16, 16]}>
            {availableCations.map((ion, index) => {
              const selected = selectedCations.find((c) => c.name === ion.name);
              const isFirstCation = index === 0 && selected;
              return (
                <Col span={24} key={ion.name}>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    padding: selected ? '8px 12px' : '4px 0',
                    background: selected ? (isFirstCation ? '#e6f7ff' : '#f6ffed') : 'transparent',
                    borderRadius: 6,
                    border: selected ? (isFirstCation ? '2px solid #1890ff' : '1px solid #b7eb8f') : 'none',
                  }}>
                    <Checkbox
                      checked={!!selected}
                      onChange={(e) => handleCationCheck(e.target.checked, ion)}
                    >
                      <strong style={{ fontSize: 14 }}>{ion.name}</strong>
                      <span style={{ marginLeft: 8, color: '#999' }}>
                        (电荷: {ion.charge > 0 ? '+' : ''}{ion.charge})
                      </span>
                      {isFirstCation && (
                        <span style={{ marginLeft: 8, color: '#1890ff', fontWeight: 500, fontSize: 12 }}>
                          ★ 溶剂比例参考此离子
                        </span>
                      )}
                    </Checkbox>
                    {selected && (
                      <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                        <span style={{ fontWeight: 600, color: '#fa541c' }}>浓度:</span>
                        <InputNumber
                          min={0.001}
                          max={10}
                          step={0.1}
                          value={selected.concentration}
                          onChange={(value) => updateCationConcentration(ion.name, value || 0)}
                          style={{ width: 130, fontWeight: 500 }}
                          addonAfter="mol/L"
                        />
                      </div>
                    )}
                  </div>
                </Col>
              );
            })}
          </Row>
          {selectedCations.length === 0 && (
            <Alert message="请至少选择一种阳离子" type="warning" showIcon style={{ marginTop: 8 }} />
          )}
        </Card>

        {/* 阴离子选择 */}
        <Divider>阴离子配置</Divider>
        <Card size="small" style={{ marginBottom: 16 }}>
          <Row gutter={[16, 16]}>
            {availableAnions.map((ion) => {
              const selected = selectedAnions.find((a) => a.name === ion.name);
              return (
                <Col span={24} key={ion.name}>
                  <div style={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    padding: selected ? '8px 12px' : '4px 0',
                    background: selected ? '#fff7e6' : 'transparent',
                    borderRadius: 6,
                    border: selected ? '1px solid #ffd591' : 'none',
                  }}>
                    <Checkbox
                      checked={!!selected}
                      onChange={(e) => handleAnionCheck(e.target.checked, ion)}
                    >
                      <strong style={{ fontSize: 14 }}>{ion.name}</strong>
                      <span style={{ marginLeft: 8, color: '#999' }}>
                        (电荷: {ion.charge})
                      </span>
                    </Checkbox>
                    {selected && (
                      <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                        <span style={{ fontWeight: 600, color: '#fa541c' }}>浓度:</span>
                        <InputNumber
                          min={0.001}
                          max={10}
                          step={0.1}
                          value={selected.concentration}
                          onChange={(value) => updateAnionConcentration(ion.name, value || 0)}
                          style={{ width: 130, fontWeight: 500 }}
                          addonAfter="mol/L"
                        />
                      </div>
                    )}
                  </div>
                </Col>
              );
            })}
          </Row>
          {selectedAnions.length === 0 && (
            <Alert message="请至少选择一种阴离子" type="warning" showIcon style={{ marginTop: 8 }} />
          )}
        </Card>

        {/* 电中性检查 */}
        {isElectricallyNeutral !== null && (
          <Alert
            message={
              isElectricallyNeutral ? (
                <Space>
                  <CheckCircleOutlined />
                  电荷平衡：体系满足电中性条件
                </Space>
              ) : (
                <Space>
                  <CloseCircleOutlined />
                  电荷不平衡：请调整离子浓度以满足电中性条件
                </Space>
              )
            }
            type={isElectricallyNeutral ? 'success' : 'error'}
            showIcon
            style={{ marginBottom: 16 }}
          />
        )}

        {/* 溶剂配置 */}
        <Divider>溶剂配置</Divider>

        {/* 重要提示 */}
        <Alert
          message={
            <span style={{ fontWeight: 600 }}>
              ⚠️ 重要提示：溶剂摩尔比是相对于第一种阳离子
              {selectedCations.length > 0 && (
                <span style={{ color: '#1890ff', marginLeft: 8 }}>
                  （当前参考：{selectedCations[0]?.name}，浓度 {selectedCations[0]?.concentration} mol/L）
                </span>
              )}
            </span>
          }
          description={
            <div style={{ marginTop: 4 }}>
              例如：若第一种阳离子浓度为 1 mol/L，溶剂摩尔比设为 3，则该溶剂浓度为 3 mol/L
            </div>
          }
          type="warning"
          showIcon
          style={{ marginBottom: 16, border: '2px solid #faad14' }}
        />

        <Form.List name="solvents" initialValue={[]}>
          {(fields, { add, remove }) => (
            <>
              {/* 分类选择器 - 快速从常用溶剂中选择 */}
              <Card size="small" style={{ marginBottom: 16, background: token.colorBgContainer }}>
                <div style={{ marginBottom: 12 }}>
                  <span style={{ fontWeight: 500, marginRight: 8 }}>从常用溶剂库选择:</span>
                  <span style={{ fontSize: 12, color: '#666' }}>点击溶剂名称即可添加</span>
                </div>
                {SOLVENT_CATEGORIES.map((cat, catIdx) => (
                  <div key={catIdx} style={{ marginBottom: 12 }}>
                    <div style={{ fontSize: 12, color: '#1890ff', marginBottom: 4 }}>
                      <strong>{cat.category}</strong>
                      <span style={{ color: '#999', marginLeft: 8 }}>{cat.description}</span>
                    </div>
                    <Space wrap size={[4, 4]}>
                      {cat.solvents.map((solvent, idx) => (
                        <Button
                          key={idx}
                          size="small"
                          type="default"
                          onClick={() => add({ name: solvent.name, smiles: solvent.smiles, molar_ratio: 1.0 })}
                          style={{ fontSize: 12 }}
                        >
                          {solvent.name}
                        </Button>
                      ))}
                    </Space>
                  </div>
                ))}
              </Card>

              {/* 已添加的溶剂列表 */}
              {fields.map((field, index) => {
                // 处理溶剂选择变化，自动填充SMILES
                const handleSolventSelect = (value: string) => {
                  const selectedSolvent = ALL_SOLVENTS.find(s => s.name === value);
                  if (selectedSolvent) {
                    const currentSolvents = form.getFieldValue('solvents') || [];
                    currentSolvents[field.name] = {
                      ...currentSolvents[field.name],
                      name: selectedSolvent.name,
                      smiles: selectedSolvent.smiles,
                    };
                    form.setFieldsValue({ solvents: currentSolvents });
                  }
                };

                return (
                  <Card
                    key={field.key}
                    size="small"
                    title={`溶剂 ${index + 1}`}
                    extra={
                      <Button
                        type="link"
                        danger
                        icon={<MinusCircleOutlined />}
                        onClick={() => remove(field.name)}
                      >
                        删除
                      </Button>
                    }
                    style={{ marginBottom: 16 }}
                  >
                    <Row gutter={16}>
                      <Col span={8}>
                        <Form.Item
                          {...field}
                          name={[field.name, 'name']}
                          label="溶剂名称"
                          rules={[{ required: true, message: '请选择或输入溶剂名称' }]}
                          tooltip="可搜索常用溶剂，或直接输入自定义名称"
                        >
                          <AutoComplete
                            placeholder="搜索或输入溶剂名称"
                            options={SOLVENT_CATEGORIES.map(cat => ({
                              label: <span style={{ fontWeight: 500, color: '#1890ff' }}>{cat.category}</span>,
                              options: cat.solvents.map(s => ({
                                value: s.name,
                                label: s.label,
                              })),
                            }))}
                            onSelect={handleSolventSelect}
                            filterOption={(inputValue, option) => {
                              // 只搜索子选项，不搜索分组标签
                              if (!option || 'options' in option) {
                                return false; // 这是分组项或空值，不直接匹配
                              }
                              const optionAny = option as any;
                              const label = optionAny?.label?.toString().toLowerCase() || '';
                              const value = optionAny?.value?.toString().toLowerCase() || '';
                              return label.includes(inputValue.toLowerCase()) || value.includes(inputValue.toLowerCase());
                            }}
                          />
                        </Form.Item>
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          {...field}
                          name={[field.name, 'smiles']}
                          label="SMILES"
                          rules={[{ required: true, message: '请输入 SMILES' }]}
                          tooltip="选择常用溶剂后自动填充，也可手动修改"
                        >
                          <Input placeholder="例如: C1COC(=O)O1" />
                        </Form.Item>
                      </Col>
                      <Col span={8}>
                        <Form.Item
                          {...field}
                          name={[field.name, 'molar_ratio']}
                          label={
                            <span style={{ color: '#fa541c', fontWeight: 600 }}>
                              摩尔比
                              {selectedCations.length > 0 && (
                                <span style={{ fontWeight: 400, color: '#666', marginLeft: 4 }}>
                                  (相对 {selectedCations[0]?.name})
                                </span>
                              )}
                            </span>
                          }
                          rules={[{ required: true, message: '请输入摩尔比' }]}
                          initialValue={1.0}
                          tooltip={`摩尔比 = 该溶剂浓度 / 第一种阳离子浓度${selectedCations.length > 0 ? `（${selectedCations[0]?.name}: ${selectedCations[0]?.concentration} mol/L）` : ''}`}
                        >
                          <InputNumber
                            min={0.1}
                            max={100}
                            step={0.1}
                            style={{ width: '100%', fontWeight: 500 }}
                          />
                        </Form.Item>
                      </Col>
                    </Row>
                  </Card>
                );
              })}

              {/* 自定义添加按钮 */}
              <Button type="dashed" onClick={() => add()} block icon={<PlusOutlined />} style={{ marginBottom: 12 }}>
                自定义添加溶剂
              </Button>

              {/* 常用溶剂组合快速添加 */}
              <div style={{ padding: '8px 12px', background: '#f5f5f5', borderRadius: 6 }}>
                <span style={{ fontSize: 12, color: '#666', marginRight: 8 }}>常用组合一键添加:</span>
                <Space wrap size="small">
                  <Button size="small" onClick={() => {
                    add({ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 });
                    add({ name: 'DMC', smiles: 'COC(=O)OC', molar_ratio: 1.0 });
                  }}>EC/DMC (1:1)</Button>
                  <Button size="small" onClick={() => {
                    add({ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 });
                    add({ name: 'DEC', smiles: 'CCOC(=O)OCC', molar_ratio: 1.0 });
                  }}>EC/DEC (1:1)</Button>
                  <Button size="small" onClick={() => {
                    add({ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 });
                    add({ name: 'EMC', smiles: 'CCOC(=O)OC', molar_ratio: 1.0 });
                  }}>EC/EMC (1:1)</Button>
                  <Button size="small" onClick={() => {
                    add({ name: 'DME', smiles: 'COCCOC', molar_ratio: 1.0 });
                    add({ name: 'DOL', smiles: 'C1COCO1', molar_ratio: 1.0 });
                  }}>DME/DOL (1:1)</Button>
                  <Button size="small" onClick={() => {
                    add({ name: 'EC', smiles: 'C1COC(=O)O1', molar_ratio: 1.0 });
                    add({ name: 'DMC', smiles: 'COC(=O)OC', molar_ratio: 1.0 });
                    add({ name: 'DEC', smiles: 'CCOC(=O)OCC', molar_ratio: 1.0 });
                  }}>EC/DMC/DEC</Button>
                  <Button size="small" onClick={() => {
                    add({ name: 'Water', smiles: 'O', molar_ratio: 10.0 });
                  }}>水溶液</Button>
                </Space>
              </div>
            </>
          )}
        </Form.List>

        {/* 盒子尺寸配置 */}
        <Divider>模拟盒子尺寸</Divider>
        <Card size="small" style={{ marginBottom: 16 }}>
          <Form.Item
            name="box_type"
            label="盒子类型"
            initialValue="cubic"
          >
            <Select>
              <Select.Option value="cubic">立方体（一个参数）</Select.Option>
              <Select.Option value="rectangular">长方体（三个参数）</Select.Option>
            </Select>
          </Form.Item>

          <Form.Item noStyle shouldUpdate={(prev, curr) => prev.box_type !== curr.box_type}>
            {({ getFieldValue }) => {
              const boxType = getFieldValue('box_type');
              if (boxType === 'cubic' || boxType === undefined) {
                return (
                  <Form.Item
                    name="box_size"
                    label="边长 (Å)"
                    rules={[{ required: boxType === 'cubic' || boxType === undefined, message: '请输入边长' }]}
                    initialValue={40}
                    extra="建议范围: 30-50 Å"
                  >
                    <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                  </Form.Item>
                );
              } else {
                return (
                  <Row gutter={16}>
                    <Col span={8}>
                      <Form.Item
                        name={['box_dimensions', 0]}
                        label="长 (Å)"
                        rules={[{ required: boxType === 'rectangular', message: '请输入长度' }]}
                        initialValue={40}
                      >
                        <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                      </Form.Item>
                    </Col>
                    <Col span={8}>
                      <Form.Item
                        name={['box_dimensions', 1]}
                        label="宽 (Å)"
                        rules={[{ required: boxType === 'rectangular', message: '请输入宽度' }]}
                        initialValue={40}
                      >
                        <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                      </Form.Item>
                    </Col>
                    <Col span={8}>
                      <Form.Item
                        name={['box_dimensions', 2]}
                        label="高 (Å)"
                        rules={[{ required: boxType === 'rectangular', message: '请输入高度' }]}
                        initialValue={40}
                      >
                        <InputNumber min={10} max={200} step={1} style={{ width: '100%' }} />
                      </Form.Item>
                    </Col>
                  </Row>
                );
              }
            }}
          </Form.Item>
        </Card>

        <Form.Item name="description" label="备注">
          <Input.TextArea rows={3} placeholder="可选的配方说明" />
        </Form.Item>
      </Form>
    </div>
  );
}

