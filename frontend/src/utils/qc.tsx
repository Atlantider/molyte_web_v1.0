/**
 * QC 相关工具函数
 */
import { Tag } from 'antd';

/**
 * 渲染溶剂模型标签
 */
export function renderSolventModel(config: any) {
  const solventConfig = config?.solvent_config;
  
  if (!solventConfig) {
    return <Tag color="default">气相</Tag>;
  }
  
  const model = solventConfig.model || 'gas';
  const solvent = solventConfig.solvent_name;
  
  if (model === 'gas') {
    return <Tag color="default">气相</Tag>;
  } else if (model === 'pcm') {
    return (
      <Tag color="blue">
        PCM{solvent ? ` (${solvent})` : ''}
      </Tag>
    );
  } else if (model === 'smd') {
    return (
      <Tag color="green">
        SMD{solvent ? ` (${solvent})` : ''}
      </Tag>
    );
  }
  
  return <Tag>{model}</Tag>;
}

/**
 * 获取溶剂模型的文本描述
 */
export function getSolventModelText(config: any): string {
  const solventConfig = config?.solvent_config;
  
  if (!solventConfig) {
    return '气相';
  }
  
  const model = solventConfig.model || 'gas';
  const solvent = solventConfig.solvent_name;
  
  if (model === 'gas') {
    return '气相';
  } else if (model === 'pcm') {
    return `PCM${solvent ? ` (${solvent})` : ''}`;
  } else if (model === 'smd') {
    return `SMD${solvent ? ` (${solvent})` : ''}`;
  }
  
  return model;
}

