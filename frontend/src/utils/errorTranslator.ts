/**
 * MD 计算错误信息翻译器
 * 将技术性错误信息翻译成用户友好的中文提示
 */

export interface TranslatedError {
  title: string;           // 错误标题
  description: string;     // 错误描述
  suggestion: string;      // 建议操作
  severity: 'error' | 'warning';  // 严重程度
  originalError?: string;  // 原始错误信息（供技术人员查看）
}

// 错误模式匹配规则
const errorPatterns: Array<{
  pattern: RegExp;
  translate: (match: RegExpMatchArray, original: string) => TranslatedError;
}> = [
  // Packmol 失败
  {
    pattern: /Failed to run Packmol|packmol.*fail|Packmol error/i,
    translate: (_, original) => ({
      title: '分子装配失败',
      description: '系统在将分子放入模拟盒子时遇到问题，通常是因为分子数量过多或盒子太小。',
      suggestion: '请尝试：1) 减少分子数量；2) 增大盒子尺寸；3) 检查分子结构是否正确',
      severity: 'error',
      originalError: original,
    }),
  },
  // LAMMPS Pair coeff 错误
  {
    pattern: /Pair coeff for hybrid has invalid style|pair_hybrid/i,
    translate: (_, original) => ({
      title: '力场参数错误',
      description: '分子间相互作用参数配置有误，可能是某些分子类型的力场参数缺失。',
      suggestion: '请检查配方中的分子是否都有对应的力场参数，或尝试使用不同的力场',
      severity: 'error',
      originalError: original,
    }),
  },
  // LAMMPS 一般错误
  {
    pattern: /LAMMPS.*ERROR[:\s]+(.+)/i,
    translate: (match, original) => ({
      title: 'MD 模拟运行错误',
      description: `模拟过程中遇到问题：${match[1] || '未知错误'}`,
      suggestion: '请检查模拟参数设置，或联系技术支持获取帮助',
      severity: 'error',
      originalError: original,
    }),
  },
  // Slurm 超时
  {
    pattern: /TIMEOUT|time.*limit|超时/i,
    translate: (_, original) => ({
      title: '计算超时',
      description: '模拟任务在规定时间内未能完成。',
      suggestion: '请尝试：1) 减少模拟步数；2) 申请更长的计算时间；3) 使用更多计算资源',
      severity: 'error',
      originalError: original,
    }),
  },
  // Slurm 内存不足
  {
    pattern: /OUT.*OF.*MEMORY|oom|内存不足|memory/i,
    translate: (_, original) => ({
      title: '内存不足',
      description: '计算过程中内存使用超出限制。',
      suggestion: '请尝试：1) 减少分子数量；2) 减小盒子尺寸；3) 申请更多内存资源',
      severity: 'error',
      originalError: original,
    }),
  },
  // Slurm 取消
  {
    pattern: /CANCELLED|用户取消/i,
    translate: (_, original) => ({
      title: '任务已取消',
      description: '该计算任务已被取消。',
      suggestion: '如需重新计算，请复制任务配置后重新提交',
      severity: 'warning',
      originalError: original,
    }),
  },
  // Slurm 节点故障
  {
    pattern: /NODE.*FAIL|节点故障/i,
    translate: (_, original) => ({
      title: '计算节点故障',
      description: '分配的计算节点发生故障。',
      suggestion: '请重新提交任务，系统会自动分配其他可用节点',
      severity: 'error',
      originalError: original,
    }),
  },
  // RESP 电荷计算失败
  {
    pattern: /RESP.*fail|resp.*error|电荷计算.*失败/i,
    translate: (_, original) => ({
      title: 'RESP 电荷计算失败',
      description: '高精度电荷计算过程中遇到问题。',
      suggestion: '请尝试：1) 检查分子结构是否正确；2) 使用 LigParGen 快速电荷计算方法',
      severity: 'error',
      originalError: original,
    }),
  },
  // 文件读写错误
  {
    pattern: /file.*not.*found|cannot.*open|IO.*error|读取.*失败|写入.*失败/i,
    translate: (_, original) => ({
      title: '文件读写错误',
      description: '计算过程中无法读取或写入必要文件。',
      suggestion: '这可能是系统临时问题，请稍后重试或联系技术支持',
      severity: 'error',
      originalError: original,
    }),
  },
  // 分子结构错误
  {
    pattern: /invalid.*smiles|smiles.*error|分子结构.*错误|structure.*error/i,
    translate: (_, original) => ({
      title: '分子结构错误',
      description: '输入的分子结构（SMILES）格式不正确或无法解析。',
      suggestion: '请检查分子的 SMILES 格式是否正确，可使用在线工具验证',
      severity: 'error',
      originalError: original,
    }),
  },
  // Slurm 一般失败（带退出码）
  {
    pattern: /Slurm状态:\s*FAILED.*退出码:\s*(\d+):(\d+)/i,
    translate: (match, original) => ({
      title: '计算任务失败',
      description: `任务运行过程中发生错误（错误码：${match[1]}）。`,
      suggestion: '请检查任务配置是否正确，或查看详细日志获取更多信息',
      severity: 'error',
      originalError: original,
    }),
  },
  // 原子重叠/距离太近
  {
    pattern: /atoms.*too.*close|overlap|原子.*重叠|距离.*太近/i,
    translate: (_, original) => ({
      title: '原子距离过近',
      description: '初始结构中存在原子过于接近的情况，可能导致能量计算异常。',
      suggestion: '请尝试：1) 增大盒子尺寸；2) 减少分子数量；3) 重新生成初始结构',
      severity: 'error',
      originalError: original,
    }),
  },
  // 能量发散
  {
    pattern: /energy.*diverge|NaN|Inf|能量.*发散|数值.*异常/i,
    translate: (_, original) => ({
      title: '模拟数值不稳定',
      description: '模拟过程中能量计算出现异常值，可能是参数设置不当或初始结构问题。',
      suggestion: '请尝试：1) 减小时间步长；2) 检查初始结构；3) 降低温度后逐步升温',
      severity: 'error',
      originalError: original,
    }),
  },
];

/**
 * 翻译错误信息
 */
export function translateError(errorMessage: string | null | undefined): TranslatedError | null {
  if (!errorMessage) {
    return null;
  }

  // 尝试匹配已知错误模式
  for (const { pattern, translate } of errorPatterns) {
    const match = errorMessage.match(pattern);
    if (match) {
      return translate(match, errorMessage);
    }
  }

  // 未匹配到已知模式，返回通用错误
  return {
    title: '计算任务失败',
    description: '任务运行过程中发生错误。',
    suggestion: '请查看详细日志或联系技术支持获取帮助',
    severity: 'error',
    originalError: errorMessage,
  };
}

/**
 * 获取错误的简短描述（用于卡片显示）
 */
export function getErrorSummary(errorMessage: string | null | undefined): string {
  const translated = translateError(errorMessage);
  if (!translated) {
    return '';
  }
  return translated.title;
}

