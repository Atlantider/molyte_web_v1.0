/**
 * MD 计算错误信息翻译器
 * 将技术性错误信息翻译成用户友好的中文提示
 */

export interface TranslatedError {
  title: string;
  description: string;
  suggestion: string;
  severity: 'error' | 'warning';
  originalError?: string;
}

const errorPatterns: Array<{
  pattern: RegExp;
  translate: (match: RegExpMatchArray, original: string) => TranslatedError;
}> = [
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
  {
    pattern: /TIMEOUT|time.*limit|超时/i,
    translate: (_, original) => ({
      title: '计算超时',
      description: '模拟任务在规定时间内未能完成。',
      suggestion: '请尝试：1) 减少模拟步数；2) 申请更长的计算时间',
      severity: 'error',
      originalError: original,
    }),
  },
  {
    pattern: /OUT.*OF.*MEMORY|oom|内存不足|memory/i,
    translate: (_, original) => ({
      title: '内存不足',
      description: '计算过程中内存使用超出限制。',
      suggestion: '请尝试：1) 减少分子数量；2) 减小盒子尺寸',
      severity: 'error',
      originalError: original,
    }),
  },
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
];

export function translateError(errorMessage: string | null | undefined): TranslatedError | null {
  if (!errorMessage) return null;

  for (const { pattern, translate } of errorPatterns) {
    const match = errorMessage.match(pattern);
    if (match) return translate(match, errorMessage);
  }

  return {
    title: '计算任务失败',
    description: '任务运行过程中发生错误。',
    suggestion: '请查看详细日志或联系技术支持获取帮助',
    severity: 'error',
    originalError: errorMessage,
  };
}

