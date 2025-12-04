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
      description: '盒子太小或分子太多',
      suggestion: '请减少分子数量或增大盒子尺寸',
      severity: 'error',
      originalError: original,
    }),
  },
  {
    pattern: /Pair coeff for hybrid has invalid style|pair_hybrid/i,
    translate: (_, original) => ({
      title: '力场参数错误',
      description: '力场参数缺失',
      suggestion: '请检查分子是否有对应的力场参数',
      severity: 'error',
      originalError: original,
    }),
  },
  {
    pattern: /LAMMPS.*ERROR[:\s]+(.+)/i,
    translate: (_, original) => ({
      title: 'MD 模拟错误',
      description: '模拟运行出错',
      suggestion: '请检查参数设置或联系管理员',
      severity: 'error',
      originalError: original,
    }),
  },
  {
    pattern: /TIMEOUT|time.*limit|超时/i,
    translate: (_, original) => ({
      title: '计算超时',
      description: '计算时间不足',
      suggestion: '请减少模拟步数或增加计算时间',
      severity: 'error',
      originalError: original,
    }),
  },
  {
    pattern: /OUT.*OF.*MEMORY|oom|内存不足|memory/i,
    translate: (_, original) => ({
      title: '内存不足',
      description: '内存超限',
      suggestion: '请减少分子数量或盒子尺寸',
      severity: 'error',
      originalError: original,
    }),
  },
  {
    pattern: /CANCELLED|用户取消/i,
    translate: (_, original) => ({
      title: '任务已取消',
      description: '任务被取消',
      suggestion: '如需重新计算请复制配置后重新提交',
      severity: 'warning',
      originalError: original,
    }),
  },
  {
    pattern: /Slurm状态:\s*FAILED.*退出码:\s*(\d+):(\d+)/i,
    translate: (_, original) => ({
      title: '计算失败',
      description: '任务运行出错',
      suggestion: '请检查配置或查看日志',
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
    title: '计算失败',
    description: '任务运行出错',
    suggestion: '请查看日志或联系管理员',
    severity: 'error',
    originalError: errorMessage,
  };
}

