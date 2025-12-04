/**
 * Cluster 高级计算结果展示面板
 */
import React, { useState, useEffect, useCallback } from 'react';
import {
  Box,
  Card,
  CardContent,
  Typography,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Chip,
  LinearProgress,
  Alert,
  Button,
  Collapse,
  IconButton,
  Tooltip,
  Grid,
  Paper,
} from '@mui/material';
import {
  ExpandMore as ExpandMoreIcon,
  ExpandLess as ExpandLessIcon,
  Refresh as RefreshIcon,
  CheckCircle as CheckCircleIcon,
  Error as ErrorIcon,
  HourglassEmpty as PendingIcon,
} from '@mui/icons-material';
import {
  getClusterAnalysisJob,
  getClusterAnalysisResults,
  getClusterAnalysisQCStatus,
  CALC_TYPE_INFO,
  AdvancedClusterJob,
  ClusterAnalysisResults,
  QCStatus,
  ClusterCalcType,
} from '../api/clusterAnalysis';

interface Props {
  jobId: number;
  onBack?: () => void;
}

const STATUS_COLORS: Record<string, 'default' | 'primary' | 'success' | 'error' | 'warning'> = {
  CREATED: 'default',
  SUBMITTED: 'primary',
  RUNNING: 'primary',
  WAITING_QC: 'warning',
  CALCULATING: 'primary',
  COMPLETED: 'success',
  FAILED: 'error',
  CANCELLED: 'default',
};

const STATUS_LABELS: Record<string, string> = {
  CREATED: '已创建',
  SUBMITTED: '已提交',
  RUNNING: '运行中',
  WAITING_QC: '等待 QC',
  CALCULATING: '计算中',
  COMPLETED: '已完成',
  FAILED: '失败',
  CANCELLED: '已取消',
};

export default function ClusterAnalysisResultsPanel({ jobId, onBack }: Props) {
  const [job, setJob] = useState<AdvancedClusterJob | null>(null);
  const [results, setResults] = useState<ClusterAnalysisResults | null>(null);
  const [qcStatus, setQcStatus] = useState<QCStatus | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [expandedSections, setExpandedSections] = useState<Record<string, boolean>>({});

  const fetchData = useCallback(async () => {
    try {
      setLoading(true);
      const [jobData, resultsData, qcData] = await Promise.all([
        getClusterAnalysisJob(jobId),
        getClusterAnalysisResults(jobId),
        getClusterAnalysisQCStatus(jobId),
      ]);
      setJob(jobData);
      setResults(resultsData);
      setQcStatus(qcData);
      setError(null);
    } catch (err) {
      setError((err as Error).message || '获取数据失败');
    } finally {
      setLoading(false);
    }
  }, [jobId]);

  useEffect(() => {
    fetchData();
    // 如果任务还在进行中，定时刷新
    const interval = setInterval(() => {
      if (job && !['COMPLETED', 'FAILED', 'CANCELLED'].includes(job.status)) {
        fetchData();
      }
    }, 10000);
    return () => clearInterval(interval);
  }, [fetchData, job]);

  const toggleSection = (section: string) => {
    setExpandedSections((prev) => ({ ...prev, [section]: !prev[section] }));
  };

  if (loading && !job) {
    return (
      <Box sx={{ p: 3 }}>
        <LinearProgress />
        <Typography sx={{ mt: 2 }}>加载中...</Typography>
      </Box>
    );
  }

  if (error) {
    return (
      <Alert severity="error" sx={{ m: 2 }}>
        {error}
        <Button onClick={fetchData} sx={{ ml: 2 }}>重试</Button>
      </Alert>
    );
  }

  if (!job) return null;

  return (
    <Box sx={{ p: 2 }}>
      {/* 顶部信息 */}
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Box>
          <Typography variant="h5" gutterBottom>
            Cluster 高级计算结果 #{jobId}
          </Typography>
          <Box sx={{ display: 'flex', gap: 1, alignItems: 'center' }}>
            <Chip
              label={STATUS_LABELS[job.status] || job.status}
              color={STATUS_COLORS[job.status] || 'default'}
              size="small"
            />
            <Typography variant="body2" color="text.secondary">
              进度: {job.progress.toFixed(0)}%
            </Typography>
          </Box>
        </Box>
        <Box>
          <Tooltip title="刷新">
            <IconButton onClick={fetchData} disabled={loading}>
              <RefreshIcon />
            </IconButton>
          </Tooltip>
          {onBack && <Button onClick={onBack} sx={{ ml: 1 }}>返回</Button>}
        </Box>
      </Box>

      {/* 进度条 */}
      {job.status !== 'COMPLETED' && job.status !== 'FAILED' && (
        <LinearProgress variant="determinate" value={job.progress} sx={{ mb: 2, height: 8, borderRadius: 4 }} />
      )}

      {/* QC 任务状态 */}
      {qcStatus && qcStatus.total_qc_jobs > 0 && (
        <Card sx={{ mb: 2 }}>
          <CardContent>
            <Typography variant="h6" gutterBottom>QC 任务状态</Typography>
            <Grid container spacing={2}>
              <Grid item xs={3}>
                <Paper sx={{ p: 1, textAlign: 'center', bgcolor: 'success.light' }}>
                  <Typography variant="h4">{qcStatus.completed}</Typography>
                  <Typography variant="body2">已完成</Typography>
                </Paper>
              </Grid>
              <Grid item xs={3}>
                <Paper sx={{ p: 1, textAlign: 'center', bgcolor: 'primary.light' }}>
                  <Typography variant="h4">{qcStatus.running}</Typography>
                  <Typography variant="body2">运行中</Typography>
                </Paper>
              </Grid>
              <Grid item xs={3}>
                <Paper sx={{ p: 1, textAlign: 'center', bgcolor: 'warning.light' }}>
                  <Typography variant="h4">{qcStatus.pending}</Typography>
                  <Typography variant="body2">等待中</Typography>
                </Paper>
              </Grid>
              <Grid item xs={3}>
                <Paper sx={{ p: 1, textAlign: 'center', bgcolor: qcStatus.failed > 0 ? 'error.light' : 'grey.200' }}>
                  <Typography variant="h4">{qcStatus.failed}</Typography>
                  <Typography variant="body2">失败</Typography>
                </Paper>
              </Grid>
            </Grid>
          </CardContent>
        </Card>
      )}

      {/* 计算类型和结果 */}
      {job.calc_types.map((calcType) => {
        const info = CALC_TYPE_INFO[calcType as ClusterCalcType];
        const calcResult = results?.results?.[calcType] as Record<string, unknown>;
        const isExpanded = expandedSections[calcType] ?? true;

        return (
          <Card key={calcType} sx={{ mb: 2 }}>
            <CardContent>
              <Box
                sx={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}
                onClick={() => toggleSection(calcType)}
              >
                <Typography variant="h6" sx={{ flexGrow: 1 }}>
                  {info?.icon} {info?.label || calcType}
                </Typography>
                {calcResult?.error ? (
                  <Chip label="失败" color="error" size="small" icon={<ErrorIcon />} />
                ) : calcResult && Object.keys(calcResult).length > 0 ? (
                  <Chip label="完成" color="success" size="small" icon={<CheckCircleIcon />} />
                ) : (
                  <Chip label="等待" color="default" size="small" icon={<PendingIcon />} />
                )}
                <IconButton size="small">
                  {isExpanded ? <ExpandLessIcon /> : <ExpandMoreIcon />}
                </IconButton>
              </Box>

              <Collapse in={isExpanded}>
                <Typography variant="body2" color="text.secondary" sx={{ my: 1 }}>
                  {info?.description}
                </Typography>
                <Typography variant="caption" color="text.secondary" sx={{ fontFamily: 'monospace' }}>
                  {info?.formula}
                </Typography>

                {/* 结果展示 */}
                {calcResult && !calcResult.error && (
                  <Box sx={{ mt: 2 }}>
                    {renderCalcTypeResult(calcType as ClusterCalcType, calcResult)}
                  </Box>
                )}

                {calcResult?.error && (
                  <Alert severity="error" sx={{ mt: 2 }}>
                    {calcResult.error as string}
                  </Alert>
                )}
              </Collapse>
            </CardContent>
          </Card>
        );
      })}

      {/* 错误信息 */}
      {job.error_message && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {job.error_message}
        </Alert>
      )}
    </Box>
  );
}

// ============================================================================
// 结果渲染函数
// ============================================================================

function renderCalcTypeResult(calcType: ClusterCalcType, result: Record<string, unknown>): React.ReactNode {
  switch (calcType) {
    case 'BINDING_TOTAL':
    case 'DESOLVATION_FULL':
      return renderBindingTotalResult(result);
    case 'BINDING_PAIRWISE':
      return renderBindingPairwiseResult(result);
    case 'DESOLVATION_STEPWISE':
      return renderDesolvationStepwiseResult(result);
    case 'REDOX':
      return renderRedoxResult(result);
    case 'REORGANIZATION':
      return renderReorganizationResult(result);
    default:
      return <pre>{JSON.stringify(result, null, 2)}</pre>;
  }
}

function renderBindingTotalResult(result: Record<string, unknown>) {
  return (
    <Box>
      <Grid container spacing={2}>
        <Grid item xs={4}>
          <Paper sx={{ p: 2, textAlign: 'center' }}>
            <Typography variant="h5" color="primary">
              {((result.e_bind_kcal_mol as number) || 0).toFixed(2)}
            </Typography>
            <Typography variant="body2">kcal/mol</Typography>
          </Paper>
        </Grid>
        <Grid item xs={4}>
          <Paper sx={{ p: 2, textAlign: 'center' }}>
            <Typography variant="h5" color="secondary">
              {((result.e_bind_ev as number) || 0).toFixed(4)}
            </Typography>
            <Typography variant="body2">eV</Typography>
          </Paper>
        </Grid>
        <Grid item xs={4}>
          <Paper sx={{ p: 2, textAlign: 'center' }}>
            <Typography variant="h5">
              {((result.e_bind_au as number) || 0).toFixed(6)}
            </Typography>
            <Typography variant="body2">Hartree</Typography>
          </Paper>
        </Grid>
      </Grid>
      <Typography variant="caption" color="text.secondary" sx={{ mt: 1, display: 'block' }}>
        E(cluster) = {((result.e_cluster_au as number) || 0).toFixed(6)} a.u.,
        E(ion) = {((result.e_ion_au as number) || 0).toFixed(6)} a.u.
      </Typography>
    </Box>
  );
}

function renderBindingPairwiseResult(result: Record<string, unknown>) {
  const bindings = (result.pairwise_bindings as Array<Record<string, unknown>>) || [];
  return (
    <Table size="small">
      <TableHead>
        <TableRow>
          <TableCell>配体</TableCell>
          <TableCell align="right">E(dimer) / a.u.</TableCell>
          <TableCell align="right">E(ligand) / a.u.</TableCell>
          <TableCell align="right">E_bind / kcal·mol⁻¹</TableCell>
          <TableCell align="right">E_bind / eV</TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {bindings.map((b, i) => (
          <TableRow key={i}>
            <TableCell>{b.ligand as string}</TableCell>
            <TableCell align="right">{((b.e_dimer_au as number) || 0).toFixed(6)}</TableCell>
            <TableCell align="right">{((b.e_ligand_au as number) || 0).toFixed(6)}</TableCell>
            <TableCell align="right" sx={{ fontWeight: 'bold' }}>
              {((b.e_bind_kcal_mol as number) || 0).toFixed(2)}
            </TableCell>
            <TableCell align="right">{((b.e_bind_ev as number) || 0).toFixed(4)}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}

function renderDesolvationStepwiseResult(result: Record<string, unknown>) {
  const steps = (result.stepwise_desolvation as Array<Record<string, unknown>>) || [];
  return (
    <Table size="small">
      <TableHead>
        <TableRow>
          <TableCell>移除配体</TableCell>
          <TableCell align="right">E(cluster) / a.u.</TableCell>
          <TableCell align="right">E(minus) / a.u.</TableCell>
          <TableCell align="right">ΔE / kcal·mol⁻¹</TableCell>
          <TableCell align="right">ΔE / eV</TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {steps.map((s, i) => (
          <TableRow key={i}>
            <TableCell>{s.ligand as string}</TableCell>
            <TableCell align="right">{((s.e_cluster_au as number) || 0).toFixed(6)}</TableCell>
            <TableCell align="right">{((s.e_minus_au as number) || 0).toFixed(6)}</TableCell>
            <TableCell align="right" sx={{ fontWeight: 'bold', color: (s.delta_e_kcal_mol as number) > 0 ? 'error.main' : 'success.main' }}>
              {((s.delta_e_kcal_mol as number) || 0).toFixed(2)}
            </TableCell>
            <TableCell align="right">{((s.delta_e_ev as number) || 0).toFixed(4)}</TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}

function renderRedoxResult(result: Record<string, unknown>) {
  const potentials = (result.redox_potentials as Array<Record<string, unknown>>) || [];
  return (
    <Table size="small">
      <TableHead>
        <TableRow>
          <TableCell>物种</TableCell>
          <TableCell align="right">E(neutral,gas) / a.u.</TableCell>
          <TableCell align="right">E(charged,sol) / a.u.</TableCell>
          <TableCell align="right">ΔG(sol) / eV</TableCell>
          <TableCell align="right">E° (vs SHE) / V</TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {potentials.map((p, i) => (
          <TableRow key={i}>
            <TableCell sx={{ fontFamily: 'monospace' }}>{p.smiles as string}</TableCell>
            <TableCell align="right">{((p.e_neutral_gas_au as number) || 0).toFixed(6)}</TableCell>
            <TableCell align="right">{((p.e_charged_sol_au as number) || 0).toFixed(6)}</TableCell>
            <TableCell align="right">{((p.delta_g_sol_ev as number) || 0).toFixed(4)}</TableCell>
            <TableCell align="right" sx={{ fontWeight: 'bold', color: 'primary.main' }}>
              {((p.oxidation_potential_v as number) || 0).toFixed(3)}
            </TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}

function renderReorganizationResult(result: Record<string, unknown>) {
  if (result.status === 'not_implemented') {
    return (
      <Alert severity="info">
        {result.message as string}
      </Alert>
    );
  }
  return <pre>{JSON.stringify(result, null, 2)}</pre>;
}

