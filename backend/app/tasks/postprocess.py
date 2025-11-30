"""
后处理任务模块

当 MD 任务完成时，自动执行后处理分析：
- RDF（径向分布函数）计算
- MSD（均方位移）计算
- 结果存储到数据库
"""

import logging
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime

from celery import Task
from sqlalchemy.orm import Session

from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.result import RDFResult
from app.core.config import settings

logger = logging.getLogger(__name__)


class DatabaseTask(Task):
    """自动管理数据库会话的任务基类"""
    _db: Session = None

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()

    @property
    def db(self) -> Session:
        if self._db is None:
            self._db = SessionLocal()
        return self._db


@celery_app.task(
    bind=True,
    base=DatabaseTask,
    name="app.tasks.postprocess.postprocess_md_job_task",
    max_retries=3,
    default_retry_delay=300,  # 5分钟后重试
)
def postprocess_md_job_task(self, job_id: int) -> Dict[str, Any]:
    """
    MD 任务完成后的自动后处理
    
    Args:
        job_id: MD 任务 ID
        
    Returns:
        Dict 包含后处理结果
    """
    db = self.db
    
    try:
        logger.info(f"[Task {self.request.id}] Starting postprocessing for job {job_id}")
        
        # 1. 获取任务
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        
        if not job:
            error_msg = f"Job {job_id} not found"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            return {"success": False, "error": error_msg}
        
        # 2. 检查任务状态（允许 COMPLETED 或 POSTPROCESSING）
        if job.status not in [JobStatus.COMPLETED, JobStatus.POSTPROCESSING]:
            error_msg = f"Job {job_id} is not completed (status: {job.status})"
            logger.warning(f"[Task {self.request.id}] {error_msg}")
            return {"success": False, "error": error_msg}
        
        # 3. 更新状态为 POSTPROCESSING
        job.status = JobStatus.POSTPROCESSING
        db.commit()
        logger.info(f"[Task {self.request.id}] Job {job_id} status updated to POSTPROCESSING")
        
        # 4. 获取工作目录
        # 优先使用 job.work_dir，如果没有则从 config 中获取
        if job.work_dir:
            work_dir = Path(job.work_dir)
        elif job.config and job.config.get("work_dir"):
            work_dir = Path(job.config.get("work_dir"))
        else:
            # 如果都没有，使用任务名称构建路径
            job_name = job.config.get("job_name") if job.config else f"job_{job_id}"
            work_dir = settings.MOLYTE_WORK_BASE_PATH / job_name
        
        if not work_dir.exists():
            error_msg = f"Work directory not found: {work_dir}"
            logger.error(f"[Task {self.request.id}] {error_msg}")
            job.status = JobStatus.COMPLETED
            job.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}
        
        # 5. 获取配方信息（用于智能 RDF 生成）
        composition = None
        if job.system:
            composition = {
                'cations': job.system.cations,
                'anions': job.system.anions,
                'solvents': job.system.solvents,
            }
            logger.info(f"[Task {self.request.id}] Using composition from system: {composition}")

        # 6. 执行 RDF 计算
        rdf_results = []
        try:
            rdf_results = calculate_default_rdfs(work_dir, job_id, composition)
            logger.info(f"[Task {self.request.id}] Calculated {len(rdf_results)} RDF pairs")

            # 保存 RDF 结果到数据库
            save_rdf_results_to_db(db, job_id, rdf_results)
            logger.info(f"[Task {self.request.id}] Saved {len(rdf_results)} RDF results to database")

            # 绘制 RDF 图表
            try:
                plot_rdf_figures(work_dir, rdf_results, job_id)
                logger.info(f"[Task {self.request.id}] Generated RDF plots")
            except Exception as e:
                logger.error(f"[Task {self.request.id}] RDF plotting failed: {e}")
                import traceback
                logger.error(traceback.format_exc())

        except Exception as e:
            logger.error(f"[Task {self.request.id}] RDF calculation failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            # RDF 失败不影响整体流程，继续执行

        # 6. 执行 MSD 计算
        msd_results = []
        try:
            from app.tasks.msd_processor import process_msd_data
            msd_results = process_msd_data(db, job_id, work_dir)
            logger.info(f"[Task {self.request.id}] Processed {len(msd_results)} MSD results")
        except Exception as e:
            logger.error(f"[Task {self.request.id}] MSD processing failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            # MSD 失败不影响整体流程，继续执行

        # 7. 更新任务配置，记录后处理结果
        if "postprocess" not in job.config:
            job.config["postprocess"] = {}

        job.config["postprocess"]["rdf_count"] = len(rdf_results)
        job.config["postprocess"]["msd_count"] = len(msd_results)
        job.config["postprocess"]["completed_at"] = datetime.now().isoformat()
        job.config["postprocess"]["rdf_pairs"] = [
            f"{r['center_label']}-{r['target_label']}" for r in rdf_results
        ]
        job.config["postprocess"]["msd_species"] = [
            r.species for r in msd_results
        ]
        
        # 8. 恢复状态为 COMPLETED
        job.status = JobStatus.COMPLETED
        db.commit()

        logger.info(f"[Task {self.request.id}] Postprocessing completed for job {job_id}")

        return {
            "success": True,
            "job_id": job_id,
            "rdf_count": len(rdf_results),
            "msd_count": len(msd_results),
        }
        
    except Exception as exc:
        logger.exception(f"[Task {self.request.id}] Exception during postprocessing: {exc}")
        
        # 恢复状态为 COMPLETED（即使后处理失败，MD 任务本身是成功的）
        try:
            job = db.query(MDJob).filter(MDJob.id == job_id).first()
            if job:
                job.status = JobStatus.COMPLETED
                job.error_message = f"Postprocessing failed: {str(exc)}"
                db.commit()
        except:
            pass
        
        raise self.retry(exc=exc)


def calculate_default_rdfs(work_dir: Path, job_id: int, composition: Dict = None) -> List[Dict[str, Any]]:
    """
    读取 LAMMPS 输出的 RDF 数据

    直接从 LAMMPS 的 out_rdf.dat 文件读取 RDF 数据，
    避免重新计算，提高效率和准确性。

    Args:
        work_dir: 工作目录
        job_id: 任务 ID
        composition: 配方数据（未使用，保留用于兼容性）

    Returns:
        List of RDF results
    """
    from app.workers.lammps_rdf_reader import LAMMPSRDFReader

    results = []

    try:
        # 使用 LAMMPS RDF 读取器
        reader = LAMMPSRDFReader(work_dir)
        rdf_data = reader.read_rdf_file()

        if not rdf_data:
            logger.warning(f"No RDF data found for job {job_id}")
            return results

        # 转换为标准格式
        for data in rdf_data:
            result = {
                'center_label': data['center_label'],
                'target_label': data['target_label'],
                'description': f"{data['center_label']} → {data['target_label']}",
                'r': data['r'],
                'g_r': data['g_r'],
                'coordination_number': data['coordination_number']
            }
            results.append(result)

        logger.info(f"Read {len(results)} RDF pairs from LAMMPS output for job {job_id}")

    except Exception as e:
        logger.error(f"Failed to read RDF data: {e}")
        import traceback
        traceback.print_exc()

    return results


def plot_rdf_figures(work_dir: Path, rdf_results: List[Dict[str, Any]], job_id: int):
    """
    绘制 RDF 图表

    Args:
        work_dir: 工作目录
        rdf_results: RDF 结果列表
        job_id: 任务 ID
    """
    from app.workers.rdf_plotter import RDFPlotter

    if not rdf_results:
        logger.warning(f"No RDF results to plot for job {job_id}")
        return

    plotter = RDFPlotter()

    # 创建图表输出目录
    plots_dir = work_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    # 1. 绘制组合图（所有 RDF 在一张图上）
    combined_file = plots_dir / "rdf_combined.png"
    plotter.plot_rdf_combined(rdf_results, combined_file, title=f"RDF Analysis - Job {job_id}")

    # 2. 绘制分类图（典型 RDF 和其他 RDF 分开）
    categorized_file = plots_dir / "rdf_categorized.png"
    plotter.plot_rdf_categorized(rdf_results, categorized_file, title=f"RDF Analysis - Job {job_id}")

    # 3. 绘制单独的图（每个 RDF 对一张图）
    individual_dir = plots_dir / "individual"
    plotter.plot_rdf_individual(rdf_results, individual_dir, prefix="rdf")

    logger.info(f"Generated RDF plots in {plots_dir}")


def save_rdf_results_to_db(db: Session, job_id: int, rdf_results: List[Dict[str, Any]]):
    """
    保存 RDF 结果到数据库

    Args:
        db: 数据库会话
        job_id: 任务 ID
        rdf_results: RDF 计算结果列表
    """
    import numpy as np
    from collections import defaultdict

    # 删除旧的 RDF 结果
    db.query(RDFResult).filter(RDFResult.md_job_id == job_id).delete()

    # 按 (center_label, target_label) 分组合并
    grouped_results = defaultdict(list)
    for result in rdf_results:
        key = (result['center_label'], result['target_label'])
        grouped_results[key].append(result)

    # 合并相同的 RDF 对
    merged_results = []
    for (center_label, target_label), group in grouped_results.items():
        if len(group) == 1:
            merged_results.append(group[0])
        else:
            # 多个相同对，合并 g(r) 和 CN
            r_values = group[0]['r']  # r 值相同
            g_r_arrays = [np.array(g['g_r']) for g in group]
            g_r_merged = np.mean(g_r_arrays, axis=0).tolist()

            cn_arrays = [np.array(g['coordination_number']) for g in group if g.get('coordination_number') is not None]
            cn_merged = np.mean(cn_arrays, axis=0).tolist() if cn_arrays else None

            merged_results.append({
                'center_label': center_label,
                'target_label': target_label,
                'r': r_values,
                'g_r': g_r_merged,
                'coordination_number': cn_merged,
            })
            logger.info(f"Merged {len(group)} RDF pairs for {center_label} -> {target_label}")

    for result in merged_results:
        # 分析 RDF 数据
        r_values = np.array(result['r'])
        g_r_values = np.array(result['g_r'])

        # 找到第一个峰的位置和高度
        first_peak_position = None
        first_peak_height = None
        coordination_number = None

        try:
            # 找到第一个峰（g(r) > 1.0 的第一个最大值）
            # 跳过 r < 1.5 Å 的部分（对于 Li-O 等配位，第一峰通常在 1.8-2.0 Å）
            valid_mask = r_values >= 1.5
            valid_r = r_values[valid_mask]
            valid_g_r = g_r_values[valid_mask]

            peaks_mask = valid_g_r > 1.0
            if peaks_mask.any():
                first_peak_idx = np.argmax(valid_g_r[peaks_mask])
                actual_idx = np.where(peaks_mask)[0][first_peak_idx]
                first_peak_position = float(valid_r[actual_idx])
                first_peak_height = float(valid_g_r[actual_idx])

                # 计算配位数（积分到第一个最小值）
                # 配位数 CN = 4πρ ∫[0 to r_min] r² g(r) dr
                # 其中 ρ = N_target / V 是目标原子的数密度
                if actual_idx + 1 < len(valid_g_r):
                    # 找到峰后的第一个最小值
                    after_peak = valid_g_r[actual_idx:]
                    min_after_peak_idx = np.argmin(after_peak[:min(20, len(after_peak))])
                    cutoff_idx_in_valid = actual_idx + min_after_peak_idx

                    # 获取必要的参数
                    n_target = result.get('target_atom_count', 1)
                    box_volume = result.get('box_volume', None)

                    if box_volume and box_volume > 0:
                        # 计算目标原子的数密度（atoms/Å³）
                        rho = n_target / box_volume

                        # 计算积分：∫r²g(r)dr（从 r=0 到第一个最小值）
                        # 需要使用完整的 r 和 g_r 数组
                        dr = r_values[1] - r_values[0] if len(r_values) > 1 else 0.05

                        # 找到 cutoff_r 在完整数组中的位置
                        cutoff_r = valid_r[cutoff_idx_in_valid]
                        cutoff_idx_full = np.searchsorted(r_values, cutoff_r)

                        # 积分到 cutoff
                        integral = np.sum(r_values[:cutoff_idx_full]**2 * g_r_values[:cutoff_idx_full] * dr)

                        # 配位数 = 4πρ * integral
                        coordination_number = float(4 * np.pi * rho * integral)

                        logger.info(f"CN calculation: rho={rho:.6f} atoms/Å³, cutoff_r={cutoff_r:.3f} Å, integral={integral:.3f}, CN={coordination_number:.2f}")
                    else:
                        logger.warning(f"Box volume not available, cannot calculate coordination number accurately")
        except Exception as e:
            logger.warning(f"Failed to analyze RDF peak: {e}")

        # 获取配位数数组
        coordination_number_values = result.get('coordination_number', None)

        # 如果配位数数组存在，获取最终配位数
        final_coordination_number = coordination_number
        if coordination_number_values is not None and len(coordination_number_values) > 0:
            final_coordination_number = float(coordination_number_values[-1])

        # 创建数据库记录
        rdf_record = RDFResult(
            md_job_id=job_id,
            center_species=result['center_label'],
            shell_species=result['target_label'],
            r_values=result['r'],
            g_r_values=result['g_r'],
            coordination_number_values=coordination_number_values,  # 保存完整数组
            first_peak_position=first_peak_position,
            first_peak_height=first_peak_height,
            coordination_number=final_coordination_number,  # 保存最终值
        )

        db.add(rdf_record)

    db.commit()
    logger.info(f"Saved {len(rdf_results)} RDF results to database for job {job_id}")

