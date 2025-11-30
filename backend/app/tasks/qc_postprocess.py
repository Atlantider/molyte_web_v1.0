"""
QC 后处理 Celery Worker

负责处理 QC 计算完成后的结果提取和ESP可视化
"""

import logging
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

from celery import Task
from sqlalchemy.orm import Session

from app.celery_app import celery_app
from app.database import SessionLocal
from app.models.qc import QCJob, QCResult, MoleculeQCCache, QCJobStatus

logger = logging.getLogger(__name__)

# Hartree to eV conversion factor
HARTREE_TO_EV = 27.2114


class DatabaseTask(Task):
    """自动管理数据库会话的任务基类"""
    _db: Session = None

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()
            self._db = None

    @property
    def db(self) -> Session:
        if self._db is None:
            self._db = SessionLocal()
        return self._db


def extract_gaussian_results(log_file: str) -> Dict[str, Any]:
    """
    从Gaussian输出文件提取计算结果
    
    Args:
        log_file: Gaussian输出文件路径
        
    Returns:
        Dict containing energy, HOMO, LUMO values
    """
    results = {
        "energy_au": None,
        "homo": None,
        "lumo": None,
    }
    
    if not os.path.exists(log_file):
        logger.error(f"Log file not found: {log_file}")
        return results
    
    # 正则表达式模式
    energy_pattern = re.compile(r"SCF Done:\s+E\([A-Z0-9]+\)\s+=\s+([-\d.]+)")
    alpha_occ_pattern = re.compile(r"Alpha\s+occ\.\seigenvalues\s+--\s+([-\d.]+(?:\s+[-\d.]+)*)")
    alpha_virt_pattern = re.compile(r"Alpha\s+virt\.\seigenvalues\s+--\s+([-\d.]+(?:\s+[-\d.]+)*)")
    
    last_energy = None
    last_homo = None
    last_lumo = None
    
    # 使用errors='replace'来处理可能的编码问题（Gaussian输出可能包含非UTF-8字符）
    with open(log_file, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            # 匹配SCF能量
            match_energy = energy_pattern.search(lines[i])
            if match_energy:
                last_energy = float(match_energy.group(1))
            
            # 匹配HOMO和LUMO
            match_occ = alpha_occ_pattern.search(lines[i])
            match_virt = alpha_virt_pattern.search(lines[i + 1]) if i + 1 < len(lines) else None
            
            if match_occ and match_virt:
                occ_values = match_occ.group(1).split()
                virt_values = match_virt.group(1).split()
                
                if occ_values:
                    last_homo = float(occ_values[-1])
                if virt_values:
                    last_lumo = float(virt_values[0])
    
    results["energy_au"] = last_energy
    results["homo"] = last_homo
    results["lumo"] = last_lumo
    
    return results


def extract_esp_values(surfanalysis_file: str) -> Tuple[Optional[float], Optional[float]]:
    """
    从Multiwfn的surfanalysis.txt提取ESP极值
    
    Returns:
        Tuple of (ESP_min, ESP_max) in kcal/mol
    """
    if not os.path.exists(surfanalysis_file):
        logger.warning(f"Surface analysis file not found: {surfanalysis_file}")
        return None, None
    
    min_pattern = r"\*\s*(\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)"
    max_pattern = r"\s*(\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)"
    
    with open(surfanalysis_file, 'r') as f:
        data = f.read()
    
    min_matches = re.findall(min_pattern, data)
    max_matches = re.findall(max_pattern, data)
    
    esp_min = None
    esp_max = None
    
    if min_matches:
        min_values = [float(match[3]) for match in min_matches]
        esp_min = min(min_values)
    
    if max_matches:
        max_values = [float(match[3]) for match in max_matches]
        esp_max = max(max_values)
    
    return esp_min, esp_max


def generate_esp_visualization(work_dir: Path, molecule_name: str, fchk_file: str) -> Optional[str]:
    """
    生成ESP可视化图像
    
    Args:
        work_dir: 工作目录
        molecule_name: 分子名称
        fchk_file: fchk文件路径
        
    Returns:
        ESP图像路径，失败返回None
    """
    try:
        # 创建Multiwfn输入文件
        txt_content = '''12
0
1
-1
-1
5
1
3
2
0
5
12
1
2'''
        txt_path = work_dir / "ESPiso.txt"
        with open(txt_path, 'w') as f:
            f.write(txt_content)
        
        # 创建VMD脚本 - 使用Tachyon渲染器用于无头模式
        # 自动计算分子大小并调整缩放比例
        vmd_content = '''# ESP可视化脚本 - 使用Tachyon渲染器（自适应缩放）
color scale method BWR
color Display Background white
axes location Off
display depthcue off
display projection Orthographic
display nearclip set 0.01
light 2 on
light 3 on

# 由于无头模式不支持EdgyGlass材质，使用Transparent
material change opacity Transparent 0.70

set nsystem 1
set colorlow -0.03
set colorhigh 0.03

for {set i 1} {$i<=$nsystem} {incr i} {
set id [expr $i-1]
mol new density$i.cub
mol addfile ESP$i.cub
mol modstyle 0 $id CPK 1.000000 0.300000 22.000000 22.000000
mol addrep $id
mol modstyle 1 $id Isosurface 0.001000 0 0 0 1 1
mol modmaterial 1 $id Transparent
mol modcolor 1 $id Volume 1
mol scaleminmax $id 1 $colorlow $colorhigh
}

# 自适应视角调整 - 根据分子大小自动缩放
display resetview
mol top 0

# 获取分子的边界框
set sel [atomselect top all]
set minmax [measure minmax $sel]
set min_coords [lindex $minmax 0]
set max_coords [lindex $minmax 1]
$sel delete

# 计算分子在三个方向的尺寸
set dx [expr {[lindex $max_coords 0] - [lindex $min_coords 0]}]
set dy [expr {[lindex $max_coords 1] - [lindex $min_coords 1]}]
set dz [expr {[lindex $max_coords 2] - [lindex $min_coords 2]}]

# 获取最大维度（考虑等密度面会比原子位置大约3-4埃）
set iso_padding 6.0
set max_dim [expr {max($dx, max($dy, $dz)) + $iso_padding}]

# 计算自适应缩放比例
# 默认VMD视野大约能显示15埃的分子，我们希望分子占据画面的50%以确保完整显示
set target_view 15.0
set fill_ratio 0.50
set auto_scale [expr {($target_view * $fill_ratio) / $max_dim}]

# 限制缩放范围，避免极端情况
if {$auto_scale > 1.0} {
    set auto_scale 1.0
}
if {$auto_scale < 0.3} {
    set auto_scale 0.3
}

puts "Molecule dimensions: $dx x $dy x $dz Angstrom"
puts "Max dimension with padding: $max_dim Angstrom"
puts "Auto scale factor: $auto_scale"

# 旋转到较好的观察角度
rotate x by 15
rotate y by 25
rotate z by 5

# 应用自适应缩放
scale by $auto_scale

# 使用Tachyon渲染器生成图片（支持无头模式）
render Tachyon ESP.dat
exit'''
        vmd_path = work_dir / "ESPiso.vmd"
        with open(vmd_path, 'w') as f:
            f.write(vmd_content)

        # 运行Multiwfn生成cube文件
        fchk_basename = os.path.basename(fchk_file)
        esp_script = f"""#!/bin/bash
cd {work_dir}
Multiwfn {fchk_basename} -ESPrhoiso 0.001 < ESPiso.txt
mv -f density.cub density1.cub
mv -f totesp.cub ESP1.cub
"""
        script_path = work_dir / "run_esp.sh"
        with open(script_path, 'w') as f:
            f.write(esp_script)
        os.chmod(script_path, 0o755)

        # 执行Multiwfn
        result = subprocess.run(
            ["bash", str(script_path)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=600
        )

        if result.returncode != 0:
            logger.warning(f"Multiwfn failed: {result.stderr}")
            return None

        # 运行VMD生成Tachyon数据文件
        result = subprocess.run(
            ["vmd", "-dispdev", "text", "-e", "ESPiso.vmd"],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=120
        )

        # 使用Tachyon渲染器将.dat转换为.tga
        tachyon_dat = work_dir / "ESP.dat"
        tachyon_tga = work_dir / "ESP.tga"
        esp_image = work_dir / "ESP.png"

        if tachyon_dat.exists():
            # 查找tachyon可执行文件路径
            tachyon_cmd = "tachyon_LINUXAMD64"  # VMD自带的tachyon
            tachyon_path = "/usr/local/lib/vmd/tachyon_LINUXAMD64"
            if os.path.exists(tachyon_path):
                tachyon_cmd = tachyon_path

            # 运行Tachyon渲染
            result = subprocess.run(
                [tachyon_cmd, "-aasamples", "12", str(tachyon_dat),
                 "-format", "TGA", "-res", "1024", "768", "-o", str(tachyon_tga)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=120
            )

            # 使用ImageMagick将TGA转换为PNG
            if tachyon_tga.exists():
                try:
                    result = subprocess.run(
                        ["convert", str(tachyon_tga), str(esp_image)],
                        cwd=str(work_dir),
                        capture_output=True,
                        text=True,
                        timeout=60
                    )
                except Exception as e:
                    logger.warning(f"ImageMagick convert failed: {e}, trying Python PIL")
                    # 尝试用PIL转换
                    try:
                        from PIL import Image as PILImage
                        img = PILImage.open(str(tachyon_tga))
                        img.save(str(esp_image))
                    except Exception as e2:
                        logger.warning(f"PIL convert failed: {e2}")

        if esp_image.exists() and esp_image.stat().st_size > 0:
            return str(esp_image)

        return None

    except Exception as e:
        logger.error(f"ESP visualization failed: {e}")
        return None


def generate_orbital_visualization(work_dir: Path, fchk_file: str, orbital_type: str) -> Optional[str]:
    """
    生成HOMO或LUMO轨道可视化图像

    Args:
        work_dir: 工作目录
        fchk_file: fchk文件路径
        orbital_type: "HOMO" 或 "LUMO"

    Returns:
        图像路径，失败返回None
    """
    try:
        orbital_lower = orbital_type.lower()
        cube_file = work_dir / f"{orbital_lower}.cub"
        image_file = work_dir / f"{orbital_type}.png"

        # 先检查是否已有cube文件
        if not cube_file.exists():
            # 使用Multiwfn生成轨道cube文件
            # 5: 在空间区域输出和绘制特定属性
            # 4: 指定单个轨道
            # h: HOMO, l: LUMO
            # 2: 中等质量网格
            # 2: 导出为cube文件
            orbital_key = "h" if orbital_type == "HOMO" else "l"

            txt_content = f'''5
4
{orbital_key}
2
2
'''
            txt_path = work_dir / f"{orbital_lower}_gen.txt"
            with open(txt_path, 'w') as f:
                f.write(txt_content)

            # 运行Multiwfn生成cube
            multiwfn_cmd = "/home/iei/share/software/Multiwfn_3.8_dev_bin_Linux/Multiwfn"
            if os.path.exists(multiwfn_cmd):
                result = subprocess.run(
                    [multiwfn_cmd, fchk_file],
                    stdin=open(txt_path, 'r'),
                    cwd=str(work_dir),
                    capture_output=True,
                    text=True,
                    timeout=300
                )

                # Multiwfn输出为 MOvalue.cub，需要重命名
                mo_cube = work_dir / "MOvalue.cub"
                if mo_cube.exists():
                    mo_cube.rename(cube_file)
                    logger.info(f"Generated {orbital_type} cube: {cube_file}")

        if not cube_file.exists():
            logger.warning(f"{orbital_type} cube file not found")
            return None

        # 使用VMD渲染轨道图像
        # 轨道用两种颜色显示正负相位
        vmd_script = f'''# {orbital_type}轨道可视化脚本 - 使用Tachyon渲染器（自适应缩放）
color Display Background white
axes location Off
display depthcue off
display projection Orthographic
display nearclip set 0.01
light 2 on
light 3 on

# 加载分子和轨道
mol new {cube_file.name}
mol modstyle 0 0 CPK 0.800000 0.300000 22.000000 22.000000

# 添加正相位等值面（蓝色）
mol addrep 0
mol modstyle 1 0 Isosurface 0.02 0 0 0 1 1
mol modcolor 1 0 ColorID 0
mol modmaterial 1 0 Transparent
material change opacity Transparent 0.75

# 添加负相位等值面（红色）
mol addrep 0
mol modstyle 2 0 Isosurface -0.02 0 0 0 1 1
mol modcolor 2 0 ColorID 1
mol modmaterial 2 0 Transparent

# 自适应视角调整
display resetview
mol top 0

# 获取分子的边界框
set sel [atomselect top all]
set minmax [measure minmax $sel]
set min_coords [lindex $minmax 0]
set max_coords [lindex $minmax 1]
$sel delete

# 计算分子在三个方向的尺寸
set dx [expr {{[lindex $max_coords 0] - [lindex $min_coords 0]}}]
set dy [expr {{[lindex $max_coords 1] - [lindex $min_coords 1]}}]
set dz [expr {{[lindex $max_coords 2] - [lindex $min_coords 2]}}]

# 获取最大维度（考虑轨道云会比原子位置大约5-6埃）
set orb_padding 8.0
set max_dim [expr {{max($dx, max($dy, $dz)) + $orb_padding}}]

# 计算自适应缩放比例
set target_view 15.0
set fill_ratio 0.50
set auto_scale [expr {{($target_view * $fill_ratio) / $max_dim}}]

if {{$auto_scale > 1.0}} {{
    set auto_scale 1.0
}}
if {{$auto_scale < 0.3}} {{
    set auto_scale 0.3
}}

puts "{orbital_type} - Molecule dimensions: $dx x $dy x $dz Angstrom"
puts "{orbital_type} - Auto scale factor: $auto_scale"

# 旋转到较好的观察角度
rotate x by 15
rotate y by 25
rotate z by 5
scale by $auto_scale

# 使用Tachyon渲染器生成图片
render Tachyon {orbital_type}.dat
exit
'''

        vmd_path = work_dir / f"{orbital_lower}_render.vmd"
        with open(vmd_path, 'w') as f:
            f.write(vmd_script)

        # 运行VMD
        result = subprocess.run(
            ["vmd", "-dispdev", "text", "-e", str(vmd_path)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=120
        )

        # 使用Tachyon渲染
        tachyon_dat = work_dir / f"{orbital_type}.dat"
        tachyon_tga = work_dir / f"{orbital_type}.tga"

        if tachyon_dat.exists():
            tachyon_cmd = "/usr/local/lib/vmd/tachyon_LINUXAMD64"
            result = subprocess.run(
                [tachyon_cmd, "-aasamples", "12", str(tachyon_dat),
                 "-format", "TGA", "-res", "1024", "768", "-o", str(tachyon_tga)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=120
            )

            # 转换TGA到PNG
            if tachyon_tga.exists():
                try:
                    from PIL import Image as PILImage
                    img = PILImage.open(str(tachyon_tga))
                    img.save(str(image_file))
                except Exception as e:
                    logger.warning(f"PIL convert failed: {e}, trying ImageMagick")
                    subprocess.run(
                        ["convert", str(tachyon_tga), str(image_file)],
                        cwd=str(work_dir),
                        capture_output=True,
                        timeout=60
                    )

        if image_file.exists() and image_file.stat().st_size > 0:
            logger.info(f"Generated {orbital_type} image: {image_file}")
            return str(image_file)

        return None

    except Exception as e:
        logger.error(f"{orbital_type} visualization failed: {e}")
        return None


def update_molecule_cache(db: Session, job: QCJob, result: QCResult):
    """更新分子QC缓存"""
    cache = db.query(MoleculeQCCache).filter(
        MoleculeQCCache.smiles == job.smiles
    ).first()

    homo_ev = result.homo * HARTREE_TO_EV if result.homo else None
    lumo_ev = result.lumo * HARTREE_TO_EV if result.lumo else None
    gap_ev = (lumo_ev - homo_ev) if (homo_ev and lumo_ev) else None

    if cache:
        # 更新现有缓存
        cache.molecule_name = job.molecule_name
        cache.basis_set = job.basis_set
        cache.functional = job.functional
        cache.energy_au = result.energy_au
        cache.homo_ev = homo_ev
        cache.lumo_ev = lumo_ev
        cache.homo_lumo_gap_ev = gap_ev
        cache.esp_min_kcal = result.esp_min_kcal
        cache.esp_max_kcal = result.esp_max_kcal
        cache.esp_image_path = result.esp_image_path
        cache.preferred_qc_result_id = result.id
        cache.calculation_count += 1
        cache.updated_at = datetime.now()
    else:
        # 创建新缓存
        cache = MoleculeQCCache(
            smiles=job.smiles,
            molecule_name=job.molecule_name,
            basis_set=job.basis_set,
            functional=job.functional,
            energy_au=result.energy_au,
            homo_ev=homo_ev,
            lumo_ev=lumo_ev,
            homo_lumo_gap_ev=gap_ev,
            esp_min_kcal=result.esp_min_kcal,
            esp_max_kcal=result.esp_max_kcal,
            esp_image_path=result.esp_image_path,
            preferred_qc_result_id=result.id,
            calculation_count=1
        )
        db.add(cache)


@celery_app.task(
    bind=True,
    base=DatabaseTask,
    name="app.tasks.qc_postprocess.postprocess_qc_job",
    max_retries=2,
    default_retry_delay=120,
)
def postprocess_qc_job(self, job_id: int) -> Dict[str, Any]:
    """
    QC任务后处理

    提取计算结果、生成ESP可视化、更新缓存
    """
    db = self.db

    try:
        logger.info(f"[Task {self.request.id}] Starting postprocessing for QC job {job_id}")

        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            return {"success": False, "error": f"QC Job {job_id} not found"}

        work_dir = Path(job.work_dir)
        if not work_dir.exists():
            job.status = QCJobStatus.FAILED
            job.error_message = "Work directory not found"
            db.commit()
            return {"success": False, "error": "Work directory not found"}

        # 1. 查找输出文件 - 智能搜索
        config = job.config or {}
        safe_name = config.get("safe_name")

        # 尝试多种可能的文件名
        possible_names = []
        if safe_name:
            possible_names.append(safe_name)
        possible_names.append(job.molecule_name)

        log_file = None
        fchk_file = None
        actual_name = None

        for name in possible_names:
            test_log = work_dir / f"{name}_out.log"
            if test_log.exists():
                log_file = test_log
                fchk_file = work_dir / f"{name}.fchk"
                actual_name = name
                break

        # 如果还没找到，在目录中搜索 *_out.log 文件
        if not log_file:
            import glob
            log_files = list(work_dir.glob("*_out.log"))
            # 排除 qc_out.log（这是Slurm的输出文件，不是Gaussian的）
            log_files = [f for f in log_files if f.name not in ("qc_out.log", "qc_err.log")]
            if log_files:
                log_file = log_files[0]
                # 从文件名提取实际名称
                actual_name = log_file.stem.replace("_out", "")
                fchk_file = work_dir / f"{actual_name}.fchk"
                logger.info(f"Found output file by glob: {log_file}")

        if not log_file or not log_file.exists():
            job.status = QCJobStatus.FAILED
            job.error_message = "Gaussian output file not found"
            db.commit()
            return {"success": False, "error": "Output file not found"}

        # 2. 提取Gaussian结果
        gaussian_results = extract_gaussian_results(str(log_file))
        logger.info(f"Extracted Gaussian results: {gaussian_results}")

        # 3. 生成可视化图像（如果fchk存在）
        esp_image_path = None
        homo_image_path = None
        lumo_image_path = None
        esp_min, esp_max = None, None

        if fchk_file and fchk_file.exists():
            # 使用实际找到的文件名或原始分子名称
            esp_name = actual_name or job.molecule_name

            # 3.1 生成ESP可视化
            esp_image_path = generate_esp_visualization(work_dir, esp_name, str(fchk_file))
            logger.info(f"ESP visualization: {esp_image_path}")

            # 提取ESP值
            surfanalysis_file = work_dir / "surfanalysis.txt"
            esp_min, esp_max = extract_esp_values(str(surfanalysis_file))
            logger.info(f"ESP values: min={esp_min}, max={esp_max}")

            # 3.2 生成HOMO轨道图像
            homo_image_path = generate_orbital_visualization(work_dir, str(fchk_file), "HOMO")
            logger.info(f"HOMO visualization: {homo_image_path}")

            # 3.3 生成LUMO轨道图像
            lumo_image_path = generate_orbital_visualization(work_dir, str(fchk_file), "LUMO")
            logger.info(f"LUMO visualization: {lumo_image_path}")

        # 4. 计算HOMO-LUMO gap
        homo_lumo_gap = None
        if gaussian_results["homo"] and gaussian_results["lumo"]:
            homo_lumo_gap = (gaussian_results["lumo"] - gaussian_results["homo"]) * HARTREE_TO_EV

        # 5. 创建QC结果记录
        qc_result = QCResult(
            qc_job_id=job.id,
            smiles=job.smiles,
            energy_au=gaussian_results["energy_au"],
            homo=gaussian_results["homo"],
            lumo=gaussian_results["lumo"],
            homo_lumo_gap=homo_lumo_gap,
            esp_min_kcal=esp_min,
            esp_max_kcal=esp_max,
            esp_image_path=esp_image_path,
            homo_image_path=homo_image_path,
            lumo_image_path=lumo_image_path,
            fchk_file_path=str(fchk_file) if fchk_file.exists() else None,
            log_file_path=str(log_file),
        )
        db.add(qc_result)
        db.flush()  # 获取result.id

        # 6. 更新分子缓存
        update_molecule_cache(db, job, qc_result)

        # 7. 更新任务状态
        job.status = QCJobStatus.COMPLETED
        job.progress = 100.0
        job.finished_at = datetime.now()
        job.log_file = str(log_file)

        # 8. 自动设置为公开（QC数据必须公开）
        job.visibility = "PUBLIC"

        db.commit()

        logger.info(f"[Task {self.request.id}] QC job {job_id} postprocessing completed")

        return {
            "success": True,
            "job_id": job_id,
            "result_id": qc_result.id,
            "energy_au": gaussian_results["energy_au"],
            "homo": gaussian_results["homo"],
            "lumo": gaussian_results["lumo"],
            "esp_min": esp_min,
            "esp_max": esp_max,
        }

    except Exception as exc:
        logger.exception(f"[Task {self.request.id}] Postprocessing failed: {exc}")

        try:
            job = db.query(QCJob).filter(QCJob.id == job_id).first()
            if job:
                job.status = QCJobStatus.FAILED
                job.error_message = f"Postprocessing failed: {str(exc)}"
                db.commit()
        except Exception as db_exc:
            logger.error(f"Failed to update job status: {db_exc}")

        raise self.retry(exc=exc, countdown=120)

