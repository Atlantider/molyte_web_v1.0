"""
批量导入API端点
"""
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, Form
from fastapi.responses import StreamingResponse
from sqlalchemy.orm import Session
from typing import Optional, List
import pandas as pd
import io
import json
from datetime import datetime
from urllib.parse import quote

from app.database import get_db
from app.models.user import User
from app.models.project import Project
from app.models.electrolyte import ElectrolyteSystem
from app.models.qc import QCJob
from app.models.job import MDJob, JobStatus
from app.schemas.batch_import import (
    BatchImportRequest, BatchImportResult,
    BatchElectrolyteRow, BatchQCRow,
    TemplateDownloadRequest
)
from app.schemas.electrolyte import MoleculeSpec, ElectrolyteCreate
from app.schemas.qc import QCJobCreate
from app.dependencies import get_current_user
from app.api.v1.electrolytes import create_electrolyte
from app.api.v1.qc import create_qc_job
from app.api.v1.jobs import create_md_job
from app.utils.electrolyte_converter import strip_ion_charge
from app.utils.hash import calculate_system_hash

router = APIRouter()


@router.get("/template/download")
async def download_template(
    template_type: str = "electrolyte",  # electrolyte, qc, combined
    include_examples: bool = True,
    current_user: User = Depends(get_current_user)
):
    """
    下载批量导入模板
    template_type: electrolyte(配方), qc(QC计算), combined(组合)
    """
    if template_type == "electrolyte":
        # 配方模板（简化版，移除队列和资源配置）
        template_data = {
            "配方名称": ["LiPF6-EC-DMC", "NaPF6-EC-EMC"] if include_examples else [],
            "阳离子名称": ["Li+", "Na+"] if include_examples else [],
            "阳离子SMILES": ["[Li+]", "[Na+]"] if include_examples else [],
            "阳离子数量": [50, 50] if include_examples else [],
            "阴离子名称": ["PF6-", "PF6-"] if include_examples else [],
            "阴离子SMILES": ["F[P-](F)(F)(F)(F)F", "F[P-](F)(F)(F)(F)F"] if include_examples else [],
            "阴离子数量": [50, 50] if include_examples else [],
            "溶剂名称": ["EC;DMC", "EC;EMC"] if include_examples else [],
            "溶剂SMILES": ["C1COC(=O)O1;COC(=O)OC", "C1COC(=O)O1;CCOC(=O)OC"] if include_examples else [],
            "溶剂数量": ["100;100", "100;100"] if include_examples else [],
            "温度(K)": [298.15, 298.15] if include_examples else [],
            "压强(atm)": [1.0, 1.0] if include_examples else [],
            "盒子尺寸(Å)": [40.0, 45.0] if include_examples else [],
        }
        filename = "配方批量导入模板.xlsx"

    elif template_type == "qc":
        # QC计算模板（简化版，移除队列和资源配置）
        template_data = {
            "分子名称": ["EC", "DMC", "Li+", "PF6-"] if include_examples else [],
            "SMILES": ["C1COC(=O)O1", "COC(=O)OC", "[Li+]", "F[P-](F)(F)(F)(F)F"] if include_examples else [],
            "分子类型": ["solvent", "solvent", "cation", "anion"] if include_examples else [],
            "泛函": ["B3LYP", "B3LYP", "B3LYP", "B3LYP"] if include_examples else [],
            "基组": ["6-31++G(d,p)", "6-31++G(d,p)", "6-31++G(d,p)", "6-31++G(d,p)"] if include_examples else [],
            "溶剂模型": ["pcm", "pcm", "pcm", "pcm"] if include_examples else [],
            "溶剂名称": ["Water", "Water", "Water", "Water"] if include_examples else [],
        }
        filename = "QC计算批量导入模板.xlsx"

    else:
        raise HTTPException(status_code=400, detail="不支持的模板类型")

    # 创建DataFrame
    df = pd.DataFrame(template_data)

    # 写入Excel
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        sheet_name = '配方' if template_type == 'electrolyte' else 'QC计算'
        df.to_excel(writer, index=False, sheet_name=sheet_name)

        # 调整列宽
        worksheet = writer.sheets[sheet_name]
        for column in worksheet.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = min(max_length + 2, 50)
            worksheet.column_dimensions[column_letter].width = adjusted_width

        # 添加说明sheet
        if template_type == 'electrolyte':
            instructions = pd.DataFrame({
                "字段名": [
                    "配方名称", "阳离子名称", "阳离子SMILES", "阳离子数量",
                    "阴离子名称", "阴离子SMILES", "阴离子数量",
                    "溶剂名称", "溶剂SMILES", "溶剂数量",
                    "温度(K)", "压强(atm)", "盒子尺寸(Å)",
                    "NPT步数", "NVT步数", "时间步长(fs)",
                    "提交MD计算", "队列", "节点数", "任务数", "每任务CPU数", "最大运行时间(分钟)"
                ],
                "是否必填": [
                    "必填", "必填", "必填", "必填",
                    "必填", "必填", "必填",
                    "可选", "可选", "可选",
                    "可选", "可选", "必填",
                    "可选", "可选", "可选",
                    "可选", "可选", "可选", "可选", "可选", "可选"
                ],
                "说明": [
                    "配方的描述性名称，系统会自动生成规范名称（EL-日期-序号-描述）",
                    "阳离子名称，多个用分号分隔，如：Li+;Na+",
                    "阳离子SMILES，多个用分号分隔，如：[Li+];[Na+]",
                    "阳离子数量，多个用分号分隔，如：50;50",
                    "阴离子名称，多个用分号分隔，如：PF6-;TFSI-",
                    "阴离子SMILES，多个用分号分隔",
                    "阴离子数量，多个用分号分隔",
                    "溶剂名称，多个用分号分隔，如：EC;DMC",
                    "溶剂SMILES，多个用分号分隔",
                    "溶剂数量，多个用分号分隔，如：100;100",
                    "模拟温度，默认298.15K",
                    "模拟压强，默认1.0atm",
                    "模拟盒子边长（立方盒子），单位Å，必须填写！",
                    "NPT平衡步数，默认5000000",
                    "NVT生产步数，默认10000000",
                    "时间步长，默认1.0fs",
                    "是否立即提交MD计算，TRUE或FALSE",
                    "Slurm队列名称，默认cpu",
                    "计算节点数，默认1",
                    "任务数，默认8",
                    "每任务CPU数，默认8",
                    "最大运行时间（分钟），默认7200"
                ]
            })
        else:  # QC
            instructions = pd.DataFrame({
                "字段名": [
                    "分子名称", "SMILES", "分子类型",
                    "泛函", "基组",
                    "溶剂模型", "溶剂名称"
                ],
                "是否必填": [
                    "必填", "必填", "必填",
                    "可选", "可选",
                    "可选", "可选"
                ],
                "说明": [
                    "分子名称，如：EC、Li+、PF6-",
                    "分子的SMILES表示，如：C1COC(=O)O1",
                    "分子类型：solvent（溶剂，电荷=0）、cation（阳离子，电荷=1）、anion（阴离子，电荷=-1）、custom（自定义）",
                    "DFT泛函，默认B3LYP",
                    "基组，默认6-31++G(d,p)",
                    "溶剂模型：gas（气相）、pcm、smd，默认pcm",
                    "溶剂名称，如：Water、Acetonitrile，默认Water"
                ]
            })

        instructions.to_excel(writer, index=False, sheet_name='填写说明')

        # 调整说明sheet的列宽
        inst_worksheet = writer.sheets['填写说明']
        inst_worksheet.column_dimensions['A'].width = 20
        inst_worksheet.column_dimensions['B'].width = 12
        inst_worksheet.column_dimensions['C'].width = 80

    output.seek(0)

    # 使用 RFC 5987 编码中文文件名
    encoded_filename = quote(filename)

    return StreamingResponse(
        output,
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": f"attachment; filename*=UTF-8''{encoded_filename}"}
    )



def parse_semicolon_list(value: str) -> List[str]:
    """解析分号分隔的字符串"""
    if not value or pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(';') if item.strip()]


def parse_electrolyte_row(row: pd.Series) -> BatchElectrolyteRow:
    """解析配方行数据"""
    return BatchElectrolyteRow(
        formulation_name=str(row.get('配方名称', row.get('formulation_name', ''))),
        cation_names=str(row.get('阳离子名称', row.get('cation_names', ''))),
        cation_smiles=str(row.get('阳离子SMILES', row.get('cation_smiles', ''))),
        cation_numbers=str(row.get('阳离子数量', row.get('cation_numbers', ''))),
        anion_names=str(row.get('阴离子名称', row.get('anion_names', ''))),
        anion_smiles=str(row.get('阴离子SMILES', row.get('anion_smiles', ''))),
        anion_numbers=str(row.get('阴离子数量', row.get('anion_numbers', ''))),
        solvent_names=str(row.get('溶剂名称', row.get('solvent_names', ''))) if pd.notna(row.get('溶剂名称', row.get('solvent_names'))) else None,
        solvent_smiles=str(row.get('溶剂SMILES', row.get('solvent_smiles', ''))) if pd.notna(row.get('溶剂SMILES', row.get('solvent_smiles'))) else None,
        solvent_numbers=str(row.get('溶剂数量', row.get('solvent_numbers', ''))) if pd.notna(row.get('溶剂数量', row.get('solvent_numbers'))) else None,
        temperature=float(row.get('温度(K)', row.get('temperature', 298.15))),
        pressure=float(row.get('压强(atm)', row.get('pressure', 1.0))),
        box_size=float(row.get('盒子尺寸(Å)', row.get('box_size'))) if pd.notna(row.get('盒子尺寸(Å)', row.get('box_size'))) else None,
        nsteps_npt=int(row.get('NPT步数', row.get('nsteps_npt', 5000000))),
        nsteps_nvt=int(row.get('NVT步数', row.get('nsteps_nvt', 10000000))),
        timestep=float(row.get('时间步长(fs)', row.get('timestep', 1.0))),
        submit_md=bool(row.get('提交MD计算', row.get('submit_md', True))),
        slurm_partition=str(row.get('队列', row.get('slurm_partition', 'cpu'))),
        slurm_nodes=int(row.get('节点数', row.get('slurm_nodes', 1))),
        slurm_ntasks=int(row.get('任务数', row.get('slurm_ntasks', 8))),
        slurm_cpus_per_task=int(row.get('每任务CPU数', row.get('slurm_cpus_per_task', 8))),
        slurm_time=int(row.get('最大运行时间(分钟)', row.get('slurm_time', 7200))),
    )


def parse_qc_row(row: pd.Series) -> BatchQCRow:
    """
    解析QC计算行数据
    智能推断电荷和自旋多重度：
    - solvent: 电荷=0, 自旋多重度=1
    - cation: 电荷=1, 自旋多重度=1
    - anion: 电荷=-1, 自旋多重度=1
    - custom: 使用用户指定的值或默认值
    """
    molecule_type = str(row.get('分子类型', row.get('molecule_type', 'custom'))).lower()

    # 智能推断电荷和自旋多重度
    if molecule_type == 'solvent':
        default_charge = 0
        default_spin = 1
    elif molecule_type == 'cation':
        default_charge = 1
        default_spin = 1
    elif molecule_type == 'anion':
        default_charge = -1
        default_spin = 1
    else:  # custom
        default_charge = 0
        default_spin = 1

    # 如果用户在模板中指定了电荷和自旋多重度，使用用户的值；否则使用智能推断的值
    charge = int(row.get('电荷', row.get('charge', default_charge))) if pd.notna(row.get('电荷', row.get('charge'))) else default_charge
    spin_multiplicity = int(row.get('自旋多重度', row.get('spin_multiplicity', default_spin))) if pd.notna(row.get('自旋多重度', row.get('spin_multiplicity'))) else default_spin

    return BatchQCRow(
        molecule_name=str(row.get('分子名称', row.get('molecule_name', ''))),
        smiles=str(row.get('SMILES', row.get('smiles', ''))),
        molecule_type=molecule_type,
        functional=str(row.get('泛函', row.get('functional', 'B3LYP'))),
        basis_set=str(row.get('基组', row.get('basis_set', '6-31++G(d,p)'))),
        charge=charge,
        spin_multiplicity=spin_multiplicity,
        solvent_model=str(row.get('溶剂模型', row.get('solvent_model', 'pcm'))) if pd.notna(row.get('溶剂模型', row.get('solvent_model'))) else None,
        solvent_name=str(row.get('溶剂名称', row.get('solvent_name', 'Water'))) if pd.notna(row.get('溶剂名称', row.get('solvent_name'))) else None,
        submit_qc=bool(row.get('提交QC计算', row.get('submit_qc', True))),
        # 使用默认的资源配置
        slurm_partition='cpu',
        slurm_nodes=1,
        slurm_ntasks=1,
        slurm_cpus_per_task=8,
        slurm_time=1440,
    )


@router.post("/upload", response_model=BatchImportResult)
async def batch_import_upload(
    file: UploadFile = File(...),
    project_id: Optional[int] = Form(None),
    project_name: Optional[str] = Form(None),
    project_description: Optional[str] = Form(None),
    sheet_name: Optional[str] = Form("配方"),  # Excel sheet名称
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    批量导入配方和QC计算
    支持Excel (.xlsx, .xls) 和 CSV (.csv) 格式
    """
    try:
        # 读取文件
        content = await file.read()

        # 根据文件类型解析
        if file.filename.endswith('.csv'):
            df = pd.read_csv(io.BytesIO(content))
        elif file.filename.endswith(('.xlsx', '.xls')):
            # 尝试读取指定的sheet
            try:
                df = pd.read_excel(io.BytesIO(content), sheet_name=sheet_name)
            except:
                # 如果指定sheet不存在，读取第一个sheet
                df = pd.read_excel(io.BytesIO(content), sheet_name=0)
        else:
            raise HTTPException(status_code=400, detail="不支持的文件格式，请上传CSV或Excel文件")

        # 检查是否为空
        if df.empty:
            raise HTTPException(status_code=400, detail="文件内容为空")

        # 确定项目
        if project_id:
            project = db.query(Project).filter(
                Project.id == project_id,
                Project.user_id == current_user.id
            ).first()
            if not project:
                raise HTTPException(status_code=404, detail="项目不存在")
        else:
            # 创建新项目
            if not project_name:
                project_name = f"批量导入_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

            project = Project(
                user_id=current_user.id,
                name=project_name,
                description=project_description or "通过批量导入创建"
            )
            db.add(project)
            db.commit()
            db.refresh(project)

        result = BatchImportResult(
            success=True,
            message="批量导入完成",
            project_id=project.id,
            project_name=project.name
        )

        # 检测文件类型（配方 or QC）
        columns = [col.lower() for col in df.columns]
        is_electrolyte = any('配方' in col or 'formulation' in col or '阳离子' in col or 'cation' in col for col in df.columns)
        is_qc = any('smiles' in col or '分子' in col or 'molecule' in col for col in df.columns)

        # 处理配方导入
        if is_electrolyte:
            result = await process_electrolyte_import(df, project, current_user, db, result)

        # 处理QC导入
        if is_qc and not is_electrolyte:
            result = await process_qc_import(df, project, current_user, db, result)

        return result

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"批量导入失败: {str(e)}")


async def process_electrolyte_import(
    df: pd.DataFrame,
    project: Project,
    current_user: User,
    db: Session,
    result: BatchImportResult
) -> BatchImportResult:
    """处理配方导入"""
    from datetime import datetime
    from app.utils.hash import calculate_system_hash

    result.total_electrolytes = len(df)

    # 获取当前日期和已有配方数量（用于生成序号）
    date_str = datetime.now().strftime("%Y%m%d")
    count = db.query(ElectrolyteSystem).count()

    for idx, row in df.iterrows():
        try:
            # 解析行数据
            electrolyte_data = parse_electrolyte_row(row)

            # 验证必填字段
            if not electrolyte_data.box_size:
                raise ValueError("盒子尺寸(Å)为必填字段，请填写")

            # 构建分子列表 - 清理离子名称中的电荷符号
            cation_names = [strip_ion_charge(name) for name in parse_semicolon_list(electrolyte_data.cation_names)]
            cation_smiles_list = parse_semicolon_list(electrolyte_data.cation_smiles)
            cation_numbers_list = [int(n) for n in parse_semicolon_list(electrolyte_data.cation_numbers)]

            anion_names = [strip_ion_charge(name) for name in parse_semicolon_list(electrolyte_data.anion_names)]
            anion_smiles_list = parse_semicolon_list(electrolyte_data.anion_smiles)
            anion_numbers_list = [int(n) for n in parse_semicolon_list(electrolyte_data.anion_numbers)]

            # 构建阳离子
            cations = [
                MoleculeSpec(name=name, smiles=smiles, number=number)
                for name, smiles, number in zip(cation_names, cation_smiles_list, cation_numbers_list)
            ]

            # 构建阴离子
            anions = [
                MoleculeSpec(name=name, smiles=smiles, number=number)
                for name, smiles, number in zip(anion_names, anion_smiles_list, anion_numbers_list)
            ]

            # 构建溶剂
            solvents = []
            if electrolyte_data.solvent_names:
                solvent_names = parse_semicolon_list(electrolyte_data.solvent_names)
                solvent_smiles_list = parse_semicolon_list(electrolyte_data.solvent_smiles)
                solvent_numbers_list = [int(n) for n in parse_semicolon_list(electrolyte_data.solvent_numbers)]

                solvents = [
                    MoleculeSpec(name=name, smiles=smiles, number=number)
                    for name, smiles, number in zip(solvent_names, solvent_smiles_list, solvent_numbers_list)
                ]

            # 生成符合命名规则的配方名称：EL-日期-序号-描述
            seq_number = count + idx + 1
            description = electrolyte_data.formulation_name
            electrolyte_name = f"EL-{date_str}-{seq_number:04d}-{description}"

            # 计算hash检查重复
            hash_key = calculate_system_hash(
                cations=[c.model_dump() for c in cations],
                anions=[a.model_dump() for a in anions],
                solvents=[s.model_dump() for s in solvents] if solvents else [],
                temperature=electrolyte_data.temperature,
                pressure=electrolyte_data.pressure
            )

            # 检查是否已存在相同配方
            existing = db.query(ElectrolyteSystem).filter(
                ElectrolyteSystem.hash_key == hash_key
            ).first()

            if existing:
                # 构建详细的对比信息
                import json
                existing_cations = json.dumps(existing.cations, ensure_ascii=False)
                existing_anions = json.dumps(existing.anions, ensure_ascii=False)
                existing_solvents = json.dumps(existing.solvents, ensure_ascii=False) if existing.solvents else "无"

                new_cations = json.dumps([c.model_dump() for c in cations], ensure_ascii=False)
                new_anions = json.dumps([a.model_dump() for a in anions], ensure_ascii=False)
                new_solvents = json.dumps([s.model_dump() for s in solvents], ensure_ascii=False) if solvents else "无"

                error_msg = (
                    f"配方已存在（ID: {existing.id}, 名称: {existing.name}）\n"
                    f"已存在配方：\n"
                    f"  阳离子: {existing_cations}\n"
                    f"  阴离子: {existing_anions}\n"
                    f"  溶剂: {existing_solvents}\n"
                    f"  温度: {existing.temperature}K, 压强: {existing.pressure}atm\n"
                    f"新导入配方：\n"
                    f"  阳离子: {new_cations}\n"
                    f"  阴离子: {new_anions}\n"
                    f"  溶剂: {new_solvents}\n"
                    f"  温度: {electrolyte_data.temperature}K, 压强: {electrolyte_data.pressure}atm\n"
                    f"如果数量不同但仍提示重复，请检查hash计算逻辑"
                )
                raise ValueError(error_msg)

            # 创建配方
            electrolyte_create = ElectrolyteCreate(
                project_id=project.id,
                name=electrolyte_name,  # 使用生成的规范名称
                cations=cations,
                anions=anions,
                solvents=solvents if solvents else None,
                temperature=electrolyte_data.temperature,
                pressure=electrolyte_data.pressure,
                box_size=electrolyte_data.box_size,
                nsteps_npt=electrolyte_data.nsteps_npt,
                nsteps_nvt=electrolyte_data.nsteps_nvt,
                timestep=electrolyte_data.timestep,
            )

            # 调用创建配方的API
            electrolyte = create_electrolyte(electrolyte_create, db, current_user)

            result.success_electrolytes += 1
            result.success_electrolyte_ids.append(electrolyte.id)
            result.electrolyte_results.append({
                "row": idx + 2,  # Excel行号（从2开始，因为有表头）
                "name": electrolyte.name,  # 返回实际创建的名称
                "original_name": description,  # 原始描述
                "id": electrolyte.id,
                "status": "success"
            })

            # 如果需要提交MD计算
            if electrolyte_data.submit_md:
                try:
                    from app.schemas.job import MDJobCreate

                    md_job_create = MDJobCreate(
                        system_id=electrolyte.id,
                        job_name=None,  # 让系统自动生成名称
                        submit_to_cluster=False,  # 批量导入时不立即提交，让用户后续批量提交
                        slurm_partition=electrolyte_data.slurm_partition,
                        slurm_nodes=electrolyte_data.slurm_nodes,
                        slurm_ntasks=electrolyte_data.slurm_ntasks,
                        slurm_cpus_per_task=electrolyte_data.slurm_cpus_per_task,
                        slurm_time=electrolyte_data.slurm_time,
                    )

                    md_job = create_md_job(md_job_create, db, current_user)

                    result.total_md_jobs += 1
                    result.success_md_jobs += 1
                    result.success_md_job_ids.append(md_job.id)
                    result.md_job_results.append({
                        "row": idx + 2,
                        "electrolyte_name": electrolyte.name,
                        "job_id": md_job.id,
                        "job_name": md_job.config.get('job_name', 'N/A'),
                        "status": "success"
                    })
                except HTTPException as e:
                    result.total_md_jobs += 1
                    result.failed_md_jobs += 1
                    error_detail = e.detail if isinstance(e.detail, str) else str(e.detail)
                    result.errors.append({
                        "row": idx + 2,
                        "type": "MD任务创建失败",
                        "message": error_detail
                    })
                except Exception as e:
                    result.total_md_jobs += 1
                    result.failed_md_jobs += 1
                    import traceback
                    error_msg = f"{str(e)}\n详细信息: {traceback.format_exc()}"
                    result.errors.append({
                        "row": idx + 2,
                        "type": "MD任务创建失败",
                        "message": error_msg
                    })

        except HTTPException as e:
            # 捕获HTTP异常（如重复配方）
            result.failed_electrolytes += 1
            error_detail = e.detail if isinstance(e.detail, str) else str(e.detail)
            result.errors.append({
                "row": idx + 2,
                "type": "配方创建失败",
                "message": error_detail
            })
        except Exception as e:
            result.failed_electrolytes += 1
            result.errors.append({
                "row": idx + 2,
                "type": "配方创建失败",
                "message": str(e)
            })

    return result


async def process_qc_import(
    df: pd.DataFrame,
    project: Project,
    current_user: User,
    db: Session,
    result: BatchImportResult
) -> BatchImportResult:
    """处理QC计算导入"""
    from app.schemas.qc import SolventConfig

    result.total_qc_jobs = len(df)

    for idx, row in df.iterrows():
        try:
            # 解析行数据
            qc_data = parse_qc_row(row)

            # 构建溶剂配置
            solvent_config = None
            if qc_data.solvent_model and qc_data.solvent_model.lower() != 'gas':
                solvent_config = SolventConfig(
                    model=qc_data.solvent_model.lower(),
                    solvent_name=qc_data.solvent_name
                )

            # 创建QC任务
            qc_job_create = QCJobCreate(
                molecule_name=qc_data.molecule_name,
                smiles=qc_data.smiles,
                molecule_type=qc_data.molecule_type,
                functional=qc_data.functional,
                basis_set=qc_data.basis_set,
                charge=qc_data.charge,
                spin_multiplicity=qc_data.spin_multiplicity,
                solvent_config=solvent_config,
                slurm_partition=qc_data.slurm_partition,
                slurm_cpus=qc_data.slurm_cpus_per_task,
                slurm_time=qc_data.slurm_time,
            )

            qc_job = create_qc_job(qc_job_create, db, current_user)

            result.success_qc_jobs += 1
            result.success_qc_job_ids.append(qc_job.id)
            result.qc_job_results.append({
                "row": idx + 2,
                "molecule_name": qc_data.molecule_name,
                "job_id": qc_job.id,
                "status": "success"
            })

        except HTTPException as e:
            # 捕获HTTP异常（如查重409错误）
            result.failed_qc_jobs += 1

            # 提取错误信息
            if isinstance(e.detail, dict):
                # 409冲突错误，包含详细信息
                error_msg = e.detail.get("message", "")
                existing_id = e.detail.get("existing_job_id")
                existing_status = e.detail.get("existing_job_status", "")

                if existing_id:
                    # 构建友好的错误信息
                    error_msg = f"重复计算：已存在相同参数的QC任务（ID: {existing_id}，状态: {existing_status}）"
                elif not error_msg:
                    # 如果没有message字段，尝试构建错误信息
                    error_msg = f"任务创建失败：{e.detail}"
            elif isinstance(e.detail, str):
                error_msg = e.detail
            else:
                error_msg = f"任务创建失败（状态码: {e.status_code}）"

            result.errors.append({
                "row": idx + 2,
                "type": "QC任务创建失败",
                "message": error_msg
            })

        except Exception as e:
            result.failed_qc_jobs += 1
            # 获取更详细的错误信息
            import traceback
            error_detail = str(e)
            if not error_detail or error_detail in ['PENDING', 'CREATED', 'RUNNING', 'COMPLETED', 'FAILED', 'CANCELLED']:
                # 如果错误信息是枚举值，获取完整的堆栈跟踪
                error_detail = f"任务创建失败: {traceback.format_exc()}"

            result.errors.append({
                "row": idx + 2,
                "type": "QC任务创建失败",
                "message": error_detail
            })

    return result

