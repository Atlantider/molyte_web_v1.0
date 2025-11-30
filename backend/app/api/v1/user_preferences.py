"""
用户偏好设置 API
包括自定义常用溶剂、常用离子组合等
"""
from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from sqlalchemy import Column, Integer, String, TIMESTAMP, ForeignKey, Text
from sqlalchemy.dialects.postgresql import JSONB
from datetime import datetime

from app.database import get_db
from app.models.user import User
from app.models.user_preferences import UserSolventCombination, UserIonCombination
from app.dependencies import get_current_user
from pydantic import BaseModel

router = APIRouter()


# ==================== Pydantic Schemas ====================

class SolventItem(BaseModel):
    """单个溶剂"""
    name: str
    smiles: str
    molar_ratio: float = 1.0


class CustomSolventCombination(BaseModel):
    """自定义溶剂组合"""
    name: str  # 组合名称，如 "我的EC/DMC"
    description: Optional[str] = None  # 描述
    solvents: List[SolventItem]  # 溶剂列表


class CustomSolventCombinationResponse(CustomSolventCombination):
    """返回的自定义溶剂组合（带ID）"""
    id: int
    user_id: int
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


class IonCombination(BaseModel):
    """离子组合"""
    name: str
    charge: int
    concentration: float


class CustomIonCombination(BaseModel):
    """自定义离子组合"""
    name: str  # 组合名称，如 "标准LiPF6"
    description: Optional[str] = None
    cations: List[IonCombination]
    anions: List[IonCombination]


class CustomIonCombinationResponse(CustomIonCombination):
    """返回的自定义离子组合（带ID）"""
    id: int
    user_id: int
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


# ==================== API 端点 ====================

@router.get("/solvent-combinations", response_model=List[CustomSolventCombinationResponse])
async def get_user_solvent_combinations(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """获取用户的自定义溶剂组合列表"""
    combinations = db.query(UserSolventCombination).filter(
        UserSolventCombination.user_id == current_user.id
    ).order_by(UserSolventCombination.created_at.desc()).all()

    return combinations


@router.post("/solvent-combinations", response_model=CustomSolventCombinationResponse)
async def create_solvent_combination(
    combination: CustomSolventCombination,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """创建新的自定义溶剂组合"""
    # 检查名称是否已存在
    existing = db.query(UserSolventCombination).filter(
        UserSolventCombination.user_id == current_user.id,
        UserSolventCombination.name == combination.name
    ).first()

    if existing:
        raise HTTPException(status_code=400, detail="该名称的组合已存在")

    # 创建新组合
    db_combination = UserSolventCombination(
        user_id=current_user.id,
        name=combination.name,
        description=combination.description,
        solvents=[s.dict() for s in combination.solvents]
    )

    db.add(db_combination)
    db.commit()
    db.refresh(db_combination)

    return db_combination


@router.put("/solvent-combinations/{combination_id}", response_model=CustomSolventCombinationResponse)
async def update_solvent_combination(
    combination_id: int,
    combination: CustomSolventCombination,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """更新自定义溶剂组合"""
    db_combination = db.query(UserSolventCombination).filter(
        UserSolventCombination.id == combination_id,
        UserSolventCombination.user_id == current_user.id
    ).first()

    if not db_combination:
        raise HTTPException(status_code=404, detail="组合不存在")

    # 更新字段
    db_combination.name = combination.name
    db_combination.description = combination.description
    db_combination.solvents = [s.dict() for s in combination.solvents]
    db_combination.updated_at = datetime.utcnow()

    db.commit()
    db.refresh(db_combination)

    return db_combination


@router.delete("/solvent-combinations/{combination_id}")
async def delete_solvent_combination(
    combination_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """删除自定义溶剂组合"""
    db_combination = db.query(UserSolventCombination).filter(
        UserSolventCombination.id == combination_id,
        UserSolventCombination.user_id == current_user.id
    ).first()

    if not db_combination:
        raise HTTPException(status_code=404, detail="组合不存在")

    db.delete(db_combination)
    db.commit()

    return {"message": "删除成功"}


@router.get("/ion-combinations", response_model=List[CustomIonCombinationResponse])
async def get_user_ion_combinations(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """获取用户的自定义离子组合列表"""
    combinations = db.query(UserIonCombination).filter(
        UserIonCombination.user_id == current_user.id
    ).order_by(UserIonCombination.created_at.desc()).all()

    return combinations


@router.post("/ion-combinations", response_model=CustomIonCombinationResponse)
async def create_ion_combination(
    combination: CustomIonCombination,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """创建新的自定义离子组合"""
    # 检查名称是否已存在
    existing = db.query(UserIonCombination).filter(
        UserIonCombination.user_id == current_user.id,
        UserIonCombination.name == combination.name
    ).first()
    
    if existing:
        raise HTTPException(status_code=400, detail="该名称的组合已存在")
    
    # 创建新组合
    db_combination = UserIonCombination(
        user_id=current_user.id,
        name=combination.name,
        description=combination.description,
        cations=[c.dict() for c in combination.cations],
        anions=[a.dict() for a in combination.anions]
    )
    
    db.add(db_combination)
    db.commit()
    db.refresh(db_combination)
    
    return db_combination


@router.delete("/ion-combinations/{combination_id}")
async def delete_ion_combination(
    combination_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """删除自定义离子组合"""
    db_combination = db.query(UserIonCombination).filter(
        UserIonCombination.id == combination_id,
        UserIonCombination.user_id == current_user.id
    ).first()

    if not db_combination:
        raise HTTPException(status_code=404, detail="组合不存在")

    db.delete(db_combination)
    db.commit()

    return {"message": "删除成功"}

