"""
计费和充值 API 端点
"""
from typing import List, Optional
from datetime import datetime
from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from app.database import get_db
from app.dependencies import get_current_active_user, get_current_admin_user
from app.models.user import User
from app.models.billing import RechargeOrder, QuotaTransaction, PaymentStatus, PaymentMethod
from app.services.billing import BillingService

router = APIRouter(prefix="/billing", tags=["billing"])


# ============ Schemas ============

class BalanceResponse(BaseModel):
    """用户余额响应"""
    balance: float
    frozen: float
    debt: float
    available: float
    price_per_hour: float
    
    class Config:
        from_attributes = True


class CreateOrderRequest(BaseModel):
    """创建充值订单请求"""
    amount: float = Field(..., ge=1, description="充值金额（元）")
    payment_method: str = Field(default="simulated", description="支付方式")


class OrderResponse(BaseModel):
    """订单响应"""
    id: int
    order_no: str
    amount: float
    cpu_hours: float
    price_per_hour: float
    payment_method: str
    payment_status: str
    created_at: datetime
    paid_at: Optional[datetime]
    
    class Config:
        from_attributes = True


class TransactionResponse(BaseModel):
    """交易记录响应"""
    id: int
    type: str
    amount: float
    balance_before: float
    balance_after: float
    description: Optional[str]
    created_at: datetime
    
    class Config:
        from_attributes = True


class SimulatePayRequest(BaseModel):
    """模拟支付请求"""
    order_id: int


# ============ 用户端 API ============

@router.get("/balance", response_model=BalanceResponse)
async def get_balance(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取当前用户余额信息"""
    balance_info = BillingService.get_user_balance(db, current_user.id)
    if not balance_info:
        raise HTTPException(status_code=404, detail="用户不存在")
    
    balance_info["price_per_hour"] = BillingService.get_cpu_hour_price(db)
    return balance_info


@router.get("/can-submit")
async def check_can_submit(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """检查用户是否可以提交任务"""
    can_submit, reason = BillingService.can_submit_job(db, current_user)
    return {
        "can_submit": can_submit,
        "reason": reason
    }


@router.post("/orders", response_model=OrderResponse)
async def create_order(
    request: CreateOrderRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """创建充值订单"""
    # 验证最低充值金额
    min_amount = float(BillingService.get_config(db, 'min_recharge_amount', '10'))
    if request.amount < min_amount:
        raise HTTPException(
            status_code=400, 
            detail=f"最低充值金额为 ¥{min_amount}"
        )
    
    # 映射支付方式
    method_map = {
        "wechat": PaymentMethod.WECHAT,
        "alipay": PaymentMethod.ALIPAY,
        "simulated": PaymentMethod.SIMULATED,
    }
    payment_method = method_map.get(request.payment_method, PaymentMethod.SIMULATED)
    
    order = BillingService.create_recharge_order(
        db=db,
        user_id=current_user.id,
        amount=request.amount,
        payment_method=payment_method
    )
    
    return OrderResponse(
        id=order.id,
        order_no=order.order_no,
        amount=order.amount,
        cpu_hours=order.cpu_hours,
        price_per_hour=order.price_per_hour,
        payment_method=order.payment_method,  # 已经是字符串
        payment_status=order.payment_status,  # 已经是字符串
        created_at=order.created_at,
        paid_at=order.paid_at
    )


@router.get("/orders", response_model=List[OrderResponse])
async def get_orders(
    skip: int = 0,
    limit: int = 20,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取用户充值订单列表"""
    orders = BillingService.get_orders(db, current_user.id, skip, limit)
    return [OrderResponse(
        id=o.id,
        order_no=o.order_no,
        amount=o.amount,
        cpu_hours=o.cpu_hours,
        price_per_hour=o.price_per_hour,
        payment_method=o.payment_method,  # 已经是字符串
        payment_status=o.payment_status,  # 已经是字符串
        created_at=o.created_at,
        paid_at=o.paid_at
    ) for o in orders]


@router.get("/transactions", response_model=List[TransactionResponse])
async def get_transactions(
    skip: int = 0,
    limit: int = 20,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取用户配额变更记录"""
    transactions = BillingService.get_transactions(db, current_user.id, skip, limit)
    return [TransactionResponse(
        id=t.id,
        type=t.type,  # 已经是字符串
        amount=t.amount,
        balance_before=t.balance_before,
        balance_after=t.balance_after,
        description=t.description,
        created_at=t.created_at
    ) for t in transactions]


@router.post("/simulate-pay")
async def simulate_payment(
    request: SimulatePayRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """模拟支付（测试用）"""
    # 验证订单属于当前用户
    order = db.query(RechargeOrder).filter(
        RechargeOrder.id == request.order_id,
        RechargeOrder.user_id == current_user.id
    ).first()

    if not order:
        raise HTTPException(status_code=404, detail="订单不存在")

    if order.payment_status != PaymentStatus.PENDING:
        raise HTTPException(status_code=400, detail=f"订单状态异常: {order.payment_status.value}")

    success, message = BillingService.complete_payment(db, order.id)

    if not success:
        raise HTTPException(status_code=400, detail=message)

    # 获取更新后的余额
    balance_info = BillingService.get_user_balance(db, current_user.id)

    return {
        "success": True,
        "message": message,
        "balance": balance_info
    }


# ============ 管理端 API ============

class PricingConfig(BaseModel):
    """定价配置"""
    cpu_hour_price: float = Field(..., ge=0.01, description="每核时单价（元）")
    min_recharge_amount: float = Field(..., ge=1, description="最低充值金额（元）")
    max_debt_cpu_hours: float = Field(..., ge=0, description="最大允许欠费机时")


class AdminAdjustRequest(BaseModel):
    """管理员调整余额请求"""
    user_id: int
    amount: float
    reason: str


@router.get("/admin/pricing", response_model=PricingConfig)
async def get_pricing_config(
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取定价配置（管理员）"""
    return PricingConfig(
        cpu_hour_price=float(BillingService.get_config(db, 'cpu_hour_price', '0.1')),
        min_recharge_amount=float(BillingService.get_config(db, 'min_recharge_amount', '10')),
        max_debt_cpu_hours=float(BillingService.get_config(db, 'max_debt_cpu_hours', '100'))
    )


@router.put("/admin/pricing")
async def update_pricing_config(
    config: PricingConfig,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """更新定价配置（管理员）"""
    BillingService.set_config(db, 'cpu_hour_price', str(config.cpu_hour_price),
                              '每核时单价（元）', admin.id)
    BillingService.set_config(db, 'min_recharge_amount', str(config.min_recharge_amount),
                              '最低充值金额（元）', admin.id)
    BillingService.set_config(db, 'max_debt_cpu_hours', str(config.max_debt_cpu_hours),
                              '最大允许欠费机时', admin.id)

    return {"success": True, "message": "定价配置已更新"}


@router.post("/admin/adjust")
async def admin_adjust_balance(
    request: AdminAdjustRequest,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """管理员手动调整用户余额"""
    success, message = BillingService.admin_adjust_balance(
        db=db,
        admin_id=admin.id,
        user_id=request.user_id,
        amount=request.amount,
        reason=request.reason
    )

    if not success:
        raise HTTPException(status_code=400, detail=message)

    return {"success": True, "message": message}


@router.get("/admin/orders", response_model=List[OrderResponse])
async def get_all_orders(
    skip: int = 0,
    limit: int = 50,
    status: Optional[str] = None,
    db: Session = Depends(get_db),
    admin: User = Depends(get_current_admin_user)
):
    """获取所有充值订单（管理员）"""
    query = db.query(RechargeOrder)

    if status:
        query = query.filter(RechargeOrder.payment_status == status)

    orders = query.order_by(RechargeOrder.created_at.desc()).offset(skip).limit(limit).all()

    return [OrderResponse(
        id=o.id,
        order_no=o.order_no,
        amount=o.amount,
        cpu_hours=o.cpu_hours,
        price_per_hour=o.price_per_hour,
        payment_method=o.payment_method,  # 已经是字符串
        payment_status=o.payment_status,  # 已经是字符串
        created_at=o.created_at,
        paid_at=o.paid_at
    ) for o in orders]

