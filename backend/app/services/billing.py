"""
计费服务 - 配额管理、充值、扣费等核心逻辑
"""
from datetime import datetime, timedelta
from typing import Optional, Dict, Any, Tuple
from sqlalchemy.orm import Session
from sqlalchemy import func
import uuid

from app.models.user import User
from app.models.job import MDJob, JobStatus
from app.models.billing import (
    SystemConfig, RechargeOrder, QuotaTransaction,
    PaymentMethod, PaymentStatus, TransactionType
)
from app.core.logger import logger


class BillingService:
    """计费服务类"""
    
    # 默认配置
    DEFAULT_CPU_HOUR_PRICE = 0.1  # 元/核时
    DEFAULT_MIN_RECHARGE = 10.0   # 最低充值金额
    DEFAULT_MAX_DEBT = 100.0      # 最大欠费机时
    
    @staticmethod
    def get_config(db: Session, key: str, default: str = None) -> Optional[str]:
        """获取系统配置"""
        config = db.query(SystemConfig).filter(SystemConfig.key == key).first()
        return config.value if config else default
    
    @staticmethod
    def set_config(db: Session, key: str, value: str, description: str = None, updated_by: int = None):
        """设置系统配置"""
        config = db.query(SystemConfig).filter(SystemConfig.key == key).first()
        if config:
            config.value = value
            config.updated_at = datetime.now()
            config.updated_by = updated_by
            if description:
                config.description = description
        else:
            config = SystemConfig(
                key=key,
                value=value,
                description=description,
                updated_by=updated_by
            )
            db.add(config)
        db.commit()
        return config
    
    @staticmethod
    def get_cpu_hour_price(db: Session) -> float:
        """获取当前机时单价"""
        price = BillingService.get_config(db, 'cpu_hour_price')
        return float(price) if price else BillingService.DEFAULT_CPU_HOUR_PRICE
    
    @staticmethod
    def get_user_balance(db: Session, user_id: int) -> Dict[str, float]:
        """获取用户余额信息"""
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return None
        return {
            "balance": user.balance_cpu_hours,
            "frozen": user.frozen_cpu_hours,
            "debt": user.debt_cpu_hours,
            "available": user.balance_cpu_hours - user.frozen_cpu_hours,
        }
    
    @staticmethod
    def can_submit_job(db: Session, user: User) -> Tuple[bool, str]:
        """检查用户是否可以提交任务"""
        # 管理员不限制
        if user.role.value == "ADMIN":
            return True, ""
        
        # 检查欠费
        if user.debt_cpu_hours > 0:
            return False, f"您有 {user.debt_cpu_hours:.2f} 核时欠费，请先充值还清欠费后再提交任务"
        
        # 检查余额
        if user.balance_cpu_hours <= 0:
            return False, "余额不足，请先充值后再提交任务"
        
        return True, ""
    
    @staticmethod
    def generate_order_no() -> str:
        """生成订单号"""
        now = datetime.now()
        return f"RC{now.strftime('%Y%m%d%H%M%S')}{uuid.uuid4().hex[:8].upper()}"
    
    @staticmethod
    def create_recharge_order(
        db: Session,
        user_id: int,
        amount: float,
        payment_method: PaymentMethod = PaymentMethod.SIMULATED,
        remark: str = None
    ) -> RechargeOrder:
        """创建充值订单"""
        price_per_hour = BillingService.get_cpu_hour_price(db)
        cpu_hours = amount / price_per_hour

        order = RechargeOrder(
            order_no=BillingService.generate_order_no(),
            user_id=user_id,
            amount=amount,
            cpu_hours=cpu_hours,
            price_per_hour=price_per_hour,
            payment_method=payment_method.value,  # 使用枚举值
            payment_status=PaymentStatus.PENDING.value,  # 使用枚举值
            expired_at=datetime.now() + timedelta(hours=2),
            remark=remark
        )
        db.add(order)
        db.commit()
        db.refresh(order)

        logger.info(f"Created recharge order: {order.order_no} for user {user_id}, amount={amount}")
        return order

    @staticmethod
    def complete_payment(db: Session, order_id: int, transaction_id: str = None) -> Tuple[bool, str]:
        """
        完成支付（模拟支付或回调确认）
        1. 增加用户余额
        2. 如有欠费，自动偿还
        3. 解锁被锁定的任务结果
        """
        order = db.query(RechargeOrder).filter(RechargeOrder.id == order_id).first()
        if not order:
            return False, "订单不存在"

        if order.payment_status != PaymentStatus.PENDING:
            return False, f"订单状态异常: {order.payment_status.value}"

        user = db.query(User).filter(User.id == order.user_id).first()
        if not user:
            return False, "用户不存在"

        # 更新订单状态
        order.payment_status = PaymentStatus.PAID.value  # 使用枚举值
        order.paid_at = datetime.now()
        order.transaction_id = transaction_id or f"SIM_{uuid.uuid4().hex[:12].upper()}"

        balance_before = user.balance_cpu_hours

        # 1. 如有欠费，先偿还欠费
        debt_repaid = 0.0
        if user.debt_cpu_hours > 0:
            if order.cpu_hours >= user.debt_cpu_hours:
                # 机时足够还清欠费
                debt_repaid = user.debt_cpu_hours
                remaining_hours = order.cpu_hours - user.debt_cpu_hours
                user.debt_cpu_hours = 0
                user.balance_cpu_hours += remaining_hours
            else:
                # 机时不够，部分偿还
                debt_repaid = order.cpu_hours
                user.debt_cpu_hours -= order.cpu_hours

            # 记录还债流水
            if debt_repaid > 0:
                debt_trans = QuotaTransaction(
                    user_id=user.id,
                    type=TransactionType.DEBT_REPAY.value,  # 使用枚举值
                    amount=-debt_repaid,
                    balance_before=balance_before,
                    balance_after=user.balance_cpu_hours,
                    reference_id=order.id,
                    reference_type="order",
                    description=f"偿还欠费 {debt_repaid:.2f} 核时"
                )
                db.add(debt_trans)
        else:
            # 无欠费，直接增加余额
            user.balance_cpu_hours += order.cpu_hours

        # 记录充值流水
        recharge_trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.RECHARGE.value,  # 使用枚举值
            amount=order.cpu_hours,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=order.id,
            reference_type="order",
            description=f"充值 ¥{order.amount:.2f}，获得 {order.cpu_hours:.2f} 核时"
        )
        db.add(recharge_trans)

        # 2. 解锁该用户所有被锁定的任务结果
        if user.debt_cpu_hours == 0:
            locked_jobs = db.query(MDJob).filter(
                MDJob.user_id == user.id,
                MDJob.result_locked == True
            ).all()
            for job in locked_jobs:
                job.result_locked = False
                job.locked_reason = None
            logger.info(f"Unlocked {len(locked_jobs)} jobs for user {user.id}")

        db.commit()
        logger.info(f"Payment completed for order {order.order_no}, user {user.id} balance: {user.balance_cpu_hours}")

        return True, f"支付成功，获得 {order.cpu_hours:.2f} 核时"

    @staticmethod
    def calculate_job_cpu_hours(job: MDJob) -> float:
        """计算任务实际消耗的机时"""
        if not job.started_at or not job.finished_at:
            return 0.0

        duration_seconds = (job.finished_at - job.started_at).total_seconds()
        duration_hours = duration_seconds / 3600
        cpu_cores = job.cpu_cores or 1

        return duration_hours * cpu_cores

    @staticmethod
    def settle_job(db: Session, job: MDJob) -> Tuple[bool, str]:
        """
        任务完成后结算
        1. 计算实际消耗机时
        2. 从用户余额扣除
        3. 余额不足则记录欠费并锁定结果
        """
        if job.billed:
            return True, "任务已结算"

        user = db.query(User).filter(User.id == job.user_id).first()
        if not user:
            return False, "用户不存在"

        # 管理员不扣费
        if user.role.value == "ADMIN":
            job.billed = True
            db.commit()
            return True, "管理员任务免费"

        # 计算实际机时
        actual_hours = BillingService.calculate_job_cpu_hours(job)
        job.actual_cpu_hours = actual_hours

        balance_before = user.balance_cpu_hours

        if user.balance_cpu_hours >= actual_hours:
            # 余额足够
            user.balance_cpu_hours -= actual_hours
            job.result_locked = False
        else:
            # 余额不足，产生欠费
            debt = actual_hours - user.balance_cpu_hours
            user.debt_cpu_hours += debt
            user.balance_cpu_hours = 0
            job.result_locked = True
            job.locked_reason = f"余额不足，欠费 {debt:.2f} 核时"
            logger.warning(f"Job {job.id} completed with debt: {debt:.2f} hours")

        # 记录消费流水
        trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.CONSUME.value,  # 使用枚举值
            amount=-actual_hours,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=job.id,
            reference_type="job",
            description=f"任务 #{job.id} 消耗 {actual_hours:.2f} 核时"
        )
        db.add(trans)

        job.billed = True
        db.commit()

        return True, f"结算完成，消耗 {actual_hours:.2f} 核时"

    @staticmethod
    def get_transactions(db: Session, user_id: int, skip: int = 0, limit: int = 20):
        """获取用户配额变更记录"""
        return db.query(QuotaTransaction).filter(
            QuotaTransaction.user_id == user_id
        ).order_by(QuotaTransaction.created_at.desc()).offset(skip).limit(limit).all()

    @staticmethod
    def get_orders(db: Session, user_id: int, skip: int = 0, limit: int = 20):
        """获取用户充值订单"""
        return db.query(RechargeOrder).filter(
            RechargeOrder.user_id == user_id
        ).order_by(RechargeOrder.created_at.desc()).offset(skip).limit(limit).all()

    @staticmethod
    def admin_adjust_balance(
        db: Session,
        admin_id: int,
        user_id: int,
        amount: float,
        reason: str
    ) -> Tuple[bool, str]:
        """管理员手动调整用户余额"""
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return False, "用户不存在"

        balance_before = user.balance_cpu_hours
        user.balance_cpu_hours += amount

        # 防止负数
        if user.balance_cpu_hours < 0:
            user.balance_cpu_hours = 0

        trans = QuotaTransaction(
            user_id=user.id,
            type=TransactionType.ADMIN_ADJUST.value,  # 使用枚举值
            amount=amount,
            balance_before=balance_before,
            balance_after=user.balance_cpu_hours,
            reference_id=admin_id,
            reference_type="admin",
            description=f"管理员调整: {reason}"
        )
        db.add(trans)
        db.commit()

        logger.info(f"Admin {admin_id} adjusted user {user_id} balance by {amount}: {reason}")
        return True, f"调整成功，当前余额 {user.balance_cpu_hours:.2f} 核时"

