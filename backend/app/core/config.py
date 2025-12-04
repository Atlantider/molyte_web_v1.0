"""
Configuration for Molyte Web application
"""
from pathlib import Path
from typing import Optional
import os


class Settings:
    """Application settings"""
    
    # Database
    DATABASE_URL: str = os.getenv(
        "DATABASE_URL",
        "postgresql://molyte_user:molyte_password@localhost:5432/molyte_db"
    )
    
    # JWT
    SECRET_KEY: str = os.getenv("SECRET_KEY", "your-secret-key-here-change-in-production")
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = int(os.getenv("ACCESS_TOKEN_EXPIRE_MINUTES", "1440"))  # 默认24小时
    
    # Molyte paths
    MOLYTE_WORK_BASE_PATH: Path = Path(os.getenv(
        "MOLYTE_WORK_BASE_PATH",
        "/public/home/xiaoji/molyte_web/data/md_work"
    ))

    MOLYTE_INITIAL_SALTS_PATH: Path = Path(os.getenv(
        "MOLYTE_INITIAL_SALTS_PATH",
        "/public/home/xiaoji/molyte_web/data/initial_salts"
    ))
    
    MOLYTE_LIGPARGEN_PATH: Path = Path(os.getenv(
        "MOLYTE_LIGPARGEN_PATH",
        "/public/software/anaconda3/envs/molyte/bin"
    ))
    
    MOLYTE_PACKMOL_PATH: Path = Path(os.getenv(
        "MOLYTE_PACKMOL_PATH",
        "/public/software/packmol-20.16.0/packmol"
    ))
    
    MOLYTE_LTEMPLIFY_PATH: Path = Path(os.getenv(
        "MOLYTE_LTEMPLIFY_PATH",
        "/public/software/moltemplate_2025-2-02/moltemplate/ltemplify.py"
    ))
    
    MOLYTE_MOLTEMPLATE_PATH: Path = Path(os.getenv(
        "MOLYTE_MOLTEMPLATE_PATH",
        "/public/software/moltemplate_2025-2-02/moltemplate/scripts/moltemplate.sh"
    ))
    
    MOLYTE_CHARGE_SAVE_PATH: Path = Path(os.getenv(
        "MOLYTE_CHARGE_SAVE_PATH",
        "/public/home/xiaoji/molyte_web/data/charges"
    ))

    # QC paths
    QC_WORK_BASE_PATH: Path = Path(os.getenv(
        "QC_WORK_BASE_PATH",
        "/public/home/xiaoji/molyte_web/data/qc_work"
    ))

    # Cluster (desolvation) work paths
    CLUSTER_WORK_BASE_PATH: Path = Path(os.getenv(
        "CLUSTER_WORK_BASE_PATH",
        "/public/home/xiaoji/molyte_web/data/cluster_work"
    ))

    # Slurm settings
    SLURM_PARTITION: str = os.getenv("SLURM_PARTITION", "cpu")
    SLURM_NODES: int = int(os.getenv("SLURM_NODES", "1"))
    SLURM_NTASKS: int = int(os.getenv("SLURM_NTASKS", "32"))
    SLURM_TIME: str = os.getenv("SLURM_TIME", "48:00:00")

    # Redis settings
    REDIS_HOST: str = os.getenv("REDIS_HOST", "localhost")
    REDIS_PORT: int = int(os.getenv("REDIS_PORT", "6379"))
    REDIS_PASSWORD: str = os.getenv("REDIS_PASSWORD", "")
    REDIS_DB: int = int(os.getenv("REDIS_DB", "0"))

    # Celery settings
    CELERY_BROKER_URL: str = os.getenv(
        "CELERY_BROKER_URL",
        "redis://localhost:6379/0"
    )
    CELERY_RESULT_BACKEND: str = os.getenv(
        "CELERY_RESULT_BACKEND",
        "redis://localhost:6379/1"
    )

    def __init__(self):
        """Initialize settings and create directories if needed"""
        # Create work directories if they don't exist
        self.MOLYTE_WORK_BASE_PATH.mkdir(parents=True, exist_ok=True)
        self.QC_WORK_BASE_PATH.mkdir(parents=True, exist_ok=True)
        self.CLUSTER_WORK_BASE_PATH.mkdir(parents=True, exist_ok=True)
        self.MOLYTE_CHARGE_SAVE_PATH.mkdir(parents=True, exist_ok=True)


# Global settings instance
settings = Settings()

