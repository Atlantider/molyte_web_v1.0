"""
Application configuration using Pydantic Settings
"""
from pydantic_settings import BaseSettings
from typing import List


class Settings(BaseSettings):
    """Application settings"""

    # Database
    DATABASE_URL: str = "postgresql://molyte_user:molyte2025@localhost:5432/molyte_web"

    # Security
    SECRET_KEY: str = "09d25e094faa6ca2556c818166b7a9563b93f7099f6f0f4caa6cf63b88e8d3e7"
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 1440  # 24小时 (原来是30分钟)

    # Application
    APP_NAME: str = "Molyte Web API"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = True

    # CORS - stored as comma-separated string in .env
    CORS_ORIGINS: str = "http://localhost:3000,http://localhost:5173,http://localhost:8080"

    # Logging
    LOG_LEVEL: str = "INFO"

    # Redis settings
    REDIS_HOST: str = "localhost"
    REDIS_PORT: int = 6379
    REDIS_PASSWORD: str = ""
    REDIS_DB: int = 0

    # Celery settings
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://localhost:6379/1"

    @property
    def cors_origins_list(self) -> List[str]:
        """Convert CORS_ORIGINS string to list"""
        return [origin.strip() for origin in self.CORS_ORIGINS.split(",")]

    class Config:
        env_file = ".env"
        case_sensitive = True


# Create global settings instance
settings = Settings()

