"""
Token schemas for authentication
"""
from pydantic import BaseModel
from typing import Optional


class Token(BaseModel):
    """Token response"""
    access_token: str
    token_type: str


class TokenData(BaseModel):
    """Token payload data"""
    username: Optional[str] = None

