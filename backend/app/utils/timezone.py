"""
Timezone utilities
"""
from datetime import datetime
import pytz


def ensure_timezone_aware(dt: datetime) -> datetime:
    """
    Ensure a datetime object is timezone-aware
    
    Args:
        dt: datetime object (may be naive or aware)
        
    Returns:
        datetime: timezone-aware datetime object (Asia/Shanghai)
    """
    if dt is None:
        return None
    
    # If already timezone-aware, return as-is
    if dt.tzinfo is not None and dt.tzinfo.utcoffset(dt) is not None:
        return dt
    
    # If naive, assume it's in Asia/Shanghai timezone
    tz = pytz.timezone('Asia/Shanghai')
    return tz.localize(dt)

