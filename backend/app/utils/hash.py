"""
Hash utilities for electrolyte system deduplication
"""
import hashlib
import json
from typing import List, Dict, Any


def normalize_molecule_list(molecules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Normalize molecule list by sorting
    
    Args:
        molecules: List of molecule specifications
        
    Returns:
        List: Sorted molecule list
    """
    return sorted(molecules, key=lambda x: (x.get("smiles", ""), x.get("name", "")))


def calculate_system_hash(
    cations: List[Dict[str, Any]],
    anions: List[Dict[str, Any]],
    solvents: List[Dict[str, Any]],
    temperature: float,
    pressure: float,
) -> str:
    """
    Calculate unique hash for an electrolyte system
    
    Args:
        cations: List of cation specifications
        anions: List of anion specifications
        solvents: List of solvent specifications
        temperature: Temperature in K
        pressure: Pressure in atm
        
    Returns:
        str: SHA256 hash of the system
    """
    # Normalize molecule lists
    norm_cations = normalize_molecule_list(cations)
    norm_anions = normalize_molecule_list(anions)
    norm_solvents = normalize_molecule_list(solvents) if solvents else []
    
    # Create canonical representation
    system_dict = {
        "cations": norm_cations,
        "anions": norm_anions,
        "solvents": norm_solvents,
        "temperature": round(temperature, 2),
        "pressure": round(pressure, 2),
    }
    
    # Convert to JSON string (sorted keys for consistency)
    system_str = json.dumps(system_dict, sort_keys=True)
    
    # Calculate SHA256 hash
    hash_obj = hashlib.sha256(system_str.encode("utf-8"))
    return hash_obj.hexdigest()

