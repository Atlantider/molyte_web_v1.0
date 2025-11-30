"""
Utility to parse ion information from LAMMPS .lt files
"""
import re
from pathlib import Path
from typing import Dict, List, Tuple
from app.core.logger import logger


def parse_lt_file_charge(lt_file_path: Path) -> float:
    """
    Parse a LAMMPS .lt file and calculate the total charge of the molecule
    
    Args:
        lt_file_path: Path to the .lt file
        
    Returns:
        Total charge of the molecule (rounded to nearest integer)
    """
    try:
        with open(lt_file_path, 'r') as f:
            content = f.read()
        
        # Find the "Data Atoms" section
        atoms_section_match = re.search(
            r'write\("Data Atoms"\)\s*\{(.*?)\}',
            content,
            re.DOTALL
        )
        
        if not atoms_section_match:
            logger.warning(f"No 'Data Atoms' section found in {lt_file_path}")
            return 0.0
        
        atoms_section = atoms_section_match.group(1)

        # Extract charges from atom lines
        # Format: $atom:id $mol:id @atom:type charge x y z
        # We need to match lines with atom definitions and extract the charge (4th field after @atom:)
        # The charge can be positive or negative, with optional sign
        atom_pattern = r'\$atom:\S+\s+\$mol[:\S]*\s+@atom:\S+\s+([-+]?\d+\.?\d*)'

        charges = []
        for match in re.finditer(atom_pattern, atoms_section):
            charge_str = match.group(1)
            charge = float(charge_str)
            charges.append(charge)
            logger.debug(f"Found atom with charge: {charge}")
        
        if not charges:
            logger.warning(f"No charges found in {lt_file_path}")
            return 0.0
        
        total_charge = sum(charges)
        
        # Round to nearest integer (ions should have integer charges)
        rounded_charge = round(total_charge)
        
        logger.info(f"Parsed {lt_file_path.name}: total_charge={total_charge:.3f}, rounded={rounded_charge}")
        
        return rounded_charge
        
    except Exception as e:
        logger.error(f"Error parsing {lt_file_path}: {e}")
        return 0.0


def scan_available_ions(salts_dir: Path) -> Dict[str, Dict[str, any]]:
    """
    Scan the salts directory and extract ion information
    
    Args:
        salts_dir: Path to the directory containing .lt files
        
    Returns:
        Dictionary with ion information:
        {
            "Li": {"charge": 1, "type": "cation", "file": "Li.lt"},
            "PF6": {"charge": -1, "type": "anion", "file": "PF6.lt"},
            ...
        }
    """
    ions_info = {}
    
    if not salts_dir.exists():
        logger.error(f"Salts directory not found: {salts_dir}")
        return ions_info
    
    # Find all .lt files
    lt_files = list(salts_dir.glob("*.lt"))
    
    for lt_file in lt_files:
        # Skip job.sh and other non-ion files
        if lt_file.stem in ['job', 'system']:
            continue
        
        ion_name = lt_file.stem
        charge = parse_lt_file_charge(lt_file)
        
        # Determine if cation or anion based on charge
        if charge > 0:
            ion_type = "cation"
        elif charge < 0:
            ion_type = "anion"
        else:
            logger.warning(f"Ion {ion_name} has zero charge, skipping")
            continue
        
        ions_info[ion_name] = {
            "charge": int(charge),
            "type": ion_type,
            "file": lt_file.name
        }
    
    logger.info(f"Scanned {len(ions_info)} ions from {salts_dir}")
    return ions_info


def get_cations_and_anions(ions_info: Dict[str, Dict]) -> Tuple[List[Dict], List[Dict]]:
    """
    Separate ions into cations and anions lists
    
    Args:
        ions_info: Dictionary from scan_available_ions()
        
    Returns:
        Tuple of (cations_list, anions_list)
        Each item in the list is {"name": "Li", "charge": 1}
    """
    cations = []
    anions = []
    
    for name, info in ions_info.items():
        ion_data = {"name": name, "charge": info["charge"]}
        
        if info["type"] == "cation":
            cations.append(ion_data)
        else:
            anions.append(ion_data)
    
    # Sort by name
    cations.sort(key=lambda x: x["name"])
    anions.sort(key=lambda x: x["name"])
    
    return cations, anions

