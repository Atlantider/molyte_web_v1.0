"""
Utility functions for converting between old and new electrolyte formats
"""
import math
import re
from typing import List, Dict
from app.schemas.electrolyte import (
    ElectrolyteCreateNew,
    IonSpec,
    SolventSpec,
    MoleculeSpec,
    BoxConfig,
)
from app.core.logger import logger

# Avogadro's number
AVOGADRO = 6.02214076e23


def strip_ion_charge(name: str) -> str:
    """
    从离子名称中移除电荷符号（+、-）

    只移除末尾的 + 或 - 符号，保留分子式中的数字。

    例如:
        Na+ -> Na
        Li+ -> Li
        PF6- -> PF6
        TFSI- -> TFSI
        Ca2+ -> Ca2  (保留分子式中的2，只去掉+)
        Ca++ -> Ca
        SO4-- -> SO4
        BF4- -> BF4

    Args:
        name: 可能带有电荷符号的离子名称

    Returns:
        移除电荷符号后的纯净名称
    """
    if not name:
        return name
    # 只移除末尾的 + 和 - 符号（一个或多个）
    # 不移除分子式中的数字
    cleaned = re.sub(r'[+\-]+$', '', name)
    return cleaned


def calculate_box_volume(box: BoxConfig) -> float:
    """
    Calculate box volume in cubic Angstroms
    
    Args:
        box: Box configuration
        
    Returns:
        Volume in Angstrom^3
    """
    if box.type == "cubic":
        side = box.dimensions[0]
        return side ** 3
    else:  # rectangular
        return box.dimensions[0] * box.dimensions[1] * box.dimensions[2]


def convert_concentration_to_count(
    concentration: float,  # mol/L
    volume_angstrom3: float,
) -> int:
    """
    Convert concentration (mol/L) to molecule count
    
    Args:
        concentration: Concentration in mol/L
        volume_angstrom3: Volume in Angstrom^3
        
    Returns:
        Number of molecules (integer)
    """
    # Convert volume from Angstrom^3 to L
    # 1 Angstrom = 1e-10 m
    # 1 Angstrom^3 = 1e-30 m^3 = 1e-27 L
    volume_liters = volume_angstrom3 * 1e-27
    
    # Calculate moles
    moles = concentration * volume_liters
    
    # Calculate number of molecules
    count = moles * AVOGADRO
    
    # Round to nearest integer
    return max(1, round(count))


def convert_new_to_old_format(data: ElectrolyteCreateNew) -> Dict:
    """
    Convert new electrolyte format (concentration-based) to old format (count-based)
    
    Args:
        data: New format electrolyte data
        
    Returns:
        Dictionary in old format compatible with database
    """
    # Calculate box volume
    volume = calculate_box_volume(data.box)
    logger.info(f"Box volume: {volume:.2f} Angstrom^3 = {volume * 1e-27:.2e} L")
    
    # Convert ions from concentration to count
    cations_old = []
    anions_old = []
    
    first_cation_count = None
    
    for cation in data.cations:
        count = convert_concentration_to_count(cation.concentration, volume)
        if first_cation_count is None:
            first_cation_count = count

        # 清理离子名称中的电荷符号
        clean_name = strip_ion_charge(cation.name)
        cations_old.append({
            "name": clean_name,
            "smiles": f"[{clean_name}{'+'*cation.charge}]",  # Generate SMILES for ion
            "number": count,
            "concentration": cation.concentration,  # Store original concentration for editing
            "charge": cation.charge,  # Store charge for editing
        })
        logger.info(f"Cation {cation.name} -> {clean_name}: {cation.concentration} mol/L -> {count} molecules")

    for anion in data.anions:
        count = convert_concentration_to_count(anion.concentration, volume)
        # 清理离子名称中的电荷符号
        clean_name = strip_ion_charge(anion.name)
        anions_old.append({
            "name": clean_name,
            "smiles": f"[{clean_name}{'-'*abs(anion.charge)}]",  # Generate SMILES for ion
            "number": count,
            "charge": anion.charge,  # Store charge for later adjustment and editing
            "concentration": anion.concentration,  # Store original concentration for editing
        })
        logger.info(f"Anion {anion.name} -> {clean_name}: {anion.concentration} mol/L -> {count} molecules")

    # Adjust last anion count to ensure electroneutrality
    if cations_old and anions_old:
        # Calculate total positive charge
        total_positive_charge = sum(c["number"] * data.cations[i].charge for i, c in enumerate(cations_old))

        # Calculate total negative charge (before adjustment)
        total_negative_charge = sum(a["number"] * abs(a["charge"]) for a in anions_old)

        # If charges don't balance, adjust the last anion
        if total_positive_charge != total_negative_charge:
            # Calculate how much charge we need from the last anion
            charge_deficit = total_positive_charge - sum(
                a["number"] * abs(a["charge"]) for a in anions_old[:-1]
            )

            # Calculate required count for last anion
            last_anion_charge = abs(anions_old[-1]["charge"])
            required_count = max(1, round(charge_deficit / last_anion_charge))

            old_count = anions_old[-1]["number"]
            anions_old[-1]["number"] = required_count

            logger.info(
                f"Adjusted last anion {anions_old[-1]['name']} count from {old_count} to {required_count} "
                f"to ensure electroneutrality (total cation charge: +{total_positive_charge})"
            )

    # Keep charge field in anions_old for QC calculations
    # The charge field is needed when creating QC jobs from the electrolyte system

    # Convert solvents from molar ratio to count
    solvents_old = []
    
    if first_cation_count is None:
        raise ValueError("At least one cation is required")
    
    for solvent in data.solvents:
        # Count = first_cation_count * molar_ratio
        count = max(1, round(first_cation_count * solvent.molar_ratio))
        solvents_old.append({
            "name": solvent.name,
            "smiles": solvent.smiles,
            "number": count,
        })
        logger.info(f"Solvent {solvent.name}: ratio {solvent.molar_ratio} -> {count} molecules")
    
    # Calculate total molecules
    total_molecules = (
        sum(c["number"] for c in cations_old) +
        sum(a["number"] for a in anions_old) +
        sum(s["number"] for s in solvents_old)
    )
    
    logger.info(f"Total molecules: {total_molecules}")

    # Calculate box size (cubic root of volume for compatibility)
    box_size = volume ** (1/3)

    # Generate description part of the name
    user_description = data.name if data.name and data.name.strip() else None

    if not user_description:
        # Auto-generate description: 阳离子-阴离子-溶剂
        parts = []

        # Add cations (sorted by concentration, take top 2) - 使用清理后的名称
        sorted_cations = sorted(data.cations, key=lambda x: x.concentration, reverse=True)
        top_cations = [strip_ion_charge(c.name) for c in sorted_cations[:2]]
        if top_cations:
            parts.append('-'.join(top_cations))

        # Add anions (sorted by concentration, take top 2) - 使用清理后的名称
        sorted_anions = sorted(data.anions, key=lambda x: x.concentration, reverse=True)
        top_anions = [strip_ion_charge(a.name) for a in sorted_anions[:2]]
        if top_anions:
            parts.append('-'.join(top_anions))

        # Add solvents (sorted by molar ratio, take top 3)
        sorted_solvents = sorted(data.solvents, key=lambda x: x.molar_ratio, reverse=True)
        top_solvents = [s.name for s in sorted_solvents[:3]]
        if top_solvents:
            parts.append('/'.join(top_solvents))

        user_description = '-'.join(parts) if parts else 'Electrolyte'
        logger.info(f"Auto-generated description: {user_description}")
    else:
        logger.info(f"User-provided description: {user_description}")

    # Build old format dictionary
    # Note: The full name with EL-YYYYMMDD-序号 prefix will be added in the API endpoint
    # after counting existing electrolytes for the day
    result = {
        "project_id": data.project_id,
        "name": user_description,  # Store description part only, will be prefixed in API
        "cations": cations_old,
        "anions": anions_old,
        "solvents": solvents_old,
        "temperature": data.temperature,
        "pressure": data.pressure,
        "box_size": box_size,
        "nsteps_npt": data.nsteps_npt,
        "nsteps_nvt": data.nsteps_nvt,
        "timestep": data.timestep,
        "force_field": data.force_field,
    }

    return result

def convert_old_to_new_format(
    electrolyte_system,
    box_dimensions: List[float] = None
) -> Dict:
    """
    Convert old electrolyte format (count-based) to new format (concentration-based)
    for editing purposes

    Args:
        electrolyte_system: ElectrolyteSystem database object
        box_dimensions: Optional box dimensions [x, y, z] or [size] for cubic

    Returns:
        Dictionary in new format for frontend editing
    """
    # Calculate box volume from box_size (assuming cubic)
    if box_dimensions is None:
        # Assume cubic box from box_size
        box_size = electrolyte_system.box_size
        volume = box_size ** 3
        box_config = {
            "type": "cubic",
            "dimensions": [round(box_size, 2)]  # Round to 2 decimal places
        }
    else:
        if len(box_dimensions) == 1:
            volume = box_dimensions[0] ** 3
            box_config = {
                "type": "cubic",
                "dimensions": [round(box_dimensions[0], 2)]
            }
        else:
            volume = box_dimensions[0] * box_dimensions[1] * box_dimensions[2]
            box_config = {
                "type": "rectangular",
                "dimensions": [round(d, 2) for d in box_dimensions]
            }

    volume_liters = volume * 1e-27

    # Convert cations from count to concentration
    cations_new = []
    first_cation_count = None

    # Get ions info for charge lookup
    from app.utils.ion_parser import scan_available_ions
    from pathlib import Path
    salts_dir = Path("/public/home/xiaoji/molyte_web/molyte_v1/initial_salts")
    ions_info = scan_available_ions(salts_dir)

    for cation in electrolyte_system.cations:
        count = cation["number"]
        if first_cation_count is None:
            first_cation_count = count

        # Use stored concentration if available, otherwise calculate
        if "concentration" in cation and cation["concentration"] is not None:
            concentration = cation["concentration"]
        else:
            # Calculate concentration: C = n / (V * N_A)
            moles = count / AVOGADRO
            concentration = moles / volume_liters

        # 清理离子名称中的电荷符号
        clean_name = strip_ion_charge(cation["name"])

        # Get charge from stored value or lookup
        if "charge" in cation and cation["charge"] is not None:
            charge = cation["charge"]
        else:
            charge = ions_info.get(clean_name, {}).get("charge", 1)

        cations_new.append({
            "name": clean_name,
            "charge": charge,
            "concentration": round(concentration, 6)
        })

    # Convert anions from count to concentration
    anions_new = []
    for anion in electrolyte_system.anions:
        count = anion["number"]

        # Use stored concentration if available, otherwise calculate
        if "concentration" in anion and anion["concentration"] is not None:
            concentration = anion["concentration"]
        else:
            # Calculate concentration
            moles = count / AVOGADRO
            concentration = moles / volume_liters

        # 清理离子名称中的电荷符号
        clean_name = strip_ion_charge(anion["name"])

        # Get charge from stored value or lookup
        if "charge" in anion and anion["charge"] is not None:
            charge = anion["charge"]
        else:
            charge = ions_info.get(clean_name, {}).get("charge", -1)

        anions_new.append({
            "name": clean_name,
            "charge": charge,
            "concentration": round(concentration, 6)
        })

    # Convert solvents from count to molar ratio
    solvents_new = []
    if first_cation_count and electrolyte_system.solvents:
        for solvent in electrolyte_system.solvents:
            count = solvent["number"]
            # Molar ratio = solvent_count / first_cation_count
            molar_ratio = count / first_cation_count if first_cation_count > 0 else 0

            solvents_new.append({
                "name": solvent["name"],
                "smiles": solvent["smiles"],
                "molar_ratio": round(molar_ratio, 2)
            })

    result = {
        "project_id": electrolyte_system.project_id,
        "name": electrolyte_system.name,
        "description": getattr(electrolyte_system, "description", ""),
        "cations": cations_new,
        "anions": anions_new,
        "solvents": solvents_new,
        "box": box_config,
        "temperature": electrolyte_system.temperature,
        "pressure": electrolyte_system.pressure,
        "nsteps_npt": electrolyte_system.nsteps_npt,
        "nsteps_nvt": electrolyte_system.nsteps_nvt,
        "timestep": electrolyte_system.timestep,
        "force_field": electrolyte_system.force_field,
    }

    logger.info(f"Converted electrolyte {electrolyte_system.id} from old to new format")

    return result


