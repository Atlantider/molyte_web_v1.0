"""
SMILES validation utility using PubChem API
"""
import requests
from typing import Dict, Optional
from app.core.logger import logger


def validate_smiles(smiles: str) -> Dict[str, any]:
    """
    Validate a SMILES string using PubChem API
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        Dictionary with validation result:
        {
            "valid": True/False,
            "name": "Molecule name" (if found),
            "iupac_name": "IUPAC name" (if found),
            "molecular_formula": "C2H6O" (if found),
            "molecular_weight": 46.07 (if found),
            "cid": 702 (PubChem CID, if found),
            "error": "Error message" (if invalid)
        }
    """
    try:
        # Step 1: Try to get compound info from PubChem by SMILES
        # Use PubChem PUG REST API
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        
        # First, try to get CID from SMILES
        cid_url = f"{base_url}/compound/smiles/{requests.utils.quote(smiles)}/cids/JSON"
        
        logger.info(f"Validating SMILES: {smiles}")
        
        response = requests.get(cid_url, timeout=5)
        
        if response.status_code == 200:
            data = response.json()
            
            if "IdentifierList" in data and "CID" in data["IdentifierList"]:
                cid = data["IdentifierList"]["CID"][0]
                
                # Step 2: Get compound properties
                props_url = f"{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,Title/JSON"
                props_response = requests.get(props_url, timeout=5)
                
                if props_response.status_code == 200:
                    props_data = props_response.json()
                    properties = props_data["PropertyTable"]["Properties"][0]
                    
                    result = {
                        "valid": True,
                        "cid": cid,
                        "molecular_formula": properties.get("MolecularFormula", ""),
                        "molecular_weight": properties.get("MolecularWeight", 0),
                        "iupac_name": properties.get("IUPACName", ""),
                        "name": properties.get("Title", ""),
                    }
                    
                    logger.info(f"SMILES validated: {smiles} -> {result['name']} (CID: {cid})")
                    return result
                else:
                    # CID found but couldn't get properties
                    return {
                        "valid": True,
                        "cid": cid,
                        "name": f"Compound CID {cid}",
                    }
        
        # If PubChem lookup failed, try basic validation
        # Check if SMILES has basic valid characters
        if _basic_smiles_check(smiles):
            logger.warning(f"SMILES not found in PubChem but passes basic check: {smiles}")
            return {
                "valid": True,
                "name": "Unknown compound",
                "note": "Not found in PubChem database, but SMILES format appears valid"
            }
        else:
            return {
                "valid": False,
                "error": "Invalid SMILES format"
            }
            
    except requests.exceptions.Timeout:
        logger.error(f"PubChem API timeout for SMILES: {smiles}")
        # Fallback to basic check
        if _basic_smiles_check(smiles):
            return {
                "valid": True,
                "name": "Unknown compound",
                "note": "PubChem API timeout, basic validation passed"
            }
        return {
            "valid": False,
            "error": "Validation timeout and basic check failed"
        }
        
    except Exception as e:
        logger.error(f"Error validating SMILES {smiles}: {e}")
        return {
            "valid": False,
            "error": f"Validation error: {str(e)}"
        }


def _basic_smiles_check(smiles: str) -> bool:
    """
    Basic SMILES format validation
    Checks for valid characters and basic structure
    """
    if not smiles or len(smiles) == 0:
        return False
    
    # Valid SMILES characters (simplified)
    # Atoms: C, N, O, S, P, F, Cl, Br, I, etc.
    # Bonds: -, =, #, :
    # Branches: ( )
    # Rings: digits
    # Aromatic: lowercase letters
    valid_chars = set("CNOSPFIBrClHcnops0123456789-=#:()[]@+\\/.%")
    
    # Check if all characters are valid
    for char in smiles:
        if char not in valid_chars:
            return False
    
    # Check balanced parentheses
    if smiles.count('(') != smiles.count(')'):
        return False
    
    if smiles.count('[') != smiles.count(']'):
        return False
    
    return True

