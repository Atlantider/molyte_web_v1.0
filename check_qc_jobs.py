#!/usr/bin/env python3
"""检查QC任务的详细信息"""
import sys
sys.path.insert(0, '/public/home/xiaoji/molyte_web/backend')

from app.database import SessionLocal
from app.models.qc import QCJob

db = SessionLocal()

try:
    # 查找所有任务
    all_jobs = db.query(QCJob).order_by(QCJob.id.desc()).limit(20).all()
    print(f"Latest 20 jobs:")
    for job in all_jobs:
        print(f"  ID={job.id}, Name={job.molecule_name}, Status={job.status}")

    # 查找特定任务
    jobs = db.query(QCJob).filter(QCJob.id.in_([157, 159])).all()

    print(f"\nFound {len(jobs)} jobs with ID 157 or 159")

    for job in jobs:
        print(f"\n{'='*80}")
        print(f"Job ID: {job.id}")
        print(f"Molecule Name: {job.molecule_name}")
        print(f"SMILES: {job.smiles}")
        print(f"Charge: {job.charge}")
        print(f"Spin Multiplicity: {job.spin_multiplicity}")
        print(f"Status: {job.status}")
        print(f"Work Dir: {job.work_dir}")
        print(f"Error Message: {job.error_message}")
        print(f"Functional: {job.functional}")
        print(f"Basis Set: {job.basis_set}")
        print(f"Molecule Type: {job.molecule_type}")

    # 测试从SMILES获取电荷和自旋
    print(f"\n{'='*80}")
    print("Testing SMILES parsing:")

    from rdkit import Chem

    test_smiles = [
        ("FSI-", "O=S(=O)(NS(=O)(F)=O)F"),
        ("TTE", "FC(F)(F)COC(C(F)(F)F)C(F)(F)F"),
    ]

    for name, smiles in test_smiles:
        print(f"\n{name}: {smiles}")
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            unpaired_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
            spin_multiplicity = unpaired_electrons + 1
            num_electrons = sum(atom.GetAtomicNum() - atom.GetFormalCharge() for atom in mol.GetAtoms())

            print(f"  Total Charge: {total_charge}")
            print(f"  Unpaired Electrons: {unpaired_electrons}")
            print(f"  Spin Multiplicity: {spin_multiplicity}")
            print(f"  Total Electrons: {num_electrons}")
        else:
            print(f"  Failed to parse SMILES")

finally:
    db.close()

