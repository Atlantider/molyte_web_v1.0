-- Migration: Add desolvation energy calculation support
-- Date: 2025-12-03
-- Description: Add desolvation_energy_results table for storing desolvation energy calculation results

-- Step 1: Add DESOLVATION_ENERGY to postprocesstype enum
DO $$
BEGIN
    -- Check if enum value already exists
    IF NOT EXISTS (
        SELECT 1 FROM pg_enum
        WHERE enumlabel = 'DESOLVATION_ENERGY'
        AND enumtypid = (SELECT oid FROM pg_type WHERE typname = 'postprocesstype')
    ) THEN
        -- Add new enum value
        ALTER TYPE postprocesstype ADD VALUE 'DESOLVATION_ENERGY';
        RAISE NOTICE 'Added DESOLVATION_ENERGY to postprocesstype enum';
    ELSE
        RAISE NOTICE 'DESOLVATION_ENERGY already exists in postprocesstype enum';
    END IF;
END$$;

-- Step 2: Create desolvation energy results table
CREATE TABLE IF NOT EXISTS desolvation_energy_results (
    id SERIAL PRIMARY KEY,
    postprocess_job_id INTEGER NOT NULL REFERENCES postprocess_jobs(id) ON DELETE CASCADE,
    solvation_structure_id INTEGER NOT NULL REFERENCES solvation_structures(id) ON DELETE CASCADE,
    
    -- Calculation parameters
    method_level VARCHAR(50) DEFAULT 'fast_xtb',
    basis_set VARCHAR(50),
    functional VARCHAR(50),
    
    -- Complete cluster energy
    e_cluster FLOAT,
    
    -- Per-ligand results (JSON)
    per_ligand_results JSONB,
    
    -- Per-type summary (JSON)
    per_type_summary JSONB,
    
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Step 3: Create indexes
CREATE INDEX IF NOT EXISTS idx_desolvation_postprocess_job_id ON desolvation_energy_results(postprocess_job_id);
CREATE INDEX IF NOT EXISTS idx_desolvation_solvation_structure_id ON desolvation_energy_results(solvation_structure_id);

-- Step 4: Add comments
COMMENT ON TABLE desolvation_energy_results IS '去溶剂化能计算结果表';
COMMENT ON COLUMN desolvation_energy_results.e_cluster IS 'Complete cluster energy in A.U.';
COMMENT ON COLUMN desolvation_energy_results.per_ligand_results IS 'Per-ligand desolvation energy results: [{ligand_id, ligand_type, ligand_label, e_ligand, e_cluster_minus, delta_e}]';
COMMENT ON COLUMN desolvation_energy_results.per_type_summary IS 'Per-type summary statistics: [{ligand_type, avg_delta_e, std_delta_e, count, min_delta_e, max_delta_e}]';
COMMENT ON TYPE postprocesstype IS 'Postprocess job types: RDF, MSD, SOLVATION, DESOLVATION_ENERGY';

