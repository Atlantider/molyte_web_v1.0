-- Migration: Add system_structures table for storing system structure XYZ content
-- Purpose: Store the last frame of system structure to avoid reading large trajectory files
-- Date: 2025-12-03

CREATE TABLE IF NOT EXISTS system_structures (
    id SERIAL PRIMARY KEY,
    md_job_id INTEGER NOT NULL UNIQUE REFERENCES md_jobs(id) ON DELETE CASCADE,
    
    -- Frame information
    frame_index INTEGER,
    total_frames INTEGER,
    atom_count INTEGER,
    
    -- Box information (JSON array: [lx, ly, lz])
    box JSONB,
    
    -- XYZ content (for 3D visualization)
    xyz_content TEXT,
    
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    
    CONSTRAINT fk_system_structures_md_job_id FOREIGN KEY (md_job_id) REFERENCES md_jobs(id) ON DELETE CASCADE
);

-- Create index for faster lookups
CREATE INDEX IF NOT EXISTS idx_system_structure_md_job_id ON system_structures(md_job_id);

-- Add comment
COMMENT ON TABLE system_structures IS 'System structure model - stores the last frame XYZ content for 3D visualization';
COMMENT ON COLUMN system_structures.md_job_id IS 'Reference to MD job';
COMMENT ON COLUMN system_structures.frame_index IS 'Frame index';
COMMENT ON COLUMN system_structures.total_frames IS 'Total number of frames';
COMMENT ON COLUMN system_structures.atom_count IS 'Number of atoms';
COMMENT ON COLUMN system_structures.box IS 'Box dimensions [lx, ly, lz]';
COMMENT ON COLUMN system_structures.xyz_content IS 'XYZ format content for 3D visualization';

