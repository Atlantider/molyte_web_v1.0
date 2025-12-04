"""
API v1 router
"""
from fastapi import APIRouter
from app.api.v1 import auth, users, projects, electrolytes, jobs, slurm, admin, research, billing, visibility, qc, user_preferences, batch_import, worker, desolvation, binding

api_router = APIRouter()

# Include all routers
api_router.include_router(auth.router, prefix="/auth", tags=["Authentication"])
api_router.include_router(users.router, prefix="/users", tags=["Users"])
api_router.include_router(projects.router, prefix="/projects", tags=["Projects"])
api_router.include_router(electrolytes.router, prefix="/electrolytes", tags=["Electrolytes"])
api_router.include_router(jobs.router, prefix="/jobs", tags=["Jobs"])
api_router.include_router(qc.router, prefix="/qc", tags=["Quantum Chemistry"])
api_router.include_router(slurm.router, prefix="/slurm", tags=["Slurm"])
api_router.include_router(admin.router, tags=["Admin"])
api_router.include_router(research.router, prefix="/research", tags=["Research"])
api_router.include_router(billing.router, tags=["Billing"])
api_router.include_router(visibility.router, prefix="/visibility", tags=["Visibility"])
api_router.include_router(user_preferences.router, prefix="/user-preferences", tags=["User Preferences"])
api_router.include_router(batch_import.router, prefix="/batch-import", tags=["Batch Import"])
api_router.include_router(worker.router, tags=["Worker"])
api_router.include_router(desolvation.router, prefix="/desolvation", tags=["Desolvation Energy"])
api_router.include_router(binding.router, prefix="/binding", tags=["Binding Analysis"])
