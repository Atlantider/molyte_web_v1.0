# Molyte Web Backend

FastAPI-based backend for Molyte Web electrolyte molecular dynamics simulation platform.

## 🚀 Quick Start

### 1. Install Dependencies

```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

### 2. Configure Environment

```bash
# Copy example environment file
cp .env.example .env

# Edit .env with your settings
nano .env
```

### 3. Run Development Server

```bash
# Start server with auto-reload
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000

# Or use Python directly
python -m app.main
```

### 4. Access API Documentation

- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc
- **Health Check**: http://localhost:8000/health

---

## 📁 Project Structure

```
backend/
├── app/
│   ├── __init__.py
│   ├── main.py              # FastAPI application entry
│   ├── config.py            # Configuration management
│   ├── database.py          # Database connection
│   ├── dependencies.py      # Dependency injection
│   ├── models/              # SQLAlchemy ORM models
│   │   ├── user.py
│   │   ├── project.py
│   │   ├── electrolyte.py
│   │   ├── job.py
│   │   └── result.py
│   ├── schemas/             # Pydantic schemas
│   │   ├── user.py
│   │   ├── project.py
│   │   ├── electrolyte.py
│   │   ├── job.py
│   │   └── token.py
│   ├── api/                 # API routes
│   │   └── v1/
│   │       ├── auth.py      # Authentication
│   │       ├── users.py     # User management
│   │       ├── projects.py  # Project management
│   │       ├── electrolytes.py  # Electrolyte systems
│   │       └── jobs.py      # Job management
│   ├── core/                # Core utilities
│   │   ├── security.py      # Password & JWT
│   │   └── logger.py        # Logging
│   └── utils/               # Utility functions
│       └── hash.py          # Hash calculation
├── tests/                   # Tests
├── requirements.txt         # Python dependencies
├── .env.example             # Environment template
└── README.md                # This file
```

---

## 🔑 API Endpoints

### Authentication
- `POST /api/v1/auth/register` - Register new user
- `POST /api/v1/auth/login` - Login and get token
- `GET /api/v1/auth/me` - Get current user info

### Users
- `GET /api/v1/users/` - List all users (admin)
- `GET /api/v1/users/{id}` - Get user by ID
- `PUT /api/v1/users/{id}` - Update user
- `DELETE /api/v1/users/{id}` - Delete user (admin)

### Projects
- `POST /api/v1/projects/` - Create project
- `GET /api/v1/projects/` - List projects
- `GET /api/v1/projects/{id}` - Get project
- `PUT /api/v1/projects/{id}` - Update project
- `DELETE /api/v1/projects/{id}` - Delete project

### Electrolyte Systems
- `POST /api/v1/electrolytes/` - Create electrolyte system
- `GET /api/v1/electrolytes/` - List electrolyte systems
- `GET /api/v1/electrolytes/{id}` - Get electrolyte system

### Jobs
- `POST /api/v1/jobs/` - Create MD job
- `GET /api/v1/jobs/` - List MD jobs
- `GET /api/v1/jobs/{id}` - Get MD job
- `PUT /api/v1/jobs/{id}` - Update MD job

---

## 🔧 Configuration

Environment variables (`.env`):

```env
# Database
DATABASE_URL=postgresql://molyte_user:molyte2025@localhost:5432/molyte_web

# Security
SECRET_KEY=your-secret-key-here
ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30

# Application
APP_NAME=Molyte Web API
APP_VERSION=1.0.0
DEBUG=True

# CORS
CORS_ORIGINS=http://localhost:3000,http://localhost:5173

# Logging
LOG_LEVEL=INFO
```

---

## 🧪 Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=app

# Run specific test file
pytest tests/test_api/test_auth.py
```

---

## 📚 Technology Stack

- **FastAPI** 0.104.1 - Web framework
- **SQLAlchemy** 2.0.23 - ORM
- **PostgreSQL** 15.14 - Database
- **Pydantic** 2.5.0 - Data validation
- **JWT** - Authentication
- **Uvicorn** - ASGI server

---

## 🔐 Security

- Passwords hashed with bcrypt
- JWT tokens for authentication
- Role-based access control (admin/user)
- CORS protection
- SQL injection protection (SQLAlchemy)

---

## 📝 Development

### Code Style

```bash
# Format code
black app/

# Check code style
flake8 app/
```

### Database Migrations

```bash
# Create migration
alembic revision --autogenerate -m "description"

# Apply migrations
alembic upgrade head

# Rollback
alembic downgrade -1
```

---

## 📞 Support

For issues and questions, please refer to the main project documentation.

