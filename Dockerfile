# Public NASA CEA + collisionless plume API (Grok Build / Render / Railway).
#
# NASA CEA: `pip install cea` (official package). Wheels ship libcea*.so;
# libgfortran5 is required at runtime. gfortran is installed so a source
# build still works if pip has no wheel for this platform.
#
# Verify (also run during image build):
#   python -c "import cea; from physics.cea_rocket import CEA_VERSION; print(cea.__version__, CEA_VERSION)"

FROM python:3.12-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran \
        gcc \
        g++ \
        libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements-api.txt /app/requirements-api.txt
RUN pip install --no-cache-dir -r /app/requirements-api.txt \
    && python -c "import cea; print('NASA CEA import ok', getattr(cea, '__version__', '?'))"

COPY backend /app/backend

ENV PYTHONPATH=/app/backend
ENV PORT=8765
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Physics wrapper import (needs PYTHONPATH=backend)
RUN python -c "import cea; from physics.cea_rocket import CEA_VERSION; print('cea', cea.__version__, 'wrapper', CEA_VERSION)"

EXPOSE 8765

# Render injects PORT. Bind 0.0.0.0. Same command as native Python on Render.
CMD ["sh", "-c", "PYTHONPATH=backend python3 -m uvicorn app:app --app-dir backend --host 0.0.0.0 --port ${PORT:-8765}"]
