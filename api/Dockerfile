# ── Stage 1: compile NIRNAYA from C++17 sources ───────────────────────────────
FROM ubuntu:22.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends g++ make && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy only what the compiler needs
COPY src/       ./src/
COPY include/   ./include/
COPY tools/     ./tools/

RUN mkdir -p bin data/refinery && \
    # Compile main solver CLI
    g++ -std=c++17 -O2 -fno-fast-math -Wall \
        -Iinclude -Isrc \
        src/driver/cli_main.cpp \
        src/driver/c_api.cpp \
        -pthread \
        -o bin/nirnaya && \
    # Compile industrial model generator
    g++ -std=c++17 -O2 -fno-fast-math -Wall \
        -Iinclude -Isrc \
        tools/gen_refinery_models.cpp \
        -pthread \
        -o bin/gen_refinery_models && \
    # Pre-generate all 6 MRPL preset MPS files
    ./bin/gen_refinery_models

# ── Stage 2: lean Python runtime image ────────────────────────────────────────
FROM python:3.11-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends libstdc++6 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy compiled artifacts from builder
COPY --from=builder /build/bin/nirnaya          ./bin/nirnaya
COPY --from=builder /build/data/refinery/       ./data/refinery/

# Install Python dependencies
COPY api/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy FastAPI application
COPY api/main.py .

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
