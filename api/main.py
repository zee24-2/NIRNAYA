"""
NIRNAYA Sovereign Optimization Solver — FastAPI Backend
SIH26119 | Team 151198 (Visionaries for Change) | MRPL
Deployed on Railway, consumed by Next.js frontend on Vercel.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from fastapi import FastAPI, UploadFile, File, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

app = FastAPI(
    title="NIRNAYA Solver API",
    description="Sovereign Optimization Solver — SIH26119 | Team 151198",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

BASE = Path(__file__).parent
BINARY = BASE / "bin" / "nirnaya"
PRESETS_DIR = BASE / "data" / "refinery"

PRESET_MODELS = {
    "crude_blending": {
        "name": "MRPL_Crude_Blending_LP",
        "domain": "Crude Blending",
        "problem_class": "LP",
        "description": "Maximize refinery revenue subject to crude blend quality & capacity constraints.",
        "color": "blue",
    },
    "lot_sizing": {
        "name": "MRPL_MultiPeriod_LotSizing_MILP",
        "domain": "Production Planning",
        "problem_class": "MILP",
        "description": "Multi-period lot-sizing with setup costs and inventory carrying constraints.",
        "color": "purple",
    },
    "crude_scheduling": {
        "name": "MRPL_Crude_Scheduling_MILP",
        "domain": "Refinery Scheduling",
        "problem_class": "MILP",
        "description": "Binary scheduling of crude distillation units across multiple time periods.",
        "color": "violet",
    },
    "logistics": {
        "name": "MRPL_Logistics_SupplyChain_MILP",
        "domain": "Logistics & Supply Chain",
        "problem_class": "MILP",
        "description": "Facility-location MILP for MRPL depot and distribution network design.",
        "color": "indigo",
    },
    "power_dispatch": {
        "name": "MRPL_Power_Economic_Dispatch_QP",
        "domain": "Power Economic Dispatch",
        "problem_class": "Convex QP",
        "description": "Quadratic cost minimization for economic dispatch of three generating units.",
        "color": "amber",
    },
    "unit_commitment": {
        "name": "MRPL_UnitCommitment_Dispatch_MIQP",
        "domain": "Unit Commitment + Dispatch",
        "problem_class": "Convex MIQP",
        "description": "Binary unit-commitment with quadratic dispatch cost. Branch-and-Bound + Mehrotra IPM.",
        "color": "rose",
    },
}


# ── Output parser ─────────────────────────────────────────────────────────────

def _sci(line: str, key: str) -> str:
    try:
        idx = line.index(key + "=")
        s = line[idx + len(key) + 1:].split()[0].rstrip("|").rstrip(",")
        return s
    except Exception:
        return "n/a"


def parse_output(text: str) -> dict:
    result: dict = {
        "status": "unknown",
        "objective": None,
        "best_bound": None,
        "gap_pct": None,
        "iterations": None,
        "nodes": None,
        "cuts": None,
        "time_s": None,
        "verification": None,
        "model_info": {},
        "raw_output": text,
    }
    for raw in text.split("\n"):
        line = raw.strip()
        if line.startswith("Status:"):
            result["status"] = line.split("Status:")[-1].strip()
        elif line.startswith("Objective:"):
            try:
                result["objective"] = float(line.split("Objective:")[-1].strip().split()[0])
            except Exception:
                pass
        elif "Best Bound:" in line and "Gap:" in line:
            try:
                parts = line.split("|")
                result["best_bound"] = float(parts[0].split("Best Bound:")[-1].strip())
                result["gap_pct"] = float(
                    parts[1].split("Gap:")[-1].strip().rstrip("%").strip()
                )
            except Exception:
                pass
        elif "Check (dd):" in line:
            result["verification"] = {
                "passed": "PASS=YES" in line,
                "max_row_viol": _sci(line, "MaxRowViol"),
                "max_bnd_viol": _sci(line, "MaxBndViol"),
                "max_int_viol": _sci(line, "MaxIntViol"),
            }
        elif "Iterations:" in line and "Nodes:" in line:
            try:
                p = [s.strip() for s in line.split("|")]
                result["iterations"] = int(p[0].split("Iterations:")[-1].strip())
                result["nodes"] = int(p[1].split("Nodes:")[-1].strip())
                result["cuts"] = int(p[2].split("Cuts:")[-1].strip())
                result["time_s"] = float(p[3].split("Time:")[-1].strip().replace("s", "").strip())
            except Exception:
                pass
        elif "Model:" in line and "Rows:" in line:
            try:
                parts = [s.strip() for s in line.split("|")]
                result["model_info"]["name"] = parts[0].split("Model:")[-1].strip()
                for p in parts[1:]:
                    if "Rows:" in p:
                        result["model_info"]["rows"] = int(p.split("Rows:")[-1].strip())
                    elif "Cols:" in p:
                        result["model_info"]["cols"] = int(p.split("Cols:")[-1].strip())
                    elif "NNZ:" in p and "Q-NNZ:" not in p:
                        result["model_info"]["nnz"] = int(p.split("NNZ:")[-1].strip())
                    elif "Ints:" in p:
                        result["model_info"]["integers"] = int(p.split("Ints:")[-1].strip())
                    elif "Q-NNZ:" in p:
                        result["model_info"]["q_nnz"] = int(p.split("Q-NNZ:")[-1].strip())
            except Exception:
                pass
    return result


def _call_solver(mps_path: str, method: str = "auto") -> dict:
    if not BINARY.exists():
        raise HTTPException(status_code=503, detail="Solver binary not found. Railway build may have failed.")
    try:
        cmd = [str(BINARY), mps_path, f"--method={method}"]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        raw = proc.stdout + proc.stderr
        return parse_output(raw)
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=408, detail="Solver timed out (120 s limit).")
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))


# ── Routes ────────────────────────────────────────────────────────────────────

@app.get("/")
def health():
    return {
        "solver": "NIRNAYA",
        "problem_id": "SIH26119",
        "team": "151198",
        "org": "Visionaries for Change",
        "problem_owner": "MRPL",
        "status": "ok",
    }


@app.get("/models")
def list_models():
    models = []
    for key, info in PRESET_MODELS.items():
        mps = PRESETS_DIR / f"{info['name']}.mps"
        models.append(
            {
                "id": key,
                **info,
                "available": mps.exists(),
            }
        )
    return {"models": models}


@app.post("/solve/preset/{model_id}")
def solve_preset(
    model_id: str,
    method: str = Query("auto", pattern="^(auto|simplex|ipm|pdhg)$"),
):
    if model_id not in PRESET_MODELS:
        raise HTTPException(status_code=404, detail=f"Unknown model id '{model_id}'.")
    info = PRESET_MODELS[model_id]
    mps = PRESETS_DIR / f"{info['name']}.mps"
    if not mps.exists():
        raise HTTPException(status_code=404, detail=f"{info['name']}.mps not found on server.")
    result = _call_solver(str(mps), method)
    result["model_meta"] = info
    return result


@app.post("/solve/upload")
async def solve_upload(
    file: UploadFile = File(...),
    method: str = Query("auto", pattern="^(auto|simplex|ipm|pdhg)$"),
):
    if not (file.filename or "").lower().endswith(".mps"):
        raise HTTPException(status_code=400, detail="Only .mps files are accepted.")
    content = await file.read()
    if len(content) > 10 * 1024 * 1024:
        raise HTTPException(status_code=413, detail="File too large (10 MB limit).")
    with tempfile.NamedTemporaryFile(suffix=".mps", delete=False, mode="wb") as tmp:
        tmp.write(content)
        tmp_path = tmp.name
    try:
        result = _call_solver(tmp_path, method)
        result["filename"] = file.filename
        return result
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass
