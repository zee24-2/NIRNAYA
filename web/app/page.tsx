"use client";

import { useState, useRef, useCallback } from "react";
import {
  CheckCircle,
  XCircle,
  Loader2,
  Upload,
  ChevronDown,
  Terminal,
  Cpu,
  Zap,
  Activity,
  BarChart3,
  AlertTriangle,
} from "lucide-react";

// ── Types ─────────────────────────────────────────────────────────────────────

type Method = "auto" | "simplex" | "ipm" | "pdhg";

interface Verification {
  passed: boolean;
  max_row_viol: string;
  max_bnd_viol: string;
  max_int_viol: string;
}

interface SolveResult {
  status: string;
  objective: number | null;
  best_bound: number | null;
  gap_pct: number | null;
  iterations: number | null;
  nodes: number | null;
  cuts: number | null;
  time_s: number | null;
  verification: Verification | null;
  model_info: Record<string, string | number>;
  model_meta?: Record<string, string>;
  filename?: string;
  raw_output: string;
}

// ── Preset Model Definitions ──────────────────────────────────────────────────

const MODELS = [
  {
    id: "crude_blending",
    label: "Crude Blending",
    subtitle: "MRPL_Crude_Blending_LP",
    badge: "LP",
    domain: "Crude Blending",
    desc: "Maximize refinery revenue subject to crude blend quality & capacity constraints.",
    accent: "from-blue-500/20 to-cyan-500/10 border-blue-500/30",
    badgeStyle: "bg-blue-500/20 text-blue-300 ring-1 ring-blue-500/40",
    icon: "🛢️",
  },
  {
    id: "lot_sizing",
    label: "Production Planning",
    subtitle: "MRPL_MultiPeriod_LotSizing_MILP",
    badge: "MILP",
    domain: "Production Planning",
    desc: "Multi-period lot-sizing with setup costs and inventory carrying constraints.",
    accent: "from-purple-500/20 to-violet-500/10 border-purple-500/30",
    badgeStyle: "bg-purple-500/20 text-purple-300 ring-1 ring-purple-500/40",
    icon: "🏭",
  },
  {
    id: "crude_scheduling",
    label: "Refinery Scheduling",
    subtitle: "MRPL_Crude_Scheduling_MILP",
    badge: "MILP",
    domain: "Refinery Scheduling",
    desc: "Binary scheduling of crude distillation units across multiple time periods.",
    accent: "from-violet-500/20 to-purple-500/10 border-violet-500/30",
    badgeStyle: "bg-violet-500/20 text-violet-300 ring-1 ring-violet-500/40",
    icon: "📅",
  },
  {
    id: "logistics",
    label: "Logistics & Supply Chain",
    subtitle: "MRPL_Logistics_SupplyChain_MILP",
    badge: "MILP",
    domain: "Logistics",
    desc: "Facility-location MILP for MRPL depot and distribution network design.",
    accent: "from-indigo-500/20 to-blue-500/10 border-indigo-500/30",
    badgeStyle: "bg-indigo-500/20 text-indigo-300 ring-1 ring-indigo-500/40",
    icon: "🚚",
  },
  {
    id: "power_dispatch",
    label: "Power Economic Dispatch",
    subtitle: "MRPL_Power_Economic_Dispatch_QP",
    badge: "Convex QP",
    domain: "Power Dispatch",
    desc: "Quadratic cost minimization for economic dispatch of three generating units.",
    accent: "from-amber-500/20 to-yellow-500/10 border-amber-500/30",
    badgeStyle: "bg-amber-500/20 text-amber-300 ring-1 ring-amber-500/40",
    icon: "⚡",
  },
  {
    id: "unit_commitment",
    label: "Unit Commitment MIQP",
    subtitle: "MRPL_UnitCommitment_Dispatch_MIQP",
    badge: "Convex MIQP",
    domain: "Unit Commitment",
    desc: "Binary commitment with quadratic dispatch cost. Branch-and-Bound + Mehrotra IPM.",
    accent: "from-rose-500/20 to-pink-500/10 border-rose-500/30",
    badgeStyle: "bg-rose-500/20 text-rose-300 ring-1 ring-rose-500/40",
    icon: "🔌",
  },
];

// ── Helpers ───────────────────────────────────────────────────────────────────

function fmt(n: number | null, decimals = 6): string {
  if (n === null || n === undefined) return "—";
  if (Number.isInteger(n)) return n.toString();
  return n.toFixed(decimals).replace(/\.?0+$/, "");
}

function statusStyle(s: string) {
  const l = s.toLowerCase();
  if (l === "optimal")
    return { dot: "bg-emerald-400", text: "text-emerald-400", badge: "bg-emerald-500/15 text-emerald-300 ring-1 ring-emerald-500/40" };
  if (l.includes("infeasible"))
    return { dot: "bg-red-400", text: "text-red-400", badge: "bg-red-500/15 text-red-300 ring-1 ring-red-500/40" };
  if (l.includes("unbounded"))
    return { dot: "bg-amber-400", text: "text-amber-400", badge: "bg-amber-500/15 text-amber-300 ring-1 ring-amber-500/40" };
  return { dot: "bg-slate-400", text: "text-slate-400", badge: "bg-slate-500/15 text-slate-300 ring-1 ring-slate-500/40" };
}

// ── Main Component ────────────────────────────────────────────────────────────

export default function Home() {
  const [method, setMethod] = useState<Method>("auto");
  const [loading, setLoading] = useState<string | null>(null); // model id or "upload"
  const [result, setResult] = useState<SolveResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [showRaw, setShowRaw] = useState(false);
  const [dragOver, setDragOver] = useState(false);
  const [activeTab, setActiveTab] = useState<"presets" | "upload">("presets");
  const fileRef = useRef<HTMLInputElement>(null);
  const resultsRef = useRef<HTMLDivElement>(null);

  const callApi = useCallback(
    async (path: string, init?: RequestInit): Promise<SolveResult> => {
      const res = await fetch(path, init);
      if (!res.ok) {
        const body = await res.json().catch(() => ({}));
        throw new Error(body.detail || `HTTP ${res.status}`);
      }
      return res.json();
    },
    []
  );

  const solvePreset = async (id: string) => {
    setLoading(id);
    setError(null);
    setResult(null);
    setShowRaw(false);
    try {
      const data = await callApi(`/api/solve/preset/${id}?method=${method}`);
      setResult(data);
      setTimeout(() => resultsRef.current?.scrollIntoView({ behavior: "smooth" }), 100);
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : "Unknown error");
    } finally {
      setLoading(null);
    }
  };

  const solveFile = async (file: File) => {
    setLoading("upload");
    setError(null);
    setResult(null);
    setShowRaw(false);
    const fd = new FormData();
    fd.append("file", file);
    try {
      const data = await callApi(`/api/solve/upload?method=${method}`, {
        method: "POST",
        body: fd,
      });
      setResult(data);
      setTimeout(() => resultsRef.current?.scrollIntoView({ behavior: "smooth" }), 100);
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : "Unknown error");
    } finally {
      setLoading(null);
    }
  };

  const onDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setDragOver(false);
    const file = e.dataTransfer.files?.[0];
    if (file) solveFile(file);
  };

  // ── Render ──────────────────────────────────────────────────────────────────
  return (
    <div className="min-h-screen bg-[#07070f] text-slate-200">

      {/* ── Header ─────────────────────────────────────────────────────────── */}
      <header className="border-b border-[#1c1c2e] bg-[#07070f]/90 backdrop-blur-md sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-6 py-4 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-emerald-400 to-teal-600 flex items-center justify-center shadow-lg shadow-emerald-500/20">
              <Zap className="w-4 h-4 text-white" />
            </div>
            <div>
              <h1 className="text-lg font-semibold tracking-tight bg-gradient-to-r from-emerald-400 to-teal-400 bg-clip-text text-transparent">
                NIRNAYA
              </h1>
              <p className="text-[10px] text-slate-500 leading-none font-mono">
                Sovereign Optimization Solver Core
              </p>
            </div>
          </div>
          <div className="hidden sm:flex items-center gap-6 text-xs text-slate-500">
            <span className="font-mono">SIH26119</span>
            <span className="w-px h-3 bg-slate-700" />
            <span>Team 151198 — <span className="text-slate-400">Visionaries for Change</span></span>
            <span className="w-px h-3 bg-slate-700" />
            <span className="text-emerald-500 font-medium">MRPL</span>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-6 py-10 space-y-10">

        {/* ── Method + Tabs Row ───────────────────────────────────────────── */}
        <div className="flex flex-col sm:flex-row gap-4 items-start sm:items-center justify-between">
          {/* Tabs */}
          <div className="flex gap-1 p-1 bg-[#0f0f1a] rounded-xl border border-[#1c1c2e]">
            {(["presets", "upload"] as const).map((tab) => (
              <button
                key={tab}
                onClick={() => setActiveTab(tab)}
                className={`px-4 py-2 text-sm font-medium rounded-lg transition-all ${
                  activeTab === tab
                    ? "bg-[#1c1c2e] text-white shadow"
                    : "text-slate-500 hover:text-slate-300"
                }`}
              >
                {tab === "presets" ? "🏭  MRPL Preset Models" : "📂  Upload MPS"}
              </button>
            ))}
          </div>

          {/* Method Selector */}
          <div className="flex items-center gap-2">
            <span className="text-xs text-slate-500">Algorithm:</span>
            <div className="relative">
              <select
                value={method}
                onChange={(e) => setMethod(e.target.value as Method)}
                className="appearance-none bg-[#0f0f1a] border border-[#1c1c2e] text-slate-200 text-sm rounded-lg px-3 py-2 pr-8 focus:outline-none focus:ring-1 focus:ring-emerald-500/50 cursor-pointer"
              >
                <option value="auto">⚡ Auto</option>
                <option value="simplex">📐 Dual Simplex</option>
                <option value="ipm">🔬 Mehrotra IPM</option>
                <option value="pdhg">🎯 PDHG (First-Order)</option>
              </select>
              <ChevronDown className="absolute right-2 top-2.5 w-4 h-4 text-slate-500 pointer-events-none" />
            </div>
          </div>
        </div>

        {/* ── Preset Models Grid ──────────────────────────────────────────── */}
        {activeTab === "presets" && (
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
            {MODELS.map((m) => {
              const isLoading = loading === m.id;
              return (
                <div
                  key={m.id}
                  className={`group relative rounded-2xl border bg-gradient-to-br p-5 transition-all hover:scale-[1.01] hover:shadow-xl hover:shadow-black/40 ${m.accent}`}
                >
                  <div className="flex items-start justify-between mb-3">
                    <div className="flex items-center gap-2">
                      <span className="text-2xl">{m.icon}</span>
                      <span className={`text-xs font-semibold px-2 py-0.5 rounded-full font-mono ${m.badgeStyle}`}>
                        {m.badge}
                      </span>
                    </div>
                  </div>
                  <h3 className="font-semibold text-sm text-white mb-1">{m.label}</h3>
                  <p className="text-xs text-slate-400 mb-1 font-mono truncate">{m.subtitle}</p>
                  <p className="text-xs text-slate-500 mb-4 leading-relaxed">{m.desc}</p>
                  <button
                    onClick={() => solvePreset(m.id)}
                    disabled={!!loading}
                    className="w-full flex items-center justify-center gap-2 py-2 px-4 rounded-lg text-sm font-medium bg-white/5 hover:bg-white/10 border border-white/10 hover:border-white/20 transition-all disabled:opacity-40 disabled:cursor-not-allowed"
                  >
                    {isLoading ? (
                      <>
                        <Loader2 className="w-4 h-4 loading-ring" />
                        Solving…
                      </>
                    ) : (
                      <>
                        <Activity className="w-4 h-4" />
                        Solve
                      </>
                    )}
                  </button>
                </div>
              );
            })}
          </div>
        )}

        {/* ── Upload Tab ──────────────────────────────────────────────────── */}
        {activeTab === "upload" && (
          <div
            onDragOver={(e) => { e.preventDefault(); setDragOver(true); }}
            onDragLeave={() => setDragOver(false)}
            onDrop={onDrop}
            onClick={() => fileRef.current?.click()}
            className={`relative rounded-2xl border-2 border-dashed p-16 flex flex-col items-center justify-center gap-4 cursor-pointer transition-all ${
              dragOver
                ? "border-emerald-500/60 bg-emerald-500/5"
                : "border-[#2d2d4a] bg-[#0f0f1a] hover:border-emerald-500/30 hover:bg-emerald-500/5"
            }`}
          >
            <input
              ref={fileRef}
              type="file"
              accept=".mps"
              className="hidden"
              onChange={(e) => { const f = e.target.files?.[0]; if (f) solveFile(f); }}
            />
            {loading === "upload" ? (
              <Loader2 className="w-10 h-10 text-emerald-400 loading-ring" />
            ) : (
              <Upload className={`w-10 h-10 transition-colors ${dragOver ? "text-emerald-400" : "text-slate-600"}`} />
            )}
            <div className="text-center">
              <p className="font-medium text-slate-300">
                {loading === "upload" ? "Solving your model…" : "Drop your .mps file here"}
              </p>
              <p className="text-sm text-slate-600 mt-1">or click to browse — max 10 MB</p>
            </div>
          </div>
        )}

        {/* ── Error Banner ────────────────────────────────────────────────── */}
        {error && (
          <div className="flex items-start gap-3 p-4 rounded-xl bg-red-500/10 border border-red-500/30 text-red-300">
            <AlertTriangle className="w-5 h-5 shrink-0 mt-0.5" />
            <div>
              <p className="font-semibold text-sm">Solver Error</p>
              <p className="text-xs text-red-400 mt-0.5">{error}</p>
            </div>
          </div>
        )}

        {/* ── Results Panel ───────────────────────────────────────────────── */}
        {result && (
          <div ref={resultsRef} className="rounded-2xl border border-[#1c1c2e] bg-[#0f0f1a] overflow-hidden">

            {/* Results Header */}
            <div className="flex flex-col sm:flex-row items-start sm:items-center justify-between gap-4 p-6 border-b border-[#1c1c2e]">
              <div className="flex items-center gap-3">
                <div className="w-9 h-9 rounded-xl bg-gradient-to-br from-emerald-500/20 to-teal-500/20 border border-emerald-500/20 flex items-center justify-center">
                  <BarChart3 className="w-5 h-5 text-emerald-400" />
                </div>
                <div>
                  <h2 className="font-semibold text-white">
                    {result.model_info?.name as string || result.filename || "Solver Result"}
                  </h2>
                  <p className="text-xs text-slate-500 font-mono mt-0.5">
                    {result.model_meta?.domain || result.model_info?.name as string}
                  </p>
                </div>
              </div>

              {/* Status Badge */}
              <div className={`flex items-center gap-2 px-3 py-1.5 rounded-full text-sm font-semibold ${statusStyle(result.status).badge}`}>
                <div className={`w-2 h-2 rounded-full ${statusStyle(result.status).dot} ${result.status === "optimal" ? "glow-green" : ""}`} />
                {result.status.toUpperCase()}
              </div>
            </div>

            {/* Model Info Pills */}
            {Object.keys(result.model_info).length > 1 && (
              <div className="px-6 pt-4 flex flex-wrap gap-2">
                {result.model_info.rows !== undefined && (
                  <span className="text-xs font-mono px-2 py-1 rounded-md bg-[#1c1c2e] text-slate-400">
                    {result.model_info.rows} rows
                  </span>
                )}
                {result.model_info.cols !== undefined && (
                  <span className="text-xs font-mono px-2 py-1 rounded-md bg-[#1c1c2e] text-slate-400">
                    {result.model_info.cols} cols
                  </span>
                )}
                {result.model_info.integers !== undefined && (
                  <span className="text-xs font-mono px-2 py-1 rounded-md bg-purple-500/10 text-purple-400 ring-1 ring-purple-500/20">
                    {result.model_info.integers} integers
                  </span>
                )}
                {result.model_info.q_nnz !== undefined && Number(result.model_info.q_nnz) > 0 && (
                  <span className="text-xs font-mono px-2 py-1 rounded-md bg-amber-500/10 text-amber-400 ring-1 ring-amber-500/20">
                    {result.model_info.q_nnz} Q-nnz
                  </span>
                )}
                <span className="text-xs font-mono px-2 py-1 rounded-md bg-[#1c1c2e] text-slate-400">
                  method: {method}
                </span>
              </div>
            )}

            {/* Objective */}
            <div className="p-6">
              {result.objective !== null && (
                <div className="mb-6">
                  <p className="text-xs text-slate-500 uppercase tracking-widest mb-1">Optimal Objective</p>
                  <p className="text-3xl font-semibold font-mono text-white">
                    {typeof result.objective === "number" ? result.objective.toLocaleString("en-US", { maximumFractionDigits: 10 }) : "—"}
                  </p>
                  {result.gap_pct !== null && result.gap_pct !== undefined && (
                    <p className="text-xs text-slate-500 mt-1 font-mono">
                      Gap: {result.gap_pct.toFixed(4)}% &nbsp;|&nbsp; Bound: {fmt(result.best_bound, 4)}
                    </p>
                  )}
                </div>
              )}

              {/* Stats Grid */}
              <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 mb-6">
                {[
                  { label: "Iterations", value: result.iterations?.toLocaleString() ?? "—", icon: <Activity className="w-4 h-4" /> },
                  { label: "B&B Nodes", value: result.nodes?.toLocaleString() ?? "—", icon: <Cpu className="w-4 h-4" /> },
                  { label: "Cuts Added", value: result.cuts?.toLocaleString() ?? "—", icon: <Zap className="w-4 h-4" /> },
                  { label: "Wall Time", value: result.time_s !== null && result.time_s !== undefined ? `${result.time_s.toFixed(4)} s` : "—", icon: <Terminal className="w-4 h-4" /> },
                ].map((stat) => (
                  <div key={stat.label} className="rounded-xl bg-[#1c1c2e]/60 border border-[#2d2d4a] p-4">
                    <div className="flex items-center gap-1.5 text-slate-500 mb-2">
                      {stat.icon}
                      <span className="text-[11px] uppercase tracking-wide">{stat.label}</span>
                    </div>
                    <p className="text-xl font-semibold font-mono text-slate-200">{stat.value}</p>
                  </div>
                ))}
              </div>

              {/* Verification */}
              {result.verification && (
                <div className={`rounded-xl border p-4 flex flex-col sm:flex-row items-start sm:items-center gap-4 ${
                  result.verification.passed
                    ? "border-emerald-500/30 bg-emerald-500/5"
                    : "border-red-500/30 bg-red-500/5"
                }`}>
                  <div className="flex items-center gap-2">
                    {result.verification.passed ? (
                      <CheckCircle className="w-6 h-6 text-emerald-400 shrink-0" />
                    ) : (
                      <XCircle className="w-6 h-6 text-red-400 shrink-0" />
                    )}
                    <div>
                      <p className={`font-semibold text-sm ${result.verification.passed ? "text-emerald-300" : "text-red-300"}`}>
                        {result.verification.passed ? "dd_real Verification PASS" : "Verification FAIL"}
                      </p>
                      <p className="text-xs text-slate-500">Independent ~32-digit certificate check</p>
                    </div>
                  </div>
                  <div className="flex gap-4 text-xs font-mono text-slate-500 sm:ml-auto">
                    <span>Row Viol: <span className="text-slate-300">{result.verification.max_row_viol}</span></span>
                    <span>Bnd Viol: <span className="text-slate-300">{result.verification.max_bnd_viol}</span></span>
                    <span>Int Viol: <span className="text-slate-300">{result.verification.max_int_viol}</span></span>
                  </div>
                </div>
              )}
            </div>

            {/* Raw Output Collapsible */}
            <div className="border-t border-[#1c1c2e]">
              <button
                onClick={() => setShowRaw((v) => !v)}
                className="w-full flex items-center justify-between px-6 py-3 text-xs text-slate-500 hover:text-slate-300 hover:bg-[#1c1c2e]/30 transition-colors"
              >
                <span className="flex items-center gap-2">
                  <Terminal className="w-3.5 h-3.5" />
                  Raw solver output
                </span>
                <ChevronDown className={`w-4 h-4 transition-transform ${showRaw ? "rotate-180" : ""}`} />
              </button>
              {showRaw && (
                <pre className="px-6 pb-6 text-[11px] font-mono text-slate-400 whitespace-pre-wrap leading-relaxed overflow-x-auto max-h-96 overflow-y-auto">
                  {result.raw_output}
                </pre>
              )}
            </div>
          </div>
        )}

        {/* ── Footer ──────────────────────────────────────────────────────── */}
        <footer className="text-center text-xs text-slate-600 py-4 border-t border-[#1c1c2e] space-y-1">
          <p>
            <span className="text-emerald-600 font-semibold">NIRNAYA</span> Sovereign Optimization Solver Core — Problem <span className="font-mono">SIH26119</span>
          </p>
          <p>Team 151198 &mdash; <em>Visionaries for Change</em> &mdash; Problem Owner: MRPL</p>
          <p className="text-slate-700">
            C++17 · Dual Simplex · Mehrotra Predictor-Corrector IPM · Parallel PDHG · Branch-and-Cut · dd_real ~32-digit verification · Zero external solver dependencies
          </p>
        </footer>
      </main>
    </div>
  );
}
