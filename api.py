"""
Tengri Structural Audit — FastAPI Backend
Режим 1: одиночный анализ (single)
Режим 2: сравнение (compare) — референс + образец
"""

from fastapi import FastAPI, UploadFile, File, HTTPException, Request
from fastapi.responses import JSONResponse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io, base64, csv, statistics
from datetime import datetime

app = FastAPI(title="Tengri Structural Audit API", version="2.0")

# ── Константы ─────────────────────────────────────────────────────────────────
MAX_FILE_MB    = 200
MAX_FILE_BYTES = MAX_FILE_MB * 1024 * 1024
FREE_ANALYSES  = 2       # бесплатных анализов на сессию
WEIGHTS        = {'G':1.0,'C':1.0,'A':0.5,'T':0.5,'U':0.5,'N':0.0}

# ── Простой in-memory счётчик (для продакшена заменить на Redis/DB) ────────────
usage_counter: dict[str, int] = {}

# ── FASTA парсер ──────────────────────────────────────────────────────────────
def parse_fasta(content: bytes) -> tuple[str, str]:
    text    = content.decode("utf-8", errors="ignore")
    lines   = text.strip().splitlines()
    header  = ""
    seq_parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:60]
        else:
            seq_parts.append(line.upper().replace("U","T"))
    seq = "".join(seq_parts)
    if not seq:
        raise ValueError("Последовательность пустая или файл некорректный")
    return header, seq

# ── WTA движок ───────────────────────────────────────────────────────────────
def wta_profile(seq: str, window: int = 200, step: int = 50) -> tuple[np.ndarray, np.ndarray]:
    E    = np.array([WEIGHTS.get(c, 0.0) for c in seq], dtype=np.float32)
    n    = len(E)
    kern = np.ones(window, dtype=np.float32) / window
    sm   = np.convolve(E, kern, mode='same')
    half = window // 2
    Sf   = np.empty(n, dtype=np.float32)
    Sb   = np.empty(n, dtype=np.float32)
    Sf[:n-half] = sm[half:]; Sf[n-half:] = sm[-1]
    Sb[half:]   = sm[:n-half]; Sb[:half]  = sm[half]
    D    = np.abs(Sf - Sb)
    pos, di = [], []
    for i in range(0, n - window, step):
        pos.append(i + window // 2)
        di.append(float(D[i + window // 2]))
    return np.array(pos), np.array(di)

def compute_stats(di: np.ndarray) -> dict:
    mean = float(np.mean(di))
    std  = float(np.std(di))
    return {
        "mean":      round(mean, 4),
        "std":       round(std,  4),
        "max_di":    round(float(np.max(di)), 4),
        "threshold": round(mean + 3.5 * std, 4),
        "scz_count": int(np.sum(di > mean + 3.5 * std)),
        "n_windows": len(di),
    }

def fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    plt.close(fig)
    return base64.b64encode(buf.read()).decode()

def di_to_csv_b64(pos, di) -> str:
    buf = io.StringIO()
    w   = csv.writer(buf)
    w.writerow(["window_index","nt_position","DI"])
    for i,(p,d) in enumerate(zip(pos, di)):
        w.writerow([i, int(p), round(d, 6)])
    return base64.b64encode(buf.getvalue().encode()).decode()

# ── Проверка лимитов ──────────────────────────────────────────────────────────
def check_usage(client_ip: str):
    count = usage_counter.get(client_ip, 0)
    if count >= FREE_ANALYSES:
        raise HTTPException(
            status_code=429,
            detail=f"Бесплатный лимит исчерпан ({FREE_ANALYSES} анализа). "
                   "Для продолжения свяжитесь с Tengri Lab."
        )
    usage_counter[client_ip] = count + 1

# ── ENDPOINT 1: одиночный анализ ─────────────────────────────────────────────
@app.post("/analyze")
async def analyze(request: Request, file: UploadFile = File(...)):
    client_ip = request.client.host
    check_usage(client_ip)

    raw = await file.read()
    if len(raw) > MAX_FILE_BYTES:
        raise HTTPException(400, f"Файл превышает {MAX_FILE_MB} МБ")

    try:
        header, seq = parse_fasta(raw)
    except ValueError as e:
        raise HTTPException(400, str(e))

    pos, di = wta_profile(seq)
    stats   = compute_stats(di)
    thresh  = stats["threshold"]

    # График
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
    fig.suptitle(f"Tengri TSA — {header}\n"
                 f"Длина: {len(seq):,} нт  ·  W=200 нт  ·  threshold=Mean+3.5σ",
                 fontsize=11)

    # Panel A
    ax1.plot(pos, di, color="#1565C0", lw=0.9, alpha=0.85)
    ax1.fill_between(pos, di, 0, alpha=0.1, color="#1565C0")
    ax1.axhline(stats["mean"],  color="#78909C", lw=0.8, ls=':', label=f'Mean={stats["mean"]:.4f}')
    ax1.axhline(thresh, color="#FF8F00", lw=1.2, ls='--', label=f'Threshold={thresh:.4f}')
    scz_mask = di > thresh
    if scz_mask.any():
        ax1.scatter(pos[scz_mask], di[scz_mask], c='#C62828', s=12, zorder=5,
                    label=f'SCZ ({stats["scz_count"]} zones)')
    ax1.set_title("Panel A — Full profile", loc='left', fontsize=9)
    ax1.set_xlabel("Nucleotide position (nt)")
    ax1.set_ylabel("Divergence Index (DI)")
    ax1.legend(fontsize=8)

    # Panel B — zoom на самый высокий SCZ
    peak_idx = int(np.argmax(di))
    w_size   = min(50, len(di)//4)
    lo, hi   = max(0, peak_idx-w_size), min(len(di), peak_idx+w_size)
    ax2.plot(pos[lo:hi], di[lo:hi], color="#C62828", lw=1.2)
    ax2.axhline(thresh, color="#FF8F00", lw=1.0, ls='--')
    ax2.scatter([pos[peak_idx]], [di[peak_idx]], c='#C62828', s=60, zorder=5,
                label=f'Peak DI={stats["max_di"]:.4f} @ nt {pos[peak_idx]}')
    ax2.set_title("Panel B — Peak SCZ zoom", loc='left', fontsize=9)
    ax2.set_xlabel("Nucleotide position (nt)")
    ax2.set_ylabel("DI")
    ax2.legend(fontsize=8)

    plt.tight_layout()

    return JSONResponse({
        "header":   header,
        "seq_len":  len(seq),
        "stats":    stats,
        "chart_b64": fig_to_base64(fig),
        "csv_b64":   di_to_csv_b64(pos, di),
        "doi":      "10.5281/zenodo.21203626",
        "remaining": FREE_ANALYSES - usage_counter.get(client_ip, 0),
    })


# ── ENDPOINT 2: сравнение референс vs образец ─────────────────────────────────
@app.post("/compare")
async def compare(
    request:   Request,
    reference: UploadFile = File(..., description="Референсный геном (стандарт)"),
    sample:    UploadFile = File(..., description="Ваш образец для сравнения")
):
    client_ip = request.client.host
    check_usage(client_ip)

    results = {}
    for name, upload in [("reference", reference), ("sample", sample)]:
        raw = await upload.read()
        if len(raw) > MAX_FILE_BYTES:
            raise HTTPException(400, f"{name}: файл превышает {MAX_FILE_MB} МБ")
        try:
            hdr, seq = parse_fasta(raw)
        except ValueError as e:
            raise HTTPException(400, f"{name}: {e}")
        pos, di = wta_profile(seq)
        results[name] = {"header": hdr, "seq": seq, "pos": pos, "di": di,
                         "stats": compute_stats(di)}

    ref = results["reference"]
    smp = results["sample"]

    # Нормализация по длине для совмещения на одном графике
    ref_pct = ref["pos"] / max(ref["pos"]) * 100 if len(ref["pos"]) else np.array([])
    smp_pct = smp["pos"] / max(smp["pos"]) * 100 if len(smp["pos"]) else np.array([])

    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    fig.suptitle("Tengri TSA — Comparative Analysis\n"
                 f"Reference: {ref['header'][:50]}  vs  Sample: {smp['header'][:50]}",
                 fontsize=11)

    # Panel A — референс
    ax = axes[0]
    ax.plot(ref_pct, ref["di"], color="#1565C0", lw=0.9, label="Reference (standard)")
    thresh_r = ref["stats"]["threshold"]
    ax.axhline(thresh_r, color="#FF8F00", lw=1.0, ls='--', label=f'Threshold={thresh_r:.4f}')
    ax.set_title(f"Panel A — Reference: {ref['header'][:60]}", loc='left', fontsize=9)
    ax.set_ylabel("DI")
    ax.legend(fontsize=8)

    # Panel B — образец
    ax = axes[1]
    ax.plot(smp_pct, smp["di"], color="#C62828", lw=0.9, label="Your sample")
    thresh_s = smp["stats"]["threshold"]
    ax.axhline(thresh_s, color="#FF8F00", lw=1.0, ls='--', label=f'Threshold={thresh_s:.4f}')
    scz_mask = smp["di"] > thresh_s
    if scz_mask.any():
        ax.scatter(smp_pct[scz_mask], smp["di"][scz_mask],
                   c='#C62828', s=10, zorder=5,
                   label=f'SCZ zones: {smp["stats"]["scz_count"]}')
    ax.set_title(f"Panel B — Sample: {smp['header'][:60]}", loc='left', fontsize=9)
    ax.set_ylabel("DI")
    ax.legend(fontsize=8)

    # Panel C — наложение
    ax = axes[2]
    ax.plot(ref_pct, ref["di"], color="#1565C0", lw=1.0, alpha=0.7, label="Reference")
    ax.plot(smp_pct, smp["di"], color="#C62828", lw=1.0, alpha=0.7, label="Sample")
    ax.fill_between(
        np.linspace(0,100,min(len(ref["di"]),len(smp["di"]))),
        np.interp(np.linspace(0,100,min(len(ref["di"]),len(smp["di"]))), ref_pct, ref["di"]),
        np.interp(np.linspace(0,100,min(len(ref["di"]),len(smp["di"]))), smp_pct, smp["di"]),
        alpha=0.15, color="#FF8F00", label="Divergence between profiles"
    )
    ax.set_title("Panel C — Overlay comparison", loc='left', fontsize=9)
    ax.set_xlabel("Normalized chromosome position (%)")
    ax.set_ylabel("DI")
    ax.legend(fontsize=8)

    plt.tight_layout()

    return JSONResponse({
        "reference": {"header": ref["header"], "seq_len": len(ref["seq"]),
                      "stats": ref["stats"]},
        "sample":    {"header": smp["header"], "seq_len": len(smp["seq"]),
                      "stats": smp["stats"]},
        "delta_max_di": round(float(smp["stats"]["max_di"] - ref["stats"]["max_di"]), 4),
        "delta_scz":    smp["stats"]["scz_count"] - ref["stats"]["scz_count"],
        "chart_b64":    fig_to_base64(fig),
        "csv_b64":      di_to_csv_b64(smp["pos"], smp["di"]),
        "doi":          "10.5281/zenodo.21203626",
        "remaining":    FREE_ANALYSES - usage_counter.get(client_ip, 0),
    })


@app.get("/")
def root():
    return {"service": "Tengri Structural Audit API v2",
            "docs": "/docs",
            "doi": "10.5281/zenodo.21203626"}

@app.get("/usage/{client_ip}")
def get_usage(client_ip: str):
    used = usage_counter.get(client_ip, 0)
    return {"used": used, "free_limit": FREE_ANALYSES,
            "remaining": max(0, FREE_ANALYSES - used)}
