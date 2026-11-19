#!/usr/bin/env python3
"""
tsa_cli.py — Tengri Structural Audit из командной строки (на ядре tengri_core_wta v2).
"""
import argparse
import csv
import json
import os
import re

import numpy as np

from tengri_core_wta import __version__, detect_scz, read_fasta, scz_zones, wta_di

DOI = "10.5281/zenodo.21203626"
NOTE = ("DI = 0.5·|GC_forward − GC_backward| между соседними окнами (веса G=C=1.0, A=T=0.5). "
        "SCZ — статистические выбросы этого показателя, а не доказанные функциональные элементы.")


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name.split()[0] if name.split() else name)[:60] or "record"


def decimate(x, y, max_points=150_000):
    if len(y) <= max_points:
        return x, y
    k = int(np.ceil(len(y) / max_points))
    pad = (-len(y)) % k
    yy = np.concatenate([y, np.full(pad, np.nan)]).reshape(-1, k)
    xx = np.concatenate([x, np.full(pad, x[-1])]).reshape(-1, k)
    with np.errstate(all="ignore"):
        return xx[:, 0], np.nanmax(yy, axis=1)


def plot(path, header, n, pos, di, thr, zones, a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    xs, ys = decimate(pos, di)
    ax1.plot(xs, ys, color="#1f5fa8", lw=0.7)
    ax1.axhline(thr, color="#e08a00", ls="--", lw=1.1, label=f"Mean+{a.sigma}σ = {thr:.4f}")
    if zones:
        ax1.scatter([z["critical_node_nt"] for z in zones], [z["peak_DI"] for z in zones],
                    c="#c62828", s=18, zorder=5, label=f"Critical Node ({len(zones)} zones)")
    ax1.set_title(f"TSA Profile: {header[:80]} ({n:,} nt)")
    ax1.set_xlabel("Position (nt)")
    ax1.set_ylabel("Divergence Index (DI)")
    ax1.legend(loc="upper right", fontsize=8)

    if zones:
        top = max(zones, key=lambda z: z["peak_DI"])
        lo = top["critical_node_nt"] - 20 * a.step
        hi = top["critical_node_nt"] + 20 * a.step
        m = (pos >= lo) & (pos <= hi)
        ax2.plot(pos[m], di[m], color="#1f5fa8", marker="o", ms=3, lw=1)
        ax2.axhline(thr, color="#e08a00", ls="--", lw=1.1)
        ax2.scatter([top["critical_node_nt"]], [top["peak_DI"]], c="#c62828", s=50, zorder=5)
        ax2.set_title(f"Zoom: {top['scz_id']} (peak DI = {top['peak_DI']:.4f} @ nt {top['critical_node_nt']:,})",
                      fontsize=10)
    else:
        ax2.text(0.5, 0.5, "No SCZ detected", ha="center", va="center", transform=ax2.transAxes)
    ax2.set_xlabel("Position (nt)")
    ax2.set_ylabel("DI")
    fig.text(0.5, 0.005,
             f"WTA, W={a.window} nt, step={a.step} nt, threshold=Mean+{a.sigma}σ, N: {a.n_mode} | "
             f"DI = 0.5·|ΔGC| | core v{__version__} | DOI: {DOI}", ha="center", fontsize=7.5)
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_record(header, seq, a, outdir):
    rid = safe(header)
    n = len(seq)
    arr = np.frombuffer(seq, dtype=np.uint8)
    n_frac = float(((arr | 0x20) == ord("n")).sum()) / n
    pos, di = wta_di(seq, a.window, a.step, n_mode=a.n_mode, max_ambig=a.max_ambig)
    is_scz, mean, std, thr = detect_scz(di, a.sigma)
    zones = scz_zones(pos, di, is_scz, a.window)

    scz_id = np.full(len(di), "", dtype=object)
    for z in zones:
        sel = is_scz & (pos >= z["nt_start"]) & (pos <= z["nt_end"] - a.window)
        scz_id[sel] = z["scz_id"]

    with open(os.path.join(outdir, f"{rid}_windows.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["window_index", "nt_start", "nt_end", "DI", "is_SCZ", "scz_id", "annotation"])
        for k in range(len(di)):
            w.writerow([k, int(pos[k]), int(pos[k]) + a.window,
                        "" if not np.isfinite(di[k]) else f"{di[k]:.6f}",
                        int(is_scz[k]), scz_id[k] or "—", "N/A"])

    warnings = []
    if a.n_mode == "zero" and n_frac > 0:
        warnings.append("n_mode=zero: N=0.0 даёт DI≈0.7 на каждой границе N-блока — эти пики артефакт кодировки.")
    masked = int((~np.isfinite(di)).sum())
    if masked:
        warnings.append(f"{masked} окон замаскированы из-за N/неоднозначных оснований (max_ambig={a.max_ambig}).")
    summary = {
        "record": header, "length_nt": n, "N_fraction": round(n_frac, 6),
        "window": a.window, "step": a.step, "sigma": a.sigma, "n_mode": a.n_mode,
        "n_windows": int(len(di)), "n_windows_masked": masked,
        "mean_DI": round(mean, 6), "std_DI": round(std, 6), "threshold": round(thr, 6),
        "n_scz_windows": int(is_scz.sum()), "n_scz_zones": len(zones), "zones": zones,
        "warnings": warnings, "note": NOTE, "core_version": __version__, "doi": DOI,
    }
    with open(os.path.join(outdir, f"{rid}_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    plot(os.path.join(outdir, f"{rid}_profile.png"), header, n, pos, di, thr, zones, a)
    print(f"{rid}: {n:,} nt | окон {len(di):,} (маскировано {masked:,}) | порог {thr:.4f} | "
          f"зон SCZ: {len(zones)}" + (f" | ⚠ {warnings[0]}" if warnings and a.n_mode == "zero" else ""))


def main():
    ap = argparse.ArgumentParser(description="Tengri Structural Audit (WTA)")
    ap.add_argument("fasta")
    ap.add_argument("--out", default="tsa_report")
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--step", type=int, default=50)
    ap.add_argument("--sigma", type=float, default=3.5)
    ap.add_argument("--n-mode", dest="n_mode", choices=["mask", "zero"], default="mask")
    ap.add_argument("--max-ambig", dest="max_ambig", type=float, default=0.10)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    found = False
    for header, seq in read_fasta(a.fasta):
        found = True
        if len(seq) < 2 * a.window + a.step:
            print(f"пропуск {safe(header)}: короче 2·window")
            continue
        run_record(header, seq, a, a.out)
    if not found:
        raise SystemExit("в файле нет последовательностей FASTA")


if __name__ == "__main__":
    main()