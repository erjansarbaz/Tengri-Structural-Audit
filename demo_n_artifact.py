"""
demo_n_artifact.py — может ли кодировка N=0.0 сама создавать «обогащение SCZ у N-блоков»?
"""
import argparse
import math
import numpy as np
from tengri_core_wta import wta_di, detect_scz, scz_zones, read_fasta

W, STEP, PROX, SIGMA, FLANK_MAX = 200, 50, 50, 3.5, 5000
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)


def n_runs(seq: bytes):
    isn = (np.frombuffer(seq, dtype=np.uint8) | 0x20) == ord("n")
    d = np.diff(isn.astype(np.int8), prepend=0, append=0)
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def dist_to_n(L, R, starts, ends):
    L = np.asarray(L, dtype=np.int64)
    R = np.asarray(R, dtype=np.int64)
    big = np.iinfo(np.int64).max
    if len(starts) == 0:
        return np.full(len(L), big, dtype=np.int64)
    idx = np.searchsorted(starts, L, side="right") - 1
    ip = np.clip(idx, 0, None)
    d_prev = np.where(idx >= 0, np.maximum(0, L - (ends[ip] - 1)), big)
    nxt = idx + 1
    inx = np.clip(nxt, 0, len(starts) - 1)
    d_next = np.where(nxt < len(starts), np.maximum(0, starts[inx] - (R - 1)), big)
    return np.minimum(d_prev, d_next)


def z_prop(hits, m, p0):
    if m == 0 or p0 <= 0 or p0 >= 1:
        return float("nan")
    return (hits / m - p0) / math.sqrt(p0 * (1 - p0) / m)


def analyze(seq: bytes):
    n = len(seq)
    starts, ends = n_runs(seq)
    p0 = float((ends - starts).sum()) / n

    # A) как в препринтах
    pos, di = wta_di(seq, W, STEP, n_mode="zero")
    is_scz, _, _, thr = detect_scz(di, SIGMA)
    zones = scz_zones(pos, di, is_scz, W, max_gap=3)
    cn = np.array([z["critical_node_nt"] for z in zones], dtype=np.int64)
    d = dist_to_n(cn - PROX, cn + PROX + 1, starts, ends) if len(cn) else np.array([])
    hits = int((d == 0).sum())
    max_di_n0 = float(np.nanmax(di))
    A = dict(zones=len(zones), hits=hits, p0=p0, z=z_prop(hits, len(zones), p0), max_di=max_di_n0)

    # B) контроль: N замаскирован, фланг 1…FLANK_MAX нт от N
    pos2, di2 = wta_di(seq, W, STEP, n_mode="mask", max_ambig=0.0)
    is2, _, _, _ = detect_scz(di2, SIGMA)
    zones2 = scz_zones(pos2, di2, is2, W, max_gap=3)
    ok = np.isfinite(di2)
    dd_all = dist_to_n(pos2[ok] - W, pos2[ok] + W, starts, ends)
    p_zone = float(((dd_all >= 1) & (dd_all <= FLANK_MAX)).sum()) / max(1, ok.sum())
    cn2 = np.array([z["critical_node_nt"] for z in zones2], dtype=np.int64)
    dd = dist_to_n(cn2 - W, cn2 + W, starts, ends) if len(cn2) else np.array([])
    in_zone = int(((dd >= 1) & (dd <= FLANK_MAX)).sum())
    B = dict(zones=len(zones2), in_zone=in_zone, p_zone=p_zone, z=z_prop(in_zone, len(zones2), p_zone),
             max_di=float(np.nanmax(di2)))
    return dict(n=n, n_runs=len(starts), A=A, B=B)


def make_random(n, gc, runs, seed):
    rng = np.random.default_rng(seed)
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    arr = rng.choice(ACGT, size=n, p=p).copy()
    for a, b in runs:
        arr[a:b] = ord("N")
    return arr.tobytes()


def fmt(r, label=None):
    A, B = r["A"], r["B"]
    head = f"{label:<10}" if label else ""
    return (f"{head}N-блоков={r['n_runs']:3d} | A (N=0.0): зон={A['zones']:4d} у N={A['hits']:3d} "
            f"({A['hits']/max(1, A['zones'])*100:5.1f}%) база={A['p0']*100:5.2f}% Z={A['z']:+7.1f} maxDI={A['max_di']:.2f}"
            f" | B (маска, фланг): зон={B['zones']:4d} во фланге={B['in_zone']:3d} "
            f"({B['in_zone']/max(1, B['zones'])*100:5.1f}%) ожид.={B['p_zone']*100:5.2f}% Z={B['z']:+5.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta")
    ap.add_argument("--seeds", type=int, default=5)
    a = ap.parse_args()

    if a.fasta:
        for h, s in read_fasta(a.fasta):
            if len(s) < 10 * W:
                continue
            print(fmt(analyze(s), label=h.split()[0][:10]))
        return

    n = 5_000_000
    cfgs = {
        "A: N≈0.3% (3 блока по 5 кб)": [(int(.2 * n), int(.2 * n) + 5000), (int(.5 * n), int(.5 * n) + 5000),
                                         (int(.8 * n), int(.8 * n) + 5000)],
        "B: N≈14% (400 кб + 3×100 кб)": [(0, 400_000), (int(.3 * n), int(.3 * n) + 100_000),
                                          (int(.55 * n), int(.55 * n) + 100_000), (int(.8 * n), int(.8 * n) + 100_000)],
    }
    print("Случайная ДНК (GC=41%, 5 Мнт, БИОЛОГИИ НЕТ). W=200, шаг 50, порог Mean+3.5σ.\n")
    for name, runs in cfgs.items():
        print(f"── {name}")
        for seed in range(1, a.seeds + 1):
            print("  " + fmt(analyze(make_random(n, 0.41, runs, seed)), label=f"seed {seed}"))
        print()


if __name__ == "__main__":
    main()