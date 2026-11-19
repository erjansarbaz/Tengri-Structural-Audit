"""
tengri_core_wta.py — единое ядро WTA для TSA (v2.0.0)

Что считает метод (без метафор):
    вес нуклеотида: G=C=1.0, A=T(U)=0.5   →   E = 0.5 + 0.5·[G или C]
    D(i) = |mean(E[i:i+W]) − mean(E[i−W:i])|  =  0.5 · |GC_вперёд − GC_назад|
То есть DI — это масштабированная разность GC-состава соседних окон (проверено численно, см. test_core.py).

Что исправлено относительно версий v1:
  1. Нет краевых артефактов: DI считается только там, где ОБА окна полные.
  2. N и IUPAC-коды не «вес 0.0», а пропуск: среднее считается по валидным основаниям,
     окна с долей невалидных > max_ambig маскируются (NaN).
  3. Память O(длина/step), а не O(длина): блочные суммы по кускам.
  4. Каждая запись FASTA считается отдельно — координаты всегда внутри записи.
"""
from __future__ import annotations

import numpy as np

__version__ = "2.0.0"

# Таблицы по байту символа: вес и признак «валидное основание»
_W = np.zeros(256, dtype=np.float32)
_OK = np.zeros(256, dtype=np.float32)
for _ch, _w in (("G", 1.0), ("C", 1.0), ("A", 0.5), ("T", 0.5), ("U", 0.5)):
    for _c in (_ch, _ch.lower()):
        _W[ord(_c)] = _w
        _OK[ord(_c)] = 1.0


def read_fasta(path):
    """Генератор (заголовок, байты последовательности) — по одной записи."""
    header, parts = None, []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if parts:
                    yield (header or "unnamed"), b"".join(parts)
                header, parts = line[1:].strip().decode("utf-8", "ignore"), []
            else:
                s = line.strip()
                if s:
                    parts.append(s)
    if parts:
        yield (header or "unnamed"), b"".join(parts)


def _block_stats(seq: bytes, step: int, n_as_zero: bool, chunk_blocks: int):
    """Сумма весов и число валидных оснований в непересекающихся блоках по `step` нт."""
    arr = np.frombuffer(seq, dtype=np.uint8)
    n_blocks = len(arr) // step          # хвост < step нт отбрасывается
    bs = np.empty(n_blocks, dtype=np.float64)
    bc = np.empty(n_blocks, dtype=np.float64)
    for b0 in range(0, n_blocks, chunk_blocks):
        b1 = min(n_blocks, b0 + chunk_blocks)
        a = arr[b0 * step: b1 * step]
        bs[b0:b1] = _W[a].reshape(-1, step).sum(axis=1, dtype=np.float64)
        if n_as_zero:
            bc[b0:b1] = float(step)
        else:
            bc[b0:b1] = _OK[a].reshape(-1, step).sum(axis=1, dtype=np.float64)
    return bs, bc


def wta_di(seq, window: int = 200, step: int = 50, n_mode: str = "mask",
           max_ambig: float = 0.10, chunk_blocks: int = 100_000):
    """
    Возвращает (pos, di):
        pos[k] — nt_start прямого окна (0-based), di[k] — Divergence Index.
    n_mode="mask": неоднозначные основания исключены из среднего; окно с долей
                   невалидных > max_ambig → NaN.
    n_mode="zero": N и любой нестандартный символ = вес 0.0 (как в препринтах v1).
    """
    if window % step:
        raise ValueError("window должен быть кратен step")
    if n_mode not in ("mask", "zero"):
        raise ValueError("n_mode: 'mask' или 'zero'")
    seq = seq if isinstance(seq, (bytes, bytearray)) else str(seq).encode()
    m = window // step
    bs, bc = _block_stats(bytes(seq), step, n_mode == "zero", chunk_blocks)
    n_blocks = len(bs)
    if n_blocks < 2 * m + 1:
        raise ValueError("последовательность короче 2·window")

    cs = np.concatenate(([0.0], np.cumsum(bs)))
    cc = np.concatenate(([0.0], np.cumsum(bc)))
    b = np.arange(m, n_blocks - m + 1)            # нужны блоки [b−m, b+m)
    sf, cf = cs[b + m] - cs[b], cc[b + m] - cc[b]
    sb, cb = cs[b] - cs[b - m], cc[b] - cc[b - m]
    with np.errstate(invalid="ignore", divide="ignore"):
        di = np.abs(sf / cf - sb / cb)
    bad = (cf < (1.0 - max_ambig) * window) | (cb < (1.0 - max_ambig) * window)
    di[bad] = np.nan
    return b * step, di


def detect_scz(di, sigma: float = 3.5):
    """Порог Mean + sigma·σ по валидным окнам. Возвращает (is_scz, mean, std, threshold)."""
    di = np.asarray(di)
    ok = np.isfinite(di)
    mean = float(np.mean(di[ok]))
    std = float(np.std(di[ok]))
    thr = mean + sigma * std
    is_scz = np.zeros(len(di), dtype=bool)
    is_scz[ok] = di[ok] > thr
    return is_scz, mean, std, thr


def scz_zones(pos, di, is_scz, window: int, max_gap: int = 3):
    """Склеивает соседние SCZ-окна (допуск max_gap окон) в зоны; Critical Node = окно с max DI."""
    idx = np.flatnonzero(is_scz)
    if len(idx) == 0:
        return []
    bounds, start, prev = [], idx[0], idx[0]
    for k in idx[1:]:
        if k - prev - 1 > max_gap:
            bounds.append((start, prev))
            start = k
        prev = k
    bounds.append((start, prev))
    zones = []
    for z, (a, c) in enumerate(bounds, 1):
        j = a + int(np.nanargmax(di[a:c + 1]))
        zones.append({
            "scz_id": f"Z{z}",
            "nt_start": int(pos[a]),
            "nt_end": int(pos[c] + window),
            "critical_node_nt": int(pos[j]),
            "peak_DI": round(float(di[j]), 6),
            "n_windows": int(c - a + 1),
        })
    return zones