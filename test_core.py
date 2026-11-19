"""
test_core.py — проверки ядра. Запуск:  python test_core.py   (или pytest)
"""
import os, tempfile, time
import numpy as np
from tengri_core_wta import wta_di, detect_scz, scz_zones, read_fasta

RNG = np.random.default_rng(12345)
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)


def rand_seq(n, gc=0.41, rng=RNG):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return rng.choice(ACGT, size=n, p=p).tobytes()


def direct_half_dgc(seq: bytes, window, step):
    """Независимая «наивная» реализация: 0.5·|GC_вперёд − GC_назад|, без весов."""
    a = np.frombuffer(seq, dtype=np.uint8)
    gc = ((a == ord("G")) | (a == ord("C"))).astype(np.float64)
    out = []
    for i in range(window, len(gc) - window + 1, step):
        out.append(0.5 * abs(gc[i:i + window].mean() - gc[i - window:i].mean()))
    return np.array(out)


def test_identity_di_is_half_abs_dgc():
    s = rand_seq(60_000)
    pos, di = wta_di(s, 200, 50)
    ref = direct_half_dgc(s, 200, 50)
    assert len(di) == len(ref)
    assert np.max(np.abs(di - ref)) < 1e-9, "DI обязан совпадать с 0.5·|ΔGC|"


def test_real_genome_identity_if_available():
    path = "/mnt/user-data/uploads/sequence__10_.fasta"
    if not os.path.exists(path):
        return
    (_, s), = list(read_fasta(path))
    pos, di = wta_di(s, 1000, 500)
    assert np.max(np.abs(di - direct_half_dgc(s, 1000, 500))) < 1e-9


def test_no_edge_artifact():
    """На случайной ДНК максимум не должен липнуть к краям."""
    n_seq, edge_hits = 300, 0
    for k in range(n_seq):
        s = rand_seq(30_000, rng=np.random.default_rng(1000 + k))
        pos, di = wta_di(s, 200, 50)
        j = int(np.nanargmax(di))
        if j < 2 or j > len(di) - 3:
            edge_hits += 1
        assert pos[0] == 200 and pos[-1] <= 30_000 - 200
    assert edge_hits / n_seq < 0.03, f"максимум на краю в {edge_hits}/{n_seq} случаев"


def test_n_masked_vs_zero():
    s = bytearray(rand_seq(100_000))
    s[40_000:45_000] = b"N" * 5_000
    s = bytes(s)
    pos_m, di_m = wta_di(s, 200, 50, n_mode="mask")
    pos_z, di_z = wta_di(s, 200, 50, n_mode="zero")
    isn = (np.frombuffer(s, dtype=np.uint8) == ord("N")).astype(np.int64)
    P = np.concatenate(([0], np.cumsum(isn)))
    fN_f = (P[pos_m + 200] - P[pos_m]) / 200.0
    fN_b = (P[pos_m] - P[pos_m - 200]) / 200.0
    heavy = np.maximum(fN_f, fN_b) > 0.10
    clean = (fN_f == 0) & (fN_b == 0)
    assert heavy.any() and clean.any()
    assert np.isnan(di_m[heavy]).all(), "окна с долей N > 10% обязаны быть замаскированы (NaN)"
    assert np.isfinite(di_m[clean]).all(), "окна без N обязаны считаться"
    assert np.nanmax(di_m) < 0.2, "после маски не должно быть скачков от N"
    assert np.nanmax(di_z) > 0.6, "при N=0.0 граница N-блока обязана давать DI≈0.7"
    _, di_s = wta_di(s, 200, 50, n_mode="mask", max_ambig=0.0)
    assert np.isnan(di_s[(fN_f > 0) | (fN_b > 0)]).all()


def test_chunk_invariance():
    s = rand_seq(80_000)
    _, a = wta_di(s, 200, 50, chunk_blocks=7)
    _, b = wta_di(s, 200, 50)
    assert np.array_equal(a, b)


def test_records_are_separate():
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        f.write(">r1 first\n" + rand_seq(5_000).decode() + "\n>r2 second\n" + rand_seq(7_000).decode() + "\n")
        path = f.name
    recs = list(read_fasta(path))
    os.unlink(path)
    assert [h.split()[0] for h, _ in recs] == ["r1", "r2"]
    assert [len(s) for _, s in recs] == [5_000, 7_000]


def test_scz_zones_basic():
    di = np.array([0.01] * 20 + [0.5, 0.6, 0.01, 0.55] + [0.01] * 20 + [0.7] + [0.01] * 5)
    pos = np.arange(len(di)) * 50
    is_scz = di > 0.3
    zones = scz_zones(pos, di, is_scz, 200, max_gap=3)
    assert len(zones) == 2 and zones[0]["n_windows"] == 4 and zones[1]["peak_DI"] == 0.7


def test_performance_and_memory():
    import tracemalloc
    n = 40_000_000
    s = rand_seq(n, rng=np.random.default_rng(7))
    tracemalloc.start()
    t0 = time.time()
    pos, di = wta_di(s, 200, 50)
    dt = time.time() - t0
    cur, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"    {n/1e6:.0f} Мнт: {dt:.1f} с ({n/1e6/dt:.0f} Мнт/с), пик памяти внутри вызова ≈ {peak/1024**2:.0f} МБ")
    assert dt < 30 and peak / 1024**2 < 200


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} прошли")
    raise SystemExit(1 if failed else 0)