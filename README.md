# 🧬 Tengri Structural Audit (TSA Engine v2.0.0)

**Deterministic GC-Divergence Pre-filtering for Genomic AI Pipelines**

This repository contains the open-source engine, datasets, and Atlases for the Tengri Structural Audit. 

Unlike heavy probabilistic models (e.g., AlphaFold) that choke on chromosome-scale data due to O(N²) memory constraints, TSA is a fast deterministic pre-filter. It scans massive genomes in milliseconds, isolating Structural Conflict Zones (SCZ) that represent high-risk areas for downstream alignment or AI inference.

## 📊 The Atlases
*   **[Human Structural Atlas (Chr 1, 6, 9, 19, 21)](ATLAS.md)**
*   **[Viral & Vaccine Atlas (Wild Virus vs. Moderna)](VIRAL_ATLAS.md)**

## 🚀 Quick Start (How to run TSA)
You don't need Docker or external servers. The engine runs locally in seconds.
Please see **[QUICKSTART.md](QUICKSTART.md)** for instructions.

## 🔬 Methodology
The engine computes the Divergence Index (DI = 0.5·|ΔGC|) across sliding windows (W=200 nt, step=50 nt). Structural Conflict Zones (SCZ) are identified via a strict statistical threshold (Mean + 3.5σ) with proper N-masking.

## 📖 Citation
Baynazarov, E. (2026). Tengri Protocol Dataset / Structural Audit Framework. Zenodo.
**DOI: 10.5281/zenodo.24244981** (Latest v4)

## 📄 License
Distributed under the MIT License.
Tengri Lab · Contact: ybaynazarovarchitect@gmail.com
