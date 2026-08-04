# Tengri Structural Audit (TSA)

Deterministic stability framework for identifying structural vulnerabilities in RNA viruses using Windowed Thermodynamic Approximation (WTA).

## Overview

The Tengri Structural Audit (TSA) is a computational method that shifts genomic analysis from descriptive statistics to deterministic physical modeling. By applying Windowed Thermodynamic Approximation (WTA), we evaluate the thermodynamic tension and structural symmetry of viral genomes to identify conserved functional "nodes" (structural seams).

## Key Capabilities

* **Topological Seam Detection:** Identification of structural shifts and divergences via Lyapunov stability parameters.
* **Adaptive Calibration:** Universal support for DNA/RNA/Synthetic profiles via adjustable weight matrices and modular input formats.
* **Deterministic Audit:** Replaces probabilistic noise analysis with physical stability indices.

## Validation Results

The framework supports universal sequence processing, with specific validation benchmarks targeting viral strains and model references:

| Strain / Accession | Target Focus | Efficiency Index | Status |
| :--- | :--- | :--- | :--- |
| `NC_004102.1` (HCV reference / model) | Viral RNA / Structural Seams | 0.98 | **Validated** |
| `NC_000913.3` (E. coli baseline / control) | Universal Sequence Baseline | 0.75 | **Under Analysis** |

*(Note: The pipeline natively handles arbitrary FASTA inputs, enabling both pathogen-specific audits and universal structural stress tests).*

## Quick Start

### Clone & Setup

```bash
git clone https://github.com/erjansarbaz/Tengri-Structural-Audit.git
cd Tengri-Structural-Audit
pip install numpy pandas
python run_global_audit.py --input data/target.fasta
Citation
If you use TSA in your research, please cite our technical protocol:

Preprint: Baynazarov, Y. (2026). Deterministic Structural Audit of Viral RNA using WTA. Preprints.org.

Zenodo: 10.5281/zenodo.21203626
