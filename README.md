Tengri Structural Audit (TSA) 🧬
Deterministic stability framework for identifying structural vulnerabilities in RNA viruses using Windowed Thermodynamic Approximation (WTA).
📌 Overview
The Tengri Structural Audit (TSA) is a computational method that shifts genomic analysis from descriptive statistics to deterministic physical modeling. By applying Windowed Thermodynamic Approximation (WTA), we evaluate the thermodynamic tension and structural symmetry of viral genomes to identify conserved functional "nodes" (structural seams).
🚀 Key Capabilities
Topological Seam Detection: Identification of artificial insertions via Lyapunov Seam divergence[cite: 15].

Adaptive Calibration: Support for DNA/RNA/Synthetic profiles via adjustable weight matrices.

Deterministic Audit: Replaces probabilistic noise analysis with physical stability indices.
📊 Validation Results (Reference: E. coli  HCV) StrainTargeting EfficiencyStatusNC_00410298.21 ValidatedNC_00982575.00 Under Analysis
🛠 Quick Start
Clone & Setup:

Bash
git clone https://github.com/erjansarbaz/Tengri-Structural-Audit.git
cd Tengri-Structural-Audit
pip install numpy pandas
Run Audit:

Bash
python run_global_audit.py --input data/target.fasta
📖 Citation
If you use TSA in your research, please cite our technical protocol:

Preprint: Baynazarov, Y. (2026). Deterministic Structural Audit of Viral RNA using WTA. Preprints.org.

Zenodo: 10.5281/zenodo.21203626
