# 🧬 Tengri Structural Audit (TSA)

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21203626.svg)](https://doi.org/10.5281/zenodo.21203626)

**Deterministic stability framework for identifying structural vulnerabilities in RNA viruses and genomic sequences using Windowed Thermodynamic Approximation (WTA).**

Developed by **Erjan Baynazarov** (Tengri Lab).

---

## 📌 Overview

Tengri Structural Audit (TSA) is an advanced computational pipeline designed to detect structural conflict zones (SCZ) in genomic and viral sequences. By applying the **Windowed Thermodynamic Approximation (WTA)** algorithm, TSA maps thermodynamic divergence across nucleotide positions, providing high-resolution structural profiling without relying on heavy deep-learning inference.

### Core Features:
- **Mode 1: Single Sequence Analysis** — Scans a FASTA file, computes divergence profiles, detects anomalous structural conflict zones (SCZ) based on statistical thresholds ($\mu + 3.5\sigma$), and generates visual reports with zoom-ins on peak regions.
- **Mode 2: Comparative Analysis** — Aligns and compares a sample sequence against a standard reference genome, highlighting structural divergences and shifts in anomaly counts.
- **Dual Architecture** — Headless high-performance **FastAPI** backend coupled with an interactive **Streamlit** graphical user interface.

---

## 🛠️ Project Structure

```text
tengri-tsa/
├── api.py                  # FastAPI backend (WTA core engine)
├── gui.py                  # Streamlit interactive frontend
├── requirements_api.txt    # Backend dependencies
├── requirements_gui.txt    # Frontend dependencies
├── Dockerfile.api          # Container config for API
├── Dockerfile.gui          # Container config for GUI
└── docker-compose.yml      # Multi-container orchestration
🚀 Quick Start (Local Installation)
1. Clone the Repository
git clone [https://github.com/erjansarbaz/Tengri-Structural-Audit.git](https://github.com/erjansarbaz/Tengri-Structural-Audit.git)
cd Tengri-Structural-Audit

2. Run with Docker (Recommended)
Make sure you have Docker and Docker Compose installed:

docker-compose up --build

Streamlit GUI: http://localhost:8501

FastAPI Docs: http://localhost:8000/docs

3. Manual Installation (Virtual Environment)

# Setup backend
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements_api.txt
uvicorn api:app --host 0.0.0.0 --port 8000 &

# Run frontend (in a separate terminal)
pip install -r requirements_gui.txt
streamlit run gui.py

🔬 Methodology (WTA)The Windowed Thermodynamic Approximation evaluates nucleotide weight distributions across sliding windows (W=200\text{ nt}, step 50\text{ nt}). The divergence index (DI) quantifies local structural stress, isolating regions that deviate significantly from baseline stability thresholds.
📖 Citation & DatasetIf you use Tengri Structural Audit in your research, please cite the corresponding dataset and framework:Baynazarov, E. (2026). 
Tengri Protocol Dataset / Structural Audit Framework. 
Zenodo.DOI: 10.5281/zenodo.21203626
📄 LicenseDistributed under the MIT License. See LICENSE for more information.Tengri Lab · 
Contact: ybaynazarovarchitect@gmail.com