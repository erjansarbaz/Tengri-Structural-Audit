# 🧬 Tengri Structural Audit (GUI Client)

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21203626.svg)](https://doi.org/10.5281/zenodo.21203626)

**Interactive graphical interface for the Tengri Structural Audit (TSA) framework, designed for identifying structural vulnerabilities in RNA viruses and genomic sequences using Windowed Thermodynamic Approximation (WTA).**

Developed by **Erjan Baynazarov** (Tengri Lab).

---

## 📌 Overview

This repository contains the **public client-side GUI** for Tengri Structural Audit (TSA). The core computational engine (FastAPI backend and proprietary WTA algorithms) operates on a secure remote server. 

The GUI connects to the backend API to provide high-resolution structural profiling of genomic FASTA sequences without exposing the underlying intellectual property or core algorithms.

### Core Features:
- **Mode 1: Single Sequence Analysis** — Upload FASTA files, compute divergence profiles, detect anomalous structural conflict zones (SCZ) based on statistical thresholds (\mu + 3.5\sigma), and download visual PNG reports or CSV data.
- **Mode 2: Comparative Analysis** — Align and compare a sample sequence against a standard reference genome to highlight structural shifts and divergence metrics.
- **Client-Server Architecture** — Lightweight Streamlit interface communicating with a high-performance backend API

## 🛠️ Repository Structure

tengri-tsa-gui/
├── gui.py                  # Streamlit interactive frontend
├── requirements_gui.txt    # Frontend dependencies
├── Dockerfile.gui          # Container configuration for GUI
└── .gitignore              # Excluded private configurations
🚀 Quick Start (Local Deployment)1. Clone the RepositoryBashgit clone [https://github.com/erjansarbaz/Tengri-Structural-Audit.git](https://github.com/erjansarbaz/Tengri-Structural-Audit.git)
cd Tengri-Structural-Audit
2. Run with DockerMake sure you have Docker installed:Bashdocker build -t tengri-gui -f Dockerfile.gui .
docker run -p 8501:8501 -e API_URL="[http://92.38.48.32:8000](http://92.38.48.32:8000)" tengri-gui
Open Streamlit GUI: http://localhost:85013. Manual Installation (Virtual Environment)Bashpython -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements_gui.txt

# Set the backend API endpoint before running
# On Windows PowerShell:
env:API_URL="[http://92.38.48.32:8000](http://92.38.48.32:8000)"

streamlit run gui.py
🔬 Methodology (WTA)The Windowed Thermodynamic Approximation evaluates nucleotide weight distributions across sliding windows (W=200\text{ nt}, step 50\text{ nt}). The divergence index (DI) quantifies local structural stress, isolating regions that deviate significantly from baseline stability thresholds.
📖 Citation & DatasetIf you use Tengri Structural Audit in your research, please cite the corresponding dataset and framework:Baynazarov, E. (2026). Tengri Protocol Dataset / Structural Audit Framework. Zenodo.DOI: 10.5281/zenodo.21203626
📄 LicenseDistributed under the MIT License. See LICENSE for more information.Tengri Lab ·
Contact: ybaynazarovarchitect@gmail.com
