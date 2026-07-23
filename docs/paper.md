---
title: 'Tengri Structural Audit: A Deterministic Framework for Genomic Stability Analysis Using Windowed Thermodynamic Approximation'
tags:
  - Python
  - bioinformatics
  - structural analysis
  - thermodynamic approximation
  - viral genomics
authors:
  - name: Yerzhan Baynazarov
    orcid: 0009-0002-0600-0266
    affiliation: 1
affiliations:
  - name: Independent Researcher, Kazakhstan
    index: 1
date: 23 July 2026
bibliography: paper.bib
---

# Summary

The Tengri Structural Audit (TSA) is a computational software package designed to identify structural vulnerabilities and high-entropy conflict zones in nucleotide sequences (viral genomes, synthetic constructs, and data streams) without relying on resource-intensive ab initio force-field simulations or probabilistic machine learning models. By applying a scale-invariant **Windowed Thermodynamic Approximation (WTA)** pipeline, the software maps sequences directly into a high-resolution thermodynamic landscape, isolating critical structural seams and anomalies via absolute local moving-average divergence.

# Statement of Need

Traditional bioinformatics pipelines for RNA secondary structure prediction (such as deep-learning-based folding or heavy thermodynamic free-energy minimizers) often suffer from extreme computational overhead and are heavily anchored to biological annotations. In high-load infrastructure, viral surveillance, or rapid screening scenarios, researchers require a fast, deterministic heuristic to detect structural stress points, artificial insertions, or mutation pressures. 

`Tengri-Structural-Audit` addresses this gap by offering a lightweight, reproducible Python toolkit that replaces probabilistic noise analysis with deterministic stability indices. It targets computational biologists, security engineers, and researchers needing high-speed structural triage of genomic files.

# Mathematical Formulation

The WTA pipeline evaluates sequences through standardized base stability mapping:
- G = 1.0, C = 1.0 (High-bond-energy baselines)
- A = 0.5, T/U = 0.5 (Lower-energy baselines)
- N = 0.0 (Neutralized residues)

Forward (S_f) and backward (S_b) local moving averages are calculated over a defined window length L, and the divergence index D(i) isolates structural conflict zones where values exceed a strict statistical threshold (> \mu + k \cdot \sigma).

# Software Architecture

The package includes:
- `run_global_audit.py`: Main execution script for processing FASTA inputs.
- `analysis_logic.py`: Core WTA mathematical operators and sliding-window engine.
- `batch_tengri.py`: Utilities for multi-file processing and batch validations.

All underlying framework versions and validation datasets are permanently archived via Zenodo (DOI: 10.5281/zenodo.21203626).
