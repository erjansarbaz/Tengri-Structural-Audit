engri Structural Audit (TSA)
Deterministic stability framework for identifying structural vulnerabilities in RNA viruses.

Tengri Structural Audit (TSA) is a computational method for analyzing the physical architecture of viral RNA through Windowed Thermodynamic Approximation (WTA). Unlike traditional sequence-alignment tools that focus solely on nucleotide identity, TSA evaluates thermodynamic tension and secondary structure stability to identify conserved functional "nodes."
Key Discovery: The HCV Vulnerability Node
During the analysis of Hepatitis C Virus (HCV) genomes, we identified a critical structural node at position 9080-9135.
Data analysis reveals a fundamental architectural divergence: HIV-1 employs a 'Cluster Isolation Strategy' to protect critical genomic nodes, whereas the T4 Phage utilizes 'Distributed Load' to maintain complex mechanical functionality. These findings validate TSA as a robust tool for viral structural classification.
Stability: This domain maintains a consistent physical architecture despite 47+ point mutations.

Validation: Our method demonstrates high recognition accuracy (>98% in reference strains).
Quick Start
Requirements
Ensure you have Python 3.9+ installed.

Bash
pip install numpy pandas
Setup
Clone this repository.

Place your target .fasta files into the repository folder.
Execution
Run the global audit:

Bash
python run_global_audit.py
Output
The script generates batch_analysis_summary.csv, containing the structural risk assessment for each provided genome.
Performance DataStrain IDTargeting Efficiency (%)NC_004102 (Reference)98.21%NC_00982575.00%NC_00982469.64%
Citation
Please cite the corresponding research paper (https://doi.org/10.5281/zenodo.21203626) when using this tool in your work.
