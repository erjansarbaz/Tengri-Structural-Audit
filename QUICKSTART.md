markdown
# 🚀 QUICKSTART: How to run TSA in 2 minutes

**Are you here from LinkedIn or Vecura? Welcome!** 

This repository contains the engine, the data, and the Atlas. If you just want to see how it works, follow these 3 steps.

### Step 1: Install dependencies
You need Python 3.9+ and numpy/pandas.
```bash
pip install numpy pandas matplotlib
Step 2: Run the engine on a test file
We included a small test file (sequence (55).fasta - Wild SARS-CoV-2). Run this command:

bash
python tsa_cli.py "sequence (55).fasta" --out my_first_audit
Step 3: Check your results
The engine will create a folder my_first_audit/ containing:

NC_045512.2_profile.png — The graph (look for the red dots).

NC_045512.2_summary.json — The numbers (SCZ zones, DI peaks).

NC_045512.2_windows.csv — Raw data for each window.

Want to run it on your own genome?
Just replace the FASTA file:

bash
python tsa_cli.py your_genome.fasta --out my_custom_audit
Where is the Atlas?
Human Chromosomes: ATLAS.md

Virus vs Vaccine: VIRAL_ATLAS.md

Found a bug?
Please open an Issue. We are in active development (v5 is coming soon with corrected N-masking and ClinVar validation).
