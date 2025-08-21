# Boltz Polaris Potency Evaluation

This repository contains evaluation scripts and results for the Boltz model's performance on the Polaris blind challenge datasets, specifically focusing on antiviral potency, ligand poses, and ADMET predictions.

## Project Structure

```
├── boltz_affinity_predictions.csv          # Affinity prediction results
├── boltz_polaris_evaluation.ipynb          # Main evaluation notebook
├── code_to_create_yamls.py                 # YAML configuration generator
├── env.yml                                 # Conda environment specification
├── evaluation/                             # Evaluation modules
│   ├── admet.py                           # ADMET evaluation functions
│   ├── bootstrapping.py                   # Statistical bootstrapping utilities
│   ├── cld.py                             # CLD (Chemical Library Design) evaluation
│   ├── ligand_poses.py                    # Ligand pose evaluation
│   ├── potency.py                         # Potency prediction evaluation
│   ├── utils.py                           # Common utilities
│   └── data/                              # Evaluation datasets
└── leaderboards/                          # Challenge leaderboard submissions
    ├── antiviral-admet-2025/              # ADMET challenge results
    ├── antiviral-ligand-poses-2025/       # Ligand poses challenge results
    └── antiviral-potency-2025/            # Potency challenge results
```

## Setup

1. Create the conda environment:
```bash
conda env create -f env.yml
conda activate boltz-polaris
```
python code_to_create_yamls.py
```
<!--  Boltz things not in this repo -->
boltz predict ./yaml_files/MOl_MERS.yaml --use_msa_server --use_potentials --no kernel

2. Run the evaluation notebook:
```bash
jupyter notebook boltz_polaris_evaluation.ipynb
```



## Results from Boltz 2
With detailed affinity predictions in [boltz_affinity_predictions.csv](boltz_affinity_predictions.csv) from Boltz 2 outputs.


## Usage

Generate YAML configurations to input into Boltz 2 for predictions:
```bash

```

## Evaluation Components

- **Potency Evaluation**: Assessment of binding affinity predictions using [`evaluation.potency`](evaluation/potency.py)
- **Ligand Poses**: Evaluation of predicted ligand conformations via [`evaluation.ligand_poses`](evaluation/ligand_poses.py)  
- **ADMET Properties**: Drug-like property predictions evaluated in [`evaluation.admet`](evaluation/admet.py)
- **Statistical Analysis**: Bootstrapping confidence intervals implemented in [`evaluation.bootstrapping`](evaluation/bootstrapping.py)

## Attribution

The evaluation framework is adapted from the [ASAP Polaris Blind Challenge Examples](https://github.com/asapdiscovery/asap-polaris-blind-challenge-examples) repository (accessed August 21st, 2025).
