import polaris as po
import pandas as pd
from rdkit import Chem

import os
from pathlib import Path
def convert_cxsmiles_to_smiles(cxsmiles):
    """Convert CXSMILES to SMILES using RDKit."""
    mol = Chem.MolFromSmiles(cxsmiles)
    if mol:
        return Chem.MolToSmiles(mol)
    else:
        return None
        print(f"Failed to convert CXSMILES: {cxsmiles}")


# load the dataset from the Hub
dataset = po.load_dataset("asap-discovery/antiviral-potency-2025-unblinded")
# convert to dataframe
df = pd.DataFrame(dataset[:])

# >8R5J_1|Chains A, B|Non-structural protein 11|Middle East respiratory syndrome-related coronavirus (1335626)

MERS_AA = "SGLVKMSHPSGDVEACMVQVTCGSMTLNGLWLDNTVWCPRHVMCPADQLSDPNYDALLISMTNHSFSVQKHIGAPANLRVVGHAMQGTLLKLTVDVANPSTPAYTFTTVKPGAAFSVLACYNGRPTGTFTVVMRPNYTIKGSFLCGSCGSVGYTKEGSVINFCYMHQMELANGTHTGSAFDGTMYGAFMDKQVHQVQLTDKYCSVNVVAWLYAAILNGCAWFVKPNRTSVVSFNEWALANQFTEFVGTQSVDMLAVKTGVAIEQLLYAIQQLYTGFQGKQILGSTMLEDEFTPEDVNMQIMGVVMQ"
# >7CAM_1|Chains A, B|3C-like proteinase|Severe acute respiratory syndrome coronavirus 2 (2697049)
SARS_AA = "GSGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ"
df_test = df[df['Set'] == 'Test']  # Display the first 5 rows of the test set


df_test['SMILES'] = df_test['CXSMILES'].apply(convert_cxsmiles_to_smiles)


def create_boltz_yaml(protein_name, protein_seq, molecule_name, smiles, output_dir="./yaml_files"):
    """
    Create a Boltz-compatible YAML file for protein-ligand prediction
    
    Args:
        protein_name: "MERS" or "SARS"
        protein_seq: Protein amino acid sequence
        molecule_name: Name of the small molecule
        smiles: SMILES string of the molecule
        output_dir: Directory to save YAML files
    """
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(exist_ok=True)
    
    # Create filename
    filename = f"{molecule_name}_{protein_name}.yaml"
    filepath = Path(output_dir) / filename
    
    # YAML content
    yaml_content = f"""version: 1
sequences:
  - protein:
      id: [A, B]
      sequence: {protein_seq}
      msa: ../msa/ASAP-0000175_{protein_name}_0.csv
  - ligand:
      id: C
      smiles: '{smiles}'

properties:
  - affinity:
      binder: C
"""
    
    # Write YAML file
    with open(filepath, 'w') as f:
        f.write(yaml_content)
    
    print(f"Created: {filepath}")
    return filepath

# Example usage with your sequences
# Replace 'your_molecule_smiles' with actual SMILES strings

# turn df_test into a dictionary of molecule names and SMILES
example_molecules = {row['Molecule Name']: row['SMILES'] for _, row in df_test.iterrows()}


# Create YAML files for each combination
for molecule_name, smiles in example_molecules.items():
    # Create MERS combinations
    create_boltz_yaml("MERS", MERS_AA, molecule_name, smiles)
    
    # Create SARS combinations  
    create_boltz_yaml("SARS", SARS_AA, molecule_name, smiles)

print("\nTo use these YAML files with Boltz, run:")
print("boltz predict <yaml_file> --use_msa_server --use_potentials")
print("\nExample:")
print("boltz predict ./yaml_files/MOl_MERS.yaml --use_msa_server --use_potentials --no kernel")