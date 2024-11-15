import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

def calculate_drug_properties(df, smiles_column='Ligand SMILES'):
    """
    Calculates drug-likeness properties for each compound in the dataset.
    
    Parameters:
    df (DataFrame): DataFrame containing the compounds with SMILES strings.
    smiles_column (str): The column name that contains the SMILES strings.

    Returns:
    DataFrame: Original DataFrame with additional columns for each drug-likeness property.
    """
    
    molecular_weights = []
    log_p_values = []
    h_donors = []
    h_acceptors = []
    rotatable_bonds = []
    psa_values = []
    qed_scores = []
    lipinski_pass = []
    veber_pass = []

    for smiles in df[smiles_column]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            molecular_weights.append(None)
            log_p_values.append(None)
            h_donors.append(None)
            h_acceptors.append(None)
            rotatable_bonds.append(None)
            psa_values.append(None)
            qed_scores.append(None)
            lipinski_pass.append(False)
            veber_pass.append(False)
        else:
            molecular_weights.append(Descriptors.MolWt(mol))
            log_p_values.append(Descriptors.MolLogP(mol))
            h_donors.append(Descriptors.NumHDonors(mol))
            h_acceptors.append(Descriptors.NumHAcceptors(mol))
            rotatable_bonds.append(Descriptors.NumRotatableBonds(mol))
            psa_values.append(Descriptors.TPSA(mol))
            qed_scores.append(QED.qed(mol))
            
            # Check Lipinski's Rule of Five
            lipinski_pass.append(
                molecular_weights[-1] < 500 and
                log_p_values[-1] < 5 and
                h_donors[-1] <= 5 and
                h_acceptors[-1] <= 10
            )
            
            # Check Veber's Rule
            veber_pass.append(
                rotatable_bonds[-1] <= 10 and psa_values[-1] <= 140
            )

    df['Molecular Weight'] = molecular_weights
    df['Log P'] = log_p_values
    df['H Donors'] = h_donors
    df['H Acceptors'] = h_acceptors
    df['Rotatable Bonds'] = rotatable_bonds
    df['PSA'] = psa_values
    df['QED Score'] = qed_scores
    df['Lipinski Pass'] = lipinski_pass
    df['Veber Pass'] = veber_pass

    return df
    
    
    
def calculate_drug_properties(df, smiles_column='Ligand SMILES'):
    """
    Calculates drug-likeness properties for each compound in the dataset.
    
    Parameters:
    df (DataFrame): DataFrame containing the compounds with SMILES strings.
    smiles_column (str): The column name that contains the SMILES strings.

    Returns:
    DataFrame: Original DataFrame with additional columns for each drug-likeness property.
    """
    
    molecular_weights = []
    log_p_values = []
    h_donors = []
    h_acceptors = []
    rotatable_bonds = []
    psa_values = []
    qed_scores = []
    lipinski_pass = []
    veber_pass = []

    for smiles in df[smiles_column]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            molecular_weights.append(None)
            log_p_values.append(None)
            h_donors.append(None)
            h_acceptors.append(None)
            rotatable_bonds.append(None)
            psa_values.append(None)
            qed_scores.append(None)
            lipinski_pass.append(False)
            veber_pass.append(False)
        else:
            molecular_weights.append(Descriptors.MolWt(mol))
            log_p_values.append(Descriptors.MolLogP(mol))
            h_donors.append(Descriptors.NumHDonors(mol))
            h_acceptors.append(Descriptors.NumHAcceptors(mol))
            rotatable_bonds.append(Descriptors.NumRotatableBonds(mol))
            psa_values.append(Descriptors.TPSA(mol))
            qed_scores.append(QED.qed(mol))
            
            # Check Lipinski's Rule of Five
            lipinski_pass.append(
                molecular_weights[-1] < 500 and
                log_p_values[-1] < 5 and
                h_donors[-1] <= 5 and
                h_acceptors[-1] <= 10
            )
            
            # Check Veber's Rule
            veber_pass.append(
                rotatable_bonds[-1] <= 10 and psa_values[-1] <= 140
            )

    df['Molecular Weight'] = molecular_weights
    df['Log P'] = log_p_values
    df['H Donors'] = h_donors
    df['H Acceptors'] = h_acceptors
    df['Rotatable Bonds'] = rotatable_bonds
    df['PSA'] = psa_values
    df['QED Score'] = qed_scores
    df['Lipinski Pass'] = lipinski_pass
    df['Veber Pass'] = veber_pass

    return df

def check_smiles_validity(smiles_list):
    """
    Checks the validity of a list of SMILES strings.

    Parameters:
    smiles_list (list): List of SMILES strings.

    Returns:
    DataFrame: DataFrame with SMILES strings and their validity.
    """
    results = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({"SMILES": smiles, "Valid": False})
        else:
            results.append({"SMILES": smiles, "Valid": True})
    return pd.DataFrame(results)

