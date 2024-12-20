import pandas as pd
import numpy as np
from transformers import AutoTokenizer, AutoModel
from tqdm import tqdm
import torch
import xml.etree.ElementTree as ET
from rdkit.Chem import Draw, AllChem
from rdkit import Chem
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset

"""
The data provided by DrugBank is given in the format of XML so we need to first read and parse the .xml file to extract the 
relevant information to DDI.
"""

def process_drugbank_to_dataframe(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()
    namespace = {'ns': 'http://www.drugbank.ca'}
    data = []
    for drug in root.findall('ns:drug', namespace):
        primary_id = None
        for drugbank_id in drug.findall('ns:drugbank-id', namespace):
            # We take primary drugbank-id as our 1st Drug ID because it is in the correct format to be mapped to BindingDB
            if drugbank_id.attrib.get('primary') == 'true': 
                primary_id = drugbank_id.text
                break

        # If no primary ID is found, skip this drug entry
        if not primary_id:
            continue

        # Extract name and description
        name = drug.find('ns:name', namespace)
        description = drug.find('ns:description', namespace)
        primary_name = name.text if name is not None else 'N/A'
        primary_description = description.text if description is not None else 'N/A'

        # For each drug there is a list of drug-interaction pairings each with different description
        interactions = drug.find('ns:drug-interactions', namespace)

        # We skip drug entries that have no DDI records
        if interactions:
            for interaction in interactions.findall('ns:drug-interaction', namespace):
                interacting_drug_id = interaction.find('ns:drugbank-id', namespace).text
                interacting_drug_name = interaction.find('ns:name', namespace).text
                interaction_description = interaction.find('ns:description', namespace).text

                data.append({
                    'primary_id': primary_id,
                    'primary_name': primary_name,
                    'primary_description': primary_description,
                    'interacting_drug_id': interacting_drug_id,
                    'interacting_drug_name': interacting_drug_name,
                    'interaction_description': interaction_description
                })

    df = pd.DataFrame(data)
    multi_index_df = df.set_index(['primary_id', 'interacting_drug_id'])

    return multi_index_df

def extract_smiles_from_drugbank(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()
    namespace = {'ns': 'http://www.drugbank.ca'}
    data = []
    
    for drug in root.findall('ns:drug', namespace):
        drugbank_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        if drugbank_id is None:
            continue
        primary_id = drugbank_id.text
        
        smiles_element = drug.find('ns:calculated-properties/ns:property[ns:kind="SMILES"]/ns:value', namespace)
        smiles = smiles_element.text if smiles_element is not None else 'N/A'
        
        data.append({
            'DrugBank ID': primary_id,
            'SMILES': smiles
        })
    
    df = pd.DataFrame(data)
    return df


"""
Experimental approach to cluster the interaction description of DDI.
We use BioBert (a pretrained BERT model on biological texts) to extract the context appropriate embeddings for out descriptions.
Then we apply K-Means to cluster the embeddings into three clusters which later are interpretated as 'Major', 'Moderate' or 'Minor' 
sensitivity.
"""

def get_biobert_embeddings(ddi, device):
    tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-v1.1")
    model = AutoModel.from_pretrained("dmis-lab/biobert-v1.1").to(device)
    texts = ddi['interaction_description'].to_list()
    def get_embeddings_in_batches(texts, batch_size=32):
        embeddings = []
        # Proceed in batches to minimize time of execution.
        for i in tqdm(range(0, len(texts), batch_size), desc="Processing batches"):
            batch_texts = texts[i:i+batch_size]
            tokens = tokenizer(batch_texts, return_tensors="pt", padding=True, truncation=True, max_length=512).to(device)
            with torch.no_grad():
                outputs = model(**tokens)
            batch_embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()
            embeddings.append(batch_embeddings)
        return np.vstack(embeddings)
    return get_embeddings_in_batches(texts, 1028)


def create_fingerprints(smiles_list, prefix):
    molecules = [Chem.MolFromSmiles(smile) for smile in smiles_list]
    fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024) for mol in molecules]
    fingerprint_arrays = [np.array(fingerprint) for fingerprint in fingerprints]
    res_df = pd.DataFrame({
        f'{prefix}_smiles': smiles_list,
        f'{prefix}_fingerprint': fingerprint_arrays
    })
    return res_df

def get_loaders(ddi):
    features_primary = np.memmap('features_primary.dat', dtype='float16', mode='w+', shape=(len(ddi), 1024))
    features_interaction = np.memmap('features_interaction.dat', dtype='float16', mode='w+', shape=(len(ddi), 1024))
    labels = np.memmap('labels.dat', dtype='float16', mode='w+', shape=(len(ddi),))

    for i, (fp1, fp2) in enumerate(zip(ddi['primary_fingerprint'], ddi['interaction_drug_fingerprint'])):
        features_primary[i] = fp1
        features_interaction[i] = fp2
    labels[:] = ddi['labels']

    X1_train, X1_test, X2_train, X2_test, y_train, y_test = train_test_split(
        features_primary, features_interaction, labels, test_size=0.2, random_state=42, stratify=labels
    )

    X1_val, X1_test_final, X2_val, X2_test_final, y_val, y_test_final = train_test_split(
        X1_test, X2_test, y_test, test_size=0.5, random_state=42, stratify=y_test
    )
    X1_train_tensor = torch.tensor(X1_train, dtype=torch.float32)
    X2_train_tensor = torch.tensor(X2_train, dtype=torch.float32)
    y_train_tensor = torch.tensor(y_train, dtype=torch.long)

    X1_val_tensor = torch.tensor(X1_val, dtype=torch.float32)
    X2_val_tensor = torch.tensor(X2_val, dtype=torch.float32)
    y_val_tensor = torch.tensor(y_val, dtype=torch.long)

    X1_test_final_tensor = torch.tensor(X1_test_final, dtype=torch.float32)
    X2_test_final_tensor = torch.tensor(X2_test_final, dtype=torch.float32)
    y_test_final_tensor = torch.tensor(y_test_final, dtype=torch.long)

    train_dataset = TensorDataset(X1_train_tensor, X2_train_tensor, y_train_tensor)
    val_dataset = TensorDataset(X1_val_tensor, X2_val_tensor, y_val_tensor)
    test_dataset = TensorDataset(X1_test_final_tensor, X2_test_final_tensor, y_test_final_tensor)

    batch_size = 2048

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    return train_loader, val_loader, test_loader

def train_model(model, optimizer, criterion, num_epochs, train_loader, val_loader, device):
    train_losses = []
    val_accuracies = []
    early_stop_threshold = 0.005
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for batch_X1, batch_X2, batch_y in train_loader:
            batch_X1, batch_X2, batch_y = batch_X1.to(device), batch_X2.to(device), batch_y.to(device)
            outputs = model(batch_X1, batch_X2)
            loss = criterion(outputs, batch_y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item()

        model.eval()
        val_loss = 0.0
        correct = 0
        total = 0
        with torch.no_grad():
            for batch_X1, batch_X2, batch_y in val_loader:
                batch_X1, batch_X2, batch_y = batch_X1.to(device), batch_X2.to(device), batch_y.to(device)
                outputs = model(batch_X1, batch_X2)
                loss = criterion(outputs, batch_y)
                val_loss += loss.item()
                _, predicted = torch.max(outputs, 1)
                total += batch_y.size(0)
                correct += (predicted == batch_y).sum().item()

        epoch_loss = running_loss / len(train_loader)
        train_losses.append(epoch_loss)
        epoch_accuracy = 100 * correct / total
        val_accuracies.append(epoch_accuracy)


        print(f"Epoch [{epoch+1}/{num_epochs}], "
            f"Train Loss: {epoch_loss:.4f}, "
            f"Val Loss: {val_loss/len(val_loader):.4f}, "
            f"Val Accuracy: {epoch_accuracy:.2f}%")
        
        if epoch > 0 and abs(train_losses[-1] - train_losses[-2]) < early_stop_threshold:
            print(f"Training Loss converged (Δ < {early_stop_threshold}). Stopping early at epoch {epoch+1}.")
            break
    return model, train_losses, val_accuracies, epoch