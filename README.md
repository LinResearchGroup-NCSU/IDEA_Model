# IDEA_Model
**Interpretable protein-DNA Energy Associative (IDEA) Model**

This repository provides a cleaned-up implementation of the IDEA model. Training and testing data are included using three human MAX transcription factors (TFs) as examples. The results shown here correspond to Figure S1 in our IDEA manuscript.

---

## 🚀 Clone the Repository

```bash
git clone https://github.com/LinResearchGroup-NCSU/IDEA_Model
cd IDEA_Model
```

---

## 🛠️ Prerequisites

Please ensure you have the following software installed:

### 1. Python 3
You can download Python 3 from the official website:  
👉 https://www.python.org/downloads/

### 2. Modeller
Modeller is used **only** to extract sequences from PDB files—not for structure modeling.

- 🔗 [Modeller Download Page](https://salilab.org/modeller/download_installation.html)
- Required script:  
  `3max_case/optimization/for_bindingE/template/sequences/buildseq.py`

### 3. Python Packages

#### Using `pip`:
```bash
pip install biopython numpy mdtraj
```

#### Or using `conda`:
```bash
conda install -c conda-forge biopython numpy mdtraj
```

---

## ⚙️ Max Protein Complexes Energy Model Training

This guide describes how to train energy models using 3 MAX TFs (PDB IDs: 1HLO, 1NLW, 1NKP).

### Step 1: Prepare Input Files

- **Protein list**:
  Create `3max_case/proteinList.txt`:
  ```
  1hlo
  1nlw
  1nkp
  ```

- **PDB files**:
  Place modified PDB files into `3max_case/PDBs`, named as `{PDB_ID}_modified.pdb`.  
  Ensure:  
  - Protein chain → A  
  - DNA chains → B and C

### Step 2: Start the Training

```bash
bash train.sh
```

### Step 3: Configuration Notes

- Modify `get_interaction_atom(residue)` in `common_functions.py` to define the atoms (e.g., CA for protein, C5 for DNA).
- In `optimize_gamma.py`, set `cutoff_mode` to filter noise (keeps top 70 eigenvalues).

### Step 4: Output Description

Output files will be found in:  
`3max_case/optimization/for_training_gamma/gammas/randomized_decoy`

- `*_gamma_filtered`: Final trained energy model  

---

## 🔬 Testing Protein-DNA Binders

To generate testing `phi` and predicted binding energies:

### Step 1: Prepare Testing Files

- Place testing PDB files into: `1hlo_phi255/PDBs`  
- DNA sequence file: `1hlo_phi255/sequences`

### Step 2: Generate Sequences (Optional)

Use `dna_half.seq` with:
```bash
python reverse_complement.py
python merge.py
```

### Step 3: Generate Testing `phi`

```bash
bash phi255.sh 1hlo
```

Outputs:
- `phi_pairwise_contact_well_native_decoys_*` (255 lines)  
- `phi_pairwise_contact_well_native_native_*` (native only)

### Step 4: Calculate Predicted Energy

```bash
cd testing_energy
python calculate_testing_energy.py
```

Results saved in: `Energy_mg.txt`

---

## 📂 Supplementary Materials

### 🔬 Trained Energy Models  
Trained gamma files used in the manuscript:  
👉 [View folder](https://github.com/LinResearchGroup-NCSU/IDEA_Model/tree/main/supplementary_materials/IDEA_trained_energy_models)

### 📊 Raw Data  
Raw inputs and results for all main figures:  
📦 [Download raw_data.zip](https://github.com/LinResearchGroup-NCSU/IDEA_Model/blob/main/supplementary_materials/raw_data.zip)

### 📁 Other Processed Published Models  
Code and data used to process published models:  
- **DBD-hunter** ([NAR, 2008](https://academic.oup.com/nar/article/36/12/3978/1130135))  
- **rCLAMPS** ([Genome Research, 2022](https://genome.cshlp.org/content/32/9/1776))  
📂 [View folder](https://github.com/LinResearchGroup-NCSU/IDEA_Model/tree/main/supplementary_materials/other_published_models)

---

## 📚 References

- Zhang, Y., Silvernail, I., Lin, Z., Lin, X. (2024).  
  *Interpretable Protein-DNA Interactions Captured by Structure-based Optimization.*  
  **bioRxiv**. [10.1101/2024.05.26.595895](https://www.biorxiv.org/content/10.1101/2024.05.26.595895v1)

- Maerkl, S. J., & Quake, S. R. (2007).  
  *A Systems Approach to Measuring the Binding Energy Landscapes of Transcription Factors.*  
  **Science**, 315(5809), 233–237. [10.1126/science.1131007](https://doi.org/10.1126/science.1131007)
