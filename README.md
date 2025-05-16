# IDEA Model: Interpretable Protein-DNA Energy Associative Model

The **Interpretable Protein-DNA Energy Associative (IDEA) Model** is a computational framework that learns protein-DNA physicochemical interactions by fusing available crystal structures and their associated sequences into an optimized energy model. We show that the model can be used to accurately predict the sequence-specific DNA binding affinities of DNA-binding proteins and is transferable across the same protein superfamily. This repository provides a clean implementation of the IDEA model, with training and testing data for three human MAX transcription factors (PDB IDs: 1hlo, 1nlw, 1nkp). Results are used in **Figure S1** of the IDEA manuscript.

## Features

- Train and visualize IDEA energy models for protein-DNA complexes.
- Generate phi values and calculate binding energies for testing binders.
- Supplementary materials.

## Installation

- **Python 3**: Download from [python.org](https://www.python.org/downloads/).
- **Python Packages**:
  ```bash
  # Using pip
  pip install biopython numpy mdtraj joblib

  # Or using conda
  conda install -c conda-forge biopython numpy mdtraj joblib
  ```

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/LinResearchGroup-NCSU/IDEA_Model
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   cd IDEA_Model
   conda env create -f IDEA.yaml -n IDEA
   conda activate IDEA
   ```
2. Ensure scripts are executable:
   ```bash
   chmod +x *.sh
   ```
  **Note**: Our old version used a single script `buildseq.py`, from the Modeller package to extract the protein and DNA sequences from the given PDB structure. To avoid further confusion and reduce the environment dependencies required to run the code, we rewrote it using basic Python code. The most time-consuming step of the IDEA is generating the phi values for decoy/testing binders. We rewrote the `template_evaluate_phi.py` by adopting parallel computation from joblib and multiprocessing, which is especially useful when we want to do genome-wide prediction. We also simplified the redundant parts of the code and added informative messages during model training and testing to make it more user-friendly.

## Usage

### Training the Energy Model

Train an energy model for three MAX protein complexes (1hlo, 1nlw, 1nkp).

1. **Prepare Input Files**:
   - Create `training/proteinList.txt` with PDB IDs:
     ```
     1hlo
     1nlw
     1nkp
     ```
   - Place PDB files in `training/PDBs`, named `{PDB ID}_modified.pdb` (e.g., `1hlo_modified.pdb`). Rename protein chains to **A** and DNA chains to **B** and **C**, respectively.

2. **Run Training**:
   ```bash
   bash train.sh
   ```

3. **Configure Settings**:
   - **Interaction Atoms**: Edit `get_interaction_atom` in `common_function.py` to implement a **coarse-grained interaction scheme**, where DNA bases are represented by the **C5** atom (or **P** if backbone-level resolution is preferred), and protein residues are represented by the **Cα (CA)** atom, although **side-chain atoms** may also be used in cases where more detailed interactions are of interest.
   - **Decoy Number**: In `training/optimization/for_bindingE/template/sequences`, the scripts generate_decoy_seq_prot.py and generate_decoy_seq_DNA.py are used to generate protein and DNA decoy sequences for training. The default sizes are 10,000 and 1,000, respectively, which we found to be robust across all tests in our manuscript. Please adjust these values based on your specific needs.
   - **Cutoff Mode**: In `training/optimization/for_training_gamma/optimize_gamma.py`, set cutoff_mode = 70 to retain the first 70 eigenvalues, replacing all others with the 70th eigenvalue. This choice typically depends on the lambda values in `training/optimization/for_training_gamma/gammas/randomized_decoy/native_trainSetFiles_phi_pairwise_contact_well-8.0_8.0_0.7_10_lamb`.

4. **Output**:
   - Results are saved in `training/optimization/for_training_gamma/gammas/randomized_decoy`.
   - Key file: `native_trainSetFiles_phi_pairwise_contact_well-8.0_8.0_0.7_10_gamma_filtered` (filtered energy model).

5. **Visualization**:
   - Run:
     ```bash
     python training/optimization/for_training_gamma/visualize.py
     ```
   - Plots are saved in `training/optimization/for_training_gamma/visualize`.

### Predicting Protein-DNA Binding Energies

Generate phi values and calculate binding energies for given testing binders (e.g., Max 255 mutated binders testing dataset in Maerkl, S. J et al.).

1. **Prepare Input Files**:
   - Place PDB files in `testing/PDBs`.
   - DNA (testing) sequences are in `testing/sequences`.

2. **Generate Testing Sequences**:
   - Use `dna_half.seq` as input. Optionally, run `reverse_complement.py` and `merge.py` to generate full double-stranded DNA sequences.
   - Ensure the length of sequences matches the native PDB structure.

3. **Generate Testing Phi**:
   ```bash
   bash test.sh 1hlo
   ```
   - Output in `testing/phis`:
     - `phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-8.0_8.0_0.7_10` (255 lines for testing sequences).
     - `phi_pairwise_contact_well_native_native_-8.0_8.0_0.7_10` (1 line for native structure).

4. **Calculate Binding Energy**:
   - Navigate to `energy_calculation` and run:
     ```bash
     python calculate_testing_energy.py
     ```
   - Output: `Energy_mg.txt` (predicted energies using **E = γΦ**).

## Supplementary Materials

- **Trained Energy Models**: [IDEA_trained_energy_models](https://github.com/LinResearchGroup-NCSU/IDEA_Model/tree/main/supplementary_materials/IDEA_trained_energy_models)
- **Raw Data**: [raw_data.zip](https://github.com/LinResearchGroup-NCSU/IDEA_Model/blob/main/supplementary_materials/raw_data.zip)
- **Processed Published Models**: Scripts for comparing with DBD-hunter and rCLAMPS: [other_published_models](https://github.com/LinResearchGroup-NCSU/IDEA_Model/tree/main/supplementary_materials/other_published_models)

## References

- Zhang, Y., Silvernail, I., Lin, Z., Lin, X. (2024). Interpretable Protein-DNA Interactions Captured by Structure-based Optimization. *bioRxiv*. [DOI:10.1101/2024.05.26.595895](https://www.biorxiv.org/content/10.1101/2024.05.26.595895v1)
- Maerkl, S. J., & Quake, S. R. (2007). A Systems Approach to Measuring the Binding Energy Landscapes of Transcription Factors. *Science*, 315(5809), 233–237. [DOI:10.1126/science.1131007](https://doi.org/10.1126/science.1131007)

## Contact

For questions or support, contact the [Lin Research Group](https://github.com/LinResearchGroup-NCSU) or open an issue on GitHub.
