# IDEA_Model
Interpretable protein-DNA Energy Associative (IDEA) Model

This is a cleaned-up version of the IDEA model
Training and testing data are provided, using the 3 human MAX TF as an example

Dataset：https://drive.google.com/drive/folders/10f85feDBfOjVS0YkCJaNOo8JDob6qWe0?usp=drive_link

## Max Protein Complexes Energy Model Training Guide

This guide walks you through the steps to train an energy model for three human Max proteins (PDB IDs: 1hlo, 1nlw, and 1nkp).

### Step 1: Prepare Input Files

1. **Protein List**:
   - Create a file named `3max_case/proteinList.txt`.
   - List the PDB IDs of the proteins you want to train, one per line. For example:
     ```
     1hlo
     1nlw
     1nkp
     ```

2. **PDB Files**:
   - Place the PDB files in the directory `3max_case/PDBs`.
   - Each PDB file should be named as `{PDB ID}_modified.pdb`. For example, `1hlo_modified.pdb`.
   - Rename all protein chains in the PDB files to chain A, and the DNA chains to B and C. This ensures consistency across all files.

### Step 2: Start the Training

Open a terminal and run the following command to start the training process:
```bash
bash train.sh
```

### Step 3: Configuration Details

1. **Choosing Interaction Atoms**:
   - Go to the file `3max_case/common_functions/common_function.py`.
   - Find the function `get_interaction_atom(residue)`.
   - This function lets you choose which atoms to use for training. For this guide, we select C5 atoms from DNA and CA atoms from protein.

2. **Setting the Cutoff Mode**:
   - Go to the file `3max_case/optimization/for_training_gamma/optimize_gamma.py`.
   - This file contains a setting called `cutoff_mode`.
   - The `cutoff_mode` helps filter out noise by keeping the first 70 eigenvalues and replacing the rest with the 70th eigenvalue. This helps improve the accuracy of the model.

### Step 4: Understanding the Output

After the training is complete, the results will be saved in the directory `3max_case/optimization/for_training_gamma/gammas/randomized_decoy`.

Key output files include:
- **Filtered Energy Model**: This file is named `native_trainSetFiles_phi_pairwise_contact_well-8.0_8.0_0.7_10_gamma_filtered`. It contains the trained energy model after filtering.
- **Other Important Files**: You will also find additional files such as `phi`, `A`, and `B` in the same folder. These files are part of the output from the training process.

By following these steps, you can successfully train an energy model for the specified Max protein complexes.

---
