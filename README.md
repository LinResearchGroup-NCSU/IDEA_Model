# IDEA_Model
Interpretable protein-DNA Energy Associative (IDEA) Model

This is a cleaned-up version of the IDEA model. Training and testing data are provided, using 3 human MAX transcription factors (TFs) as an example. The final result here is used in Figure S1 in our IDEA manuscript.

## Max Protein Complexes Energy Model Training Guide

This guide will help you train an energy model for three human Max proteins (PDB IDs: 1hlo, 1nlw, and 1nkp).

### Step 1: Prepare Input Files

1. **Protein List**:
   Create a file named `3max_case/proteinList.txt` and list the PDB IDs of the proteins you want to train, one per line. For example:
   ```
   1hlo
   1nlw
   1nkp
   ```

2. **PDB Files**:
   Place the PDB files in the directory `3max_case/PDBs`. Each PDB file should be named as `{PDB ID}_modified.pdb`. For example, `1hlo_modified.pdb`. Rename all protein chains in the PDB files to chain A, and the DNA chains to B and C. This ensures consistency across all files.

### Step 2: Start the Training

Open a terminal and run the following command to start the training process:
```bash
bash train.sh
```

### Step 3: Configuration Details

1. **Choosing Interaction Atoms**:
   Go to the file `3max_case/common_functions/common_function.py` and find the function `get_interaction_atom(residue)`. This function lets you choose which atoms to use for training. For this guide, we select C5 atoms from DNA and CA atoms from protein.

2. **Setting the Cutoff Mode**:
   Go to the file `3max_case/optimization/for_training_gamma/optimize_gamma.py`. This file contains a setting called `cutoff_mode`. The `cutoff_mode` helps filter out noise by keeping the first 70 eigenvalues and replacing the rest with the 70th eigenvalue. This helps improve the accuracy of the model.

### Step 4: Understanding the Output

After the training is complete, the results will be saved in the directory `3max_case/optimization/for_training_gamma/gammas/randomized_decoy`.

Key output files include:
- **Filtered Energy Model**: This file is named `native_trainSetFiles_phi_pairwise_contact_well-8.0_8.0_0.7_10_gamma_filtered`. It contains the trained energy model after filtering.
- **Other Important Files**: You will also find additional files such as `phi`, `A`, and `B` in the same folder. These files are part of the output from the training process.

## Testing Protein-DNA Binders and Calculating Binding Energy

Once you have the trained energy model, you may want to generate the `phi` for testing protein-DNA binders and calculate the associated binding energy. The `1hlo_phi255` directory provides an example of generating the `phi` for 255 testing binders corresponding to the Max Mitomi data reported in Maerkl and Quake (2007).

### Step 1: Prepare Testing Input Files

1. **PDB Files**:
   Place the PDB files in the directory `1hlo_phi255/PDBs`.

2. **DNA Sequence File**:
   The DNA sequence file is `1hlo_phi255/sequences`.

### Step 2: Generate Testing Sequences

The input for this program starts from `dna_half.seq`. The scripts `reverse_complement.py` and `merge.py` help generate the full testing sequence. These steps are not necessary, but please ensure that your testing sequences are consistent with the DNA sequence in the native structure file (PDB).

### Step 3: Generate the Testing Phi

Run the following command to generate the testing `phi`:
```bash
bash 1hlo_phi255.sh
```
The generated testing `phi` will be in `1hlo_phi255/phis`, named `phi_pairwise_contact_well_native_decoys_CPLEX_randomization_-8.0_8.0_0.7_10`. It should have 255 lines as we have 255 testing sequences. The `phi_pairwise_contact_well_native_native_-8.0_8.0_0.7_10` file is the `phi` for the native structure (only one line as we have one structure as input).

### Step 4: Calculate Predicted Energy

Calculate the predicted energy based on the gamma file generated in `3max_case` and the testing `phi` file using the formula  E = γΦ.

1. **Run the Python Script**:
   Navigate to the `testing_energy` directory and run the script `calculate_testing_energy.py` to calculate the energy.
   ```bash
   python calculate_testing_energy.py
   ```
   The generated energy will be in `Energy_mg.txt`. 

## References

Maerkl, S. J., & Quake, S. R. (2007). A Systems Approach to Measuring the Binding Energy Landscapes of Transcription Factors. *Science*, 315(5809), 233-237. DOI: [10.1126/science.1131007](https://doi.org/10.1126/science.1131007)
```
