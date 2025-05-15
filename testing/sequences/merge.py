def merge_files(file1, file2, output_file):
    with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output_file, 'w') as out:
        for line1, line2 in zip(f1, f2):
            merged_line = line1.strip() + line2.strip()
            out.write(merged_line + '\n')

merge_files('dna_half.seq', 'dna_half_complement.seq', 'dna.seq')

