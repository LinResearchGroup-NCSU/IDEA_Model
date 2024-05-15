def replace_and_generate_file(native_file, modeller_file, output_file):
    with open(native_file, 'r') as f:
        native_content = f.readline().strip()

    import re
    lowercase_match = re.search(r'[a-z]+', native_content)
    if not lowercase_match:
        print("No lowercase sequence found in native.seq!")
        return
    lowercase_sequence = lowercase_match.group(0)

    with open(modeller_file, 'r') as f:
        modeller_lines = f.readlines()

    with open(output_file, 'w') as out_f:
        for modeller_line in modeller_lines:
            modeller_line = modeller_line.strip()
            new_content = native_content.replace(lowercase_sequence, modeller_line)
            out_f.write(new_content + "\n")


replace_and_generate_file("native.seq", "dna_modeller.seq", "native.decoys")

