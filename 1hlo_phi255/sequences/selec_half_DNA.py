# read
with open('native.seq', 'r') as file:
    content = file.read()

# select lowercase
lowercase_chars = ''.join([char for char in content if char.islower()])

# select the first half
half_length = len(lowercase_chars) // 2
first_half = lowercase_chars[:half_length]

# change letter
translation = first_half.maketrans({'e': 'A', 'l': 'G', 'j': 'C', 't': 'T'})
modified_sequence = first_half.translate(translation)

# save
with open('native_dna_half.seq', 'w') as output_file:
    output_file.write(modified_sequence)
