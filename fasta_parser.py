def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary of sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: A dictionary with sequence IDs as keys and sequences as values.
    """
    sequences = {}
    with open(file_path, 'r') as file:
        current_id = None
        current_seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]  # Skip the ">"
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    return sequences
