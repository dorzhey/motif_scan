from collections import Counter

def calculate_gc_content(sequence):
    """
    Calculates the GC content of a DNA sequence.

    Args:
        sequence (str): A DNA sequence.

    Returns:
        float: GC content as a percentage of the total sequence length.
    """
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0.0

def calculate_kmer_distribution(sequence, k=3):
    """
    Calculates the k-mer distribution of a DNA sequence.

    Args:
        sequence (str): A DNA sequence.
        k (int): Length of the k-mers.

    Returns:
        dict: A dictionary with k-mers as keys and their counts as values.
    """
    kmer_counts = Counter(sequence[i:i+k] for i in range(len(sequence) - k + 1))
    return dict(kmer_counts)

def compute_tf_idf_kmers(sequences, k=3):
    """
    Computes the TF-IDF scores for k-mers across multiple sequences.

    Args:
        sequences (list): A list of DNA sequences.
        k (int): Length of the k-mers.

    Returns:
        dict: A dictionary of k-mers with their TF-IDF scores.
    """
    from math import log
    kmer_doc_count = Counter()
    kmer_counts_per_seq = []

    # Count k-mer occurrences per sequence and globally
    for seq in sequences:
        kmers = calculate_kmer_distribution(seq, k)
        kmer_counts_per_seq.append(kmers)
        kmer_doc_count.update(kmers.keys())

    total_docs = len(sequences)
    tf_idf_scores = {}

    # Compute TF-IDF for each sequence
    for doc_kmers in kmer_counts_per_seq:
        doc_scores = {}
        total_kmers_in_doc = sum(doc_kmers.values())
        for kmer, count in doc_kmers.items():
            tf = count / total_kmers_in_doc  # Term Frequency
            idf = log(1 + (total_docs / (1 + kmer_doc_count[kmer])))  # Adjusted IDF
            doc_scores[kmer] = tf * idf
        tf_idf_scores.update(doc_scores)

    return tf_idf_scores

