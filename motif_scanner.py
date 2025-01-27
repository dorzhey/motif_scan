import pandas as pd
import requests
import numpy as np
from sequence_analysis import calculate_gc_content, calculate_kmer_distribution, compute_tf_idf_kmers

class MotifScanner:
    def __init__(self, species):
        """
        Initializes the MotifScanner and fetches all motif IDs for the specified species.

        Args:
            species (str): Either 'human' or 'mouse'.
        """
        self.species = species.lower()
        self.motif_ids = self.fetch_motif_ids()

    def fetch_motif_ids(self, arity="mono"):
        """
        Fetches all MotifIDs for the given species and arity from the HOCOMOCO API.

        Args:
            arity (str): 'mono' for mononucleotide models (default).

        Returns:
            list: A list of MotifIDs.
        """
        url = f"https://hocomoco11.autosome.org/{self.species}/{arity}.json"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            raise Exception(f"Failed to fetch motif IDs: {response.status_code} {response.text}")

    def fetch_pwm(self, motif_id):
        """
        Fetches the PWM for a given MotifID.

        Args:
            motif_id (str): The MotifID for which to fetch the PWM.

        Returns:
            numpy.ndarray: PWM as a 2D array.
        """
        url = f"https://hocomoco11.autosome.org/motif/{motif_id}/pwm.json"
        response = requests.get(url)
        if response.status_code == 200:
            pwm = response.json()
            return np.array(pwm)
        else:
            print(f"Failed to fetch PWM for {motif_id}: {response.status_code} {response.text}")
            return None

    def calculate_pwm_score(self, sequence, pwm):
        """
        Calculates the PWM score for a given sequence.

        Args:
            sequence (str): DNA sequence to score.
            pwm (numpy.ndarray): PWM matrix.

        Returns:
            float: The PWM score for the sequence.
        """
        if len(sequence) != pwm.shape[0]:
            raise ValueError(f"Sequence length must match PWM length ({pwm.shape[0]}).")

        nucleotide_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}
        score = 0
        for i, nucleotide in enumerate(sequence):
            score += pwm[i, nucleotide_to_index[nucleotide]]
        return score

    def find_motif_hits(self, sequences):
        """
        Finds motif hits in DNA sequences using PWMs.

        Args:
            sequences (dict): A dictionary of {sequence_id: DNA_sequence}.

        Returns:
            dict: {sequence_id: {motif_id: score}}
        """
        motif_hits = {}
        for motif_id in self.motif_ids[:10]:
            pwm = self.fetch_pwm(motif_id)
            if pwm is None:
                continue  # Skip if PWM fetch fails

            for seq_id, sequence in sequences.items():
                if len(sequence) < pwm.shape[0]:
                    continue  # Skip sequences shorter than the PWM length

                # Sliding window across the sequence
                max_score = None
                for i in range(len(sequence) - pwm.shape[0] + 1):
                    window = sequence[i:i+pwm.shape[0]]
                    score = self.calculate_pwm_score(window, pwm)
                    if max_score is None or score > max_score:
                        max_score = score

                if max_score is not None:
                    if seq_id not in motif_hits:
                        motif_hits[seq_id] = {}
                    motif_hits[seq_id][motif_id] = max_score

        # Retain only top 5 motifs per sequence
        for seq_id in motif_hits:
            motif_hits[seq_id] = dict(sorted(motif_hits[seq_id].items(), key=lambda x: x[1], reverse=True)[:5])

        return motif_hits

    def scan(self, fasta_file):
        """
        Scans DNA sequences for motif hits and calculates sequence features.

        Args:
            fasta_file (str): Path to the FASTA file containing sequences.

        Returns:
            pd.DataFrame: Table of sequence features and top 5 motif hits.
        """
        from fasta_parser import parse_fasta
        sequences = parse_fasta(fasta_file)

        results = []
        all_tf_idf_scores = compute_tf_idf_kmers(list(sequences.values()), k=3)

        # Find motif hits
        motif_hits = self.find_motif_hits(sequences)

        for seq_id, sequence in sequences.items():
            gc_content = calculate_gc_content(sequence)
            kmer_distribution = calculate_kmer_distribution(sequence, k=3)

            # Collect all features for this sequence
            row = {
                "sequence_id": seq_id,
                "gc_content": gc_content,
                "kmer_distribution": kmer_distribution,
                "tf_idf_scores": all_tf_idf_scores,
                "top_5_motifs": motif_hits.get(seq_id, {}),
            }
            results.append(row)

        return pd.DataFrame(results)
