# MotifScanner: DNA Motif Discovery and Analysis

MotifScanner is a Python-based tool for identifying transcription factor binding sites in DNA sequences. It combines traditional PWM-based motif scanning with advanced sequence analysis features, such as GC content calculation, k-mer distribution, and TF-IDF scoring. The system dynamically retrieves motifs from the HOCOMOCO database using API calls, making it efficient and user-friendly. It also provides detailed outputs and visualizations to enhance the interpretability of results.

---

## Features
- **Motif Scanning**:
  - Identify transcription factor binding motifs using Position Weight Matrices (PWMs).
  - Dynamic retrieval of motifs via HOCOMOCO API.
  - Support for species-specific queries (*Homo sapiens* and *Mus musculus*).

- **Sequence Analysis**:
  - GC content computation for DNA sequences.
  - K-mer distribution and frequency analysis.
  - TF-IDF scoring for k-mer importance.

- **Visualization**:
  - Bar plots for k-mer distributions.
  - Heatmaps for TF-IDF scores.
  - Top motifs displayed with bar charts for each sequence.

---

## Repository Structure

- `api_integration.py`: Functions for interacting with the HOCOMOCO API.
- `explore.ipynb`: Jupyter Notebook for interactive exploration and testing.
- `fasta_parser.py`: A parser for reading and handling FASTA files.
- `motif_scanner.py`: Core logic for motif scanning and feature extraction.
- `sequence_analysis.py`: Functions for analyzing sequences, including GC content, k-mers, and TF-IDF.
- `visualize_results.py`: Visualization utilities for presenting analysis results.
