import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def visualize_results(results):
    """
    Visualizes the results of the MotifScanner.

    Args:
        results (pd.DataFrame): The output DataFrame from MotifScanner.scan().
    """
    # GC Content Bar Chart
    plt.figure(figsize=(10, 5))
    plt.bar(results['sequence_id'], results['gc_content'], color='skyblue')
    plt.title('GC Content per Sequence')
    plt.xlabel('Sequence ID')
    plt.ylabel('GC Content (%)')
    plt.show()

    # K-mer Distribution Bar Charts
    for index, row in results.iterrows():
        kmer_dist = row['kmer_distribution']
        plt.figure(figsize=(12, 6))
        plt.bar(kmer_dist.keys(), kmer_dist.values(), color='lightgreen')
        plt.title(f'K-mer Distribution for {row["sequence_id"]}')
        plt.xlabel('K-mers')
        plt.ylabel('Frequency')
        plt.xticks(rotation=45)
        plt.show()

    # TF-IDF Scores Heatmap
    # Flatten the TF-IDF scores into a DataFrame
    tf_idf_flat = pd.DataFrame(results['tf_idf_scores'].to_list(), index=results['sequence_id']).fillna(0)
    plt.figure(figsize=(12, 8))
    sns.heatmap(tf_idf_flat, cmap='viridis', annot=False)
    plt.title('TF-IDF Scores Heatmap')
    plt.xlabel('K-mers')
    plt.ylabel('Sequence ID')
    plt.show()

    # Top 5 Motifs Bar Charts
    for index, row in results.iterrows():
        top_motifs = row['top_5_motifs']
        plt.figure(figsize=(12, 6))
        plt.bar(top_motifs.keys(), top_motifs.values(), color='orange')
        plt.title(f'Top 5 Motifs for {row["sequence_id"]}')
        plt.xlabel('Motif IDs')
        plt.ylabel('Scores')
        plt.xticks(rotation=45)
        plt.show()