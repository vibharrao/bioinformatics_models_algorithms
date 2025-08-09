import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from seaborn import clustermap

def load_interaction_matrix(protein_links_path):
    # Load the interaction data
    data = pd.read_csv(protein_links_path, sep='\t')
    
    # Extract unique proteins
    proteins = np.unique(data[['protein1', 'protein2']].values)
    protein_index = {protein: i for i, protein in enumerate(proteins)}
    
    # Initialize an empty symmetric matrix for interactions
    n_proteins = len(proteins)
    interaction_matrix = np.zeros((n_proteins, n_proteins))
    
    # Populate the matrix with interaction scores
    for _, row in data.iterrows():
        i, j = protein_index[row['protein1']], protein_index[row['protein2']]
        interaction_matrix[i, j] = row['combined_score']
        interaction_matrix[j, i] = row['combined_score']  # Symmetric assumption
    
    return interaction_matrix, proteins

def apply_nmf_to_protein_data(interaction_matrix, n_components=10):
    # Apply NMF with 10 components
    nmf_model = NMF(n_components=n_components, init='random', random_state=42)
    W = nmf_model.fit_transform(interaction_matrix)
    H = nmf_model.components_
    return W, H

def find_cluster_with_fewest_proteins(W):
    # Assign each protein to the cluster with the highest value in W
    cluster_assignments = np.argmax(W, axis=1)
    
    # Find the cluster with the fewest assigned proteins
    unique, counts = np.unique(cluster_assignments, return_counts=True)
    smallest_cluster = unique[np.argmin(counts)]
    
    # Get indices of proteins in the smallest cluster
    proteins_in_smallest_cluster = np.where(cluster_assignments == smallest_cluster)[0]
    
    return proteins_in_smallest_cluster, smallest_cluster

def visualize_protein_cluster(interaction_matrix, proteins_in_smallest_cluster, output_filename='Protein_Cluster_Heatmap.png'):
    # Extract the submatrix for proteins in the selected cluster
    cluster_interactions = interaction_matrix[np.ix_(proteins_in_smallest_cluster, proteins_in_smallest_cluster)]
    
    # Perform hierarchical clustering on the selected cluster
    linkage_matrix = linkage(pdist(cluster_interactions), method='average')
    
    # Plot the heatmap with dendrogram
    sns.clustermap(cluster_interactions, row_linkage=linkage_matrix, col_linkage=linkage_matrix, cmap='viridis')
    plt.savefig(output_filename)
    plt.close()

if __name__ == "__main__":
    import sys
    protein_links_path = sys.argv[1]
    
    # Load data and interaction matrix
    interaction_matrix, proteins = load_interaction_matrix(protein_links_path)
    
    # Apply NMF
    W, H = apply_nmf_to_protein_data(interaction_matrix, n_components=10)
    
    # Find the cluster with the fewest proteins
    proteins_in_smallest_cluster, smallest_cluster = find_cluster_with_fewest_proteins(W)
    
    # Visualize the selected cluster
    visualize_protein_cluster(interaction_matrix, proteins_in_smallest_cluster, output_filename='Protein_Cluster_Heatmap1.png')
    print(f"Protein interaction heatmap saved as 'Protein_Cluster_Heatmap1.png' for cluster {smallest_cluster + 1}.")
