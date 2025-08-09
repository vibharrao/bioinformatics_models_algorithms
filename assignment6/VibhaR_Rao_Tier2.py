import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

def load_interaction_matrix(protein_links_path):
    # Load the protein interaction data
    data = pd.read_csv(protein_links_path, sep='\t')
    
    # Create a symmetric interaction matrix
    proteins = sorted(set(data['protein1']).union(set(data['protein2'])))
    protein_index = {protein: idx for idx, protein in enumerate(proteins)}
    n_proteins = len(proteins)
    interaction_matrix = np.zeros((n_proteins, n_proteins))

    for _, row in data.iterrows():
        i = protein_index[row['protein1']]
        j = protein_index[row['protein2']]
        interaction_matrix[i, j] = row['combined_score']
        interaction_matrix[j, i] = row['combined_score']  # Make it symmetric
    
    return interaction_matrix, proteins

def apply_nmf(interaction_matrix, n_components=10):
    # Apply NMF to decompose the interaction matrix into W and H matrices
    nmf_model = NMF(n_components=n_components, init='random', random_state=42)
    W = nmf_model.fit_transform(interaction_matrix)  # Shape: (n_proteins, 10)
    H = nmf_model.components_                        # Shape: (10, n_proteins)
    return W, H

def select_least_populated_cluster(W):
    # Assign each protein to the cluster with the highest value in W
    cluster_assignments = np.argmax(W, axis=1)
    
    # Count the number of proteins in each cluster
    cluster_counts = np.bincount(cluster_assignments)
    
    # Find the cluster with the lowest protein count
    least_populated_cluster = np.argmin(cluster_counts)
    selected_proteins = np.where(cluster_assignments == least_populated_cluster)[0]
    
    return selected_proteins

def plot_cluster_heatmap(interaction_matrix, selected_proteins, proteins, output_filename='Protein_Cluster_Heatmap.png'):
    # Subset the interaction matrix for the selected proteins
    selected_interaction_matrix = interaction_matrix[np.ix_(selected_proteins, selected_proteins)]
    
    # Get the names of the selected proteins
    selected_protein_names = [proteins[i] for i in selected_proteins]
    
    # Perform hierarchical clustering on the selected subset
    row_linkage = linkage(pdist(selected_interaction_matrix), method='ward')
    
    # Create a heatmap with dendrogram using seaborn
    plt.figure(figsize=(16, 16))  # Increase figure size for better readability
    clustermap = sns.clustermap(
    selected_interaction_matrix, 
    row_linkage=row_linkage, 
    col_linkage=row_linkage,
    cmap="viridis", 
    xticklabels=selected_protein_names, 
    yticklabels=selected_protein_names,
    figsize=(16, 16), 
    cbar_kws={"label": "Interaction Strength"}
)

    # Set font size for x and y tick labels
    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), fontsize=6)
    clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_yticklabels(), fontsize=6)

    # Save the plot
    plt.savefig(output_filename, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # Load command-line arguments
    protein_links_path = sys.argv[1]
    
    # Load the interaction matrix
    interaction_matrix, proteins = load_interaction_matrix(protein_links_path)
    
    # Apply NMF
    W, H = apply_nmf(interaction_matrix, n_components=10)
    
    # Select the cluster with the fewest proteins
    selected_proteins = select_least_populated_cluster(W)
    
    # Plot the heatmap with dendrogram
    plot_cluster_heatmap(interaction_matrix, selected_proteins, proteins, output_filename='Protein_Cluster_Heatmap.png')
