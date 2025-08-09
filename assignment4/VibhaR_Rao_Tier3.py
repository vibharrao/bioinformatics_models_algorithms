import numpy as np
from sklearn.cluster import AgglomerativeClustering
from collections import Counter
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# Load Dogs data
dogs_X = np.load('dogs_X.npy')
dogs_clades = np.load('dogs_clades.npy', allow_pickle=True)

# Perform hierarchical clustering with average linkage
hierarchical_clustering_dogs = AgglomerativeClustering(n_clusters=30, linkage='average')
cluster_labels_dogs = hierarchical_clustering_dogs.fit_predict(dogs_X)

# Function to calculate clustering error
def calculate_clustering_error(labels, true_labels, k):
    total_error = 0
    for i in range(k):
        cluster_indices = np.where(labels == i)[0]
        true_labels_in_cluster = true_labels[cluster_indices]
        most_common_label = Counter(true_labels_in_cluster).most_common(1)[0][0]
        error = np.sum(true_labels_in_cluster != most_common_label)
        total_error += error
    return total_error

# Compute the clustering error
error_hierarchy_dogs = calculate_clustering_error(cluster_labels_dogs, dogs_clades, 30)
print(f'Hierarchical Clustering Error for Dogs Dataset: {error_hierarchy_dogs}')

# Find the most common clade for each cluster
most_common_clades = []
for i in range(30):
    cluster_indices = np.where(cluster_labels_dogs == i)[0]
    true_clades_in_cluster = dogs_clades[cluster_indices]
    if len(true_clades_in_cluster) > 0:  # Check if the cluster is not empty
        most_common_clade = Counter(true_clades_in_cluster).most_common(1)[0][0]
        most_common_clades.append(most_common_clade)
    else:
        most_common_clades.append('None')  # Handle empty clusters

# Generate the linkage matrix
Z_dogs = linkage(dogs_X, 'average')

# Plot the dendrogram
plt.figure(figsize=(12, 8))
dendro_dogs = dendrogram(Z_dogs, truncate_mode='lastp', p=30, show_leaf_counts=True)

# Get the current x-ticks from the dendrogram
ax = plt.gca()  # Get the current axis
x_tick_positions = ax.get_xticks()
ax.get_yaxis().set_visible(False)

# Ensure the length of tick labels matches the number of leaves (30 in this case)
if len(x_tick_positions) == len(most_common_clades):
    ax.set_xticks(x_tick_positions)
    ax.set_xticklabels(most_common_clades, rotation=45, ha='right')  # Set the x-tick labels to the most common clades
else:
    print("Warning: The number of x-tick positions does not match the number of clusters.")

# Set the title and save the figure
plt.title('Hierarchical Clustering Dendrogram for Dogs Dataset (Average Linkage)')
plt.savefig('Dogs_dendrogram.png')
plt.show(block=False)  # Non-blocking mode, the script won't wait for you to close the plot
plt.close('all')
