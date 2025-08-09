from sklearn.cluster import AgglomerativeClustering
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import Counter

# Load MNIST data
MNIST_X = np.load('MNIST_X_subset.npy')
MNIST_y = np.load('MNIST_y_subset.npy')

# Perform hierarchical clustering with average linkage
hierarchical_clustering = AgglomerativeClustering(n_clusters=10, linkage='average')
cluster_labels = hierarchical_clustering.fit_predict(MNIST_X)

# Calculate clustering error
def calculate_clustering_error_hierarchy(labels, true_labels, k):
    total_error = 0
    for i in range(k):
        cluster_indices = np.where(labels == i)[0]
        true_labels_in_cluster = true_labels[cluster_indices]
        if len(true_labels_in_cluster) > 0:  # Ensure the cluster is not empty
            most_common_label = Counter(true_labels_in_cluster).most_common(1)[0][0]
            error = np.sum(true_labels_in_cluster != most_common_label)
            total_error += error
    return total_error

error_hierarchy = calculate_clustering_error_hierarchy(cluster_labels, MNIST_y, 10)
print(f'Hierarchical Clustering Error = {error_hierarchy}')

# Get the most common ground truth class for each of the 10 clusters from clustering labels
mostcommonlabel = []
for i in range(10):
    cluster_indices = np.where(cluster_labels == i)[0]
    true_labels_in_cluster = MNIST_y[cluster_indices]
    if len(true_labels_in_cluster) > 0:  # Check if the cluster has any data points
        most_common_label = Counter(true_labels_in_cluster).most_common(1)[0][0]
        mostcommonlabel.append(most_common_label)

# Generate the linkage matrix for dendrogram visualization
Z = linkage(MNIST_X[:1000], 'average')  # Use a subset to reduce computational load

# Plot the dendrogram first to get the correct number of leaf nodes
plt.figure(figsize=(10, 7))
dendro = dendrogram(Z, truncate_mode='lastp', p=10, show_leaf_counts=True)

# Get the current x-ticks from the dendrogram
ax = plt.gca()  # Get the current axis
ax.get_yaxis().set_visible(False)  # Turn off the y-axis labels and ticks

# Ensure the length of tick labels matches the number of leaves
if len(ax.get_xticks()) == len(mostcommonlabel):
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(mostcommonlabel, rotation=45, ha='right')  # Set the x-tick labels to the most common digits
else:
    print("Warning: The number of x-tick positions does not match the number of clusters.")

plt.title('Hierarchical Clustering Dendrogram (Average Linkage)')
plt.savefig('MNIST_dendrogram.png')
plt.show(block=False)  # Non-blocking mode, the script won't wait for you to close the plot
plt.close('all')
