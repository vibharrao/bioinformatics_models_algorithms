import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter

# Load the dataset
MNIST_X = np.load('MNIST_X_subset.npy')
MNIST_y = np.load('MNIST_y_subset.npy')

# Reshape and visualize the first example (index 0)
first_image = MNIST_X[0].reshape(28, 28)
plt.imshow(first_image, cmap='gray')
plt.title(f'Label: {MNIST_y[0]}')
plt.savefig('first_image.png')  # Save the image to a file
plt.close()  # This will close the figure


# Perform K-Means clustering for K=10 and K=11
kmeans_10 = KMeans(n_clusters=10, random_state=42).fit(MNIST_X)
kmeans_11 = KMeans(n_clusters=11, random_state=42).fit(MNIST_X)

def visualize_centroids(kmeans, k, filename):
    centroids = kmeans.cluster_centers_
    fig, axes = plt.subplots(2, 5 if k == 10 else 6, figsize=(10, 5))
    for i, ax in enumerate(axes.flat):
        if i < k:
            ax.imshow(centroids[i].reshape(28, 28), cmap='gray')
        ax.axis('off')
    plt.savefig(filename)
    plt.close()


visualize_centroids(kmeans_10, 10, 'centroids_k10.png')
visualize_centroids(kmeans_11, 11, 'centroids_k11.png')


def calculate_clustering_error(kmeans, labels, k):
    total_error = 0
    for i in range(k):
        cluster_indices = np.where(kmeans.labels_ == i)[0]
        true_labels = labels[cluster_indices]
        most_common_label = Counter(true_labels).most_common(1)[0][0]
        error = np.sum(true_labels != most_common_label)
        total_error += error
    return total_error

error_10 = calculate_clustering_error(kmeans_10, MNIST_y, 10)
error_11 = calculate_clustering_error(kmeans_11, MNIST_y, 11)

print(f'K=10 Error={error_10}')
print(f'K=11 Error={error_11}')

