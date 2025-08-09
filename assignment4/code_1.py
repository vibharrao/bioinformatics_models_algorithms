import numpy as np
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt

# Load data
X = np.load("MNIST_X_subset.npy")
y = np.load("MNIST_y_subset.npy")


def visualize_digit(image):
  """Reshapes and plots a single MNIST digit image."""
  plt.imshow(image.reshape(28, 28), cmap="gray")
  plt.axis("off")
  plt.show()


# Visualize the first digit
visualize_digit(X[0])


def kmeans_clustering(X, k, filename):
  """Performs K-Means clustering and visualizes centroids."""
  kmeans = KMeans(n_clusters=k, random_state=42)
  kmeans.fit(X)

  # Reshape centroids and plot
  centroids = kmeans.cluster_centers_
  centroid_images = centroids.reshape(-1, 28, 28)

  # Adjust the number of subplots to match the number of centroids
  num_rows = int(np.ceil(k / 5))
  fig, axes = plt.subplots(num_rows, 5, figsize=(10, 5 * num_rows))

  for i, ax in enumerate(axes.flat):
    if i < k:
      ax.imshow(centroid_images[i], cmap="gray")
      ax.axis("off")
    else:
      ax.axis("off")  # Hide empty subplots

  fig.suptitle(f"K-Means Centroids (K={k})", fontsize=12)
  plt.tight_layout()
  plt.savefig(filename)


def calculate_clustering_error(X, y, labels):
  """Calculates the total number of misclassified samples."""
  error = 0
  for cluster in np.unique(labels):
    cluster_data = X[labels == cluster]
    cluster_labels = y[labels == cluster]
    most_common_label = np.bincount(cluster_labels).argmax()
    error += np.sum(cluster_labels != most_common_label)
  return error


# Perform K-Means clustering with K=10 and K=11
kmeans_clustering(X, 10, "centroids_k10.png")
kmeans_clustering(X, 11, "centroids_k11.png")

# Calculate and print clustering errors
labels_k10 = KMeans(n_clusters=10, random_state=42).fit_predict(X)
error_k10 = calculate_clustering_error(X, y, labels_k10)
print(f"K=10 Error={error_k10}")

labels_k11 = KMeans(n_clusters=11, random_state=42).fit_predict(X)
error_k11 = calculate_clustering_error(X, y, labels_k11)
print(f"K=11 Error={error_k11}")