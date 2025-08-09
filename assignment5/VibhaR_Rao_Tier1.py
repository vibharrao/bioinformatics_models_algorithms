import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.preprocessing import LabelEncoder
from mpl_toolkits.mplot3d import Axes3D


# Set random state for reproducibility
random_state = 42


# Part 1: PCA on MNIST Dataset
# Load MNIST dataset
mnist_X = np.load('MNIST_X_subset.npy')
mnist_y = np.load('MNIST_y_subset.npy')


# Perform PCA to reduce dimensions to 2
pca_mnist = PCA(n_components=2, random_state=random_state)
mnist_X_pca_2D = pca_mnist.fit_transform(mnist_X)


# Save 2D PCA-transformed plot
plt.figure(figsize=(8, 6))
scatter = plt.scatter(mnist_X_pca_2D[:, 0], mnist_X_pca_2D[:, 1], c=mnist_y, cmap='tab10', s=10)
plt.colorbar(scatter, ticks=range(10), label="Digit Labels")
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('MNIST PCA 2D Projection')
plt.savefig('MNIST_PCA_2D.png')
plt.close()


# Reconstruct an image using the first 2 principal components
mnist_img_reduced = pca_mnist.transform(mnist_X[0].reshape(1, -1))
mnist_img_reconstructed = pca_mnist.inverse_transform(mnist_img_reduced).reshape(28, 28)


# Save reconstructed image and original image
plt.imsave('MNIST_reconstructed_2PC.png', mnist_img_reconstructed, cmap='gray')
plt.imsave('MNIST_original.png', mnist_X[0].reshape(28, 28), cmap='gray')


# Manually choose a coordinate to resemble digit "1" in PCA-transformed space and reconstruct
chosen_coord = np.array([-800, 600])  # coordinates based on 2D PCA scatter plot for digit "1"
mnist_img_from_coord = pca_mnist.inverse_transform(chosen_coord).reshape(28, 28)
plt.imsave('MNIST_reconstructed_1_from_coord.png', mnist_img_from_coord, cmap='gray')


# Part 2: PCA on Dogs SNP Dataset
# Load Dogs SNP dataset
dogs_X = np.load('dogs_X.npy')
dogs_clades = np.load('dogs_clades.npy', allow_pickle=True)  # Use allow_pickle=True for object arrays


# Convert categorical clade labels to numerical codes
label_encoder = LabelEncoder()
dogs_clades_numerical = label_encoder.fit_transform(dogs_clades)


# Perform PCA to reduce dimensions to 2
pca_dogs = PCA(n_components=2, random_state=random_state)
dogs_X_pca_2D = pca_dogs.fit_transform(dogs_X)


# Save 2D PCA-transformed plot
plt.figure(figsize=(8, 6))
scatter = plt.scatter(dogs_X_pca_2D[:, 0], dogs_X_pca_2D[:, 1], c=dogs_clades_numerical, cmap='viridis', s=10)
plt.colorbar(scatter, label="Clade Labels (Encoded)")
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('Dogs SNP PCA 2D Projection')
plt.savefig('Dogs_PCA_2D.png')
plt.close()


# Part 3: MDS on Molecular Distance Matrix
# Load molecular distance matrix, dropping unnecessary columns
molecule_data = pd.read_csv('molecule_distances.tsv', sep='\t')
molecule_distances = molecule_data.drop(columns=['Atom Index', 'Element'])


# Perform MDS to reconstruct 3D coordinates
mds = MDS(n_components=3, dissimilarity='euclidean', random_state=random_state, metric = True)
molecule_3D_coords = mds.fit_transform(molecule_distances)


# Save the coordinates to CSV
molecule_coords_df = pd.DataFrame(molecule_3D_coords, columns=['X', 'Y', 'Z'])
molecule_coords_df['Atom'] = molecule_data['Atom Index']
molecule_coords_df['Element'] = molecule_data['Element'] 
molecule_coords_df.to_csv('molecule_coordinates.csv', index=False)

