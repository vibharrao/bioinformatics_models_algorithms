import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF

def load_data(dog_X_path, dog_clades_path):
    # Load the Dogs dataset
    dog_X = np.load(dog_X_path)
    dog_clades = np.load(dog_clades_path, allow_pickle=True)   
    return dog_X, dog_clades

def apply_nmf(dog_X, n_components=5):
    # Apply NMF to decompose the dataset into W and H matrices with increased max_iter
    nmf_model = NMF(n_components=n_components, init='random', random_state=42, max_iter=500)
    W = nmf_model.fit_transform(dog_X)  # W contains the contributions of each component for each dog. shape: [1355, 5]
    H = nmf_model.components_          # H contains the learned feature representation for each component. shape: [5, 784]
    return W, H

def normalize_W(W):
    # Normalize the W matrix so each row (dog) sums to 1 (proportion)
    return W / W.sum(axis=1, keepdims=True) #Divides each row in W by the sum of its elements, ensuring the row sums to 1.

def plot_dominant_clusters(W_norm, dog_clades, output_filename='NMF_Dogs.png'):
    # Determine the dominant cluster for each dog and sort by it
    dominant_clusters = np.argmax(W_norm, axis=1) # finds the highest-value cluster for each dog
    dominant_values = np.max(W_norm, axis=1) # retrieves the corresponding values
    
    # Sort dogs by dominant cluster and proportion (descending)
    sorted_indices = np.lexsort((-dominant_values, dominant_clusters))
    W_sorted = W_norm[sorted_indices]
    
    # Create a stacked bar plot
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.stackplot(range(W_sorted.shape[0]), W_sorted.T, labels=[f'Cluster {i+1}' for i in range(W_sorted.shape[1])])
    ax.set_xlabel('Dog Samples (Sorted by Dominant Cluster)')
    
    # Set the same x-limits for both axes
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())  # Corrected line to use ax.get_xlim()

    # Identify the unique clades and their indices
    unique_clades, clade_indices = np.unique(dog_clades[sorted_indices], return_index=True)
    # Mark the clade intervals with vertical lines
    for i in range(0, len(clade_indices)):
        ax2.axvline(clade_indices[i], color='gray', linestyle='--', linewidth=0.5)

    # Set x-axis ticks at the position of each new clade
    ax2.set_xticks(clade_indices)
    ax2.set_xticklabels(unique_clades, rotation=90, fontsize=8)
    ax2.set_xlabel('Dog Clades (Sorted by Dominant Cluster)', fontsize=12)
    ax.set_ylabel('Proportion')
    ax.set_title('NMF Clustering Proportion of Dogs Dataset')
    ax.legend(loc='upper right')
    plt.savefig(output_filename)


def identify_dominant_component(W_norm, dog_clades, breed='Basenji'):
    # Remove '**' prefix and normalize the breed names in dog_clades by converting to lowercase and stripping whitespace
    normalized_clades = np.array([str(clade).replace("**", "").strip().lower() for clade in dog_clades])
    
    # Normalize the target breed name
    target_breed = breed.strip().lower()
    
    # Identify the indices of the specified breed
    breed_indices = np.where(normalized_clades == target_breed)[0]
    
    # Check if the breed exists in the dataset
    if breed_indices.size == 0:
        print(f"No samples found for breed '{breed}' in the dataset.")
        return None
    
    dominant_components = np.argmax(W_norm[breed_indices], axis=1) # For the specified breed, finds the most common dominant component by counting occurrences.
    # breed_indices is an array of indices for dogs of that breed
    #dominant_components is an array where each element represents 
    #the index of the most dominant component for a specific dog in the breed
    

    # Find the most common dominant component
    unique, counts = np.unique(dominant_components, return_counts=True)
    most_common_component = unique[np.argmax(counts)]
    print(f"The dominant component for {breed} samples is Cluster {most_common_component + 1}")

    return most_common_component

if __name__ == "__main__":
    # Load command-line arguments
    dog_X_path = sys.argv[1]
    dog_clades_path = sys.argv[2]
    
    # Load data
    dog_X, dog_clades = load_data(dog_X_path, dog_clades_path)
    
    # Apply NMF
    W, H = apply_nmf(dog_X, n_components=5)
    
    # Normalize the W matrix
    W_norm = normalize_W(W)
    
    # Plot the dominant clusters as a stacked bar plot
    plot_dominant_clusters(W_norm, dog_clades, output_filename='NMF_Dogs.png')
    
    # Identify the dominant component for Basenji and Wolf samples
    basenji_dominant_component = identify_dominant_component(W_norm, dog_clades, breed='Basenji')
    wolf_dominant_component = identify_dominant_component(W_norm, dog_clades, breed='Wolf')
    
    if basenji_dominant_component is not None and basenji_dominant_component == wolf_dominant_component:
        print(f"The dominant component for Basenji and Wolf samples is the same: Cluster {basenji_dominant_component + 1}")