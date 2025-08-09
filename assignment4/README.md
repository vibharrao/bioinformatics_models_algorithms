# Assignment 4

## Description

This assignment involves performing unsupervised learning on a subset of the MNIST dataset and, optionally, the Dogs SNP dataset. 
The primary objective is to apply K-Means and hierarchical clustering to these datasets, evaluate the clustering performance, and explore variations to improve results.

### Background on MNIST

The MNIST dataset is a widely used image dataset in the engineering and computer science communities. It is publicly available at [The MNIST Database](https://yann.lecun.com/exdb/mnist/).

Each image in the MNIST dataset is a grayscale image of size **28x28 pixels**, resulting in **784** pixels when vectorized into a 1D array. The dataset contains images of handwritten digits from 0 to 9, encompassing 10 distinct classes. Due to its simplicity and the ease of visualizing digits, MNIST is widely used for testing and demonstrating image processing and machine learning techniques.

Although MNIST consists of image data, the unsupervised learning techniques applied in this assignment, such as K-Means and hierarchical clustering, are general-purpose methods. They can be applied to any dataset, including biological datasets, making the skills learned here transferable to many domains.

### Background on Dogs SNP Dataset

In this assignment, the Dogs SNP Dataset provides an optional component where you will work with **genetic data** from dogs. The dataset consists of **784 SNP features** for each of 1355 dog samples, representing variations in their DNA at specific positions in the genome.

The samples are grouped into **clades** (genetically related groups of breeds), provided in `dogs_clades.npy`, which serve as the true labels for clustering evaluation.

---

### Assignment (Tier 1)

1. **Load a subset of the MNIST dataset and visualize the first example**:
   - Load the provided NumPy arrays `MNIST_X_subset.npy` and `MNIST_y_subset.npy`.
   - This subset consists of **6000 images** from the MNIST dataset, with **600 images for each digit** (0-9).
   - Reshape the first example (index 0) and visualize it.

2. **Perform K-Means clustering** on the dataset with K means. Your script should perform clustering with K=10 AND K=11.

3. **Visualize the centroids** of the clusters:
   - Reshape each centroid into the original image format (28x28 pixels).
   - Save the plot of the centroids as an image file.
   - Your program should produce two centroid image files: `centroids_k10.png` and `centroids_k11.png`

4. **Compute the clustering error**:
   - For each cluster:
     - Determine the majority actual digit label.
     - Count the number of samples in the cluster that do **not** match the majority label.
   - Sum these counts across all clusters to obtain the **total clustering error**.
   - Report the clustering error in terms of the number of misclassified samples.

5. **Output the clustering error**:
   - Print your clustering error to stdout in the following format: 
      
      ```
      K=10 Error=123
      K=11 Error=123
      ```

      Where 123 is the actual error you calculated.


**Execution Command**:

    ```
    python Firstname_Lastname_Tier1.py
    ```

**Output**:
 - Centroid image files: `centroids_k10.png` and `centroids_k11.png`
 - K and error printed to stdout.

---

### Extra Credit 1 (Tier 2)

*Complete this if all of the above is finished in less than 4 hours.*

**Hierarchical Clustering on MNIST** 

   1. Perform hierarchical clustering on the MNIST dataset with **k=10** with average linkage.
   2. Output the clustering error.
   3. Visualize the hierarchy of clusters as a dendrogram (tree), label each of the 10 terminal clusters with the most common ground truth class digit at that location, and save it as `MNIST_dendrogram.png`.
   4. **Provide a SHORT (2-3 sentences) explanation** of the hierarchy of clusters that you observe (or some part of it), and discuss whether it aligns with your expectations for the digits.

**Execution Command**:

    ```
    python Firstname_Lastname_Tier2.py
    ```

**Output**:

- Hierarchical clustering plot (`MNIST_dendrogram.png`).
- Explanation saved in `MNIST_paragraph.txt`.



### Extra Credit 2 (Tier 3)

*Do this if you finish all of the above in less than 4 hours.*

**Hierarchical Clustering on Dogs Dataset**

1. **Load the Dogs dataset**:

   - `dogs_X.npy`: SNP data of 1355 dog samples, each with 784 SNP features.
   - `dogs_clades.npy`: Clade information for each dog sample.

2. **Perform hierarchical clustering** on the Dogs dataset.

3. **Compute the clustering error**:

   - Use the clade information as the true labels.
   - Follow the same error computation method as before.

4. **Visualize the dendrogram** of the hierarchical clustering (save as `Dogs_dendrogram.png`), labeling each of the 30 terminal nodes with the most common clade.

**Execution Command**:

    ```
    python Firstname_Lastname_Tier3.py
    ```

**Output**:

- Hierarchical clustering plot (`Dogs_dendrogram.png`).

---

## Assignment Requirements

### Data

- **MNIST Dataset**:

  - Provided as `MNIST_X_subset.npy` and `MNIST_y_subset.npy`.
  - Contains 6000 images (600 per digit), each of size 28x28 pixels (784 when vectorized).

- **Dogs Dataset** (Tier 3):

  - Provided as `dogs_X.npy` and `dogs_clades.npy`.
  - Contains SNP data for 1355 dog samples with 784 SNP features.

### Execution

- The scripts should be executed using the specified commands.
- Outputs may vary due to different initializations, but should be within the expected range for validation.
- Ensure that all code is well-documented and generates the expected outputs.

---

Submission Requirements:

**Note**: the scripts should output all expected plots and clustering errors mentioned before. You also have to independently submit the plots that the script generates.

1. You **must** submit your solution and outputs for Tier 1 with the following naming scheme:

        Firstname_Lastname_Tier1.py, centroids_k10.png, centroids_k11.png

2. If you attempted Tier 2 for extra credit, **you must submit all of the above, plus the following**:

        Firstname_Lastname_Tier2.py, MNIST_dendrogram.png, MNIST_paragraph.txt

3. If you attempted Tier 3, **you must submit all of the above, plus the following**:

        Firstname_Lastname_Tier3.py, Dogs_dendrogram.png


## Hints

1. **Libraries**:

   - Use `numpy`, `scikit-learn`, and `matplotlib` for computations and visualizations.
   - For K-Means clustering, use `KMeans` from `sklearn.cluster`.
   - For hierarchical clustering, use `scipy.cluster.hierarchy` or `AgglomerativeClustering` from `sklearn.cluster`.

2. **Data Preprocessing**:

   - Flatten the 2D images into 1D arrays if necessary (for MNIST data).

3. **Visualization**:

   - Reshape centroids back into 28x28 pixel images for visualization.
   - Use `matplotlib` to display and save images and plots.
   - For dendrograms, use `scikit-learn` or `scipy.cluster.hierarchy.dendrogram`.

4. **Clustering Error Calculation**:

   - For each cluster:
     - Identify the most common true label.
     - Count misclassified samples within the cluster.
   - Sum the misclassifications across all clusters for the total error.

5. **Hierarchical Clustering**:

   - Optionally, experiment with different linkage methods (e.g., 'ward', 'complete') and distance metrics.

---

## Notes

- Fix random seeds (e.g., use `random_state` in `KMeans`) for reproducibility.
- Comment your code to explain steps and logic.
- Use `np.load(path/to/file, allow_pickle=True)` when loading `dogs_clades.npy`.
