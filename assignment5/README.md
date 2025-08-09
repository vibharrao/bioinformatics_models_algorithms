# Assignment 5

## Description

This assignment focuses on applying dimensionality reduction techniques, specifically Principal Component Analysis (PCA) and Multidimensional Scaling (MDS), to visualize and reconstruct data from high-dimensional datasets. You will work with the MNIST dataset, the Dogs SNP dataset, and a molecular distance matrix. The primary objectives are to visualize high-dimensional data in two dimensions using PCA, reconstruct images using a reduced number of principal components, and reconstruct a three-dimensional molecular structure using MDS.

### Background on PCA and MDS

**Principal Component Analysis (PCA)** is a statistical technique used to reduce the dimensionality of datasets while preserving as much variance as possible. It achieves this by transforming the original variables into a new set of uncorrelated variables called principal components, ordered by the amount of variance they capture from the data.

**Multidimensional Scaling (MDS)** is a means of visualizing the level of similarity or dissimilarity of individual cases of a dataset. It positions each data point in N-dimensional space such that the between-point distances are preserved as well as possible. This is particularly useful when you have a distance matrix and wish to reconstruct the spatial configuration of the data points.

---

### Assignment (Tier 1)

#### Part 1: PCA on MNIST dataset

1. **Load the MNIST subset**:

   - Use the same MNIST subset provided in the previous assignment (`MNIST_X_subset.npy` and `MNIST_y_subset.npy`).
   - This subset consists of **6000 images** from the MNIST dataset, with **600 images for each digit** (0-9).

2. **Perform PCA to reduce dimensions to 2**:

   - Apply PCA to reduce the data from **784 dimensions** to **2 dimensions**.
   - Store the transformed data for visualization.

3. **Visualize the 2D PCA-transformed data**:

   - Create a scatter plot of the transformed data.
   - Color-code the points based on their actual digit labels.
   - Save the plot as `MNIST_PCA_2D.png`.

4. **Reconstruct an image using the first 2 principal components**:

   - Take the first example (index 0) from the dataset.
   - Project it onto the first 2 principal components to obtain its reduced representation.
   - Reconstruct the image from this reduced representation back to the original **784-dimensional** space.
   - Visualize and save the reconstructed image as `MNIST_reconstructed_2PC.png`.
   - Also, visualize and save the original image as `MNIST_original.png` for comparison.

5. **Reconstruct an image from a selected 2D point**:

   - Observing the 2D plot obtained in the previous space, manually choose the 2D coordinates of a location (not a data point) in the PCA-transformed space to visualize an image that resembles the digit 1.
   - Use the PCA model to reconstruct this point back to the original **784-dimensional** space.
   - Visualize and save the reconstructed image as `MNIST_reconstructed_1_from_coord.png`. You have just generated a new handwritten-like 1. This is a form of generative AI.

#### Part 2: PCA on Dogs SNP dataset

1. **Load the Dogs SNP dataset**:

   - Use the provided `dogs_X.npy` and `dogs_clades.npy` files.

2. **Perform PCA to reduce dimensions to 2**:

   - Apply PCA to reduce the SNP data from **784 dimensions** to **2 dimensions**.
   - Store the transformed data for visualization.

3. **Visualize the 2D PCA-transformed data**:

   - Create a scatter plot of the transformed data.
   - Color-code the points based on their clade labels from `dogs_clades.npy`.
   - Save the plot as `Dogs_PCA_2D.png`.

#### Part 3: MDS on molecular distance matrix

1. **Load the molecular distance matrix and atom types**:

   - A CSV file `molecule_distances.tsv` is provided, containing the pairwise distances between atoms in a molecule and the element label for each atom. This type of data (element and atomic distances) is produced by NMR on a chemical sample.

2. **Perform MDS to reconstruct 3D coordinates**:

   - Apply MDS to the distance matrix to obtain 3D coordinates of the atoms.
   - Use metric MDS with Euclidean distances.
   - Output the coordinates in a CSV file `molecule_coordinates.csv`.

**Execution command**:

    ```
    python Firstname_Lastname_Tier1.py
    ```

**Note**: The script should generate all the files mentioned above.

---


### Extra Credit 1 (Tier 2)

*Complete this if all of the above is finished in less than 4 hours.*

**3D Visualization of the molecular structure**

1. **Visualize the molecule in 3D**:

   - Using the coordinates obtained from MDS, create a 3D scatter plot of the atoms.
   - Color-code the atoms using the [CPK coloring convention](https://en.wikipedia.org/wiki/CPK_coloring).
   - Adjust the size of the markers (i.e., larger for heavier atoms).

2. **Connect atoms based on distance**:

   - For each pair of atoms, if their distance in the distance matrix is $< 1.6$, draw a line (edge) between them to represent a bond.
   - Use a 3D plotting library that supports line segments between points, such as `matplotlib`'s 3D toolkit or `Plotly`.

3. **Save the 3D plot**:

   - Save the visualization as `molecule_3D_plot.png`.

**Execution command**:

    ```
    python Firstname_Lastname_Tier2.py
    ```

---

### Extra Credit 2 (Tier 3)

*Do this if you finish all of the above in less than 4 hours.*

**Identify the molecule and bring a leaf**

1. **Identify the molecule**:

   - Using the reconstructed 3D structure and atom types, identify the molecule.
   - Research its chemical structure and properties.
   - Find out in which tree this molecule can be found.

   - In `Firstname_Lastname_Tier3.txt`, include the name of the molecule, and the type of tree on campus where it can be found.

2. **Bring a leaf**:

   - On Thursday, bring a leaf from this tree to class.

---

## Assignment Requirements

### Data

- **MNIST Dataset**:

  - Provided as `MNIST_X_subset.npy` and `MNIST_y_subset.npy` .
  - Contains 6000 images (600 per digit), each of size 28x28 pixels (784 when vectorized).

- **Dogs SNP Dataset**:

  - Provided as `dogs_X.npy` and `dogs_clades.npy`.
  - Contains SNP data for 1355 dog samples with 784 SNP features.

- **Molecular Data**:

  - `molecule_distances.csv`: Contains the pairwise distances between atoms and the element numbers.

### Execution

- The scripts should be executed using the specified commands.
- Outputs may vary due to different initializations but should be within the expected range for validation.
- Ensure that all code is well-documented and generates the expected outputs.

---

## Submission Requirements

**Note**: The scripts should output all expected plots and data files mentioned before. You also have to independently submit the plots that the scripts generate.

1. You **must** submit your solution for Tier 1 with the following naming scheme:

        Firstname_Lastname_Tier1.py
        MNIST_PCA_2D.png
        MNIST_original.png
        MNIST_reconstructed_2PC.png
        MNIST_reconstructed_1_from_coord.png
        Dogs_PCA_2D.png
        molecule_coordinates.csv

2. If you attempted Tier 2 for extra credit, **you must submit all of the above, plus the following**:

        Firstname_Lastname_Tier2.py
        molecule_3D_plot.png

3. If you attempted Tier 3, **you must submit all of the above, plus the following**:

        Firstname_Lastname_Tier3.txt
        A leaf (bring to class)

---

## Hints

1. **Libraries**:

   - Use `numpy`, `pandas`, `scikit-learn`, and `matplotlib` for computations and visualizations.
   - For PCA, use `PCA` from `sklearn.decomposition`.
   - For MDS, use `MDS` from `sklearn.manifold`.
   - For 3D plotting, consider using `matplotlib's` `mpl_toolkits.mplot3d`.

2. **Visualization**:

   - Use `matplotlib.pyplot.scatter` for 2D scatter plots.
   - When reconstructing images, reshape the vectors back to 28x28 pixel grids for visualization using `imshow`.
   - For 3D plots, label axes and include legends to indicate atom types.

3. **Image Reconstruction**:

   - When reconstructing images using a reduced number of principal components, use the `inverse_transform` method of the PCA object.
   - The reconstructed images will be an approximation and may appear fuzzy.

6. **Identifying the Molecule**:

   - Use the elements and 3D shape to deduce the molecule.
   - Consult chemical databases like PubChem or ChemSpider.
   - The molecule is commonly found in a specific type of tree (plenty around here).

---

## Notes

- Set random seeds (e.g., `random_state` parameter) for reproducibility where applicable.
- Comment your code to explain the logic and steps clearly.
- Ensure all file paths are correctly specified and files are properly saved.
- For the molecule visualization, adjust viewing angles or use interactive plots for better insight.
