# Assignment 6

## Description

In this assignment, you will explore eigenvalues and eigenvectors in the context of the Fibonacci sequence and apply Non-Negative Matrix Factorization (NMF) to the Dogs dataset to uncover underlying genetic ancestry patterns. The primary objectives are to derive the explicit formula for the Fibonacci sequence using matrix diagonalization and to apply NMF to high-dimensional data for dimensionality reduction and pattern recognition.

### Background on NMF

**Non-Negative Matrix Factorization (NMF)** is a dimensionality reduction technique that factors a non-negative matrix into two non-negative matrices. Unlike PCA, NMF results in additive, non-negative components, which can lead to more interpretable parts-based representations.

---

### Assignment (Tier 1)

#### Part 1: Fibonacci Sequence Using Eigenvalues and Eigenvectors

1. **Derive the Explicit Formula for the Fibonacci Sequence**:

   - Consider the 2x2 matrix:

      [ 1  1 ]

      [ 1  0 ]

   - This matrix can generate the Fibonacci sequence through matrix multiplication.
   - Your task is to find the eigenvalues and eigenvectors of matrix \( A \).
   - Use these to derive the explicit formula for the \( n \)-th term of the Fibonacci sequence.
   - The formula should not use recursion.
   - Write a function that calculates the \( n \)-th Fibonacci number using this explicit formula, which must be from the eigenvectors.

2. **Code Implementation**:

   - Implement the function in Python, taking as input two starting digits.
   - The code will be tested with different starting vectors of two digits.
   - Ensure that your code does not use recursion.

#### Part 2: Non-Negative Matrix Factorization on Dogs Dataset

1. **Load the Dogs dataset**:

   - Use the same Dogs dataset provided in previous assignments (`dogs_X.npy` and `dogs_clades.npy`).
   - This dataset consists of 784 SNPs for 1355 dogs. And there are 30 different clades.

2. **Apply Non-Negative Matrix Factorization (NMF)**:

   - Apply NMF to decompose the dataset into two non-negative matrices \( W \) and \( H \), where \( X \approx WH \).
   - Set the number of components \( n\_components = 5 \).
   - Use the `NMF` class from `sklearn.decomposition` and use `n_components=5`, `init='random'`, `random_state=42`.

3. **Create an admixture plot**:

   - Normalize the W matrix to represent proportions, (so the components for each individual sum to 1.)
   - Determine the dominant (largest value) cluster for each dog.
   - Sort the dogs by dominant cluster and proportion (highest to lowest).
   - Create a stacked plot (use `stackplot` from `matplotlib`) where each bar represents a dog, and the y-axis represents the proportion to each cluster. 
   - Visualize the admixture plot
   - Save the plot with filenames `NMF_Dogs.png`.

<!--

4. **Identify Basenji and Wolf samples**:

   - Identify the cluster containing the Bisenji (an ancient dog breed from Africa).
   - How many of the Wolf samples are in the same cluster?

-->

**Execution Command**:

    ```bash
    python Firstname_Lastname_Tier1.py
    ```

---

### Extra Credit 1 (Tier 2)

*Complete this if all of the above is finished in less than 4 hours.*

**NMF on Protein-Protein Interaction Dataset**

1. **Load the Protein-Protein Interaction Dataset**:

   - A dataset `protein_links.tsv` is provided.
   - This dataset contains interaction strengths between different proteins.
   - Create a symmetric interaction matrix, assuming undirected interactions, setting to zero the scores for the missing interactions.


2. **Apply NMF to the Dataset**:

   - Apply NMF to decompose the interaction matrix into \( W \) and \( H \) with 10 components.
   - Use scikit-learn's implementation of NMF and use `n_components=10`, `init='random'`, `random_state=42`.

3. **Filter proteins in selected cluster**:
   - Assign each protein to the cluster with the highest value in W – W has shape (n_proteins, n_clusters).
   - Select the cluster with the lowest protein count.

4. **Visualization**:

   - Perform hierarchical clustering within the selected cluster to group similar proteins together. Use `linkage` from `scipy.cluster.hierarchy`.
   - Create a heatmap with dendrogram of the interaction strengths for the proteins in the selected cluster. Use `clustermap` from `seaborn`.
   - Save the heatmap as `Protein_Cluster_Heatmap.png`.

**Execution Command**:

    ```bash
    python Firstname_Lastname_Tier2.py
    ```

---

<!--

### Extra Credit 2 (Tier 3)

*Do this if you finish all of the above in less than 4 hours.*

**Biological Interpretation**

   - Investigate the biological processes involved in the selected cluster.
   - Look up in the provided `protein_info.tsv` the names and annotations of the proteins in the cluster.
   - Use online databases or literature to provide insights into the functions of these proteins.
   - Discuss its potential involvement in specific biological pathways or diseases.
   - Reference at least two scientific articles that support your analysis.
   - Write a brief report summarizing your findings and save it as `Protein_Cluster_Report.txt`.

---
-->

## Assignment Requirements

### Data

- **Dogs Dataset**:

  - Provided with previous assignments as `dogs_X.npy` and `dogs_clades.npy`.

- **Protein-Protein Interaction Dataset** (for Tier 2 and Tier 3):

  - Provided as `protein_links.tsv`.
  - Contains interaction strengths between proteins.
  <!-- - Supplementary `protein_info.tsv` with protein names and detailed annotations. -->

### Execution

- The scripts will be executed using the specified commands.
- Ensure that all code is well-documented and generates the expected outputs.

---

## Submission Requirements

**Note**: The scripts should output all expected files and plots mentioned before. You also have to independently submit the plots and text files that the scripts generate.

1. You **must** submit your solution for Tier 1 with the following naming scheme:

        Firstname_Lastname_Tier1.py
        NMF_Dogs.png
        

2. If you attempted Tier 2 for extra credit, **you must submit all of the above, plus the following**:

        Firstname_Lastname_Tier2.py
        Protein_Cluster_Heatmap.png
<!--
3. If you attempted Tier 3, **you must submit all of the above, plus the following**:

        Protein_Cluster_Report.txt
-->
---

## Hints

1. **Libraries**:

   - Use `numpy`, `scipy`, `pandas`, `scikit-learn`, and `matplotlib` for computations and visualizations.
   - For eigenvalues and eigenvectors, use functions from `numpy.linalg` or `scipy.linalg`.
   - For NMF, use `NMF` from `sklearn.decomposition`.

3. **Fibonacci Sequence**:

   - The explicit formula involves the golden ratio \( \phi = \frac{1 + \sqrt{5}}{2} \) and its conjugate.
   - Be careful with floating-point precision when dealing with irrational numbers.

4. **Non-Negative Matrix Factorization**:

   - Ensure that the input data is non-negative.
   - The components in NMF can be interpreted as parts of the original data.
   - When visualizing components, consider normalizing them for better contrast.

5. **Protein-Protein Interaction Analysis**:

   - Use the fixed `random_state` to ensure reproducibility.
   - Clusters can be interpreted using biological databases like Gene Ontology, KEGG, or Reactome.
   - For visualization, use the `seaborn` library to create heatmaps.

---

## Notes

- Set random seeds (e.g., `random_state=42`) for reproducibility where applicable.
- Comment your code to explain the logic and steps clearly.
- Ensure all file paths are correctly specified and files are properly saved.
- When dealing with large datasets, ensure that your code is optimized for efficiency.
