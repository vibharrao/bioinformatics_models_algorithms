import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations
import matplotlib.patches as mpatches

# Load the MDS coordinates and atom information
molecule_coords_df = pd.read_csv('molecule_coordinates.csv')  # Contains X, Y, Z coordinates and Atom Index, Element
molecule_distances = pd.read_csv('molecule_distances.tsv', sep='\t').drop(columns=['Atom Index', 'Element'])

# CPK Coloring Convention with Hexadecimal Web Colors
cpk_colors_hex = {
    1: '#FFFFFF',     # Hydrogen
    6: '#909090',     # Carbon
    8: '#FF0D0D',     # Oxygen   
}

# Default color for elements not specified
default_color = '#FF1493'  # pink (Hex: FF1493)


# Determine atom colors based on atomic number in 'Element' column
colors = [cpk_colors_hex.get(atomic_number, default_color) for atomic_number in molecule_coords_df['Element']]

# Set atom sizes using atomic numbers with a scaling factor for visibility
scaling_factor = 50  # Adjust to control the overall size of atoms in the plot
sizes = [scaling_factor * (atomic_number ** 0.5) for atomic_number in molecule_coords_df['Element']]

# Create the 3D scatter plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot atoms with colors and sizes
sc = ax.scatter(
    molecule_coords_df['X'], 
    molecule_coords_df['Y'], 
    molecule_coords_df['Z'], 
    c=colors, 
    s=sizes, 
    edgecolor='k', 
    alpha=0.8
)

# Add a custom legend based on atomic numbers and CPK colors
unique_atomic_numbers = molecule_coords_df['Element'].unique()
legend_handles = []

for atomic_number in unique_atomic_numbers:
    color = cpk_colors_hex.get(atomic_number, default_color)
    # Create a Circle patch with a black edge for each unique atomic number
    legend_handles.append(
        mpatches.Circle((0, 0), radius=0.5, facecolor=color, edgecolor='black', label=f"Atom {atomic_number}")
    )

# Display the legend with custom handles
ax.legend(handles=legend_handles, title="Atoms", loc="upper right", markerscale=1.5, fontsize=10)

# Draw bonds if distance < 1.6
for i, j in combinations(range(len(molecule_coords_df)), 2):
    distance = molecule_distances.iloc[i, j]
    if distance < 1.6:
        x_values = [molecule_coords_df.iloc[i]['X'], molecule_coords_df.iloc[j]['X']]
        y_values = [molecule_coords_df.iloc[i]['Y'], molecule_coords_df.iloc[j]['Y']]
        z_values = [molecule_coords_df.iloc[i]['Z'], molecule_coords_df.iloc[j]['Z']]
        ax.plot(x_values, y_values, z_values, color='gray', linewidth=0.5)

# Set axis labels and plot title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title("3D Visualization of Molecule with CPK Coloring and Bonds")

ax.view_init(elev=30, azim=45)
# Save the plot
plt.savefig('molecule_3D_plot.png', dpi=300)  # Save with 300 dpi for high quality

# Display the interactive plot
plt.show()
