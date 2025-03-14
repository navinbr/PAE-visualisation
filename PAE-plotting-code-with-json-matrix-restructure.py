import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Load JSON data for WT and mutant
wt_file = "fold_fbxo22_human_q8nez5_full_data_0.json"  # Replace with actual filename
mut_file = "fold_fbxo22_human_q8nez5_val222del_full_data_0.json"  # Replace with actual filename

with open(wt_file) as f:
    wt_data = json.load(f)
with open(mut_file) as f:
    mut_data = json.load(f)

# Extract PAE matrices from the JSON data
wt_pae_matrix = np.array(wt_data["pae"])
mut_pae_matrix = np.array(mut_data["pae"])

# Function to introduce a gap at residue 222 in the mutant PAE matrix
def insert_gap(matrix, gap_position):
    size = matrix.shape[0] + 1
    new_matrix = np.full((size, size), np.nan)  # Initialize with NaNs
    new_matrix[:gap_position, :gap_position] = matrix[:gap_position, :gap_position]
    new_matrix[gap_position + 1:, :gap_position] = matrix[gap_position:, :gap_position]
    new_matrix[:gap_position, gap_position + 1:] = matrix[:gap_position, gap_position:]
    new_matrix[gap_position + 1:, gap_position + 1:] = matrix[gap_position:, gap_position:]
    return new_matrix

# Apply gap insertion for the mutant PAE matrix
mut_pae_matrix_with_gap = insert_gap(mut_pae_matrix, 222)

# Compute the difference matrix (Mutant - WT)
wt_pae_trimmed = wt_pae_matrix[:mut_pae_matrix_with_gap.shape[0], :mut_pae_matrix_with_gap.shape[1]]
difference_matrix = mut_pae_matrix_with_gap - wt_pae_trimmed

# Function to plot PAE matrices with correct Y-axis labels
def plot_pae_fix_y_top_to_bottom(matrix, title, ax, total_residues=403):
    cmap = plt.cm.viridis if "Difference" not in title else plt.cm.coolwarm
    cmap.set_bad(color='white')  # Represent gaps with white
    im = ax.imshow(matrix, cmap=cmap, origin='upper')  # Top-left is (1,1), bottom-left is (1,403)
    ax.set_title(title)
    ax.set_xlabel("Scored Residue")
    ax.set_ylabel("Aligned Residue")
    ax.set_xticks(range(0, total_residues, 50))
    ax.set_xticklabels(range(1, total_residues + 1, 50))  # X-axis from 1 to 403
    ax.set_yticks(range(0, total_residues, 50))
    ax.set_yticklabels(range(1, total_residues + 1, 50))  # Y-axis from 1 to 403 (top to bottom)
    plt.colorbar(im, ax=ax, label="Expected Position Error (Å)" if "Difference" not in title else "Difference in Expected Position Error (Å)")

# Create single-page plots for all three matrices
fig, axes = plt.subplots(1, 3, figsize=(24, 8))

# WT PAE
plot_pae_fix_y_top_to_bottom(wt_pae_matrix, "WT PAE", axes[0])

# Mutant PAE
plot_pae_fix_y_top_to_bottom(mut_pae_matrix_with_gap, "Mutant PAE (Val222del with Gap)", axes[1])

# Difference Plot
plot_pae_fix_y_top_to_bottom(difference_matrix, "Difference (Mutant - WT)", axes[2])

# Save as PDF and PNG
pdf_path = "PAE_Plots_One_Page.pdf"
png_path = "PAE_Plots_One_Page.png"

fig.savefig(png_path, dpi=300, bbox_inches='tight')
with PdfPages(pdf_path) as pdf:
    pdf.savefig(fig)
plt.close(fig)

print(f"Saved PDF: {pdf_path}")
print(f"Saved PNG: {png_path}")