from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# File paths
# -----------------------------
input_file = Path(
    "Results/GSE130727_Multimodel/GSE130727_pydeseq2_candidate_results.tsv"
)

output_figure = Path(
    "Figures/Multimodel/Figure7_GSE130727_pydeseq2_heatmap.png"
)

output_matrix = Path(
    "Results/GSE130727_Multimodel/GSE130727_pydeseq2_heatmap_matrix.tsv"
)

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(input_file, sep="\t")

print("Columns found:")
print(df.columns.tolist())

# -----------------------------
# Keep only needed columns
# -----------------------------
required_columns = ["Comparison", "Gene", "log2FoldChange", "padj"]
missing = [col for col in required_columns if col not in df.columns]

if missing:
    raise ValueError(f"Missing required columns: {missing}")

plot_df = df[required_columns].copy()

# -----------------------------
# Desired order
# -----------------------------
comparison_order = [
    "HAEC_IR",
    "HUVEC_IR",
    "IMR90_IR",
    "IMR90_Replicative",
    "WI38_Doxorubicin",
    "WI38_HRAS",
    "WI38_IR",
    "WI38_Replicative",
]

gene_order = ["FOS", "IL1A", "CXCL8", "CDKN2B", "ICAM1"]

# -----------------------------
# Build heatmap matrix
# -----------------------------
heatmap_matrix = plot_df.pivot(
    index="Comparison",
    columns="Gene",
    values="log2FoldChange"
)

heatmap_matrix = heatmap_matrix.reindex(index=comparison_order, columns=gene_order)

# significance matrix
sig_matrix = plot_df.pivot(
    index="Comparison",
    columns="Gene",
    values="padj"
)
sig_matrix = sig_matrix.reindex(index=comparison_order, columns=gene_order)

# Save matrix
heatmap_matrix.to_csv(output_matrix, sep="\t")

print("\nHeatmap matrix:")
print(heatmap_matrix)

# -----------------------------
# Prepare plot
# -----------------------------
data = heatmap_matrix.values.astype(float)

# symmetric color scale around zero
max_abs = np.nanmax(np.abs(data))
vmin = -max_abs
vmax = max_abs

fig, ax = plt.subplots(figsize=(10, 7))

im = ax.imshow(data, cmap="coolwarm", vmin=vmin, vmax=vmax, aspect="auto")

# Axis labels
ax.set_xticks(np.arange(len(gene_order)))
ax.set_yticks(np.arange(len(comparison_order)))
ax.set_xticklabels(gene_order, fontsize=12, fontweight="bold")
ax.set_yticklabels(comparison_order, fontsize=11)

ax.set_xlabel("Candidate gene", fontsize=12)
ax.set_ylabel("Senescence model", fontsize=12)

ax.set_title(
    "Figure 7. GSE130727 multimodel formal validation heatmap\n"
    "PyDESeq2 log2FC across 8 senescence models for 5 candidate genes",
    fontsize=15,
    fontweight="bold"
)

# -----------------------------
# Cell annotations
# -----------------------------
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        value = data[i, j]

        if pd.isna(value):
            text_str = "NA"
        else:
            text_str = f"{value:.2f}"

        # add star if padj < 0.05
        padj_value = sig_matrix.iloc[i, j]
        if pd.notna(padj_value) and padj_value < 0.05:
            text_str = text_str + "  *"

        # choose text color for readability
        if pd.isna(value):
            text_color = "black"
        elif abs(value) > max_abs * 0.45:
            text_color = "white"
        else:
            text_color = "black"

        ax.text(
            j, i, text_str,
            ha="center", va="center",
            color=text_color,
            fontsize=10
        )

# colorbar
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("PyDESeq2 log2FC (Senescent vs Control)", fontsize=11)

# footnote
fig.text(
    0.5, 0.02,
    "* padj < 0.05",
    ha="center",
    fontsize=10,
    style="italic"
)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig(output_figure, dpi=300, bbox_inches="tight")
plt.show()

print("\nFigure saved as:")
print(output_figure)

print("\nMatrix saved as:")
print(output_matrix)