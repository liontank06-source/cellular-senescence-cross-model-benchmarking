from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

input_file = Path(
    "Results/MetaAnalysis/"
    "GSE130727_candidate_reliability_metrics.tsv"
)

figure_directory = Path(
    "Figures/MetaAnalysis"
)

png_output = figure_directory / (
    "Figure10_GSE130727_candidate_reliability_landscape.png"
)

pdf_output = figure_directory / (
    "Figure10_GSE130727_candidate_reliability_landscape.pdf"
)


# ---------------------------------------------------------
# Confirm input and output locations
# ---------------------------------------------------------

if not input_file.exists():
    raise FileNotFoundError(
        f"Could not find:\n{input_file.resolve()}"
    )

figure_directory.mkdir(
    parents=True,
    exist_ok=True
)


# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

data = pd.read_csv(
    input_file,
    sep="\t"
)


# ---------------------------------------------------------
# Validate required columns
# ---------------------------------------------------------

required_columns = [
    "Gene",
    "Evidence_coverage_percent",
    "Positive_direction_percent",
    "All_LOO_same_direction_as_full",
    "Pooled_log2FoldChange",
    "HKSJ_pvalue",
    "I2_percent"
]

missing_columns = [
    column
    for column in required_columns
    if column not in data.columns
]

if missing_columns:
    raise ValueError(
        "Missing required columns: "
        + ", ".join(missing_columns)
    )


# ---------------------------------------------------------
# Convert numeric columns
# ---------------------------------------------------------

numeric_columns = [
    "Evidence_coverage_percent",
    "Positive_direction_percent",
    "Pooled_log2FoldChange",
    "HKSJ_pvalue",
    "I2_percent"
]

for column in numeric_columns:

    data[column] = pd.to_numeric(
        data[column],
        errors="coerce"
    )

if data[numeric_columns].isna().any().any():
    raise ValueError(
        "Missing or non-numeric reliability metrics were found."
    )


# ---------------------------------------------------------
# Convert leave-one-out stability to Boolean
# ---------------------------------------------------------

data["All_LOO_same_direction_as_full"] = (
    data["All_LOO_same_direction_as_full"]
    .astype(str)
    .str.strip()
    .str.lower()
    .map(
        {
            "true": True,
            "false": False
        }
    )
)

if data[
    "All_LOO_same_direction_as_full"
].isna().any():
    raise ValueError(
        "Could not interpret the leave-one-out "
        "direction-stability column."
    )


# ---------------------------------------------------------
# Arrange genes
# ---------------------------------------------------------

gene_order = [
    "CDKN2B",
    "IL1A",
    "CXCL8",
    "ICAM1",
    "FOS"
]

data["Gene"] = pd.Categorical(
    data["Gene"],
    categories=gene_order,
    ordered=True
)

data = (
    data
    .sort_values("Gene")
    .reset_index(drop=True)
)


# ---------------------------------------------------------
# Separate directionally stable and unstable candidates
# ---------------------------------------------------------

stable_data = data[
    data["All_LOO_same_direction_as_full"]
].copy()

unstable_data = data[
    ~data["All_LOO_same_direction_as_full"]
].copy()


# ---------------------------------------------------------
# Create figure
# ---------------------------------------------------------

figure, axis = plt.subplots(
    figsize=(11.5, 7.5)
)


# Candidates whose pooled direction changed
# in at least one leave-one-model-out analysis
axis.scatter(
    unstable_data["Evidence_coverage_percent"],
    unstable_data["Positive_direction_percent"],
    s=180,
    marker="o",
    label="Direction-sensitive\nin leave-one-model-out",
    zorder=3
)


# Candidates whose pooled direction stayed the same
# in every leave-one-model-out analysis
#
# Hollow diamonds allow overlapping points to remain visible.
axis.scatter(
    stable_data["Evidence_coverage_percent"],
    stable_data["Positive_direction_percent"],
    s=260,
    marker="D",
    facecolors="none",
    edgecolors="black",
    linewidths=2,
    label="Direction-stable\nin leave-one-model-out",
    zorder=4
)


# ---------------------------------------------------------
# Reference line
# ---------------------------------------------------------

axis.axhline(
    y=50,
    linestyle="--",
    linewidth=1.2
)


# ---------------------------------------------------------
# Gene-label positions
# ---------------------------------------------------------

label_offsets = {
    "CDKN2B": (12, 10),
    "IL1A": (-12, -18),
    "CXCL8": (14, 14),
    "ICAM1": (14, -24),
    "FOS": (-12, 14)
}


for _, row in data.iterrows():

    gene = str(row["Gene"])

    x_value = row[
        "Evidence_coverage_percent"
    ]

    y_value = row[
        "Positive_direction_percent"
    ]

    x_offset, y_offset = label_offsets[gene]

    if gene in [
        "IL1A",
        "FOS"
    ]:
        horizontal_alignment = "right"
    else:
        horizontal_alignment = "left"

    axis.annotate(
        gene,
        xy=(
            x_value,
            y_value
        ),
        xytext=(
            x_offset,
            y_offset
        ),
        textcoords="offset points",
        ha=horizontal_alignment,
        va="center",
        fontsize=11,
        fontweight="bold",
        annotation_clip=False
    )


# ---------------------------------------------------------
# Titles and labels
# ---------------------------------------------------------

figure.suptitle(
    "Senescence Marker Reliability Landscape",
    fontsize=18,
    fontweight="bold",
    y=0.97
)

figure.text(
    0.5,
    0.922,
    (
        "Cross-model evidence frequency and positive "
        "directionality across eight GSE130727 senescence models"
    ),
    ha="center",
    fontsize=11,
    fontstyle="italic"
)

axis.set_xlabel(
    "Evidence coverage "
    "(percentage of models with adjusted p < 0.05)",
    fontsize=11
)

axis.set_ylabel(
    "Positive direction among significant models (%)",
    fontsize=11
)


# ---------------------------------------------------------
# Axis limits and grid
# ---------------------------------------------------------

axis.set_xlim(
    55,
    95
)

axis.set_ylim(
    42,
    90
)

axis.grid(
    linestyle=":",
    alpha=0.5
)


# ---------------------------------------------------------
# Explanation inside the figure
# ---------------------------------------------------------

axis.text(
    56,
    44,
    (
        "50% line = equal numbers of significant "
        "increases and decreases"
    ),
    fontsize=8,
    fontstyle="italic",
    va="bottom"
)


# ---------------------------------------------------------
# Legend
# ---------------------------------------------------------

axis.legend(
    loc="upper left",
    bbox_to_anchor=(1.01, 1.00),
    borderaxespad=0.0,
    frameon=True,
    fontsize=8.5,
    labelspacing=0.8,
    handletextpad=0.8,
    borderpad=0.5
)


# ---------------------------------------------------------
# Figure note
# ---------------------------------------------------------

figure.text(
    0.5,
    0.02,
    (
        "Evidence coverage = significant models / 8. "
        "Positive direction = significant increases / "
        "all significant models. Hollow diamonds indicate "
        "that the pooled-effect direction remained unchanged "
        "after omitting any single model."
    ),
    ha="center",
    fontsize=8,
    fontstyle="italic"
)


figure.tight_layout(
    rect=[
        0.04,
        0.07,
        0.90,
        0.89
    ]
)


# ---------------------------------------------------------
# Save figure
# ---------------------------------------------------------

figure.savefig(
    png_output,
    dpi=300,
    bbox_inches="tight"
)

figure.savefig(
    pdf_output,
    bbox_inches="tight"
)

plt.close(
    figure
)


# ---------------------------------------------------------
# Terminal output
# ---------------------------------------------------------

print(
    "\nRELIABILITY LANDSCAPE CREATED"
)

print(
    "-----------------------------"
)

print(png_output)
print(pdf_output)