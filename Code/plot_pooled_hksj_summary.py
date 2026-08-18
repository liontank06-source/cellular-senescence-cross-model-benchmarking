from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

input_file = Path(
    "Results/MetaAnalysis/"
    "GSE130727_candidate_HKSJ_meta_analysis_summary.tsv"
)

figure_directory = Path(
    "Figures/MetaAnalysis"
)

png_output = figure_directory / (
    "Figure9_GSE130727_pooled_HKSJ_summary.png"
)

pdf_output = figure_directory / (
    "Figure9_GSE130727_pooled_HKSJ_summary.pdf"
)


# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

if not input_file.exists():
    raise FileNotFoundError(
        f"Could not find:\n{input_file.resolve()}"
    )

figure_directory.mkdir(
    parents=True,
    exist_ok=True
)

data = pd.read_csv(
    input_file,
    sep="\t"
)


# ---------------------------------------------------------
# Validate columns
# ---------------------------------------------------------

required_columns = [
    "Gene",
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
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
# Arrange genes
# ---------------------------------------------------------

gene_order = [
    "CDKN2B",
    "ICAM1",
    "CXCL8",
    "IL1A",
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
# Extract values
# ---------------------------------------------------------

effects = data[
    "Pooled_log2FoldChange"
].to_numpy(dtype=float)

lower_limits = data[
    "HKSJ_CI95_lower"
].to_numpy(dtype=float)

upper_limits = data[
    "HKSJ_CI95_upper"
].to_numpy(dtype=float)

pvalues = data[
    "HKSJ_pvalue"
].to_numpy(dtype=float)

i2_values = data[
    "I2_percent"
].to_numpy(dtype=float)

lower_errors = effects - lower_limits
upper_errors = upper_limits - effects

y_positions = np.arange(
    len(data),
    0,
    -1
)


# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

figure, axis = plt.subplots(
    figsize=(11, 6.5)
)

axis.errorbar(
    effects,
    y_positions,
    xerr=[
        lower_errors,
        upper_errors
    ],
    fmt="D",
    markersize=8,
    capsize=5,
    linewidth=2
)

axis.axvline(
    x=0,
    linestyle="--",
    linewidth=1.2
)

axis.set_yticks(
    y_positions
)

axis.set_yticklabels(
    data["Gene"].astype(str),
    fontsize=12,
    fontweight="bold"
)

axis.tick_params(
    axis="y",
    pad=120
)
# Move gene names farther left to create space for the effect estimates
axis.tick_params(
    axis="y",
    pad=120
)

axis.set_xlabel(
    "Pooled log2 fold change "
    "(Senescent versus Control)",
    fontsize=11
)

figure.suptitle(
    "Cross-Model Pooled Effects for Candidate Senescence Markers",
    fontsize=16,
    fontweight="bold",
    y=0.97
)

axis.grid(
    axis="x",
    linestyle=":",
    alpha=0.5
)


# ---------------------------------------------------------
# Add statistical annotations
# ---------------------------------------------------------

# Add p-value and I-squared in a separate right column
for index, y_position in enumerate(
    y_positions
):

    right_annotation = (
        f"{pvalues[index]:.3f}; "
        f"{i2_values[index]:.1f}%"
    )

    axis.text(
        1.025,
        y_position,
        right_annotation,
        transform=axis.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=9,
        clip_on=False
    )


# Header for the numerical column on the right
axis.text(
    1.025,
    1.025,
    "HKSJ p; I²",
    transform=axis.transAxes,
    ha="left",
    va="bottom",
    fontsize=9,
    fontweight="bold",
    clip_on=False
)


# Add pooled effect and HKSJ confidence interval on the left
for index, y_position in enumerate(
    y_positions
):

    left_annotation = (
        f"{effects[index]:.2f} "
        f"[{lower_limits[index]:.2f}, "
        f"{upper_limits[index]:.2f}]"
    )

    axis.text(
        -0.025,
        y_position,
        left_annotation,
        transform=axis.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=9
    )


# Header for the numerical column on the left
axis.text(
    -0.025,
    1.025,
    "Pooled log2FC [HKSJ 95% CI]",
    transform=axis.transAxes,
    ha="right",
    va="bottom",
    fontsize=9,
    fontweight="bold"
)


# ---------------------------------------------------------
# Adjust horizontal limits
# ---------------------------------------------------------

x_min = min(
    lower_limits
) - 0.5

x_max = max(
    upper_limits
) + 0.5


# ---------------------------------------------------------
# Figure notes
# ---------------------------------------------------------

figure.text(
    0.5,
    0.025,
    (
        "Pooled estimates use Paule–Mandel random effects "
        "with Hartung–Knapp–Sidik–Jonkman-adjusted "
        "95% confidence intervals."
    ),
    ha="center",
    fontsize=8,
    fontstyle="italic"
)

figure.subplots_adjust(
    left=0.34,
    right=0.80,
    bottom=0.16,
    top=0.82
)


# ---------------------------------------------------------
# Save
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


print(
    "\nPOOLED HKSJ SUMMARY FIGURE CREATED"
)

print(
    "----------------------------------"
)

print(png_output)
print(pdf_output)