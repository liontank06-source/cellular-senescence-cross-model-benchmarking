from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# Input and output files
input_file = Path(
    "Results/GSE175533_WI38/"
    "GSE175533_validation_summary.tsv"
)

output_file = Path(
    "Figures/Validation/"
    "Figure5_GSE175533_integrated_validation_summary.png"
)

if not input_file.exists():
    raise FileNotFoundError(
        f"The validation summary was not found: {input_file}"
    )

# Read the integrated validation results
data = pd.read_csv(
    input_file,
    sep="\t"
)

required_columns = [
    "gene",
    "Endpoint_log2FC",
    "Endpoint_adjusted_p_value",
    "RS_time_log2FC",
    "RS_time_padj",
    "Validation_class"
]

missing_columns = [
    column
    for column in required_columns
    if column not in data.columns
]

if missing_columns:
    raise ValueError(
        f"Missing required columns: {missing_columns}"
    )


# Preferred visual order
gene_order = [
    "IL1A",
    "CDKN2B",
    "CXCL8",
    "ICAM1",
    "FOS"
]

data["gene"] = pd.Categorical(
    data["gene"],
    categories=gene_order,
    ordered=True
)

data = data.sort_values("gene").reset_index(drop=True)


def format_p_value(value):
    """Format adjusted p-values for the final figure."""

    if pd.isna(value):
        return "NA"

    if value < 0.001:
        return "<0.001"

    if value < 0.01:
        return f"{value:.4f}"

    return f"{value:.3f}"


def format_validation_class(value):
    """Wrap the longest validation label."""

    if value == "Not validated as increasing marker":
        return "Not validated as\nincreasing marker"

    return value


# Build the rows displayed in the figure
table_rows = []

for _, row in data.iterrows():

    table_rows.append(
        [
            row["gene"],
            f'{row["Endpoint_log2FC"]:.3f}',
            format_p_value(
                row["Endpoint_adjusted_p_value"]
            ),
            f'{row["RS_time_log2FC"]:.3f}',
            format_p_value(
                row["RS_time_padj"]
            ),
            format_validation_class(
                row["Validation_class"]
            )
        ]
    )


column_headers = [
    "Gene",
    "Endpoint log2FC\n(PDL53 vs PDL20)",
    "Endpoint\nadjusted p-value",
    "RS time-course\nlog2FC",
    "RS time-course\nadjusted p-value",
    "Validation class"
]


# Colors for validation classes
class_colors = {
    "Strong validation": "#D9EAD3",
    "Partial validation": "#FFF2CC",
    "Not validated as increasing marker": "#F4CCCC"
}

class_text_colors = {
    "Strong validation": "#1B5E20",
    "Partial validation": "#7F6000",
    "Not validated as increasing marker": "#990000"
}


# Create the figure
fig, ax = plt.subplots(
    figsize=(14, 8.5)
)

ax.axis("off")

fig.suptitle(
    "Figure 5. Integrated Validation Summary in GSE175533 (WI-38)",
    fontsize=19,
    fontweight="bold",
    y=0.97
)

ax.text(
    0.5,
    0.91,
    (
        "Endpoint PDL20 vs PDL53 analysis integrated with "
        "RS_deseq2 time-course results"
    ),
    ha="center",
    va="center",
    fontsize=12,
    style="italic",
    transform=ax.transAxes
)


# Create the table
table = ax.table(
    cellText=table_rows,
    colLabels=column_headers,
    cellLoc="center",
    colLoc="center",
    colWidths=[
        0.11,
        0.21,
        0.16,
        0.16,
        0.17,
        0.22
    ],
    bbox=[
        0.02,
        0.20,
        0.96,
        0.65
    ]
)

table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1, 1.7)


# Format the header
for column_index in range(len(column_headers)):

    header_cell = table[0, column_index]

    header_cell.set_facecolor("#EDEDED")
    header_cell.set_edgecolor("black")
    header_cell.set_linewidth(1.1)

    header_cell.get_text().set_fontweight("bold")
    header_cell.get_text().set_fontsize(10.5)


# Format all data cells
for row_index in range(1, len(table_rows) + 1):

    # Gene names
    gene_cell = table[row_index, 0]
    gene_cell.get_text().set_fontweight("bold")
    gene_cell.get_text().set_fontsize(12)

    # Add borders to every cell
    for column_index in range(len(column_headers)):

        cell = table[row_index, column_index]

        cell.set_edgecolor("black")
        cell.set_linewidth(0.9)

    # Color the validation-class cell
    validation_class = data.loc[
        row_index - 1,
        "Validation_class"
    ]

    validation_cell = table[row_index, 5]

    validation_cell.set_facecolor(
        class_colors[validation_class]
    )

    validation_cell.get_text().set_color(
        class_text_colors[validation_class]
    )

    validation_cell.get_text().set_fontweight("bold")
    validation_cell.get_text().set_wrap(True)


# Explanatory legend
ax.text(
    0.03,
    0.135,
    "Strong",
    color=class_text_colors["Strong validation"],
    fontsize=11,
    fontweight="bold",
    transform=ax.transAxes
)

ax.text(
    0.115,
    0.135,
    "= endpoint and time-course support a senescence-associated increase",
    fontsize=10.5,
    transform=ax.transAxes
)

ax.text(
    0.03,
    0.090,
    "Partial",
    color=class_text_colors["Partial validation"],
    fontsize=11,
    fontweight="bold",
    transform=ax.transAxes
)

ax.text(
    0.115,
    0.090,
    "= endpoint support, but time-course support is limited",
    fontsize=10.5,
    transform=ax.transAxes
)

ax.text(
    0.03,
    0.045,
    "Not validated",
    color=class_text_colors[
        "Not validated as increasing marker"
    ],
    fontsize=11,
    fontweight="bold",
    transform=ax.transAxes
)

ax.text(
    0.145,
    0.045,
    "= does not support an increasing senescence marker in this dataset",
    fontsize=10.5,
    transform=ax.transAxes
)


# Save without blocking the terminal
plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("\nFigure saved as:")
print(output_file)