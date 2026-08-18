from pathlib import Path

import numpy as np
import pandas as pd


input_folder = Path(
    "Data/GSE130727_Multimodel"
)

output_folder = Path(
    "Results/GSE130727_Multimodel"
)

matrix_file = (
    input_folder
    / "GSE130727_gene_count_matrix.tsv.gz"
)

mapping_file = (
    input_folder
    / "GSE130727_candidate_ensembl_mapping.tsv"
)

design_file = (
    input_folder
    / "GSE130727_comparison_design.tsv"
)

qc_file = (
    input_folder
    / "GSE130727_count_matrix_QC.tsv"
)

long_output = (
    output_folder
    / "GSE130727_candidate_CPM_long.tsv"
)

summary_output = (
    output_folder
    / "GSE130727_candidate_model_summary.tsv"
)


# Preferred orders
gene_order = [
    "FOS",
    "IL1A",
    "CXCL8",
    "CDKN2B",
    "ICAM1"
]

comparison_order = [
    "HAEC_IR",
    "HUVEC_IR",
    "IMR90_IR",
    "IMR90_Replicative",
    "WI38_Doxorubicin",
    "WI38_HRAS",
    "WI38_IR",
    "WI38_Replicative"
]


# Read input files
matrix = pd.read_csv(
    matrix_file,
    sep="\t",
    compression="gzip"
)

mapping = pd.read_csv(
    mapping_file,
    sep="\t"
)

design = pd.read_csv(
    design_file,
    sep="\t"
)

qc = pd.read_csv(
    qc_file,
    sep="\t"
)


# Clean Ensembl IDs
matrix["Ensembl_ID"] = (
    matrix["gene"]
    .astype(str)
    .str.split(".")
    .str[0]
)

mapping["Ensembl_ID"] = (
    mapping["Ensembl_ID"]
    .astype(str)
    .str.split(".")
    .str[0]
)


# Keep the five candidate genes
candidate_ids = set(
    mapping["Ensembl_ID"]
)

candidate_matrix = matrix[
    matrix["Ensembl_ID"].isin(candidate_ids)
].copy()


if len(candidate_matrix) != 5:
    raise ValueError(
        "Expected 5 candidate genes, but found "
        f"{len(candidate_matrix)}."
    )


# Connect Ensembl IDs to gene symbols
id_to_symbol = dict(
    zip(
        mapping["Ensembl_ID"],
        mapping["Gene"]
    )
)

candidate_matrix["Gene"] = (
    candidate_matrix["Ensembl_ID"]
    .map(id_to_symbol)
)


# Identify sample columns
sample_columns = [
    column
    for column in matrix.columns
    if column.startswith("GSM")
]


if len(sample_columns) != 37:
    raise ValueError(
        "Expected 37 sample columns, but found "
        f"{len(sample_columns)}."
    )


# Convert the candidate matrix to long format
candidate_long = candidate_matrix.melt(
    id_vars=[
        "Ensembl_ID",
        "Gene"
    ],
    value_vars=sample_columns,
    var_name="GSM",
    value_name="Raw_count"
)


# Add total library counts for normalization
candidate_long = candidate_long.merge(
    qc[
        [
            "GSM",
            "Total_counts"
        ]
    ],
    on="GSM",
    how="left"
)


if candidate_long["Total_counts"].isna().any():
    raise ValueError(
        "Some samples are missing total-count information."
    )


# Counts per million
candidate_long["CPM"] = (
    candidate_long["Raw_count"]
    / candidate_long["Total_counts"]
    * 1_000_000
)


# Add each sample to its matched comparison
# IMR-90 PDL15 controls correctly appear in two comparisons
candidate_long = design[
    [
        "GSM",
        "Title",
        "cell_line",
        "Comparison",
        "Status"
    ]
].merge(
    candidate_long,
    on="GSM",
    how="left"
)


if candidate_long["Gene"].isna().any():
    raise ValueError(
        "Some comparison samples could not be connected "
        "to candidate-gene counts."
    )


# Save the sample-level normalized table
candidate_long.to_csv(
    long_output,
    sep="\t",
    index=False
)


# Create one summary row per gene and model
summary_rows = []

for (
    comparison,
    gene
), group in candidate_long.groupby(
    [
        "Comparison",
        "Gene"
    ]
):

    control_values = group.loc[
        group["Status"] == "Control",
        "CPM"
    ].astype(float)

    senescent_values = group.loc[
        group["Status"] == "Senescent",
        "CPM"
    ].astype(float)

    if control_values.empty or senescent_values.empty:
        raise ValueError(
            f"{comparison} / {gene} is missing a group."
        )

    control_mean = control_values.mean()
    senescent_mean = senescent_values.mean()

    descriptive_log2fc = np.log2(
        (senescent_mean + 1)
        /
        (control_mean + 1)
    )

    summary_rows.append(
        {
            "Comparison": comparison,
            "Gene": gene,
            "Control_n": len(control_values),
            "Senescent_n": len(senescent_values),
            "Control_mean_CPM": control_mean,
            "Senescent_mean_CPM": senescent_mean,
            "Control_SD_CPM": control_values.std(),
            "Senescent_SD_CPM": senescent_values.std(),
            "Descriptive_log2FC": descriptive_log2fc,
            "Direction": (
                "Higher in Senescent"
                if descriptive_log2fc > 0
                else "Higher in Control"
            )
        }
    )


summary = pd.DataFrame(
    summary_rows
)


# Preserve the preferred order
summary["Comparison"] = pd.Categorical(
    summary["Comparison"],
    categories=comparison_order,
    ordered=True
)

summary["Gene"] = pd.Categorical(
    summary["Gene"],
    categories=gene_order,
    ordered=True
)

summary = summary.sort_values(
    [
        "Comparison",
        "Gene"
    ]
)


summary.to_csv(
    summary_output,
    sep="\t",
    index=False
)


print("\nCandidate genes normalized successfully.")

print("\nCandidate genes:")
print(
    summary["Gene"]
    .astype(str)
    .unique()
)

print("\nNumber of comparisons:")
print(
    summary["Comparison"]
    .nunique()
)

print("\nSummary shape:")
print(summary.shape)

print("\nDirection counts:")
print(
    summary["Direction"]
    .value_counts()
)

print("\nFirst fifteen summary rows:")
print(
    summary.head(15)
    .round(4)
    .to_string(index=False)
)

print("\nSample-level CPM table saved as:")
print(long_output)

print("\nModel summary saved as:")
print(summary_output)