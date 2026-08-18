from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

formal_results_file = Path(
    "Results/GSE130727_Multimodel/"
    "GSE130727_pydeseq2_candidate_results.tsv"
)

hksj_results_file = Path(
    "Results/MetaAnalysis/"
    "GSE130727_candidate_HKSJ_meta_analysis_summary.tsv"
)

leave_one_out_file = Path(
    "Results/MetaAnalysis/"
    "GSE130727_candidate_leave_one_model_out_summary.tsv"
)

output_directory = Path(
    "Results/MetaAnalysis"
)

output_file = output_directory / (
    "GSE130727_candidate_reliability_metrics.tsv"
)


# ---------------------------------------------------------
# Confirm input files exist
# ---------------------------------------------------------

for file_path in [
    formal_results_file,
    hksj_results_file,
    leave_one_out_file
]:
    if not file_path.exists():
        raise FileNotFoundError(
            f"Could not find:\n{file_path.resolve()}"
        )

output_directory.mkdir(
    parents=True,
    exist_ok=True
)


# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

formal_data = pd.read_csv(
    formal_results_file,
    sep="\t"
)

hksj_data = pd.read_csv(
    hksj_results_file,
    sep="\t"
)

leave_one_out_data = pd.read_csv(
    leave_one_out_file,
    sep="\t"
)


# ---------------------------------------------------------
# Validate required columns
# ---------------------------------------------------------

formal_required_columns = [
    "Gene",
    "Comparison",
    "log2FoldChange",
    "padj"
]

hksj_required_columns = [
    "Gene",
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent"
]

leave_one_out_required_columns = [
    "Gene",
    "All_LOO_same_direction_as_full",
    "Minimum_LOO_pooled_log2FoldChange",
    "Maximum_LOO_pooled_log2FoldChange",
    "Most_influential_omitted_model"
]


def check_columns(dataframe, required_columns, table_name):
    """Confirm that all required columns are present."""

    missing_columns = [
        column
        for column in required_columns
        if column not in dataframe.columns
    ]

    if missing_columns:
        raise ValueError(
            f"{table_name} is missing columns: "
            + ", ".join(missing_columns)
        )


check_columns(
    formal_data,
    formal_required_columns,
    "Formal PyDESeq2 results"
)

check_columns(
    hksj_data,
    hksj_required_columns,
    "HKSJ results"
)

check_columns(
    leave_one_out_data,
    leave_one_out_required_columns,
    "Leave-one-model-out results"
)


# ---------------------------------------------------------
# Convert numeric columns
# ---------------------------------------------------------

for column in [
    "log2FoldChange",
    "padj"
]:
    formal_data[column] = pd.to_numeric(
        formal_data[column],
        errors="coerce"
    )

for column in [
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent"
]:
    hksj_data[column] = pd.to_numeric(
        hksj_data[column],
        errors="coerce"
    )

for column in [
    "Minimum_LOO_pooled_log2FoldChange",
    "Maximum_LOO_pooled_log2FoldChange"
]:
    leave_one_out_data[column] = pd.to_numeric(
        leave_one_out_data[column],
        errors="coerce"
    )

if formal_data[
    ["log2FoldChange", "padj"]
].isna().any().any():
    raise ValueError(
        "Missing or non-numeric formal results were found."
    )


# ---------------------------------------------------------
# Convert leave-one-out direction column to Boolean
# ---------------------------------------------------------

leave_one_out_data[
    "All_LOO_same_direction_as_full"
] = (
    leave_one_out_data[
        "All_LOO_same_direction_as_full"
    ]
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

if leave_one_out_data[
    "All_LOO_same_direction_as_full"
].isna().any():
    raise ValueError(
        "Could not interpret the leave-one-out "
        "direction-stability column."
    )


# ---------------------------------------------------------
# Gene order
# ---------------------------------------------------------

gene_order = [
    "CDKN2B",
    "IL1A",
    "CXCL8",
    "ICAM1",
    "FOS"
]


# ---------------------------------------------------------
# Calculate transparent reliability metrics
# ---------------------------------------------------------

reliability_rows = []

for gene in gene_order:

    gene_data = formal_data[
        formal_data["Gene"] == gene
    ].copy()

    total_models = (
        gene_data["Comparison"].nunique()
    )

    if total_models != 8:
        raise ValueError(
            f"{gene} contains {total_models} models "
            "instead of 8."
        )

    significant_mask = (
        gene_data["padj"] < 0.05
    )

    increase_mask = (
        significant_mask
        & (
            gene_data["log2FoldChange"] > 0
        )
    )

    decrease_mask = (
        significant_mask
        & (
            gene_data["log2FoldChange"] < 0
        )
    )

    significant_models = int(
        significant_mask.sum()
    )

    significant_increases = int(
        increase_mask.sum()
    )

    significant_decreases = int(
        decrease_mask.sum()
    )

    not_significant = (
        total_models
        - significant_models
    )

    evidence_coverage_percent = (
        significant_models
        / total_models
        * 100
    )

    if significant_models > 0:

        positive_direction_percent = (
            significant_increases
            / significant_models
            * 100
        )

        signed_direction_balance_percent = (
            (
                significant_increases
                - significant_decreases
            )
            / significant_models
            * 100
        )

        absolute_directional_consistency_percent = (
            abs(
                significant_increases
                - significant_decreases
            )
            / significant_models
            * 100
        )

    else:

        positive_direction_percent = np.nan
        signed_direction_balance_percent = np.nan
        absolute_directional_consistency_percent = np.nan

    reliability_rows.append(
        {
            "Gene": gene,
            "Total_models": total_models,
            "Significant_models": significant_models,
            "Significant_increases":
                significant_increases,
            "Significant_decreases":
                significant_decreases,
            "Not_significant": not_significant,
            "Evidence_coverage_percent":
                evidence_coverage_percent,
            "Positive_direction_percent":
                positive_direction_percent,
            "Signed_direction_balance_percent":
                signed_direction_balance_percent,
            "Absolute_directional_consistency_percent":
                absolute_directional_consistency_percent
        }
    )


reliability_table = pd.DataFrame(
    reliability_rows
)


# ---------------------------------------------------------
# Keep only the complete eight-model HKSJ results
# ---------------------------------------------------------

if "Scenario" in hksj_data.columns:

    hksj_full = hksj_data[
        hksj_data["Scenario"] == "All_8_models"
    ].copy()

else:

    hksj_full = hksj_data.copy()


hksj_columns_to_keep = [
    "Gene",
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent"
]

hksj_full = hksj_full[
    hksj_columns_to_keep
].copy()


# ---------------------------------------------------------
# Keep selected leave-one-model-out metrics
# ---------------------------------------------------------

leave_one_out_columns_to_keep = [
    "Gene",
    "All_LOO_same_direction_as_full",
    "Minimum_LOO_pooled_log2FoldChange",
    "Maximum_LOO_pooled_log2FoldChange",
    "Most_influential_omitted_model"
]

leave_one_out_selected = (
    leave_one_out_data[
        leave_one_out_columns_to_keep
    ]
    .copy()
)


# ---------------------------------------------------------
# Merge all reliability evidence
# ---------------------------------------------------------

final_table = (
    reliability_table
    .merge(
        hksj_full,
        on="Gene",
        how="left",
        validate="one_to_one"
    )
    .merge(
        leave_one_out_selected,
        on="Gene",
        how="left",
        validate="one_to_one"
    )
)


# ---------------------------------------------------------
# Final validation
# ---------------------------------------------------------

if final_table.isna().any().any():
    missing_summary = (
        final_table
        .isna()
        .sum()
    )

    raise ValueError(
        "Missing values were found after merging:\n"
        + missing_summary[
            missing_summary > 0
        ].to_string()
    )


# ---------------------------------------------------------
# Save output
# ---------------------------------------------------------

final_table.to_csv(
    output_file,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------

display_columns = [
    "Gene",
    "Significant_increases",
    "Significant_decreases",
    "Not_significant",
    "Evidence_coverage_percent",
    "Positive_direction_percent",
    "All_LOO_same_direction_as_full",
    "Pooled_log2FoldChange",
    "HKSJ_pvalue",
    "I2_percent"
]

print(
    "\nCANDIDATE RELIABILITY METRICS CREATED"
)

print(
    "-------------------------------------"
)

print(
    final_table[
        display_columns
    ]
    .round(2)
    .to_string(index=False)
)

print(
    "\nFILE SAVED"
)

print(
    "----------"
)

print(output_file)