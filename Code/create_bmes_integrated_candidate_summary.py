from pathlib import Path

import pandas as pd


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

consensus_file = Path(
    "Results/Integrated/"
    "Final_candidate_consensus.tsv"
)

reliability_file = Path(
    "Results/MetaAnalysis/"
    "GSE130727_candidate_reliability_metrics.tsv"
)

output_directory = Path(
    "Results/Integrated"
)

output_file = output_directory / (
    "BMES_integrated_candidate_summary.tsv"
)


# ---------------------------------------------------------
# Confirm that input files exist
# ---------------------------------------------------------

for file_path in [
    consensus_file,
    reliability_file
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
# Load tables
# ---------------------------------------------------------

consensus_data = pd.read_csv(
    consensus_file,
    sep="\t"
)

reliability_data = pd.read_csv(
    reliability_file,
    sep="\t"
)


# ---------------------------------------------------------
# Required columns
# ---------------------------------------------------------

consensus_required_columns = [
    "Gene",
    "BJ_pathway_count",
    "BJ_log2FC_Senescent_vs_Young",
    "BJ_padj",
    "BJ_discovery_support",
    "GSE175533_validation_class",
    "GSE130727_significant_increases",
    "GSE130727_significant_decreases",
    "GSE130727_not_significant",
    "GSE130727_multimodel_class",
    "Final_consensus"
]

reliability_required_columns = [
    "Gene",
    "Total_models",
    "Significant_models",
    "Significant_increases",
    "Significant_decreases",
    "Not_significant",
    "Evidence_coverage_percent",
    "Positive_direction_percent",
    "Signed_direction_balance_percent",
    "Absolute_directional_consistency_percent",
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent",
    "All_LOO_same_direction_as_full",
    "Minimum_LOO_pooled_log2FoldChange",
    "Maximum_LOO_pooled_log2FoldChange",
    "Most_influential_omitted_model"
]


def check_columns(
    dataframe,
    required_columns,
    table_name
):
    """Confirm that all expected columns exist."""

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
    consensus_data,
    consensus_required_columns,
    "Final candidate consensus"
)

check_columns(
    reliability_data,
    reliability_required_columns,
    "Reliability metrics"
)


# ---------------------------------------------------------
# Keep relevant columns
# ---------------------------------------------------------

consensus_selected = consensus_data[
    consensus_required_columns
].copy()

reliability_selected = reliability_data[
    reliability_required_columns
].copy()


# ---------------------------------------------------------
# Confirm that old and new GSE130727 counts agree
# ---------------------------------------------------------

count_check = consensus_selected.merge(
    reliability_selected[
        [
            "Gene",
            "Significant_increases",
            "Significant_decreases",
            "Not_significant"
        ]
    ],
    on="Gene",
    how="inner",
    validate="one_to_one"
)

increase_match = (
    count_check[
        "GSE130727_significant_increases"
    ]
    == count_check[
        "Significant_increases"
    ]
)

decrease_match = (
    count_check[
        "GSE130727_significant_decreases"
    ]
    == count_check[
        "Significant_decreases"
    ]
)

not_significant_match = (
    count_check[
        "GSE130727_not_significant"
    ]
    == count_check[
        "Not_significant"
    ]
)

if not (
    increase_match.all()
    and decrease_match.all()
    and not_significant_match.all()
):
    raise ValueError(
        "The original consensus counts do not match "
        "the new reliability metrics."
    )


# ---------------------------------------------------------
# Merge original and upgraded evidence
# ---------------------------------------------------------

integrated_table = consensus_selected.merge(
    reliability_selected,
    on="Gene",
    how="inner",
    validate="one_to_one"
)


# ---------------------------------------------------------
# Convert leave-one-out stability to Boolean
# ---------------------------------------------------------

integrated_table[
    "All_LOO_same_direction_as_full"
] = (
    integrated_table[
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

if integrated_table[
    "All_LOO_same_direction_as_full"
].isna().any():
    raise ValueError(
        "Could not interpret leave-one-model-out stability."
    )


# ---------------------------------------------------------
# Transparent BMES reliability classification
# ---------------------------------------------------------

def classify_bmes_reliability(row):
    """
    Classify candidates using transparent cross-model metrics.
    This does not replace the original project consensus.
    """

    stable = row[
        "All_LOO_same_direction_as_full"
    ]

    positive_percent = row[
        "Positive_direction_percent"
    ]

    coverage = row[
        "Evidence_coverage_percent"
    ]

    pooled_effect = row[
        "Pooled_log2FoldChange"
    ]

    if (
        stable
        and positive_percent >= 75
        and pooled_effect > 0
    ):
        return (
            "Highest directional robustness"
        )

    if (
        stable
        and pooled_effect > 0
    ):
        return (
            "Direction-stable but context-dependent"
        )

    if coverage >= 80:
        return (
            "High-coverage but direction-sensitive"
        )

    if positive_percent == 50:
        return (
            "Balanced mixed-direction evidence"
        )

    return (
        "Lower-coverage direction-sensitive evidence"
    )


integrated_table[
    "BMES_reliability_class"
] = integrated_table.apply(
    classify_bmes_reliability,
    axis=1
)


# ---------------------------------------------------------
# Create concise display fields
# ---------------------------------------------------------

integrated_table[
    "GSE130727_result_pattern"
] = integrated_table.apply(
    lambda row: (
        f"{int(row['Significant_increases'])} increases, "
        f"{int(row['Significant_decreases'])} decreases, "
        f"{int(row['Not_significant'])} not significant"
    ),
    axis=1
)

integrated_table[
    "HKSJ_pooled_summary"
] = integrated_table.apply(
    lambda row: (
        f"{row['Pooled_log2FoldChange']:.2f} "
        f"[{row['HKSJ_CI95_lower']:.2f}, "
        f"{row['HKSJ_CI95_upper']:.2f}]; "
        f"p={row['HKSJ_pvalue']:.3f}"
    ),
    axis=1
)

integrated_table[
    "Leave_one_out_stability"
] = integrated_table[
    "All_LOO_same_direction_as_full"
].map(
    {
        True: "Direction-stable",
        False: "Direction-sensitive"
    }
)


# ---------------------------------------------------------
# Arrange final output
# ---------------------------------------------------------

final_columns = [
    "Gene",
    "BJ_pathway_count",
    "BJ_log2FC_Senescent_vs_Young",
    "BJ_padj",
    "BJ_discovery_support",
    "GSE175533_validation_class",
    "GSE130727_result_pattern",
    "Evidence_coverage_percent",
    "Positive_direction_percent",
    "Leave_one_out_stability",
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent",
    "Most_influential_omitted_model",
    "GSE130727_multimodel_class",
    "Final_consensus",
    "BMES_reliability_class"
]

final_table = integrated_table[
    final_columns
].copy()


# ---------------------------------------------------------
# Apply preferred gene order
# ---------------------------------------------------------

gene_order = [
    "CDKN2B",
    "IL1A",
    "CXCL8",
    "ICAM1",
    "FOS"
]

final_table["Gene"] = pd.Categorical(
    final_table["Gene"],
    categories=gene_order,
    ordered=True
)

final_table = (
    final_table
    .sort_values("Gene")
    .reset_index(drop=True)
)


# ---------------------------------------------------------
# Final validation and save
# ---------------------------------------------------------

if len(final_table) != 5:
    raise ValueError(
        f"Expected 5 genes, but found {len(final_table)}."
    )

if final_table.isna().any().any():
    missing_counts = (
        final_table
        .isna()
        .sum()
    )

    raise ValueError(
        "Missing values were found:\n"
        + missing_counts[
            missing_counts > 0
        ].to_string()
    )

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
    "GSE175533_validation_class",
    "GSE130727_result_pattern",
    "Evidence_coverage_percent",
    "Positive_direction_percent",
    "Leave_one_out_stability",
    "HKSJ_pooled_summary",
    "I2_percent",
    "BMES_reliability_class"
]

# HKSJ_pooled_summary was created before selecting columns,
# so attach it temporarily for the terminal display.
display_table = final_table.merge(
    integrated_table[
        [
            "Gene",
            "HKSJ_pooled_summary"
        ]
    ],
    on="Gene",
    how="left",
    validate="one_to_one"
)

print(
    "\nBMES INTEGRATED CANDIDATE SUMMARY CREATED"
)

print(
    "-----------------------------------------"
)

print(
    display_table[
        display_columns
    ]
    .round(2)
    .to_string(index=False)
)

print(
    "\nCOUNT CONSISTENCY CHECK: PASSED"
)

print(
    "\nFILE SAVED"
)

print(
    "----------"
)

print(output_file)