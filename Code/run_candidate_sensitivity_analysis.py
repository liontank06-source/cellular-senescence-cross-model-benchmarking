from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from statsmodels.stats.meta_analysis import combine_effects


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

input_file = Path(
    "Results/GSE130727_Multimodel/"
    "GSE130727_pydeseq2_candidate_results.tsv"
)

output_directory = Path(
    "Results/MetaAnalysis"
)

leave_one_out_output = output_directory / (
    "GSE130727_candidate_leave_one_model_out.tsv"
)

leave_one_out_summary_output = output_directory / (
    "GSE130727_candidate_leave_one_model_out_summary.tsv"
)

imr90_sensitivity_output = output_directory / (
    "GSE130727_candidate_IMR90_sensitivity.tsv"
)


# ---------------------------------------------------------
# Load and validate data
# ---------------------------------------------------------

if not input_file.exists():
    raise FileNotFoundError(
        f"Could not find:\n{input_file.resolve()}"
    )

output_directory.mkdir(
    parents=True,
    exist_ok=True
)

data = pd.read_csv(
    input_file,
    sep="\t"
)

required_columns = [
    "Gene",
    "Comparison",
    "log2FoldChange",
    "lfcSE",
    "padj"
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

meta_input = data[required_columns].copy()

for column in [
    "log2FoldChange",
    "lfcSE",
    "padj"
]:
    meta_input[column] = pd.to_numeric(
        meta_input[column],
        errors="coerce"
    )

if meta_input[
    ["log2FoldChange", "lfcSE", "padj"]
].isna().any().any():
    raise ValueError(
        "Missing or non-numeric values were found."
    )

if (meta_input["lfcSE"] <= 0).any():
    raise ValueError(
        "All lfcSE values must be greater than zero."
    )

meta_input["variance"] = (
    meta_input["lfcSE"] ** 2
)


# ---------------------------------------------------------
# Random-effects meta-analysis function
# ---------------------------------------------------------

def run_random_effects(subset):
    """
    Run a Paule-Mandel random-effects meta-analysis.
    """

    subset = (
        subset
        .sort_values("Comparison")
        .reset_index(drop=True)
    )

    effects = subset[
        "log2FoldChange"
    ].to_numpy(dtype=float)

    variances = subset[
        "variance"
    ].to_numpy(dtype=float)

    model_names = (
        subset["Comparison"]
        .astype(str)
        .tolist()
    )

    if len(effects) < 3:
        raise ValueError(
            "At least three models are required."
        )

    result = combine_effects(
        effect=effects,
        variance=variances,
        method_re="iterated",
        row_names=model_names
    )

    pooled_effect = float(
        result.mean_effect_re
    )

    pooled_se = float(
        result.sd_eff_w_re
    )

    ci_lower = (
        pooled_effect
        - 1.96 * pooled_se
    )

    ci_upper = (
        pooled_effect
        + 1.96 * pooled_se
    )

    if pooled_se > 0:
        z_statistic = (
            pooled_effect / pooled_se
        )

        pooled_pvalue = float(
            2 * norm.sf(
                abs(z_statistic)
            )
        )
    else:
        z_statistic = np.nan
        pooled_pvalue = np.nan

    q_statistic = float(
        result.q
    )

    q_df = len(effects) - 1

    q_pvalue = float(
        chi2.sf(
            q_statistic,
            q_df
        )
    )

    if q_statistic > 0:
        i2_percent = max(
            0.0,
            (
                q_statistic - q_df
            )
            / q_statistic
            * 100
        )
    else:
        i2_percent = 0.0

    tau_squared = max(
        0.0,
        float(result.tau2)
    )

    ci_excludes_zero = bool(
        ci_lower > 0
        or ci_upper < 0
    )

    if ci_lower > 0:
        overall_direction = (
            "Overall increase"
        )
    elif ci_upper < 0:
        overall_direction = (
            "Overall decrease"
        )
    else:
        overall_direction = (
            "No clear overall direction"
        )

    return {
        "Models_included": len(effects),
        "Pooled_log2FoldChange": pooled_effect,
        "Pooled_SE": pooled_se,
        "CI95_lower": ci_lower,
        "CI95_upper": ci_upper,
        "Pooled_pvalue": pooled_pvalue,
        "Q": q_statistic,
        "Q_df": q_df,
        "Q_pvalue": q_pvalue,
        "Tau_squared": tau_squared,
        "I2_percent": i2_percent,
        "CI_excludes_zero": ci_excludes_zero,
        "Overall_direction": overall_direction
    }


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
# Leave-one-model-out sensitivity analysis
# ---------------------------------------------------------

leave_one_out_rows = []
leave_one_out_summary_rows = []

for gene in gene_order:

    gene_data = (
        meta_input[
            meta_input["Gene"] == gene
        ]
        .copy()
    )

    if gene_data["Comparison"].nunique() != 8:
        raise ValueError(
            f"{gene} does not contain exactly 8 models."
        )

    full_result = run_random_effects(
        gene_data
    )

    full_effect = full_result[
        "Pooled_log2FoldChange"
    ]

    gene_leave_one_out_rows = []

    models = sorted(
        gene_data["Comparison"]
        .astype(str)
        .unique()
    )

    for omitted_model in models:

        reduced_data = gene_data[
            gene_data["Comparison"]
            != omitted_model
        ].copy()

        sensitivity_result = (
            run_random_effects(
                reduced_data
            )
        )

        pooled_effect = sensitivity_result[
            "Pooled_log2FoldChange"
        ]

        effect_change = (
            pooled_effect - full_effect
        )

        same_direction = bool(
            (pooled_effect >= 0)
            == (full_effect >= 0)
        )

        row = {
            "Gene": gene,
            "Omitted_model": omitted_model,
            "Full_pooled_log2FoldChange":
                full_effect,
            "Effect_change_from_full":
                effect_change,
            "Absolute_effect_change":
                abs(effect_change),
            "Same_direction_as_full":
                same_direction,
            **sensitivity_result
        }

        leave_one_out_rows.append(
            row
        )

        gene_leave_one_out_rows.append(
            row
        )

    gene_results = pd.DataFrame(
        gene_leave_one_out_rows
    )

    most_influential_index = (
        gene_results[
            "Absolute_effect_change"
        ]
        .idxmax()
    )

    most_influential_row = (
        gene_results.loc[
            most_influential_index
        ]
    )

    leave_one_out_summary_rows.append(
        {
            "Gene": gene,
            "Full_pooled_log2FoldChange":
                full_effect,
            "Minimum_LOO_pooled_log2FoldChange":
                gene_results[
                    "Pooled_log2FoldChange"
                ].min(),
            "Maximum_LOO_pooled_log2FoldChange":
                gene_results[
                    "Pooled_log2FoldChange"
                ].max(),
            "All_LOO_same_direction_as_full":
                gene_results[
                    "Same_direction_as_full"
                ].all(),
            "Any_LOO_CI_excludes_zero":
                gene_results[
                    "CI_excludes_zero"
                ].any(),
            "Most_influential_omitted_model":
                most_influential_row[
                    "Omitted_model"
                ],
            "Maximum_absolute_effect_change":
                most_influential_row[
                    "Absolute_effect_change"
                ],
            "Minimum_LOO_I2_percent":
                gene_results[
                    "I2_percent"
                ].min(),
            "Maximum_LOO_I2_percent":
                gene_results[
                    "I2_percent"
                ].max()
        }
    )


# ---------------------------------------------------------
# Shared-control IMR-90 sensitivity analysis
# ---------------------------------------------------------

imr90_rows = []

imr90_scenarios = {
    "All_8_models": [],
    "Exclude_IMR90_IR": [
        "IMR90_IR"
    ],
    "Exclude_IMR90_Replicative": [
        "IMR90_Replicative"
    ],
    "Exclude_both_IMR90_models": [
        "IMR90_IR",
        "IMR90_Replicative"
    ]
}

for gene in gene_order:

    gene_data = (
        meta_input[
            meta_input["Gene"] == gene
        ]
        .copy()
    )

    for scenario, excluded_models in (
        imr90_scenarios.items()
    ):

        scenario_data = gene_data[
            ~gene_data["Comparison"].isin(
                excluded_models
            )
        ].copy()

        result = run_random_effects(
            scenario_data
        )

        imr90_rows.append(
            {
                "Gene": gene,
                "Scenario": scenario,
                "Excluded_models": (
                    ", ".join(excluded_models)
                    if excluded_models
                    else "None"
                ),
                **result
            }
        )


# ---------------------------------------------------------
# Save files
# ---------------------------------------------------------

leave_one_out_table = pd.DataFrame(
    leave_one_out_rows
)

leave_one_out_summary_table = (
    pd.DataFrame(
        leave_one_out_summary_rows
    )
)

imr90_table = pd.DataFrame(
    imr90_rows
)

leave_one_out_table.to_csv(
    leave_one_out_output,
    sep="\t",
    index=False
)

leave_one_out_summary_table.to_csv(
    leave_one_out_summary_output,
    sep="\t",
    index=False
)

imr90_table.to_csv(
    imr90_sensitivity_output,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------

summary_display_columns = [
    "Gene",
    "Full_pooled_log2FoldChange",
    "Minimum_LOO_pooled_log2FoldChange",
    "Maximum_LOO_pooled_log2FoldChange",
    "All_LOO_same_direction_as_full",
    "Most_influential_omitted_model",
    "Maximum_absolute_effect_change"
]

print(
    "\nLEAVE-ONE-MODEL-OUT ANALYSIS COMPLETE"
)
print(
    "-------------------------------------"
)

print(
    leave_one_out_summary_table[
        summary_display_columns
    ]
    .round(4)
    .to_string(index=False)
)

imr90_display_columns = [
    "Gene",
    "Scenario",
    "Pooled_log2FoldChange",
    "CI95_lower",
    "CI95_upper",
    "I2_percent",
    "Overall_direction"
]

print(
    "\nIMR-90 SHARED-CONTROL SENSITIVITY"
)
print(
    "---------------------------------"
)

print(
    imr90_table[
        imr90_display_columns
    ]
    .round(4)
    .to_string(index=False)
)

print(
    "\nFILES SAVED"
)
print(
    "-----------"
)
print(
    leave_one_out_output
)
print(
    leave_one_out_summary_output
)
print(
    imr90_sensitivity_output
)