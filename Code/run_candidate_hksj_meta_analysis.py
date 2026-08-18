from pathlib import Path

import pandas as pd
from scipy.stats import t
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

full_output = output_directory / (
    "GSE130727_candidate_HKSJ_meta_analysis_summary.tsv"
)

imr90_output = output_directory / (
    "GSE130727_candidate_HKSJ_IMR90_sensitivity.tsv"
)


# ---------------------------------------------------------
# Load and validate input
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
    "lfcSE"
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
    "lfcSE"
]:
    meta_input[column] = pd.to_numeric(
        meta_input[column],
        errors="coerce"
    )

if meta_input[
    ["log2FoldChange", "lfcSE"]
].isna().any().any():
    raise ValueError(
        "Missing or non-numeric values were found."
    )

if (meta_input["lfcSE"] <= 0).any():
    raise ValueError(
        "Every lfcSE value must be greater than zero."
    )

meta_input["variance"] = (
    meta_input["lfcSE"] ** 2
)


# ---------------------------------------------------------
# Interpretation functions
# ---------------------------------------------------------

def classify_direction(ci_lower, ci_upper):
    """Classify the pooled direction using the HKSJ CI."""

    if ci_lower > 0:
        return "Overall increase"

    if ci_upper < 0:
        return "Overall decrease"

    return "No clear overall direction"


def classify_heterogeneity(i2_percent):
    """Classify I-squared descriptively."""

    if i2_percent < 25:
        return "Low"

    if i2_percent < 50:
        return "Moderate"

    if i2_percent < 75:
        return "Substantial"

    return "Considerable"


# ---------------------------------------------------------
# HKSJ random-effects function
# ---------------------------------------------------------

def run_hksj_analysis(subset):
    """
    Paule-Mandel random-effects meta-analysis with
    HKSJ standard error and t-based confidence interval.
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
        row_names=model_names,
        use_t=True,
        alpha=0.05
    )

    pooled_effect = float(
        result.mean_effect_re
    )

    # HKSJ-adjusted standard error
    hksj_se = float(
        result.sd_eff_w_re_hksj
    )

    degrees_of_freedom = int(
        result.df_resid
    )

    # conf_int returns:
    # 0 = fixed effect, standard scale
    # 1 = random effect, standard scale
    # 2 = fixed effect, HKSJ
    # 3 = random effect, HKSJ
    hksj_ci = result.conf_int(
        alpha=0.05,
        use_t=True
    )[3]

    hksj_ci_lower = float(
        hksj_ci[0]
    )

    hksj_ci_upper = float(
        hksj_ci[1]
    )

    if hksj_se > 0:
        t_statistic = (
            pooled_effect / hksj_se
        )

        hksj_pvalue = float(
            2
            * t.sf(
                abs(t_statistic),
                df=degrees_of_freedom
            )
        )
    else:
        t_statistic = float("nan")
        hksj_pvalue = float("nan")

    # Standard normal random-effects CI for comparison
    normal_ci = result.conf_int(
        alpha=0.05,
        use_t=False
    )[1]

    normal_ci_lower = float(
        normal_ci[0]
    )

    normal_ci_upper = float(
        normal_ci[1]
    )

    q_statistic = float(
        result.q
    )

    q_df = len(effects) - 1

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

    return {
        "Models_included": len(effects),
        "Degrees_of_freedom": degrees_of_freedom,
        "Pooled_log2FoldChange": pooled_effect,
        "HKSJ_SE": hksj_se,
        "HKSJ_t_statistic": t_statistic,
        "HKSJ_pvalue": hksj_pvalue,
        "HKSJ_CI95_lower": hksj_ci_lower,
        "HKSJ_CI95_upper": hksj_ci_upper,
        "Normal_CI95_lower": normal_ci_lower,
        "Normal_CI95_upper": normal_ci_upper,
        "Tau_squared": max(
            0.0,
            float(result.tau2)
        ),
        "HKSJ_scale": float(
            result.scale_hksj_re
        ),
        "I2_percent": i2_percent,
        "Heterogeneity":
            classify_heterogeneity(i2_percent),
        "HKSJ_overall_direction":
            classify_direction(
                hksj_ci_lower,
                hksj_ci_upper
            )
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
# Full eight-model HKSJ analysis
# ---------------------------------------------------------

full_rows = []

for gene in gene_order:

    gene_data = meta_input[
        meta_input["Gene"] == gene
    ].copy()

    if gene_data["Comparison"].nunique() != 8:
        raise ValueError(
            f"{gene} does not contain exactly 8 models."
        )

    result = run_hksj_analysis(
        gene_data
    )

    full_rows.append(
        {
            "Gene": gene,
            "Scenario": "All_8_models",
            **result
        }
    )


# ---------------------------------------------------------
# IMR-90 sensitivity with HKSJ intervals
# ---------------------------------------------------------

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

imr90_rows = []

for gene in gene_order:

    gene_data = meta_input[
        meta_input["Gene"] == gene
    ].copy()

    for scenario, excluded_models in (
        imr90_scenarios.items()
    ):

        scenario_data = gene_data[
            ~gene_data["Comparison"].isin(
                excluded_models
            )
        ].copy()

        result = run_hksj_analysis(
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
# Save results
# ---------------------------------------------------------

full_table = pd.DataFrame(
    full_rows
)

imr90_table = pd.DataFrame(
    imr90_rows
)

full_table.to_csv(
    full_output,
    sep="\t",
    index=False
)

imr90_table.to_csv(
    imr90_output,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Terminal output
# ---------------------------------------------------------

full_display_columns = [
    "Gene",
    "Pooled_log2FoldChange",
    "Normal_CI95_lower",
    "Normal_CI95_upper",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent",
    "HKSJ_overall_direction"
]

print(
    "\nFULL HKSJ META-ANALYSIS COMPLETE"
)
print(
    "--------------------------------"
)

print(
    full_table[
        full_display_columns
    ]
    .round(4)
    .to_string(index=False)
)

imr90_display_columns = [
    "Gene",
    "Scenario",
    "Pooled_log2FoldChange",
    "HKSJ_CI95_lower",
    "HKSJ_CI95_upper",
    "HKSJ_pvalue",
    "I2_percent",
    "HKSJ_overall_direction"
]

print(
    "\nHKSJ IMR-90 SENSITIVITY"
)
print(
    "-----------------------"
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
print(full_output)
print(imr90_output)