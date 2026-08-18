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

summary_output = output_directory / (
    "GSE130727_candidate_meta_analysis_summary.tsv"
)

weights_output = output_directory / (
    "GSE130727_candidate_meta_analysis_model_weights.tsv"
)


# ---------------------------------------------------------
# Confirm that files and folders are available
# ---------------------------------------------------------

if not input_file.exists():
    raise FileNotFoundError(
        f"Could not find the input file:\n{input_file.resolve()}"
    )

output_directory.mkdir(
    parents=True,
    exist_ok=True
)


# ---------------------------------------------------------
# Load and validate the PyDESeq2 results
# ---------------------------------------------------------

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
        "The following required columns are missing: "
        + ", ".join(missing_columns)
    )

meta_input = data[required_columns].copy()

numeric_columns = [
    "log2FoldChange",
    "lfcSE",
    "padj"
]

for column in numeric_columns:
    meta_input[column] = pd.to_numeric(
        meta_input[column],
        errors="coerce"
    )

if meta_input[numeric_columns].isna().any().any():
    raise ValueError(
        "Missing or non-numeric values were found in the "
        "meta-analysis input."
    )

if (meta_input["lfcSE"] <= 0).any():
    raise ValueError(
        "Every lfcSE value must be greater than zero."
    )

# Sampling variance of each log2 fold-change estimate
meta_input["variance"] = meta_input["lfcSE"] ** 2


# ---------------------------------------------------------
# Functions used to interpret heterogeneity
# ---------------------------------------------------------

def classify_heterogeneity(i2_percent):
    """Return a descriptive category for I-squared."""

    if i2_percent < 25:
        return "Low"
    if i2_percent < 50:
        return "Moderate"
    if i2_percent < 75:
        return "Substantial"

    return "Considerable"


def classify_overall_effect(ci_lower, ci_upper):
    """Classify the direction of the pooled effect."""

    if ci_lower > 0:
        return "Overall increase"

    if ci_upper < 0:
        return "Overall decrease"

    return "No clear overall direction"


# ---------------------------------------------------------
# Random-effects meta-analysis
# ---------------------------------------------------------

gene_order = [
    "CDKN2B",
    "IL1A",
    "CXCL8",
    "ICAM1",
    "FOS"
]

summary_rows = []
weight_rows = []

for gene in gene_order:

    gene_data = (
        meta_input[
            meta_input["Gene"] == gene
        ]
        .copy()
        .sort_values("Comparison")
        .reset_index(drop=True)
    )

    if gene_data.empty:
        raise ValueError(
            f"No results were found for {gene}."
        )

    effects = gene_data[
        "log2FoldChange"
    ].to_numpy(dtype=float)

    variances = gene_data[
        "variance"
    ].to_numpy(dtype=float)

    model_names = gene_data[
        "Comparison"
    ].astype(str).tolist()

    # Paule-Mandel random-effects estimator
    result = combine_effects(
        effect=effects,
        variance=variances,
        method_re="iterated",
        row_names=model_names
    )

    pooled_log2fc = float(
        result.mean_effect_re
    )

    pooled_se = float(
        result.sd_eff_w_re
    )

    ci_lower = pooled_log2fc - 1.96 * pooled_se
    ci_upper = pooled_log2fc + 1.96 * pooled_se

    if pooled_se > 0:
        z_statistic = pooled_log2fc / pooled_se
        pooled_pvalue = 2 * norm.sf(
            abs(z_statistic)
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

    overall_effect = classify_overall_effect(
        ci_lower,
        ci_upper
    )

    heterogeneity_category = classify_heterogeneity(
        i2_percent
    )

    summary_rows.append(
        {
            "Gene": gene,
            "Models": len(effects),
            "Estimator": "Paule-Mandel random effects",
            "Pooled_log2FoldChange": pooled_log2fc,
            "Pooled_SE": pooled_se,
            "CI95_lower": ci_lower,
            "CI95_upper": ci_upper,
            "Pooled_pvalue": pooled_pvalue,
            "Pooled_fold_change": 2 ** pooled_log2fc,
            "Q": q_statistic,
            "Q_df": q_df,
            "Q_pvalue": q_pvalue,
            "Tau_squared": tau_squared,
            "I2_percent": i2_percent,
            "Heterogeneity": heterogeneity_category,
            "Overall_effect": overall_effect
        }
    )

    random_weights = np.asarray(
        result.weights_re,
        dtype=float
    )

    weight_percentages = (
        random_weights
        / random_weights.sum()
        * 100
    )

    for index, row in gene_data.iterrows():

        weight_rows.append(
            {
                "Gene": gene,
                "Comparison": row["Comparison"],
                "log2FoldChange": row["log2FoldChange"],
                "lfcSE": row["lfcSE"],
                "padj": row["padj"],
                "Random_effects_weight_percent":
                    weight_percentages[index]
            }
        )


# ---------------------------------------------------------
# Save the results
# ---------------------------------------------------------

summary_table = pd.DataFrame(
    summary_rows
)

weights_table = pd.DataFrame(
    weight_rows
)

summary_table.to_csv(
    summary_output,
    sep="\t",
    index=False
)

weights_table.to_csv(
    weights_output,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Print a concise terminal summary
# ---------------------------------------------------------

display_columns = [
    "Gene",
    "Pooled_log2FoldChange",
    "CI95_lower",
    "CI95_upper",
    "Pooled_pvalue",
    "I2_percent",
    "Heterogeneity",
    "Overall_effect"
]

print("\nRANDOM-EFFECTS META-ANALYSIS COMPLETE")
print("-------------------------------------")
print(
    summary_table[
        display_columns
    ].round(4).to_string(index=False)
)

print("\nFILES SAVED")
print("-----------")
print(summary_output)
print(weights_output)