from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------
# File locations
# ---------------------------------------------------------

candidate_results_file = Path(
    "Results/GSE130727_Multimodel/"
    "GSE130727_pydeseq2_candidate_results.tsv"
)

hksj_summary_file = Path(
    "Results/MetaAnalysis/"
    "GSE130727_candidate_HKSJ_meta_analysis_summary.tsv"
)

figure_directory = Path(
    "Figures/MetaAnalysis"
)

results_directory = Path(
    "Results/MetaAnalysis"
)

forest_data_output = results_directory / (
    "GSE130727_candidate_forest_plot_data.tsv"
)


# ---------------------------------------------------------
# Confirm files and folders
# ---------------------------------------------------------

if not candidate_results_file.exists():
    raise FileNotFoundError(
        f"Could not find:\n{candidate_results_file.resolve()}"
    )

if not hksj_summary_file.exists():
    raise FileNotFoundError(
        f"Could not find:\n{hksj_summary_file.resolve()}"
    )

figure_directory.mkdir(
    parents=True,
    exist_ok=True
)

results_directory.mkdir(
    parents=True,
    exist_ok=True
)


# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

candidate_data = pd.read_csv(
    candidate_results_file,
    sep="\t"
)

hksj_data = pd.read_csv(
    hksj_summary_file,
    sep="\t"
)


# ---------------------------------------------------------
# Validate required columns
# ---------------------------------------------------------

candidate_required_columns = [
    "Gene",
    "Comparison",
    "log2FoldChange",
    "lfcSE",
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

missing_candidate_columns = [
    column
    for column in candidate_required_columns
    if column not in candidate_data.columns
]

missing_hksj_columns = [
    column
    for column in hksj_required_columns
    if column not in hksj_data.columns
]

if missing_candidate_columns:
    raise ValueError(
        "Missing candidate columns: "
        + ", ".join(missing_candidate_columns)
    )

if missing_hksj_columns:
    raise ValueError(
        "Missing HKSJ columns: "
        + ", ".join(missing_hksj_columns)
    )


# ---------------------------------------------------------
# Convert numeric columns
# ---------------------------------------------------------

for column in [
    "log2FoldChange",
    "lfcSE",
    "padj"
]:
    candidate_data[column] = pd.to_numeric(
        candidate_data[column],
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

if candidate_data[
    ["log2FoldChange", "lfcSE"]
].isna().any().any():
    raise ValueError(
        "Missing or non-numeric model effects were found."
    )

if hksj_data[
    [
        "Pooled_log2FoldChange",
        "HKSJ_CI95_lower",
        "HKSJ_CI95_upper"
    ]
].isna().any().any():
    raise ValueError(
        "Missing or non-numeric HKSJ results were found."
    )


# ---------------------------------------------------------
# Model and gene order
# ---------------------------------------------------------

model_order = [
    "HAEC_IR",
    "HUVEC_IR",
    "IMR90_IR",
    "IMR90_Replicative",
    "WI38_Doxorubicin",
    "WI38_HRAS",
    "WI38_IR",
    "WI38_Replicative"
]

model_label_map = {
    "HAEC_IR": "HAEC — ionizing radiation",
    "HUVEC_IR": "HUVEC — ionizing radiation",
    "IMR90_IR": "IMR-90 — ionizing radiation",
    "IMR90_Replicative": "IMR-90 — replicative exhaustion",
    "WI38_Doxorubicin": "WI-38 — doxorubicin",
    "WI38_HRAS": "WI-38 — HRAS",
    "WI38_IR": "WI-38 — ionizing radiation",
    "WI38_Replicative": "WI-38 — replicative exhaustion"
}

gene_order = [
    "CDKN2B",
    "IL1A",
    "CXCL8",
    "ICAM1",
    "FOS"
]


# ---------------------------------------------------------
# Calculate individual model confidence intervals
# ---------------------------------------------------------

candidate_data["Model_CI95_lower"] = (
    candidate_data["log2FoldChange"]
    - 1.96 * candidate_data["lfcSE"]
)

candidate_data["Model_CI95_upper"] = (
    candidate_data["log2FoldChange"]
    + 1.96 * candidate_data["lfcSE"]
)


# ---------------------------------------------------------
# Store all plotted data
# ---------------------------------------------------------

forest_data_rows = []


# ---------------------------------------------------------
# Generate one forest plot per gene
# ---------------------------------------------------------

for gene in gene_order:

    gene_data = (
        candidate_data[
            candidate_data["Gene"] == gene
        ]
        .copy()
    )

    if gene_data["Comparison"].nunique() != 8:
        raise ValueError(
            f"{gene} does not contain exactly 8 models."
        )

    gene_data["Comparison"] = pd.Categorical(
        gene_data["Comparison"],
        categories=model_order,
        ordered=True
    )

    gene_data = (
        gene_data
        .sort_values("Comparison")
        .reset_index(drop=True)
    )

    pooled_rows = hksj_data[
        hksj_data["Gene"] == gene
    ]

    if "Scenario" in hksj_data.columns:
        pooled_rows = pooled_rows[
            pooled_rows["Scenario"] == "All_8_models"
        ]

    if len(pooled_rows) != 1:
        raise ValueError(
            f"Expected one full HKSJ result for {gene}, "
            f"but found {len(pooled_rows)}."
        )

    pooled_row = pooled_rows.iloc[0]

    pooled_effect = float(
        pooled_row["Pooled_log2FoldChange"]
    )

    pooled_lower = float(
        pooled_row["HKSJ_CI95_lower"]
    )

    pooled_upper = float(
        pooled_row["HKSJ_CI95_upper"]
    )

    pooled_pvalue = float(
        pooled_row["HKSJ_pvalue"]
    )

    i2_percent = float(
        pooled_row["I2_percent"]
    )

    # Highest model appears at the top
    model_y_positions = np.arange(
        len(gene_data),
        0,
        -1
    )

    pooled_y_position = 0

    effects = gene_data[
        "log2FoldChange"
    ].to_numpy(dtype=float)

    lower_limits = gene_data[
        "Model_CI95_lower"
    ].to_numpy(dtype=float)

    upper_limits = gene_data[
        "Model_CI95_upper"
    ].to_numpy(dtype=float)

    lower_errors = effects - lower_limits
    upper_errors = upper_limits - effects

    figure, axis = plt.subplots(
        figsize=(12.5, 7)
    )

    # Individual model results
    axis.errorbar(
        effects,
        model_y_positions,
        xerr=[
            lower_errors,
            upper_errors
        ],
        fmt="o",
        markersize=7,
        capsize=4,
        linewidth=1.5,
        label="Individual senescence model"
    )

    # Pooled HKSJ result
    axis.errorbar(
        [pooled_effect],
        [pooled_y_position],
        xerr=[
            [pooled_effect - pooled_lower],
            [pooled_upper - pooled_effect]
        ],
        fmt="D",
        markersize=9,
        capsize=5,
        linewidth=2.0,
        label="Random-effects pooled estimate"
    )

    # Null-effect reference line
    axis.axvline(
        x=0,
        linestyle="--",
        linewidth=1.2
    )

    y_positions = list(
        model_y_positions
    ) + [pooled_y_position]

    y_labels = [
        model_label_map[model]
        for model in model_order
    ] + ["Pooled HKSJ estimate"]

    axis.set_yticks(
        y_positions
    )

    axis.set_yticklabels(
        y_labels,
        fontsize=10
    )

    axis.set_xlabel(
        "log2 fold change (Senescent versus Control)",
        fontsize=11
    )

    figure.suptitle(
        f"{gene}: Cross-Model Random-Effects Meta-Analysis",
        fontsize=18,
        fontweight="bold",
        y=0.965
    )

    summary_text = (
    f"Pooled log2FC = {pooled_effect:.2f} "
    f"[HKSJ 95% CI: {pooled_lower:.2f}, {pooled_upper:.2f}]"
    f"  |  HKSJ p = {pooled_pvalue:.3f}"
    f"  |  I² = {i2_percent:.1f}%"
)

    figure.text(
        0.5,
        0.922,
        summary_text,
        ha="center",
        va="center",
        fontsize=11
    )

    axis.grid(
        axis="x",
        linestyle=":",
        alpha=0.5
    )

    axis.legend(
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0,
    frameon=True,
    fontsize=9
    )

    figure.text(
        0.5,
        0.015,
        (
            "Individual-model intervals: log2FC ± 1.96 × lfcSE. "
            "Pooled interval: Paule–Mandel random effects with HKSJ adjustment."
        ),
        horizontalalignment="center",
        fontsize=8,
        fontstyle="italic"
    )

    figure.tight_layout(
        rect=[
            0.02,
            0.06,
            0.80,
            0.90
        ]
    )

    png_output = figure_directory / (
        f"ForestPlot_{gene}_GSE130727_HKSJ.png"
    )

    pdf_output = figure_directory / (
        f"ForestPlot_{gene}_GSE130727_HKSJ.pdf"
    )

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

    # Save model rows used in each forest plot
    for _, row in gene_data.iterrows():

        forest_data_rows.append(
            {
                "Gene": gene,
                "Result_type": "Individual_model",
                "Comparison": str(row["Comparison"]),
                "log2FoldChange":
                    row["log2FoldChange"],
                "CI95_lower":
                    row["Model_CI95_lower"],
                "CI95_upper":
                    row["Model_CI95_upper"],
                "lfcSE":
                    row["lfcSE"],
                "padj":
                    row["padj"]
            }
        )

    forest_data_rows.append(
        {
            "Gene": gene,
            "Result_type": "Pooled_HKSJ",
            "Comparison": "All_8_models",
            "log2FoldChange": pooled_effect,
            "CI95_lower": pooled_lower,
            "CI95_upper": pooled_upper,
            "lfcSE": pd.NA,
            "padj": pooled_pvalue
        }
    )


# ---------------------------------------------------------
# Save plotted values
# ---------------------------------------------------------

forest_data_table = pd.DataFrame(
    forest_data_rows
)

forest_data_table.to_csv(
    forest_data_output,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------

print(
    "\nFOREST PLOTS CREATED SUCCESSFULLY"
)

print(
    "---------------------------------"
)

for gene in gene_order:
    print(
        figure_directory
        / f"ForestPlot_{gene}_GSE130727_HKSJ.png"
    )

print(
    "\nPLOT DATA SAVED"
)

print(
    "---------------"
)

print(
    forest_data_output
)