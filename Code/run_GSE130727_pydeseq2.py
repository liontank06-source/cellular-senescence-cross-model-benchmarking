from contextlib import redirect_stdout
from pathlib import Path
import io

import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


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

design_file = (
    input_folder
    / "GSE130727_comparison_design.tsv"
)

mapping_file = (
    input_folder
    / "GSE130727_candidate_ensembl_mapping.tsv"
)

results_folder = (
    output_folder
    / "DESeq2_results"
)

candidate_output = (
    output_folder
    / "GSE130727_pydeseq2_candidate_results.tsv"
)

run_summary_output = (
    output_folder
    / "GSE130727_pydeseq2_run_summary.tsv"
)


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


def classify_result(log2fc, padj):
    """Classify the formal DESeq2 result."""

    if pd.isna(padj):
        return "padj unavailable"

    if padj >= 0.05:
        return "Not significant"

    if log2fc > 0:
        return "Significant increase"

    if log2fc < 0:
        return "Significant decrease"

    return "Significant, no directional change"


def main():

    results_folder.mkdir(
        parents=True,
        exist_ok=True
    )

    print("\nReading count matrix...")

    matrix = pd.read_csv(
        matrix_file,
        sep="\t",
        compression="gzip"
    )

    design = pd.read_csv(
        design_file,
        sep="\t"
    )

    mapping = pd.read_csv(
        mapping_file,
        sep="\t"
    )

    # Clean Ensembl IDs
    matrix["gene"] = (
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

    if matrix["gene"].duplicated().any():
        raise ValueError(
            "Duplicated gene IDs were found in the count matrix."
        )

    matrix = matrix.set_index("gene")

    candidate_template = (
        mapping[
            [
                "Gene",
                "Ensembl_ID"
            ]
        ]
        .drop_duplicates("Ensembl_ID")
        .set_index("Ensembl_ID")
    )

    candidate_results = []
    run_summary_rows = []

    # Use two CPU cores to avoid overwhelming the laptop
    inference = DefaultInference(
        n_cpus=2
    )

    for comparison in comparison_order:

        model_design = design[
            design["Comparison"] == comparison
        ].copy()

        # A GSM should appear only once inside each comparison
        model_design = model_design.drop_duplicates(
            subset="GSM"
        )

        control_n = (
            model_design["Status"]
            .eq("Control")
            .sum()
        )

        senescent_n = (
            model_design["Status"]
            .eq("Senescent")
            .sum()
        )

        print(
            f"\nRunning {comparison}: "
            f"{control_n} Control vs "
            f"{senescent_n} Senescent"
        )

        if control_n < 2 or senescent_n < 2:
            raise ValueError(
                f"{comparison} does not have enough replicates."
            )

        sample_ids = (
            model_design["GSM"]
            .astype(str)
            .tolist()
        )

        missing_samples = [
            gsm
            for gsm in sample_ids
            if gsm not in matrix.columns
        ]

        if missing_samples:
            raise ValueError(
                f"{comparison} is missing samples: "
                f"{missing_samples}"
            )

        # PyDESeq2 requires samples as rows and genes as columns
        counts = (
            matrix[sample_ids]
            .transpose()
            .copy()
        )

        counts = counts.apply(
            pd.to_numeric,
            errors="raise"
        )

        counts = counts.round().astype(int)

        if (counts < 0).any().any():
            raise ValueError(
                f"Negative counts found in {comparison}."
            )

        # Remove genes with extremely low total counts
        genes_to_keep = (
            counts.sum(axis=0) >= 10
        )

        counts = counts.loc[
            :,
            genes_to_keep
        ]

        print(
            "Genes retained after filtering:",
            counts.shape[1]
        )

        metadata = (
            model_design
            .set_index("GSM")
            .loc[
                counts.index,
                ["Status"]
            ]
            .copy()
        )

        metadata["Status"] = pd.Categorical(
            metadata["Status"],
            categories=[
                "Control",
                "Senescent"
            ],
            ordered=True
        )

        # Fit the negative-binomial model
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design="~Status",
            refit_cooks=False,
            inference=inference,
            quiet=True,
            low_memory=True
        )

        dds.deseq2()

        # Positive log2FC means Senescent > Control
        stats = DeseqStats(
            dds,
            contrast=[
                "Status",
                "Senescent",
                "Control"
            ],
            alpha=0.05,
            cooks_filter=True,
            independent_filter=True,
            inference=inference,
            quiet=True
        )

        # summary() creates results_df.
        # Redirect its large terminal output.
        with redirect_stdout(io.StringIO()):
            stats.summary()

        results = stats.results_df.copy()

        results.index = (
            results.index
            .astype(str)
            .str.split(".")
            .str[0]
        )

        results.index.name = "Ensembl_ID"

        results = results.reset_index()

        results["Gene"] = (
            results["Ensembl_ID"]
            .map(
                mapping.set_index(
                    "Ensembl_ID"
                )["Gene"]
            )
        )

        results["Comparison"] = comparison

        # Save the complete result for this model
        complete_output = (
            results_folder
            / (
                f"{comparison}_"
                "DESeq2_all_genes.tsv.gz"
            )
        )

        results.to_csv(
            complete_output,
            sep="\t",
            index=False,
            compression="gzip"
        )

        # Preserve all five candidates, even if one was filtered
        candidate_model = candidate_template.join(
            results.set_index(
                "Ensembl_ID"
            )[
                [
                    "baseMean",
                    "log2FoldChange",
                    "lfcSE",
                    "stat",
                    "pvalue",
                    "padj"
                ]
            ],
            how="left"
        ).reset_index()

        candidate_model["Comparison"] = comparison
        candidate_model["Control_n"] = control_n
        candidate_model["Senescent_n"] = senescent_n

        candidate_model["Tested_by_DESeq2"] = (
            candidate_model["baseMean"]
            .notna()
        )

        candidate_model["Direction"] = np.select(
            [
                candidate_model[
                    "log2FoldChange"
                ] > 0,
                candidate_model[
                    "log2FoldChange"
                ] < 0
            ],
            [
                "Higher in Senescent",
                "Higher in Control"
            ],
            default="No change or unavailable"
        )

        candidate_model["Formal_result"] = (
            candidate_model.apply(
                lambda row: classify_result(
                    row["log2FoldChange"],
                    row["padj"]
                ),
                axis=1
            )
        )

        candidate_results.append(
            candidate_model
        )

        run_summary_rows.append(
            {
                "Comparison": comparison,
                "Control_n": control_n,
                "Senescent_n": senescent_n,
                "Genes_before_filtering": (
                    matrix.shape[0]
                ),
                "Genes_tested": (
                    counts.shape[1]
                ),
                "Full_results_file": (
                    complete_output.name
                )
            }
        )

        print(
            f"{comparison} completed successfully."
        )

    # Combine candidate results from all eight models
    candidate_results = pd.concat(
        candidate_results,
        ignore_index=True
    )

    candidate_results["Comparison"] = pd.Categorical(
        candidate_results["Comparison"],
        categories=comparison_order,
        ordered=True
    )

    gene_order = [
        "FOS",
        "IL1A",
        "CXCL8",
        "CDKN2B",
        "ICAM1"
    ]

    candidate_results["Gene"] = pd.Categorical(
        candidate_results["Gene"],
        categories=gene_order,
        ordered=True
    )

    candidate_results = candidate_results.sort_values(
        [
            "Comparison",
            "Gene"
        ]
    )

    candidate_results.to_csv(
        candidate_output,
        sep="\t",
        index=False
    )

    run_summary = pd.DataFrame(
        run_summary_rows
    )

    run_summary.to_csv(
        run_summary_output,
        sep="\t",
        index=False
    )

    print(
        "\nAll eight PyDESeq2 analyses "
        "completed successfully."
    )

    print("\nCandidate result counts:")
    print(
        candidate_results[
            "Formal_result"
        ].value_counts()
    )

    print("\nCandidate results:")
    print(
        candidate_results[
            [
                "Comparison",
                "Gene",
                "baseMean",
                "log2FoldChange",
                "pvalue",
                "padj",
                "Formal_result"
            ]
        ]
        .round(5)
        .to_string(index=False)
    )

    print("\nCandidate results saved as:")
    print(candidate_output)

    print("\nRun summary saved as:")
    print(run_summary_output)

    print("\nFull results folder:")
    print(results_folder)


if __name__ == "__main__":
    main()