from pathlib import Path
import pandas as pd


# Input files
endpoint_file = Path(
    "Data/GSE175533_WI38/"
    "GSE175533_candidate_statistics.tsv"
)

time_course_file = Path(
    "Data/GSE175533_WI38/"
    "GSE175533_RS_deseq2_candidate_results.tsv"
)

# Final output file
output_file = Path(
    "Results/GSE175533_WI38/"
    "GSE175533_validation_summary.tsv"
)

candidate_order = [
    "FOS",
    "IL1A",
    "CXCL8",
    "CDKN2B",
    "ICAM1"
]

# Read both analyses
endpoint = pd.read_csv(endpoint_file, sep="\t")
time_course = pd.read_csv(time_course_file, sep="\t")

# Rename columns so their meaning is clear
endpoint = endpoint.rename(
    columns={
        "Gene": "gene",
        "PDL20_mean_TPM": "Endpoint_PDL20_mean_TPM",
        "PDL53_mean_TPM": "Endpoint_PDL53_mean_TPM",
        "log2FC_PDL53_vs_PDL20": "Endpoint_log2FC",
        "p_value": "Endpoint_p_value",
        "adjusted_p_value": "Endpoint_adjusted_p_value"
    }
)

time_course = time_course.rename(
    columns={
        "baseMean": "RS_baseMean",
        "log2FoldChange": "RS_time_log2FC",
        "pvalue": "RS_time_p_value",
        "padj": "RS_time_padj"
    }
)

# Keep only the columns needed for the final table
endpoint_columns = [
    "gene",
    "Endpoint_PDL20_mean_TPM",
    "Endpoint_PDL53_mean_TPM",
    "Endpoint_log2FC",
    "Endpoint_p_value",
    "Endpoint_adjusted_p_value"
]

time_course_columns = [
    "gene",
    "RS_baseMean",
    "RS_time_log2FC",
    "RS_time_p_value",
    "RS_time_padj"
]

endpoint = endpoint[endpoint_columns]
time_course = time_course[time_course_columns]

# Combine the two analyses using the gene name
summary = endpoint.merge(
    time_course,
    on="gene",
    how="outer"
)


def classify_validation(row):
    endpoint_increase = row["Endpoint_log2FC"] > 0
    endpoint_significant = (
        row["Endpoint_adjusted_p_value"] < 0.05
    )

    time_increase = row["RS_time_log2FC"] > 0
    time_significant = row["RS_time_padj"] < 0.05

    if (
        endpoint_increase
        and endpoint_significant
        and time_increase
        and time_significant
    ):
        return "Strong validation"

    if (
        endpoint_increase
        and endpoint_significant
        and time_increase
        and not time_significant
    ):
        return "Partial validation"

    if time_significant and not time_increase:
        return "Not validated as increasing marker"

    return "Inconclusive"


summary["Validation_class"] = summary.apply(
    classify_validation,
    axis=1
)


def biological_interpretation(row):
    gene = row["gene"]

    interpretations = {
        "FOS": (
            "No endpoint increase and significant decrease across "
            "the time course; likely context- or time-dependent."
        ),
        "IL1A": (
            "Strong increase in both analyses, although endpoint TPM "
            "remains relatively low."
        ),
        "CXCL8": (
            "Strong late-stage increase, but the change is not "
            "significant across the entire time course."
        ),
        "CDKN2B": (
            "Strong and consistent increase supporting "
            "senescence-associated cell-cycle arrest."
        ),
        "ICAM1": (
            "Clear endpoint increase with a positive but "
            "not statistically significant time-course trend."
        )
    }

    return interpretations.get(
        gene,
        "No interpretation available."
    )


summary["Interpretation"] = summary.apply(
    biological_interpretation,
    axis=1
)

# Preserve our preferred gene order
summary["gene"] = pd.Categorical(
    summary["gene"],
    categories=candidate_order,
    ordered=True
)

summary = summary.sort_values("gene")

# Save the final summary
summary.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nGSE175533 validation summary:")
print(
    summary[
        [
            "gene",
            "Endpoint_log2FC",
            "Endpoint_adjusted_p_value",
            "RS_time_log2FC",
            "RS_time_padj",
            "Validation_class"
        ]
    ].round(5).to_string(index=False)
)

print("\nSummary saved as:")
print(output_file)