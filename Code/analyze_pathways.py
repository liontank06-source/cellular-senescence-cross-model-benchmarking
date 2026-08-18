import pandas as pd

file_path = "Data/BJ_Discovery/gProfiler_BJ_Senescent_enrichment_results.csv"

results = pd.read_csv(file_path)

print("Columns:")
print(results.columns.tolist())

sasp = results[
    (results["source"] == "REAC")
    & results["term_name"].str.contains(
        "Senescence-Associated Secretory Phenotype",
        case=False,
        na=False
    )
]

print("\nSASP result:")
print(sasp.to_string(index=False))

# Extract the genes responsible for the SASP enrichment
sasp_genes = sasp.iloc[0]["intersections"].split(",")

print("\nGenes overlapping with the SASP pathway:")
for gene in sasp_genes:
    print(gene)

print("\nNumber of SASP genes:", len(sasp_genes))

# Pathways that are most relevant to our senescence project
pathways_to_check = [
    "Senescence-Associated Secretory Phenotype (SASP)",
    "Cellular Senescence",
    "DNA Damage/Telomere Stress Induced Senescence",
    "Oxidative Stress Induced Senescence",
    "Cytokine Signaling in Immune system",
    "Extracellular matrix organization"
]

pathway_summary = []
pathway_gene_rows = []

for pathway_name in pathways_to_check:

    pathway_result = results[
        (results["source"] == "REAC")
        & results["term_name"].str.contains(
            pathway_name,
            case=False,
            na=False,
            regex=False
        )
    ]

    if pathway_result.empty:
        print(f"\nPathway not found: {pathway_name}")
        continue

    row = pathway_result.iloc[0]

    genes = [
        gene.strip()
        for gene in str(row["intersections"]).split(",")
        if gene.strip()
    ]

    print(f"\nPathway: {row['term_name']}")
    print("Adjusted p-value:", row["adjusted_p_value"])
    print("Intersection size:", row["intersection_size"])
    print("Genes:")

    for gene in genes:
        print(gene)

        pathway_gene_rows.append({
            "Pathway": row["term_name"],
            "Gene": gene
        })

    pathway_summary.append({
        "Pathway": row["term_name"],
        "Adjusted_p_value": row["adjusted_p_value"],
        "Negative_log10_padj": row[
            "negative_log10_of_adjusted_p_value"
        ],
        "Pathway_size": row["term_size"],
        "Intersection_size": row["intersection_size"]
    })

# Save a summary of the selected pathways
pd.DataFrame(pathway_summary).to_csv(
    "Results/BJ_Discovery/BJ_selected_pathways_summary.tsv",
    sep="\t",
    index=False
)

# Save every gene associated with each selected pathway
pd.DataFrame(pathway_gene_rows).to_csv(
    "Results/BJ_Discovery/BJ_selected_pathway_genes.tsv",
    sep="\t",
    index=False
)

# Compare genes across the selected pathways
pathway_genes = pd.read_csv(
    "Results/BJ_Discovery/BJ_selected_pathway_genes.tsv",
    sep="\t"
)

senescent_expression = pd.read_csv(
    "Results/BJ_Discovery/BJ_higher_in_Senescent.tsv",
    sep="\t"
)

# Count how many different pathways contain each gene
gene_counts = (
    pathway_genes.groupby("Gene")["Pathway"]
    .nunique()
    .reset_index(name="Pathway_count")
)

# Put the pathway names together for each gene
gene_pathways = (
    pathway_genes.groupby("Gene")["Pathway"]
    .apply(lambda pathways: " | ".join(sorted(set(pathways))))
    .reset_index(name="Pathways")
)

# Combine pathway information
gene_summary = gene_counts.merge(
    gene_pathways,
    on="Gene"
)

# Add gene-expression results
gene_summary = gene_summary.merge(
    senescent_expression[
        ["Symbol", "log2FoldChange", "padj", "baseMean", "Description"]
    ],
    left_on="Gene",
    right_on="Symbol",
    how="left"
)

gene_summary = gene_summary.drop(columns="Symbol")

# Keep genes appearing in at least two pathways
repeated_genes = gene_summary[
    gene_summary["Pathway_count"] >= 2
].copy()

repeated_genes["abs_log2FoldChange"] = (
    repeated_genes["log2FoldChange"].abs()
)

repeated_genes = repeated_genes.sort_values(
    ["Pathway_count", "abs_log2FoldChange"],
    ascending=[False, False]
)

print("\nGenes appearing in multiple pathways:")
print(
    repeated_genes[
        [
            "Gene",
            "Pathway_count",
            "log2FoldChange",
            "baseMean",
            "Pathways"
        ]
    ].to_string(index=False)
)

repeated_genes.to_csv(
    "Results/BJ_Discovery/BJ_repeated_pathway_genes.tsv",
    sep="\t",
    index=False
)

# Remove histone genes temporarily to focus on interpretable candidates
non_histone_candidates = repeated_genes[
    ~repeated_genes["Gene"].str.match(r"^H[1-4]", na=False)
].copy()

# Prioritize genes by:
# 1. Number of relevant pathways
# 2. Magnitude of expression change
# 3. Expression level
prioritized_candidates = non_histone_candidates.sort_values(
    ["Pathway_count", "abs_log2FoldChange", "baseMean"],
    ascending=[False, False, False]
)

print("\nPrioritized non-histone candidate genes:")
print(
    prioritized_candidates[
        [
            "Gene",
            "Pathway_count",
            "log2FoldChange",
            "padj",
            "baseMean",
            "Description",
            "Pathways"
        ]
    ].head(20).to_string(index=False)
)

prioritized_candidates.to_csv(
    "Results/BJ_Discovery/BJ_prioritized_senescence_candidates.tsv",
    sep="\t",
    index=False
)
