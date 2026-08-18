import pandas as pd

file_path = "Data/BJ_Discovery/GEO2R_BJ_Young_vs_Senescent_full_table.tsv"

genes = pd.read_csv(file_path, sep="\t")

print(genes.head())
print("Number of rows and columns:", genes.shape)
print("Columns:", genes.columns.tolist())

print("\nData types:")
print(genes.dtypes)

print("\nMissing values:")
print(genes.isna().sum())

# Keep genes that are statistically significant and change strongly
significant_genes = genes[
    (genes["padj"] < 0.05)
    & (genes["log2FoldChange"].abs() > 1)
].copy()

print("\nNumber of significant genes:", len(significant_genes))

print("\nFirst 10 significant genes:")
print(
    significant_genes[
        ["Symbol", "log2FoldChange", "padj", "baseMean", "Description"]
    ].head(10)
)

# Save the filtered table
significant_genes.to_csv(
    "Results/BJ_Discovery/BJ_significant_genes.tsv",
    sep="\t",
    index=False
)
# Separate significant genes by direction of change
positive_genes = significant_genes[
    significant_genes["log2FoldChange"] > 1
].sort_values("log2FoldChange", ascending=False).copy()

negative_genes = significant_genes[
    significant_genes["log2FoldChange"] < -1
].sort_values("log2FoldChange", ascending=True).copy()

print("\nPositive log2FoldChange genes:", len(positive_genes))
print("Negative log2FoldChange genes:", len(negative_genes))

columns_to_show = [
    "Symbol",
    "log2FoldChange",
    "padj",
    "baseMean",
    "Description"
]

print("\nTop 10 positive genes:")
print(positive_genes[columns_to_show].head(10))

print("\nTop 10 negative genes:")
print(negative_genes[columns_to_show].head(10))

# Save both tables
positive_genes.to_csv(
    "Results/BJ_Discovery/BJ_positive_log2FC_genes.tsv",
    sep="\t",
    index=False
)

negative_genes.to_csv(
    "Results/BJ_Discovery/BJ_negative_log2FC_genes.tsv",
    sep="\t",
    index=False
)
# Rename the two groups using the confirmed comparison direction
higher_in_young = positive_genes.copy()
higher_in_senescent = negative_genes.copy()

print("\nGenes higher in Young:", len(higher_in_young))
print("Genes higher in Senescent:", len(higher_in_senescent))

higher_in_young.to_csv(
    "Results/BJ_Discovery/BJ_higher_in_Young.tsv",
    sep="\t",
    index=False
)

higher_in_senescent.to_csv(
    "Results/BJ_Discovery/BJ_higher_in_Senescent.tsv",
    sep="\t",
    index=False
)
# Select strong senescent candidates with reasonable expression
top_senescent_candidates = higher_in_senescent[
    higher_in_senescent["baseMean"] >= 10
].head(20).copy()

columns_to_show = [
    "Symbol",
    "log2FoldChange",
    "padj",
    "baseMean",
    "Description"
]

print("\nTop 20 genes higher in Senescent:")
print(top_senescent_candidates[columns_to_show])

top_senescent_candidates[columns_to_show].to_csv(
    "Results/BJ_Discovery/BJ_top20_Senescent_candidates.tsv",
    sep="\t",
    index=False
)
# Create a clean gene-symbol list for pathway analysis
senescent_gene_list = (
    higher_in_senescent["Symbol"]
    .dropna()
    .drop_duplicates()
)

print("\nGenes available for pathway analysis:", len(senescent_gene_list))
print(senescent_gene_list.head(20))

senescent_gene_list.to_csv(
    "Results/BJ_Discovery/BJ_Senescent_gene_symbols.txt",
    index=False,
    header=False
)
# Create the background list from all genes tested by GEO2R
background_gene_list = (
    genes["Symbol"]
    .dropna()
    .drop_duplicates()
)

print("\nGenes available as background:", len(background_gene_list))

background_gene_list.to_csv(
    "Results/BJ_Discovery/BJ_background_gene_symbols.txt",
    index=False,
    header=False
)