# %%

import pandas as pd

# Read the data files
genes_df = pd.read_csv("output/ora_compound_gene.csv")
kegg_df = pd.read_csv("output/ora_gene_pathway_kegg.csv")
kegg_df["set_name"] = kegg_df["set_name"].str.replace(" - Homo sapiens (human)", "")
wiki_df = pd.read_csv("output/ora_gene_pathway_wikipathways.csv")


# %%
# Function to format p-values and FDR
def format_value(x):
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"


# Prepare gene table (show only fdr < 1)
gene_table = genes_df[genes_df["fdr"] < 0.99].copy()
gene_table["p_value"] = gene_table["p_value"].apply(format_value)
gene_table["fdr"] = gene_table["fdr"].apply(format_value)
gene_md = gene_table[
    ["set_name", "overlap_size", "set_size", "p_value", "fdr"]
].to_markdown(
    index=False,
    headers=["Gene", "Overlap", "Set Size", "P-value", "FDR"],
)

# %%
# Prepare pathway table
# Add source column to each dataset
kegg_df["source"] = "KEGG"
wiki_df["source"] = "Wiki"

# Combine pathway data
combined_pathways = pd.concat([kegg_df, wiki_df])
pathway_table = combined_pathways[combined_pathways["fdr"] < 0.99].copy()
pathway_table["p_value"] = pathway_table["p_value"].apply(format_value)
pathway_table["fdr"] = pathway_table["fdr"].apply(format_value)

pathway_md = pathway_table[
    ["set_name", "source", "overlap_size", "set_size", "p_value", "fdr"]
].to_markdown(
    index=False,
    headers=["Pathway", "Source", "Overlap", "Set Size", "P-value", "FDR"],
)

# %%

# Print the markdown tables
print("Table 1: Gene Enrichment Analysis Results (fdr < 0.99)")
print()
print(gene_md)
print("\n\n")
print("Table 2: Pathway Enrichment Analysis Results (fdr < 0.99)")
print()
print(pathway_md)


# %%
