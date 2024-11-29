# %%
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import matplotlib.pyplot as plt
import seaborn as sns


def process_compound_data(metadata_path: str) -> tuple[dict, dict]:
    """
    Process compound metadata to create gene and SMILES mappings.

    Args:
        metadata_path: Path to TSV file with columns: broad_sample, target_list, smiles

    Returns:
        tuple of (gene_to_compounds, compound_to_smiles) mappings
    """
    df = pd.read_csv(metadata_path, sep="\t")
    df = df.dropna(subset=["target_list", "smiles"])

    # Map compounds to SMILES
    compound_to_smiles = dict(zip(df["broad_sample"], df["smiles"]))

    # Map genes to their targeting compounds
    gene_to_compounds = {}
    for _, row in df.iterrows():
        for gene in row["target_list"].split("|"):
            if gene not in gene_to_compounds:
                gene_to_compounds[gene] = []
            gene_to_compounds[gene].append(row["broad_sample"])

    return gene_to_compounds, compound_to_smiles


def calculate_tanimoto_matrix(smiles_dict: dict) -> pd.DataFrame:
    """
    Calculate Tanimoto similarity matrix for compounds using Morgan fingerprints.

    Args:
        smiles_dict: Dictionary mapping compound IDs to SMILES strings

    Returns:
        DataFrame with pairwise Tanimoto similarities
    """
    # Initialize Morgan fingerprint generator (radius 2 is standard for ECFP4)
    fpgen = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
        countSimulation=False,  # Ensure binary fingerprints for Tanimoto
    )

    # Generate fingerprints
    fps = []
    valid_compounds = []

    for cid, smiles in smiles_dict.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            fps.append(fpgen.GetFingerprintAsNumPy(mol))
            valid_compounds.append(cid)

    if not fps:
        raise ValueError("No valid compounds found")

    # Calculate Tanimoto similarity
    fp_matrix = np.stack(fps)
    intersection = np.dot(fp_matrix, fp_matrix.T)
    sum_fps = np.sum(fp_matrix, axis=1)
    union = np.add.outer(sum_fps, sum_fps) - intersection
    similarity_matrix = intersection / union

    return pd.DataFrame(
        similarity_matrix, index=valid_compounds, columns=valid_compounds
    )


def analyze_gene_diversity(
    gene_to_compounds: dict, similarity_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Analyze chemical diversity of compounds targeting each gene.

    Args:
        gene_to_compounds: Dictionary mapping genes to their targeting compounds
        similarity_df: DataFrame of pairwise compound similarities

    Returns:
        DataFrame with diversity metrics for each gene
    """
    results = []

    for gene, compounds in gene_to_compounds.items():
        # Get valid compounds (those with similarity data)
        valid_compounds = [c for c in compounds if c in similarity_df.index]

        if len(valid_compounds) >= 2:
            # Extract similarity submatrix for these compounds
            submatrix = similarity_df.loc[valid_compounds, valid_compounds]

            # Calculate mean pairwise similarity (excluding self-similarity)
            n = len(valid_compounds)
            mean_similarity = (submatrix.sum().sum() - n) / (n * (n - 1))

            results.append(
                {
                    "gene": gene,
                    "mean_similarity": mean_similarity,
                    "compound_count": len(valid_compounds),
                    "diversity_score": 1
                    - mean_similarity,  # Convert similarity to diversity
                }
            )

    return pd.DataFrame(results)


def plot_diversity_distribution(results_df: pd.DataFrame) -> None:
    """
    Plot distribution of chemical diversity scores across genes.

    Args:
        results_df: DataFrame from analyze_gene_diversity
    """
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")

    # Create distribution plot
    sns.histplot(data=results_df, x="diversity_score", bins=30)

    # Add mean and median lines
    plt.axvline(
        results_df["diversity_score"].mean(),
        color="red",
        linestyle="--",
        label=f'Mean: {results_df["diversity_score"].mean():.2f}',
    )
    plt.axvline(
        results_df["diversity_score"].median(),
        color="green",
        linestyle="--",
        label=f'Median: {results_df["diversity_score"].median():.2f}',
    )

    # Customize plot
    plt.title("Distribution of Chemical Diversity Scores", fontsize=12)
    plt.xlabel("Diversity Score (0 = identical, 1 = completely different)", fontsize=10)
    plt.ylabel("Number of Genes", fontsize=10)
    plt.legend()

    # Add summary statistics
    stats = (
        f'n = {len(results_df)}\n'
        f'mean = {results_df["diversity_score"].mean():.2f}\n'
        f'std = {results_df["diversity_score"].std():.2f}'
    )
    plt.text(
        0.95,
        0.95,
        stats,
        transform=plt.gca().transAxes,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()


def main(
    compound_metadata_path: str, gene_metadata_path: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run complete analysis pipeline.

    Args:
        compound_metadata_path: Path to compound metadata TSV file
        gene_metadata_path: Path to gene metadata TSV file

    Returns:
        tuple of (similarity_matrix, diversity_results)
    """
    # Load and process data
    gene_to_compounds, compound_to_smiles = process_compound_data(
        compound_metadata_path
    )

    # Calculate similarities
    similarity_df = calculate_tanimoto_matrix(compound_to_smiles)

    # Analyze diversity
    results_df = analyze_gene_diversity(gene_to_compounds, similarity_df)

    # Create visualization
    plot_diversity_distribution(results_df)

    # Analyze diversity but only for genes present in the gene metadata
    gene_metadata = pd.read_csv(gene_metadata_path, sep="\t")

    results_df_filtered = results_df[results_df["gene"].isin(gene_metadata["gene"])]

    plot_diversity_distribution(results_df_filtered)

    return similarity_df, results_df, results_df_filtered


if __name__ == "__main__":
    compound_metadata_path = "https://raw.githubusercontent.com/jump-cellpainting/JUMP-Target/refs/heads/master/JUMP-Target-1_compound_metadata.tsv"
    gene_metadata_path = "https://raw.githubusercontent.com/jump-cellpainting/JUMP-Target/refs/heads/master/JUMP-Target-1_orf_metadata.tsv"
    similarity_df, results_df, results_df_filtered = main(
        compound_metadata_path, gene_metadata_path
    )

    # Print top 10 most diverse gene targets
    print("\nMost chemically diverse gene targets:")
    print(
        results_df.nlargest(10, "diversity_score")[
            ["gene", "diversity_score", "compound_count"]
        ]
    )

    print("\nMost chemically diverse gene targets (filtered):")
    print(
        results_df_filtered.nlargest(10, "diversity_score")[
            ["gene", "diversity_score", "compound_count"]
        ]
    )

# %%
