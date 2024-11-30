#!/usr/bin/env python3
"""Process compound prediction data and merge with repurposing annotations.

Output CSV columns:
- pert_iname_compound: Compound ID in genemod project
- gene: Target gene for prediction (distinct from annotated targets)
- gene_compound_global_rank: Prediction rank for compound-gene pair
- disease_area, target, moa, indication, clinical_phase: Clinical annotations from https://repo-hub.broadinstitute.org/
- pert_iname_normalized: Normalized compound name for matching
"""

import pandas as pd
from pathlib import Path
import re


def normalize_string(text: str) -> str:
    """Convert string to lowercase and remove hyphens."""
    return re.sub("-", "", str(text).lower())


def main():
    script_dir = Path(__file__).parent

    # Load metrics data
    metrics = pd.read_parquet(
        script_dir.parent
        / "4.figures/output/prediction_df_transformer_crispr_orf_max.parquet"
    )

    # Filter for compound formulation
    metrics = metrics[metrics["formulation"] == "compound"]

    # Load repurposing data
    rephub = pd.read_csv(
        script_dir / "input/repurposing_drugs_20200324.txt", sep="\t", comment="!"
    )

    # Add normalized columns
    metrics["pert_iname_normalized"] = metrics["pert_iname_compound"].apply(
        normalize_string
    )
    rephub["pert_iname_normalized"] = rephub["pert_iname"].apply(normalize_string)

    # Merge datasets
    metrics_augmented = pd.merge(
        metrics, rephub, on="pert_iname_normalized", how="inner"
    )

    # Select and process columns
    output_columns = [
        "pert_iname_compound",
        "gene_compound_global_rank",
        "moa",
        "gene",
        "disease_area",
        "target",
        "indication",
        "pert_iname_normalized",
        "clinical_phase",
    ]

    # Process well-predicted compounds
    well_predicted = metrics_augmented[output_columns].sort_values(
        "gene_compound_global_rank", ascending=True
    )
    # Convert gene_compound_global_rank to integer
    well_predicted["gene_compound_global_rank"] = well_predicted[
        "gene_compound_global_rank"
    ].astype(int)

    # Save results
    Path("output").mkdir(exist_ok=True)
    well_predicted.to_csv("output/well_predicted_compounds.csv", index=False)


if __name__ == "__main__":
    main()
