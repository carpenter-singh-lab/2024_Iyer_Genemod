import pandas as pd

# Read the parquet file
prediction_df_transformer_crispr_orf_max_full = pd.read_parquet(
    "../4.figures/output/prediction_df_transformer_crispr_orf_max_full.parquet"
)

# Filter for compound formulation
prediction_df_transformer_crispr_orf_max_full_compound = (
    prediction_df_transformer_crispr_orf_max_full[
        prediction_df_transformer_crispr_orf_max_full["formulation"] == "compound"
    ]
)

# Calculate threshold
formulation_compound_threshold = prediction_df_transformer_crispr_orf_max_full_compound[
    "y_actual"
].sum()

# Transform the data
prediction_df_transformer_crispr_orf_max_full_compound = (
    prediction_df_transformer_crispr_orf_max_full_compound.assign(
        y_predicted=lambda df: df["gene_compound_global_rank"]
        <= formulation_compound_threshold,
        y_actual=lambda df: df["y_actual"] == 1.0,
    )[["pert_iname_compound", "gene", "y_actual", "y_predicted"]]
)

# Sort and write to parquet
(
    prediction_df_transformer_crispr_orf_max_full_compound.sort_values(
        ["pert_iname_compound", "gene"]
    ).to_parquet(
        "output/crispr_orf_max_lo_compound_full_predictions.parquet", index=False
    )
)
