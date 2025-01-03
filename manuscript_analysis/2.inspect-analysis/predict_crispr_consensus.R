predict_crispr_consensus <- function(prediction_df, y_pos_frac_thresh) {
  prediction_df_crispr_consensus <-
    prediction_df %>%
    filter(features %in% c("crispr", "crispr_orf")) %>%
    group_by(formulation,
             features,
             Metadata_broad_sample_Compound,
             Metadata_genes2,
             y_actual) %>%
    summarise(
      y_prob_mean = mean(y_prob),
      y_prob_max = max(y_prob),
      .groups = "keep"
    ) %>%
    ungroup() %>%
    mutate(
      y_pred_mean = 1.0 * (y_prob_mean > y_pos_frac_thresh),
      y_pred_max = 1.0 * (y_prob_max > y_pos_frac_thresh)
    ) %>%
    select(
      formulation,
      features,
      Metadata_genes_CRISPR = Metadata_genes2,
      Metadata_broad_sample_Compound,
      y_actual,
      y_prob_mean,
      y_prob_max,
      y_pred_mean,
      y_pred_max
    )

  prediction_df_crispr_consensus_mean <-
    prediction_df_crispr_consensus %>%
    select(
      formulation,
      features,
      Metadata_genes_CRISPR,
      Metadata_broad_sample_Compound,
      y_actual,
      y_prob = y_prob_mean,
      y_pred = y_pred_mean
    ) %>%
    mutate(features = str_c(features, "__consensus"))

  prediction_df_crispr_consensus_max <-
    prediction_df_crispr_consensus %>%
    select(
      formulation,
      features,
      Metadata_genes_CRISPR,
      Metadata_broad_sample_Compound,
      y_actual,
      y_prob = y_prob_max,
      y_pred = y_pred_max
    ) %>%
    mutate(features = str_c(features, "__max_consensus"))

  prediction_df_crispr_consensus <-
    bind_rows(prediction_df_crispr_consensus_mean,
              prediction_df_crispr_consensus_max)

  prediction_df_crispr_consensus
}
