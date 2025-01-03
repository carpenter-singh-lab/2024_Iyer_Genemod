output_format <- "github_document"

render_notebook <-
  function(notebook_name, output_suffix = "", ...) {
    output_file <- paste0(notebook_name, output_suffix, ".md")

    rmarkdown::render(
      glue::glue("{notebook_name}.Rmd"),
      output_file = output_file,
      output_dir = "knit_notebooks",
      output_format = output_format,
      ...
    )

    output_file_rel <- file.path("knit_notebooks", output_file)

    read_lines(output_file_rel) %>%
      str_remove_all(file.path(getwd(), "knit_notebooks/")) %>%
      write_lines(output_file_rel)

  }


#---- 0.inspect-metadata ----

render_notebook("0.inspect-metadata")

#---- 1.inspect-splits-predictions ----

formulation_features <- read_csv("input/formulation_features_transformer.csv")

render_notebook_inspect_splits_predictions <-
  function(formulation_features_df) {
    render_notebook(
      "1.inspect-splits-predictions",
      output_suffix =
        paste(
          "",
          formulation_features_df$formulation,
          formulation_features_df$features,
          sep = "__"
        ),
      params = as.list(formulation_features_df)
    )
  }

for (i in seq(nrow(formulation_features))) {
  render_notebook_inspect_splits_predictions(formulation_features[i, ])
}

splits_reports <-
  list.files("output", pattern = "splits_report__.*.csv", full.names = T) %>%
  map(read_csv, col_types = cols()) %>%
  reduce(inner_join, by = "count")

splits_reports %>%
  write_csv("output/splits_report_combined.csv")

# - `n_gene_compound_connected_FALSE_leak_TRUE`: 0 for all configurations.
# - `n_gene_compound_connected_TRUE_leak_TRUE`: Same as above
# - `n_gene_connected_FALSE_leak_TRUE`: 0 for all `gene` (i.e. Leave-out-gene) formulations
# - `n_gene_connected_TRUE_leak_TRUE`: Same as above
# - `n_compound_connected_FALSE_leak_TRUE`: 0 for all `compound` (i.e. Leave-out-compound) and `compound-together` (i.e. Leave-out-compound-together) formuations.
# - `n_compound_connected_TRUE_leak_TRUE`: Same as above

bind_cols(
  splits_reports %>% filter(
    count %in% c(
      "n_gene_compound_connected_FALSE_leak_TRUE",
      "n_gene_compound_connected_TRUE_leak_TRUE"
    )
  ) %>% select(-count) %>% pivot_longer(everything()) %>% summarize(gene_compound_check = sum(value)),
  splits_reports %>% filter(
    count %in% c(
      "n_gene_connected_FALSE_leak_TRUE",
      "n_gene_connected_TRUE_leak_TRUE"
    )
  ) %>% select(matches("^gene_")) %>% pivot_longer(everything()) %>% summarize(gene_check = sum(value)),
  splits_reports %>% filter(
    count %in% c(
      "n_compound_connected_FALSE_leak_TRUE",
      "n_compound_connected_TRUE_leak_TRUE"
    )
  ) %>% select(matches("^compound_")) %>% pivot_longer(everything()) %>% summarize(compound_check = sum(value))
) %>%
  write_csv("output/splits_report_checks.csv")


predictions_reports <-
  list.files("output", pattern = "predictions_report__.*.csv", full.names = T) %>%
  map(read_csv, col_types = cols()) %>%
  reduce(inner_join, by = "count")

predictions_reports %>%
  write_csv("output/predictions_report_combined.csv")

# this needs https://imagemagick.org/

system("montage -label %t -tile 1x -geometry +1+1 -pointsize 60 output/gene_compound_matrix__*.png output/montage_gene_compound_matrix.png")

system("montage -label %t -tile 1x -geometry +1+1 -pointsize 60 output/gene_compound_matrix_sampled__*.png output/montage_gene_compound_matrix_sampled.png")


#---- 1.1.cosine_baseline_setup.Rmd ----

render_notebook("1.1.cosine_baseline_setup")

#---- 1.2.cosine_baseline_unsupervised.Rmd ----

render_notebook("1.2.cosine_baseline_unsupervised")

#---- 1.3.cosine_baseline_supervised.Rmd ----

render_notebook(
  "1.3.cosine_baseline_supervised",
  params = list(formulation = "compound",
                formulation_column = "pert_iname_compound"),
  output_suffix = "_compound_compound"
)

render_notebook(
  "1.3.cosine_baseline_supervised",
  params = list(formulation = "gene",
                formulation_column = "gene"),
  output_suffix = "_gene_gene"
)

render_notebook(
  "1.3.cosine_baseline_supervised",
  params = list(formulation = "pair",
                formulation_column = "pert_iname_compound"),
  output_suffix = "_pair_compound"
)

render_notebook(
  "1.3.cosine_baseline_supervised",
  params = list(formulation = "pair",
                formulation_column = "gene"),
  output_suffix = "_pair_gene"
)

#---- 1.4.cosine_baseline_supervised_pair.Rmd ----

render_notebook("1.4.cosine_baseline_supervised_pair")

#---- 1.5.crispr_consensus_transformer.Rmd ----

render_notebook("1.5.crispr_consensus_transformer")

#---- 2.score-predictions ----

# render_notebook("2.score-predictions",
#                 params = list(run_bootstrap = FALSE))

render_notebook("2.score-predictions",
                params = list(run_bootstrap = TRUE,
                              n_bootstrap = 100,
                              top_k_percent = 5,
                              n_cores = 14))

system("montage -tile 1x -geometry +1+1 -pointsize 60 output/waterfall_compound__*.png output/montage_waterfall_compound.png")

system("montage -tile 1x -geometry +1+1 -pointsize 60 output/waterfall_gene__*.png output/montage_waterfall_gene.png")

#---- 3.inspect-splits-consistency.Rmd ----

render_notebook("3.inspect-splits-consistency")

#---- 4.inspect-granular.Rmd ----

render_notebook("4.inspect-granular")

#---- 5.inspect-granular-compare.Rmd ----

render_notebook("5.inspect-granular-compare")

#---- 6.counts.Rmd ----

render_notebook("6.counts")

