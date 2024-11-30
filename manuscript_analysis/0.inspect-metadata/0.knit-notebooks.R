
notebook_name <- "0.update-metadata"
output_format <- "github_document"
output_file <- file.path("knit_notebooks", paste0(notebook_name, ".md"))

rmarkdown::render(
  glue::glue("{notebook_name}.Rmd"),
  output_file = output_file,
  output_dir = "knit_notebooks",
  output_format = output_format
)

