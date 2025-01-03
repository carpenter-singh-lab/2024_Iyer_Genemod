process_metadata <- function(df) {
  if ("Metadata_broad_sample_Compound" %in% names(df)) {
    df <-
      df %>%
      mutate(Metadata_pert_id_compound = str_sub(Metadata_broad_sample_Compound, 1, 13))
  }

  df <-
    df %>%
    select(matches("Metadata_"), everything())

  names(df) <- str_remove_all(names(df), "Metadata_")

  names(df) <- str_to_lower(names(df))

  df
}

format_df <- function(df) {
  is_not_integer_ <- function(x)
    is.numeric(x) && any(x %% 1 != 0)

  is_integer_ <- function(x)
    is.numeric(x) && all(x %% 1 == 0)

  df <-
    df %>%
    ungroup() %>%
    mutate(across(where(is_not_integer_), round, digits = 4)) %>%
    mutate(across(where(is_integer_), as.integer))

  df
}