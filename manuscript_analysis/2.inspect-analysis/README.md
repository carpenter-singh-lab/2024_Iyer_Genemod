# Inspect analysis

- [knit_notebooks/0.inspect-metadata.md](knit_notebooks/0.inspect-metadata.md): Inspect metadata used in predictive models: Here, we compare the connections used for predictive modeling (“current connections”), with the expanded list of connections (“updated connections”) generated in [../0.inspect-metadata/knit_notebooks/0.update-metadata.md](../0.inspect-metadata/knit_notebooks/0.update-metadata.md) ([../0.inspect-metadata/output/JUMP-Target_compounds_crispr_orf_connections.csv](../0.inspect-metadata/output/JUMP-Target_compounds_crispr_orf_connections.csv))
- [knit_notebooks/1.inspect-splits-predictions__compound__crispr.md](knit_notebooks/1.inspect-splits-predictions__compound__crispr.md): (Each configuration has a corresponding file). For the `compound__crispr` configuration, perform various checks on the data splits and verify that the predictions are made on the test set. The notebook produces:
  - [output/splits_report__compound__crispr.csv](output/splits_report__compound__crispr.csv): a `splits_report` which is later collated as [output/splits_report_combined.csv](output/splits_report_combined.csv).
    - `n_gene_compound_connected_FALSE_leak_TRUE`: 0 for all configurations.
    - `n_gene_compound_connected_TRUE_leak_TRUE`: Same as above
    - `n_gene_connected_FALSE_leak_TRUE`: 0 for all `gene` (i.e. Leave-out-gene) formulations
    - `n_gene_connected_TRUE_leak_TRUE`: Same as above
    - `n_compound_connected_FALSE_leak_TRUE`: 0 for all `compound` (i.e. Leave-out-compound) formulations.
    - `n_compound_connected_TRUE_leak_TRUE`: Same as above
    - The output of all these checks, across configurations, is at [output/splits_report_checks.csv](output/splits_report_checks.csv) (all values should be 0).
  - [output/gene_compound_matrix__compound__crispr.png](output/gene_compound_matrix__compound__crispr.png): A gene x compound matrix indicating `TRUE`/ `FALSE` connections, and `TRAIN` / `TEST` connection.
  - [output/gene_compound_matrix__compound__crispr.tsv](output/gene_compound_matrix__compound__crispr.tsv): A text file with the same information, in a format that can be loaded in https://software.broadinstitute.org/morpheus/
  - [output/gene_compound_matrix_sampled__compound__crispr.png](output/gene_compound_matrix_sampled__compound__crispr.png): A randomly downsampled version of the same as above.
  - [output/splits__compound__crispr.csv.gz](output/splits__compound__crispr.csv.gz): This combined the information from [../1.analysis/splits/compound__crispr__test.csv](../1.analysis/splits/compound__crispr__test.csv), [../1.analysis/splits/compound__crispr__train.csv](../1.analysis/splits/compound__crispr__train.csv), and [../1.analysis/splits/compound__crispr__val.csv](../1.analysis/splits/compound__crispr__val.csv), and annotates with gene name and compound name, rather the `broad_sample` identifier.
  - [output/predictions_report__compound__crispr.csv](output/predictions_report__compound__crispr.csv): a `predictions_report` which is later collated as [output/predictions_report_combined.csv](output/predictions_report_combined.csv); this reports whether the predictions are made on the test set.


## Configurations

In all configurations we use, we ensure that
- the same gene-compound pair is never present in both train and test
- the splits across train/test/validate is 60-20-20 – keeping the classes (positive and negative) as balanced as possible

All 4 x 3 = 12 combinations of Formulations-Features below are meaningful

### Formulations

- Standard (`pair`): A gene or compound can be present in both train and test/validation.
- Leave-out-compound (`compound`): A compound is present in either in test/validation or train, but not both
    - Biological insight gained: Can I predict the target if I’ve never seen this compound before?
- Leave-out-gene (`gene`): A gene is present in either in test/validation or train, but not both
    - Biological insight: Can I predict the target if I’ve never seen this gene target before?


### Features

- CRISPR + ORF features (`crispr_orf`)
    - Biological insight: How well can we predict if both CRISPR and ORF profiles are available
    - Only one CRISPR profile is used at a time (see below for other scenarios)
- ORF features (`orf`)
    - Biological insight: How well can we predict if only ORF profiles are available
- CRISPR features (`crispr`)
    - Biological insight: How well can we predict if only CRISPR profiles are available
    - Only one CRISPR profile is used at a time (see below for other scenarios)

We make two additional sets of predictions by combining CRISPR guides targeting the same gene
- `_max` (or `_max__consensus`): take the maximum of the probabilities of the two guides
- `_avg` (or just `__consensus`): take the average of the probabilities of the two guides
