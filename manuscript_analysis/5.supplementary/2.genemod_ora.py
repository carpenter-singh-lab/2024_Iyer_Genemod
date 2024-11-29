import polars as pl
from pathlib import Path
from ora import perform_ora
from pathway import get_pathway_info

if __name__ == "__main__":
    script_dir = Path(__file__).parent
    output_dir = script_dir / "output"
    well_predicted_compounds = pl.read_csv(output_dir / "well_predicted_compounds.csv")

    # Compound-Gene analysis
    mappings = well_predicted_compounds.select(["pert_iname_compound", "gene"]).sort(
        "pert_iname_compound", "gene"
    )
    compound_hits = (
        well_predicted_compounds.filter(
            pl.col("gene_compound_global_rank") <= len(well_predicted_compounds)
        )
        .select("pert_iname_compound")
        .unique()
    )

    results = perform_ora(
        sets=mappings,
        hits=compound_hits,
        id_column="pert_iname_compound",
        set_column="gene",
        min_set_size=2,
    )
    print("\nResults of Compound-Gene ORA analysis:", results)
    pl.DataFrame(results).write_csv(output_dir / "ora_compound_gene.csv")

    # Gene-Pathway analysis
    gene_hits = (
        well_predicted_compounds.filter(
            pl.col("gene_compound_global_rank") <= len(well_predicted_compounds)
        )
        .select("gene")
        .unique()
    )

    pathway_df = get_pathway_info(mappings["gene"].unique().to_list())

    # Run ORA for both KEGG and WIKIPATHWAYS
    for source in ["KEGG", "WIKIPATHWAYS"]:
        pathway_mappings = pathway_df.filter(pl.col("source") == source).select(
            ["gene", "pathway_name"]
        )

        results = perform_ora(
            sets=pathway_mappings,
            hits=gene_hits,
            id_column="gene",
            set_column="pathway_name",
            min_set_size=2,
        )
        print(f"\nResults of {source} Pathway ORA analysis:", results)
        pl.DataFrame(results).write_csv(
            output_dir / f"ora_gene_pathway_{source.lower()}.csv"
        )
