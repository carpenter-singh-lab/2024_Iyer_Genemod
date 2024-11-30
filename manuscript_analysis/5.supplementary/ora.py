#!/usr/bin/env python3

import polars as pl
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import click
from typing import Union, List, Optional


def perform_ora(
    sets: Union[pl.DataFrame, str],
    hits: Union[List[str], pl.DataFrame, str],
    universe: Optional[Union[List[str], pl.DataFrame, str]] = None,
    id_column: str = "item_id",
    set_column: str = "set_name",
    min_set_size: int = 2,
) -> pl.DataFrame:
    """Performs overrepresentation analysis using hypergeometric test.

    Args:
        sets: Set definitions as DataFrame or path
        hits: Items to test for enrichment
        universe: Background set (default: all items in sets)
        id_column: Column name for item IDs
        set_column: Column name for set names
        min_set_size: Minimum set size to test

    Returns:
        DataFrame with overlap stats, p-values, and FDR
    """
    # Handle input types
    if isinstance(sets, str):
        sets = pl.read_csv(sets)

    if isinstance(hits, str):
        hits = pl.read_csv(hits).item().to_list()
    elif isinstance(hits, pl.DataFrame):
        hits = hits[hits.columns[0]].to_list()

    if universe is not None:
        if isinstance(universe, str):
            universe = pl.read_csv(universe).item().to_list()
        elif isinstance(universe, pl.DataFrame):
            universe = universe[universe.columns[0]].to_list()
    else:
        universe = sets[id_column].unique().to_list()

    universe_set = set(universe)
    hit_set = set(hits)

    set_sizes = (
        sets.group_by(set_column)
        .agg(pl.count(id_column).alias("size"))
        .filter(pl.col("size") >= min_set_size)
    )

    valid_sets = set_sizes[set_column].to_list()

    M = len(universe_set)  # Total number of items in universe
    N = len(hit_set)  # Number of hits

    results = []
    for current_set in valid_sets:
        set_members = set(
            sets.filter(pl.col(set_column) == current_set)[id_column].to_list()
        )

        overlap = hit_set & set_members
        x = len(overlap)  # Number of hits in current set
        n = len(set_members)  # Size of current set

        p_value = hypergeom.sf(x - 1, M, n, N)

        results.append(
            {
                "set_name": current_set,
                "overlap_size": x,
                "set_size": n,
                "hit_size": N,
                "universe_size": M,
                "p_value": p_value,
            }
        )

    results_df = pl.DataFrame(results)

    p_values = results_df["p_value"].to_numpy()
    fdr_values = multipletests(p_values, method="fdr_bh")[1]

    return results_df.with_columns(
        [pl.Series("fdr", fdr_values).round(5), pl.col("p_value").round(5)]
    ).sort("p_value")


@click.command()
@click.option(
    "--sets",
    required=True,
    type=click.Path(exists=True),
    help="CSV file containing set mappings. Must have columns for IDs and set names",
)
@click.option(
    "--hits",
    required=True,
    type=click.Path(exists=True),
    help="CSV file containing hits (single column, any header)",
)
@click.option(
    "--universe",
    type=click.Path(exists=True),
    help="Optional CSV file containing universe items (single column, any header)",
)
@click.option(
    "--id-column",
    default="item_id",
    help="Name of column containing item IDs in sets CSV",
)
@click.option(
    "--set-column",
    default="set_name",
    help="Name of column containing set names in sets CSV",
)
@click.option(
    "--min-set-size", default=2, show_default=True, help="Minimum set size to consider"
)
@click.option(
    "--output",
    default="ora_results.csv",
    show_default=True,
    help="Output CSV file path",
)
def main(sets, hits, universe, id_column, set_column, min_set_size, output):
    """Runs ORA analysis and saves results to CSV."""
    results = perform_ora(
        sets=sets,
        hits=hits,
        universe=universe,
        id_column=id_column,
        set_column=set_column,
        min_set_size=min_set_size,
    )

    results.write_csv(output)
    click.echo(f"\nResults saved to {output}")

    significant = results.filter(pl.col("fdr") < 0.05)
    if len(significant) > 0:
        click.echo("\nSignificant results (FDR < 0.05):")
        click.echo(significant)
    else:
        click.echo("\nNo significant results found at FDR < 0.05")


def test():
    """Tests ORA with minimal example data."""
    # Create test data
    sets_data = {
        "item_id": [
            # Glycolysis items
            "HK1",
            "HK2",
            "GPI",
            "PFKL",
            "ALDOA",
            # TCA cycle items
            "CS",
            "ACO1",
            "IDH1",
            "OGDH",
            "SDHA",
        ],
        "set_name": ["Glycolysis"] * 5 + ["TCA_Cycle"] * 5,
    }
    sets_df = pl.DataFrame(sets_data)

    # Define hits (simulating upregulation of glycolysis)
    hits = ["HK2", "PFKL", "ALDOA", "IDH1"]

    # Run analysis
    results = perform_ora(sets=sets_df, hits=hits)
    print("\nResults of test ORA:")
    print(results)


if __name__ == "__main__":
    # Comment/uncomment to switch between CLI and test
    # main()
    test()
