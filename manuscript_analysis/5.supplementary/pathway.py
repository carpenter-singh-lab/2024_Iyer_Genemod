#!/usr/bin/env python3
"""Maps genes to pathways using mygene.info service."""

import polars as pl
import mygene
import click
from itertools import chain
from typing import List


def process_pathway_source(hit: dict, source_key: str) -> List[dict]:
    """Processes pathway data from a specific database source.

    Args:
        hit: Single gene result from mygene query
        source_key: Dictionary key for pathway source

    Returns:
        List of pathway mappings with gene, name, id, and source
    """
    gene_symbol = hit.get("symbol", "")
    if "pathway" in hit and source_key in hit["pathway"]:
        source_pathways = hit["pathway"][source_key]
        if not isinstance(source_pathways, list):
            source_pathways = [source_pathways]

        return [
            {
                "gene": gene_symbol,
                "pathway_name": pathway.get("name", ""),
                "pathway_id": pathway.get("id", ""),
                "source": source_key.upper(),
            }
            for pathway in source_pathways
        ]
    return []


def get_pathway_info(genes: List[str]) -> pl.DataFrame:
    """Maps genes to pathways using mygene.info service.

    Args:
        genes: List of gene symbols to query

    Returns:
        DataFrame with gene-pathway mappings and sources
    """
    SCHEMA = {"gene": str, "pathway_name": str, "pathway_id": str, "source": str}

    if not genes:
        return pl.DataFrame(schema=SCHEMA)

    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        genes,
        scopes="symbol",
        fields=["symbol", "pathway.kegg", "pathway.wikipathways"],
        species="human",
        returnall=True,
    )

    sources = ("kegg", "wikipathways")
    pathway_data = list(
        chain.from_iterable(
            process_pathway_source(hit, source)
            for hit in results["out"]
            for source in sources
        )
    )

    df = pl.DataFrame(pathway_data if pathway_data else [], schema=SCHEMA)

    if results.get("missing") or results.get("notfound"):
        print(f"Warning: {len(results.get('missing', []))} genes not found")
        print(f"Warning: {len(results.get('notfound', []))} genes had no results")

    return df


@click.command()
@click.argument("genes", type=click.File("r"))
@click.option(
    "--output",
    type=click.Path(writable=True),
    default="pathways.csv",
    show_default=True,
    help="Output CSV file path",
)
def main(genes: click.File, output: str) -> None:
    """Maps genes to pathways using mygene.info.

    GENES should be a file with one gene symbol per line.

    Example:
        $ echo -e "TP53\\nBRCA1" > genes.txt
        $ python pathway.py genes.txt --output pathways.csv
    """
    gene_list = [line.strip() for line in genes if line.strip()]
    df = get_pathway_info(gene_list)
    df.write_csv(output)
    print(f"Found {len(df)} pathway mappings for {len(gene_list)} genes")


def test():
    """Simple test function for debugging."""
    genes = ["TP53", "BRCA1"]
    df = get_pathway_info(genes)
    print(df)


if __name__ == "__main__":
    # Comment/uncomment to switch between CLI and test
    # main()
    test()
