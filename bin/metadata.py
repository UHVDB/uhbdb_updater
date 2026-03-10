#!/usr/bin/env python

import argparse
import polars as pl

def parse_args(args=None):
    description = "Create a TSV file of unique fastas output by clusty."
    epilog = "Example usage: python unique_fastas.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to metadata TSV file.",
    )
    parser.add_argument(
        "-f",
        "--fast_files_tsv",
        help="TSV file containing paths to fastas.",
    )
    parser.add_argument(
        "-c",
        "--clusters",
        help="TSV file containing cluster information.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output TSV file.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    ### Read fast files TSV file (downloaded fasta files)
    fast_files_set = set(
        pl.read_csv(args.fast_files_tsv, separator="\t", has_header=False)['column_1']
    )

    ### Read input metadata TSV file
    input_df = (
        pl.read_csv(args.input, separator="\t")
        .filter(
            pl.col('genome_id').is_in(fast_files_set)
        )
    )

    ### Identify clusters
    clusters_df = (
        pl.read_csv(args.clusters, separator="\t", has_header=True)
            .rename({"object": "genome_id"})
    )

    ### Join metadata with clusters
    joined_df = input_df.join(clusters_df, on="genome_id", how="inner", coalesce=True)
    
    ### Write output TSV file
    (
        joined_df
            .sort('n50', descending=True)
            .write_csv(args.output, separator="\t", include_header=True)
    )

if __name__ == "__main__":
    main()
