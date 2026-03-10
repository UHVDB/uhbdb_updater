#!/usr/bin/env python

import argparse
import polars as pl

def parse_args(args=None):
    description = "Create a TSV file of unique fastas output by clusty."
    epilog = "Example usage: python unique_fastas.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-l",
        "--old_metadata",
        help="TSV file containing old metadata.",
    )
    parser.add_argument(
        "-n",
        "--new_metadata",
        help="Path to new metadata TSV file.",
    )
    parser.add_argument(
        "-f",
        "--new_fasta_tsv",
        help="Path to new FASTA TSV file.",
    )
    parser.add_argument(
        "-t",
        "--combined_clusters",
        help="Path to combined clusters TSV file.",
    )
    parser.add_argument(
        "-m",
        "--output_metadata",
        help="Path to output metadata TSV file.",
    )
    parser.add_argument(
        "-o",
        "--output_new_unique_fastas",
        help="Path to output TXT file containing new unique fastas.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    ### Read fast files TSV file (downloaded fasta files)
    fast_files_set = set(
        pl.read_csv(args.new_fasta_tsv, separator="\t", has_header=False)['column_1']
    )

    ### Read old metadata
    old_metadata_df = pl.read_csv(args.old_metadata, separator="\t")

    ### Read new metadata
    new_metadata_df = (
        pl.read_csv(args.new_metadata, separator="\t")
        .filter(pl.col('genome_id').is_in(fast_files_set))
    )

    ### Read combined clusters file
    clusters_df = (
        pl.read_csv(args.combined_clusters, separator="\t", has_header=True)
            .rename({"object": "genome_id", "cluster": "combined_cluster"})
    )

    ### Create output metadata
    (
        pl.concat([old_metadata_df, new_metadata_df], how="diagonal")
            .join(clusters_df, on="genome_id", how="left")
            .with_columns([
                pl.when(pl.col('combined_cluster').is_not_null())
                    .then(pl.col('combined_cluster'))
                    .otherwise(pl.col('cluster'))
                    .alias('cluster')
            ])
            .drop('combined_cluster')
            .sort('n50', descending=True)
            .write_csv(args.output_metadata, separator="\t", include_header=True)
    )

    ### Create output new unique fastas file
    (
        pl.read_csv(args.new_fasta_tsv, separator="\t", has_header=False, new_columns=['genome_id', 'path'])
            .filter(pl.col('genome_id').is_in(set(clusters_df['combined_cluster'])))
            [['path']]
            .write_csv(args.output_new_unique_fastas, separator="\t", include_header=False)
    )

if __name__ == "__main__":
    main()
