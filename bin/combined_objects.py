#!/usr/bin/env python

import argparse
import polars as pl

def parse_args(args=None):
    description = "Create a TSV file of unique fastas output by clusty."
    epilog = "Example usage: python unique_fastas.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-n",
        "--new_metadata",
        help="Path to new metadata TSV file.",
    )
    parser.add_argument(
        "-l",
        "--old_metadata",
        help="TSV file containing old metadata.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to output TXT file.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    ### Read old metadata
    old_metadata_df = (
        pl.read_csv(args.old_metadata, separator="\t")
            .filter(pl.col('genome_id') == pl.col('cluster'))
            .select(['genome_id', 'n50'])
    )

    ### Read new metadata
    new_metadata_df = pl.read_csv(args.new_metadata, separator="\t")[['genome_id', 'n50']]

    ### Combine old and new unique metadata
    combined_metadata_df = pl.concat([old_metadata_df, new_metadata_df], how="vertical")

    ### Write combined objects file sorted by n50
    (
        combined_metadata_df
            .sort('n50', descending=True)
            [['genome_id']]
            .write_csv(args.output, separator="\t", include_header=True)
    )


if __name__ == "__main__":
    main()
