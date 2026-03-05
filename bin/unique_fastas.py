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
        help="Path to TSV file input into sketchlib.",
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

    ### Read input TSV file
    input_df = pl.read_csv(args.input, separator="\t", has_header=False, new_columns=["genome_id", "path"])

    ### Identify reference fasta
    ref_id = set(
        pl.read_csv(args.clusters, separator="\t", has_header=True)
            .head(1)
            ['cluster']
    )

    ### Write out reference fasta path
    (
        input_df
            .filter(pl.col("genome_id").is_in(ref_id))
            [['path']]
            .write_csv(args.output + '.ref.txt', separator="\t", include_header=False)
    )

    ### Identify cluster reps (unique fastas)
    unique_genome_ids = set(
        pl.read_csv(args.clusters, separator="\t", has_header=True)
            .filter(~pl.col('cluster').is_in(set(ref_id)))
            ['cluster']
    )

    ### Filter input TSV to retain only unique fastas
    unique_df = input_df.filter(pl.col("genome_id").is_in(unique_genome_ids))

    ### Write output TSV file
    unique_df[['path']].write_csv(args.output + '.unique.txt', separator="\t", include_header=False)

if __name__ == "__main__":
    main()
