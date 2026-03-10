#!/usr/bin/env python

import argparse
import polars as pl

def parse_args(args=None):
    description = "Split input samplesheet by genus."
    epilog = "Example usage: python split_by_genus.py --help"

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input TSV file.",
    )
    parser.add_argument(
        "-l",
        "--large_threshold",
        type=int,
        help="Threshold for large genera.",
    )
    parser.add_argument('--version', action='version', version='1.0.0')
    return parser.parse_args(args)



def main(args=None):
    args = parse_args(args)

    ### Read input TSV file
    input_df = pl.read_csv(args.input, separator="\t")

    ### Filter input file
    # completeness >= 90
    # contamination <= 5
    # gunc_css <= 0.45
    filtered_df = input_df.filter(
        (pl.col("completeness") >= 90) &
        (pl.col("contamination") <= 5) &
        (pl.col("gunc_css") <= 0.45)
    )

    ### Sort by n50
    sorted_df = filtered_df.sort("n50", descending=True)

    ### Group by genus
    genera_lst = sorted_df["genus"].unique().to_list()
    for genus in genera_lst:
        group = sorted_df.filter(pl.col("genus") == genus)
        if group.height >= args.large_threshold:
            prefix = "large_genera/"
        else:
            prefix = "small_genera/"

        urls_path = f"{prefix}{genus}_urls.txt"
        local_paths = f"{prefix}{genus}_local_fastas.txt"
        object_file = f"{prefix}{genus}_objects.txt"
        meta_file = f"{prefix}{genus}_metadata.tsv"

        # write out aria2c file
        (
            group
                .filter((pl.col("url").is_not_null()) & (pl.col('path').is_null()))
                .with_columns([
                    # get file suffix from url (e.g. .fna.gz)
                    pl.col('url').str.extract(r'(\.fna\.gz|\.fa\.gz|\.fasta\.gz|\.fna|\.fa|.fasta)$').alias('suffix'),
                ])
                .with_columns([
                    # create new column with
                    (pl.col('url') + "\n out=" + pl.col('source_db') + "_" + pl.col('genome') + pl.col('suffix')).alias('aria2c_entry')
                ])
                [['aria2c_entry']]
                .write_csv(urls_path, include_header=False)
        )

        # write out local symlink file
        (
            group
                .filter(pl.col("path").is_not_null())
                .with_columns([
                    # get file suffix from url (e.g. .fna.gz)
                    pl.col('path').str.extract(r'(\.fna\.gz|\.fa\.gz|\.fasta\.gz|\.fna|\.fa|.fasta)$').alias('suffix'),
                ])
                .with_columns([
                    # create new column with
                    (pl.col('path') + " ./local_fastas/" + pl.col('source_db') + "_" + pl.col('genome') + pl.col('suffix')).alias('symlink_entry')
                ])
                [['symlink_entry']]
                .write_csv(local_paths, include_header=False)
        )
        
        # write out objects file (just genome names sorted by n50)
        (
            group
                .with_columns([
                    # create new column with
                    (pl.col('source_db') + "_" + pl.col('genome')).alias('object')
                ])
                [['object']]
                .write_csv(object_file, separator="\t")
        )
    
        # write out filtered metadata file
        (
            group
                .with_columns([
                    # create new column with
                    (pl.col('source_db') + "_" + pl.col('genome')).alias('genome_id')
                ])
                .drop(['path'])
                .write_csv(meta_file, separator="\t")
        )

if __name__ == "__main__":
    main()
