#!/usr/bin/env python3

# Kranti Konganti

import argparse
import glob
import inspect
import logging
import os
import pprint
import re
from collections import defaultdict


# Multiple inheritence for pretty printing of help text.
class MultiArgFormatClasses(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


# Main
def main() -> None:
    """
    This script will take the final taxonomic classification files and create a
    global relative abundance type file in the current working directory. The
    relative abundance type files should be in CSV or TSV format and should have
    the lineage or taxonomy in first column and samples in the subsequent columns.
    """
    # Set logging.
    logging.basicConfig(
        format="\n"
        + "=" * 55
        + "\n%(asctime)s - %(levelname)s\n"
        + "=" * 55
        + "\n%(message)s\n\n",
        level=logging.DEBUG,
    )

    # Debug print.
    ppp = pprint.PrettyPrinter(width=55)
    prog_name = inspect.stack()[0].filename

    parser = argparse.ArgumentParser(
        prog=prog_name, description=main.__doc__, formatter_class=MultiArgFormatClasses
    )

    required = parser.add_argument_group("required arguments")

    required.add_argument(
        "-abn",
        dest="rel_abn_dir",
        default=False,
        required=True,
        help="Absolute UNIX path to the parent directory that contains the\n"
        + "abundance type files.",
    )
    parser.add_argument(
        "-op",
        dest="out_prefix",
        default="nowayout.tblsum",
        required=False,
        help="Set the output file(s) prefix for output(s) generated\nby this program.",
    )
    parser.add_argument(
        "-header",
        dest="header",
        action="store_true",
        default=True,
        required=False,
        help="Do the relative abundance files have a header.",
    )
    parser.add_argument(
        "-filepat",
        dest="file_pat",
        default="*.lineage_summary.tsv",
        required=False,
        help="Files will be searched by this suffix for merged output generation\nby this program.",
    )
    parser.add_argument(
        "-failedfilepat",
        dest="failed_file_pat",
        default="*FAILED.txt",
        required=False,
        help="Files will be searched by this suffix for merged output generation\nby this program.",
    )
    parser.add_argument(
        "-delim",
        dest="delim",
        default="\t",
        required=False,
        help="The delimitor by which the fields are separated in the file.",
    )

    args = parser.parse_args()
    rel_abn_dir = args.rel_abn_dir
    is_header = args.header
    out_prefix = args.out_prefix
    file_pat = args.file_pat
    failed_file_pat = args.failed_file_pat
    delim = args.delim
    suffix = re.sub(r"^\*", "", file_pat)
    rel_abn_comb = os.path.join(os.getcwd(), out_prefix + ".txt")
    rel_abn_files = glob.glob(os.path.join(rel_abn_dir, file_pat))
    failed_rel_abn_files = glob.glob(os.path.join(rel_abn_dir, failed_file_pat))
    empty_results = "Relative abundance results did not pass thresholds"
    sample2lineage, seen_lineage = (defaultdict(defaultdict), defaultdict(int))

    if len(rel_abn_files) == 0:
        logging.info(
            "Unable to find any files with .tsv extentsion.\nNow trying .csv extension."
        )
        rel_abn_files = glob.glob(os.path.join(rel_abn_dir, "*.csv"))
        delim = ","

    if len(failed_rel_abn_files) == 0:
        logging.info(
            f"Unable to find any files with patttern {failed_file_pat}.\n"
            + "The failed samples will not appear in the final aggregate file."
        )

    if rel_abn_dir:
        if not os.path.isdir(rel_abn_dir):
            logging.error("UNIX path\n" + f"{rel_abn_dir}\n" + "does not exist!")
            exit(1)
        if len(rel_abn_files) <= 0:
            with open(rel_abn_comb, "w") as rel_abn_comb_fh:
                rel_abn_comb_fh.write(f"Sample\n{empty_results} in any samples\n")
            rel_abn_comb_fh.close()
            exit(0)

        for failed_rel_abn in failed_rel_abn_files:
            with open(failed_rel_abn, "r") as failed_fh:
                sample2lineage[failed_fh.readline().strip()].setdefault(
                    "unclassified", []
                ).append(float("1.0"))
            failed_fh.close()

        for rel_abn_file in rel_abn_files:
            sample_name = re.match(r"(^.+?)\..*$", os.path.basename(rel_abn_file))[1]

            with open(rel_abn_file, "r") as rel_abn_fh:
                if is_header:
                    sample_names = rel_abn_fh.readline().strip().split(delim)[1:]
                    if len(sample_names) > 2:
                        logging.error(
                            "The individual relative abundance file has more "
                            + "\nthan 1 sample. This is rare in the context of running the "
                            + "\n nowayout Nextflow workflow."
                        )
                        exit(1)
                    elif len(sample_names) < 2:
                        sample_name = re.sub(suffix, "", os.path.basename(rel_abn_file))
                        logging.info(
                            "Seems like there is no sample name in the lineage summary file."
                            + f"\nTherefore, sample name has been extracted from file name: {sample_name}."
                        )
                    else:
                        sample_name = sample_names[0]

                for line in rel_abn_fh.readlines():
                    cols = line.strip().split(delim)
                    lineage = cols[0]
                    abn = cols[1]
                    sample2lineage[sample_name].setdefault(lineage, []).append(
                        float(abn)
                    )
                    seen_lineage[lineage] = 1

        with open(rel_abn_comb, "w") as rel_abn_comb_fh:
            samples = sorted(sample2lineage.keys())
            rel_abn_comb_fh.write(f"Lineage{delim}" + delim.join(samples) + "\n")

            for lineage in sorted(seen_lineage.keys()):
                rel_abn_comb_fh.write(lineage)
                for sample in samples:
                    if lineage in sample2lineage[sample].keys():
                        rel_abn_comb_fh.write(
                            delim
                            + "".join(
                                [str(abn) for abn in sample2lineage[sample][lineage]]
                            )
                        )
                    else:
                        rel_abn_comb_fh.write(f"{delim}0.0")
                rel_abn_comb_fh.write("\n")
        rel_abn_comb_fh.close()


if __name__ == "__main__":
    main()
