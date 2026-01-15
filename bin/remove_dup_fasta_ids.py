#!/usr/bin/env python3

import argparse
import gzip
import inspect
import logging
import os
import pprint
import shutil
from typing import BinaryIO, TextIO, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genericpath import isdir


# Multiple inheritence for pretty printing of help text.
class MultiArgFormatClasses(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


def write_fasta(seq: str, id: str, fh: Union[TextIO, BinaryIO]) -> None:
    """
    Write sequence with no description to specified file.
    """
    SeqIO.write(
        SeqRecord(Seq(seq), id=id, description=str()),
        fh,
        "fasta",
    )


# Main
def main() -> None:
    """
    This script takes:
        1. A FASTA file in gzip or non-gzip (ASCII TXT) format and

    and then generates a new FASTA file with duplicate FASTA IDs replaced
    with a unique ID.
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
    prog_name = os.path.basename(inspect.stack()[0].filename)

    parser = argparse.ArgumentParser(
        prog=prog_name, description=main.__doc__, formatter_class=MultiArgFormatClasses
    )

    required = parser.add_argument_group("required arguments")

    required.add_argument(
        "-fna",
        dest="fna",
        default=False,
        required=True,
        help="Absolute UNIX path to .fna or .fna.gz file.",
    )
    parser.add_argument(
        "-lin",
        dest="lineages",
        default=False,
        required=False,
        help="Absolute UNIX path to lineages.csv file for which the"
        + "\nthe duplicate IDs will be made unique corresponding to"
        + "\nthe FASTA IDs",
    )
    parser.add_argument(
        "-outdir",
        dest="out_folder",
        default=os.getcwd(),
        required=False,
        help="By default, the output is written to this\nfolder.",
    )
    parser.add_argument(
        "-f",
        dest="force_write_out",
        default=False,
        action="store_true",
        required=False,
        help="Force overwrite the output file.",
    )
    parser.add_argument(
        "--fna-suffix",
        dest="fna_suffix",
        default=".fna",
        required=False,
        help="Suffix of the output FASTA file.",
    )

    # Parse defaults
    args = parser.parse_args()
    fna = args.fna
    lineages = args.lineages
    outdir = args.out_folder
    overwrite = args.force_write_out
    fna_suffix = args.fna_suffix
    new_fna = os.path.join(
        outdir, os.path.basename(fna).split(".")[0] + "_dedup_ids" + fna_suffix
    )
    lin_header = False
    new_lin = False
    seen_ids = dict()
    seen_lineages = dict()

    # Basic checks
    if not overwrite and os.path.exists(new_fna):
        logging.warning(
            f"Output destination [{os.path.basename(new_fna)}] already exists!"
            + "\nPlease use -f to delete and overwrite."
        )
    elif overwrite and os.path.exists(new_fna):
        logging.info(f"Overwrite requested. Deleting {os.path.basename(new_fna)}...")
        if os.path.isdir(new_fna):
            shutil.rmtree(new_fna)
        else:
            os.remove(new_fna)

    # Prepare for writing
    new_fna_fh = open(new_fna, "+at")

    # If lineages file is mentioned, index it.
    if lineages and os.path.exists(lineages) and os.path.getsize(lineages) > 0:
        new_lin = os.path.join(os.getcwd(), os.path.basename(lineages) + "_dedup.csv")
        new_lin_fh = open(new_lin, "w")
        with open(lineages, "r") as l_fh:
            lin_header = l_fh.readline()
            for line in l_fh:
                cols = line.strip().split(",")
                if len(cols) < 9:
                    logging.error(
                        f"The row in the lineages file {os.path.basename(lineages)}"
                        + f"\ndoes not have 9 required columns: {len(cols)}"
                        + f"\n\n{lin_header.strip()}\n{line.strip()}"
                    )
                    exit(1)
                elif len(cols) > 9:
                    logging.info(
                        f"The row in the lineages file {os.path.basename(lineages)}"
                        + f"\nhas more than 9 required columns: {len(cols)}"
                        + f"\nRetaining only 9 columns of the following 10 columns."
                        + f"\n\n{lin_header.strip()}\n{line.strip()}"
                    )

                if cols[0] not in seen_lineages.keys():
                    seen_lineages[cols[0]] = ",".join(cols[1:9])

        new_lin_fh.write(lin_header)
        l_fh.close()

    # Read FASTA and create unique FASTA IDs.
    logging.info(f"Creating new FASTA with unique IDs.")
    try:
        fna_fh = gzip.open(fna, "rt")
        _ = fna_fh.readline()
    except gzip.BadGzipFile:
        logging.info(
            f"Input FASTA file {os.path.basename(fna)} is not in\nGZIP format."
            + "\nAttempting text parsing."
        )
        fna_fh = open(fna, "r")

    for record in SeqIO.parse(fna_fh, format="fasta"):
        seq_id = record.id

        if record.id not in seen_ids.keys():
            seen_ids[record.id] = 1
        else:
            seen_ids[record.id] += 1

        if seen_ids[seq_id] > 1:
            seq_id = str(record.id) + str(seen_ids[record.id])

        if new_lin:
            new_lin_fh.write(",".join([seq_id, seen_lineages[record.id]]) + "\n")

        write_fasta(record.seq, seq_id, new_fna_fh)

    if new_lin:
        new_lin_fh.close()

    logging.info("Done!")


if __name__ == "__main__":

    main()
