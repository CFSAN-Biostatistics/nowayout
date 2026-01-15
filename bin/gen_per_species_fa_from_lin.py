#!/usr/bin/env python3

import argparse
import gzip
import inspect
import logging
import os
import pprint
import re
import shutil
from collections import defaultdict
from typing import BinaryIO, TextIO, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Multiple inheritence for pretty printing of help text.
class MultiArgFormatClasses(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


def get_lineages(csv: os.PathLike) -> defaultdict:
    """
    Parse the lineages.csv file and store a list of
    accessions.
    """
    lineages = dict()
    if csv == None or not (os.path.exists(csv) or os.path.getsize(csv) > 0):
        logging.error(
            f"The CSV file [{os.path.basename(csv)}] is empty or does not exist!"
        )
        exit(1)

    logging.info(f"Indexing {os.path.basename(csv)}...")

    with open(csv, "r") as csv_fh:
        _ = csv_fh.readline().strip().split(",")
        for line in csv_fh:
            cols = line.strip().split(",")

            if len(cols) < 9:
                logging.error(
                    f"The CSV file {os.path.basename(csv)} should have a mandatory 9 columns."
                    + "\n\nEx: identifiers,superkingdom,phylum,class,order,family,genus,species,strain"
                    + "\nAB211151.1,Eukaryota,Arthropoda,Malacostraca,Decapoda,Majidae,Chionoecetes,Chionoecetes opilio,"
                    + f"\n\nGot:\n{line}"
                )
                exit(1)

            lineages[cols[0]] = re.sub(r"\W+", "-", "_".join(cols[7].split(" ")))

    csv_fh.close()
    return lineages


def get_unique_dir(file_num: int, k=3) -> str:
    """
    Return a unique directory name to manage a
    hundred's of thousands of files.
    """
    dir_name = list()
    for i in range(k - 1, -1, -1):
        letter = chr(ord("A") + (file_num // (26**i)) % 26)  # Start from ASCII 65 (A).
        dir_name.append(letter)
    return "".join(dir_name)


def write_fasta(
    recs: list, basedir: os.PathLike, name: str, suffix: str, name_re: re
) -> None:
    """
    Write sequence with no description to a specified file.
    """
    sanitized_name = name_re.sub("", name)

    if not os.path.exists(basedir):
        print(f"Writing into {os.path.basename(basedir)}")
        os.makedirs(basedir)

    SeqIO.write(
        recs,
        os.path.join(basedir, sanitized_name + suffix),
        "fasta",
    )


def parse_fasta(fh: Union[TextIO, BinaryIO], sp2accs: dict) -> list:
    """
    Parse the sequences and create per species FASTA record.
    """
    records = defaultdict()

    for record in SeqIO.parse(fh, "fasta"):

        id = record.id
        seq = record.seq

        if id in sp2accs.keys():
            records.setdefault(sp2accs[id], []).append(
                SeqRecord(Seq(seq), id=id, description=str())
            )
        else:
            print(f"Lineage row does not exist for accession: {id}")

    logging.info(f"Collected FASTA records for {len(records.keys())} species'.")
    fh.close()
    return records


# Main
def main() -> None:
    """
    This script takes:
        1. The FASTA file and,
        2. Takes the corresponding lineages.csv file and,

    then generates a folder containing individual FASTA sequence files
    per species.
    """

    # Set logging.
    logging.basicConfig(
        format="\n"
        + "=" * 55
        + "\n%(asctime)s - %(levelname)s\n"
        + "=" * 55
        + "\n%(message)s\r\r",
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
        "-fa",
        dest="fna",
        default=False,
        required=True,
        help="Absolute UNIX path to the FASTA file that corresponds"
        + "\nto the lineages.csv file.",
    )
    required.add_argument(
        "-csv",
        dest="csv",
        default=False,
        required=True,
        help="Absolute UNIX path to lineages.csv which has a guaranteed 9 "
        + "\ncolumns with the first being an accession.",
    )
    parser.add_argument(
        "-out",
        dest="out_folder",
        default=os.path.join(os.getcwd(), "species"),
        required=False,
        help="By default, the output is written to this\nfolder.",
    )
    parser.add_argument(
        "-f",
        dest="force_write_out",
        default=False,
        action="store_true",
        required=False,
        help="Force overwrite output directory contents.",
    )
    parser.add_argument(
        "-suffix",
        dest="fna_suffix",
        default=".fna",
        required=False,
        help="Suffix of the individual species FASTA files\nthat will be saved.",
    )
    parser.add_argument(
        "-dlen",
        dest="dir_name_len",
        default=int(3),
        required=False,
        help="Name of the unique directory\nthat will be generated to manage FASTA files.",
    )
    parser.add_argument(
        "-nfiles",
        dest="num_files_per_dir",
        default=int(1000),
        required=False,
        help="Number of FASTA files per unique directory.",
    )

    # Parse defaults
    args = parser.parse_args()
    csv = args.csv
    fna = args.fna
    out = args.out_folder
    overwrite = args.force_write_out
    fna_suffix = args.fna_suffix
    dir_name_len = args.dir_name_len
    num_files_per_dir = args.num_files_per_dir
    name_re = re.compile(r"^\W+")

    # Basic checks
    if not overwrite and os.path.exists(out):
        logging.warning(
            f"Output destination [{os.path.basename(out)}] already exists!"
            + "\nPlease use -f to delete and overwrite."
        )
    elif overwrite and os.path.exists(out):
        logging.info(f"Overwrite requested. Deleting {os.path.basename(out)}...")
        shutil.rmtree(out)

    # Get  taxonomy from ncbitax2lin
    lineages = get_lineages(csv)

    logging.info(f"Creating new squences per species...")

    if not os.path.exists(out):
        os.makedirs(out)

    try:
        gz_fh = gzip.open(fna, "rt")
        fa_recs = parse_fasta(gz_fh, lineages)
    except gzip.BadGzipFile:
        logging.info(
            f"Input FASTA file {os.path.basename(fna)} is not in\nGZIP format."
        )
        txt_fh = open(fna, "r")
        fa_recs = parse_fasta(txt_fh, lineages)
    finally:
        logging.info("Assigned FASTA records per species...")

    logging.info("Writing FASTA records per species...")

    species_list = list(fa_recs.keys())
    for i in range(0, len(species_list)):
        sp = species_list[i]
        unique_out_dir = os.path.join(
            out, get_unique_dir(i // num_files_per_dir, k=dir_name_len)
        )
        write_fasta(fa_recs[sp], unique_out_dir, sp, fna_suffix, name_re)


if __name__ == "__main__":
    main()

# ~/apps/nowayout/bin/gen_per_species_fa_from_bold.py -tsv BOLD_Public.05-Feb-2024.tsv -csv ../tax.csv                                  ─╯

# =======================================================
# 2024-02-08 21:37:28,541 - INFO
# =======================================================
# Indexing tax.csv...

# =======================================================
# 2024-02-08 21:38:06,567 - INFO
# =======================================================
# Creating new squences per species...

# =======================================================
# 2024-02-08 21:38:06,572 - INFO
# =======================================================
# Input TSV file BOLD_Public.05-Feb-2024.tsv is not in
# GZIP format.

# =======================================================
# 2024-02-08 22:01:04,554 - INFO
# =======================================================
# Collected FASTA records for 497421 species'.

# =======================================================
# 2024-02-08 22:24:35,000 - INFO
# =======================================================
# No. of raw records present in `ncbitax2lin` [tax.csv]: 2550767
# No. of valid records collected from `ncbitax2lin` [tax.csv]: 2134980
# No. of raw records in TSV [BOLD_Public.05-Feb-2024.tsv]: 9735210
# No. of valid records in TSV [BOLD_Public.05-Feb-2024.tsv]: 4988323
# No. of FASTA records for which new lineages were created: 4069202
# No. of FASTA records for which only genus, species and/or strain information were created: 919121

# =======================================================
# 2024-02-08 22:24:35,001 - INFO
# =======================================================
# Succesfully created lineages and FASTA records! Done!!
