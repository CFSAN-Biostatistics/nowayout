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


def get_lineages(csv: os.PathLike, cols: list) -> list:
    """
    Parse the output from `ncbitax2lin` tool and
    return a dict of lineages where the key is
    genusspeciesstrain.
    """
    lineages = dict()
    if csv == None or not (os.path.exists(csv) or os.path.getsize(csv) > 0):
        logging.error(
            f"The CSV file [{os.path.basename(csv)}] is empty or does not exist!"
        )
        exit(1)

    logging.info(f"Indexing {os.path.basename(csv)}...")

    with open(csv, "r") as csv_fh:
        header_cols = csv_fh.readline().strip().split(",")
        user_req_cols = [
            tcol_i for tcol_i, tcol in enumerate(header_cols) if tcol in cols
        ]
        cols_not_found = [tcol for tcol in cols if tcol not in header_cols]
        raw_recs = 0

        if len(cols_not_found) > 0:
            logging.error(
                f"The following columns do not exist in the"
                + f"\nCSV file [ {os.path.basename(csv)} ]:\n"
                + "".join(cols_not_found)
            )
            exit(1)
        elif len(user_req_cols) > 9:
            logging.error(
                f"Only a total of 9 columns are needed!"
                + "\ntax_id,kindom,phylum,class,order,family,genus,species,strain"
            )
            exit(1)

        for tax in csv_fh:
            raw_recs += 1
            lcols = tax.strip().split(",")

            if bool(lcols[user_req_cols[8]]):
                lineages[lcols[user_req_cols[8]]] = ",".join(
                    [lcols[l] for l in user_req_cols[1:]]
                )
            elif bool(lcols[user_req_cols[7]]):
                lineages[lcols[user_req_cols[7]]] = ",".join(
                    [lcols[l] for l in user_req_cols[1:8]] + [str()]
                )

    csv_fh.close()
    return lineages, raw_recs


def write_fasta(recs: list, basedir: os.PathLike, name: str, suffix: str) -> None:
    """
    Write sequence with no description to a specified file.
    """
    SeqIO.write(
        recs,
        os.path.join(basedir, name + suffix),
        "fasta",
    )


def check_and_get_cols(pat: re, cols: str, delim: str) -> list:
    """
    Check if header column matches the pattern and return
    columns.
    """
    if not pat.match(cols):
        logging.error(
            f"Supplied columns' names {cols} should only have words"
            f"\n(alphanumeric) separated by: {delim}."
        )
        exit(1)
    else:
        cols = re.sub("\n", "", cols).split(delim)

    return cols


def parse_tsv(fh: Union[TextIO, BinaryIO], tcols: list, delim: str) -> list:
    """
    Parse the TSV file and produce the required per
    species FASTA's.
    """
    records, sp2accs = (defaultdict(list), defaultdict(list))
    header = fh.readline().strip().split(delim)
    raw_recs = 0

    if not all(col in header for col in tcols):
        logging.error(
            "The following columns were not found in the"
            + f"\nheader row of file {os.path.basename(fh.name)}\n"
            + "\n".join([ele for ele in tcols if ele not in header])
        )

    id_i, genus_i, species_i, strain_i, seq_i = [
        i for i, ele in enumerate(header) if ele in tcols
    ]

    for record in fh:
        raw_recs += 1

        id = record.strip().split(delim)[id_i]
        genus = record.strip().split(delim)[genus_i]
        species = re.sub(r"[\/\\]+", "-", record.strip().split(delim)[species_i])
        strain = record.strip().split(delim)[strain_i]
        seq = re.sub(r"[^ATGC]+", "", record.strip().split(delim)[seq_i], re.IGNORECASE)

        if re.match(r"None|Null", species, re.IGNORECASE):
            continue

        # print(id)
        # print(genus)
        # print(species)
        # print(strain)
        # print(seq)

        records.setdefault(species, []).append(
            SeqRecord(Seq(seq), id=id, description=str())
        )
        sp2accs.setdefault(species, []).append(id)

    logging.info(f"Collected FASTA records for {len(records.keys())} species'.")
    fh.close()
    return records, sp2accs, raw_recs


# Main
def main() -> None:
    """
    This script takes:
        1. The TSV file from BOLD systems,
        2. Takes as input a .csv file generated by `ncbitax2lin`.

    and then generates a folder containing individual FASTA sequence files
    per species. This is only possible if the full taxonomy of the barcode
    sequence is present in the FASTA header.
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
        "-tsv",
        dest="tsv",
        default=False,
        required=True,
        help="Absolute UNIX path to the TSV file from BOLD systems"
        + "\nin uncompressed TXT format.",
    )
    required.add_argument(
        "-csv",
        dest="csv",
        default=False,
        required=True,
        help="Absolute UNIX path to .csv or .csv.gz file which is generated "
        + "\nby the `ncbitax2lin` tool.",
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
        "-ccols",
        dest="csv_cols",
        default="tax_id,superkingdom,phylum,class,order,family,genus,species,strain",
        required=False,
        help="Taxonomic lineage will be built using these columns from the output of"
        + "\n`ncbitax2lin`\ntool.",
    )
    parser.add_argument(
        "-ccols-sep",
        dest="csv_delim",
        default=",",
        required=False,
        help="The delimitor of the fields in the CSV file.",
    )
    parser.add_argument(
        "-tcols",
        dest="tsv_cols",
        default="processid\tgenus\tspecies\tsubspecies\tnucraw",
        required=False,
        help="For each species, the nucletide sequences will be\naggregated.",
    )
    parser.add_argument(
        "-tcols-sep",
        dest="tsv_delim",
        default="\t",
        required=False,
        help="The delimitor of the fields in the TSV file.",
    )

    # Parse defaults
    args = parser.parse_args()
    tsv = args.tsv
    csv = args.csv
    csep = args.csv_delim
    tsep = args.tsv_delim
    csv_cols = args.csv_cols
    tsv_cols = args.tsv_cols
    out = args.out_folder
    overwrite = args.force_write_out
    fna_suffix = args.fna_suffix
    ccols_pat = re.compile(f"^[\w\{csep}]+?\w$")
    tcols_pat = re.compile(f"^[\w\{tsep}]+?\w$")
    final_lineages = os.path.join(out, "lineages.csv")
    lineages_not_found = os.path.join(out, "lineages_not_found.csv")
    base_fasta_dir = os.path.join(out, "fasta")

    # Basic checks
    if not overwrite and os.path.exists(out):
        logging.warning(
            f"Output destination [{os.path.basename(out)}] already exists!"
            + "\nPlease use -f to delete and overwrite."
        )
    elif overwrite and os.path.exists(out):
        logging.info(f"Overwrite requested. Deleting {os.path.basename(out)}...")
        shutil.rmtree(out)

    # Validate user requested columns
    passed_ccols = check_and_get_cols(ccols_pat, csv_cols, csep)
    passed_tcols = check_and_get_cols(tcols_pat, tsv_cols, tsep)

    # Get  taxonomy from ncbitax2lin
    lineages, raw_recs = get_lineages(csv, passed_ccols)

    # Finally, read BOLD tsv if lineage exists.
    logging.info(f"Creating new squences per species...")

    if not os.path.exists(out):
        os.makedirs(out)

    try:
        gz_fh = gzip.open(tsv, "rt")
        records, sp2accs, traw_recs = parse_tsv(gz_fh, passed_tcols, tsep)
    except gzip.BadGzipFile:
        logging.info(f"Input TSV file {os.path.basename(tsv)} is not in\nGZIP format.")
        txt_fh = open(tsv, "r")
        records, sp2accs, traw_recs = parse_tsv(txt_fh, passed_tcols, tsep)

    passed_tax_check = 0
    failed_tax_check = 0
    fasta_recs_written = 0
    l_fh = open(final_lineages, "w")
    ln_fh = open(lineages_not_found, "w")
    l_fh.write(
        "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
    )
    ln_fh.write("fna_id,parsed_org\n")

    if not os.path.exists(base_fasta_dir):
        os.makedirs(base_fasta_dir)

    for genus_species in records.keys():
        fasta_recs_written += len(records[genus_species])
        write_fasta(
            records[genus_species],
            base_fasta_dir,
            "_".join(genus_species.split(" ")),
            fna_suffix,
        )
        org_words = genus_species.split(" ")

        for id in sp2accs[genus_species]:
            if genus_species in lineages.keys():
                this_line = ",".join([id, lineages[genus_species]]) + "\n"

                if len(org_words) > 2:
                    this_line = (
                        ",".join(
                            [id, lineages[genus_species].rstrip(","), genus_species]
                        )
                        + "\n"
                    )

                l_fh.write(this_line)
                passed_tax_check += 1
            else:
                this_line = (
                    ",".join(
                        [
                            id,
                            "",
                            "",
                            "",
                            "",
                            "",
                            org_words[0],
                            genus_species,
                            "",
                        ]
                    )
                    + "\n"
                )
                if len(org_words) > 2:
                    this_line = (
                        ",".join(
                            [
                                id,
                                "",
                                "",
                                "",
                                "",
                                "",
                                org_words[0],
                                org_words[0] + " " + org_words[1],
                                genus_species,
                            ]
                        )
                        + "\n"
                    )
                l_fh.write(this_line)
                ln_fh.write(",".join([id, genus_species]) + "\n")
                failed_tax_check += 1

    logging.info(
        f"No. of raw records present in `ncbitax2lin` [{os.path.basename(csv)}]: {raw_recs}"
        + f"\nNo. of valid records collected from `ncbitax2lin` [{os.path.basename(csv)}]: {len(lineages.keys())}"
        + f"\nNo. of raw records in TSV [{os.path.basename(tsv)}]: {traw_recs}"
        + f"\nNo. of valid records in TSV [{os.path.basename(tsv)}]: {passed_tax_check + failed_tax_check}"
        + f"\nNo. of FASTA records for which new lineages were created: {passed_tax_check}"
        + f"\nNo. of FASTA records for which only genus, species and/or strain information were created: {failed_tax_check}"
    )

    if (passed_tax_check + failed_tax_check) != fasta_recs_written:
        logging.error(
            f"The number of input FASTA records [{fasta_recs_written}]"
            + f"\nis not equal to number of lineages created [{passed_tax_check + failed_tax_check}]!"
        )
        exit(1)
    else:
        logging.info("Succesfully created lineages and FASTA records! Done!!")


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
