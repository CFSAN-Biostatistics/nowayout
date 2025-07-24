#!/usr/bin/env python3

import argparse
import gzip
import inspect
import logging
import os
import pprint
import re
import shutil
import ssl
import tempfile
from html.parser import HTMLParser
from urllib.request import urlopen

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Multiple inheritence for pretty printing of help text.
class MultiArgFormatClasses(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


# HTMLParser override class to get fna.gz and gbff.gz
class NCBIHTMLParser(HTMLParser):
    def __init__(self, *, convert_charrefs: bool = ...) -> None:
        super().__init__(convert_charrefs=convert_charrefs)
        self.reset()
        self.href_data = list()

    def handle_data(self, data):
        self.href_data.append(data)


# Download organelle FASTA and GenBank file.
def dl_mito_seqs_and_flat_files(url: str, suffix: re, out: os.PathLike) -> os.PathLike:
    """
    Method to save .fna.gz and .gbff.gz files for the
    RefSeq mitochondrion release.
    """
    contxt = ssl.create_default_context()
    contxt.check_hostname = False
    contxt.verify_mode = ssl.CERT_NONE

    if url == None:
        logging.error(
            "Please provide the base URL where .fna.gz and .gbff.gz"
            + "\nfiles for RefSeq mitochondrion can be found."
        )
        exit(1)

    if os.path.exists(out):
        for file in os.listdir(out):
            file_path = os.path.join(out, file)

            if suffix.match(file_path) and os.path.getsize(file_path) > 0:
                logging.info(
                    f"The required mitochondrion file(s)\n[{os.path.basename(file_path)}]"
                    + " already exists.\nSkipping download from NCBI..."
                    + "\nPlease use -f to delete and overwrite."
                )
                return file_path
    else:
        os.makedirs(out)

    html_parser = NCBIHTMLParser()
    logging.info(f"Finding latest NCBI RefSeq mitochondrion release at:\n{url}")

    with urlopen(url, context=contxt) as response:
        with tempfile.NamedTemporaryFile(delete=False) as tmp_html_file:
            shutil.copyfileobj(response, tmp_html_file)

    with open(tmp_html_file.name, "r") as html:
        html_parser.feed("".join(html.readlines()))

    file = suffix.search("".join(html_parser.href_data)).group(0)
    file_url = "/".join([url, file + ".gz"])
    file_at = os.path.join(out, file)

    logging.info(f"Found NCBI RefSeq mitochondrian file(s):\n{file_url}")

    logging.info(f"Saving to:\n{file_at}")

    with tempfile.NamedTemporaryFile(delete=False) as tmp_gz:
        with urlopen(file_url, context=contxt) as response:
            tmp_gz.write(response.read())

    with open(file_at, "w") as fh:
        with gzip.open(tmp_gz.name, "rb") as web_gz:
            fh.write(web_gz.read().decode("utf-8"))

    html.close()
    tmp_gz.close()
    tmp_html_file.close()
    os.unlink(tmp_gz.name)
    os.unlink(tmp_html_file.name)
    fh.close()
    web_gz.close()
    response.close()

    return file_at


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
                    [re.sub(r"[\,\"]", "", lcols[l]) for l in user_req_cols[1:]]
                )
            elif bool(lcols[user_req_cols[7]]):
                lineages[lcols[user_req_cols[7]]] = ",".join(
                    [re.sub(r"[\,\"]", "", lcols[l]) for l in user_req_cols[1:8]]
                    + [str()]
                )

            if lcols[7] == "Rondeletia bicolor Goode & Bean":
                print(lineages[lcols[8]])
                exit(0)

    csv_fh.close()
    return lineages, raw_recs


def from_genbank(gbk: os.PathLike, min_len: int) -> dict:
    """
    Method to parse GenBank file and return
    organism to latest accession mapping.
    """
    accs2orgs = dict()
    sanitize_pat = re.compile(r"[\,\"]")

    if not (os.path.exists(gbk) or os.path.getsize(gbk) > 0):
        logging.info(
            f"The GenBank file [{os.path.basename(gbk)}] does not exist"
            + "\nor is of size 0."
        )
        exit(1)

    logging.info(f"Indexing {os.path.basename(gbk)}...")

    # a = open("./_accs", "w")
    try:
        for record in SeqIO.parse(gbk, "genbank"):
            if len(record.seq) < min_len:
                continue
            else:
                # a.write(f"{record.id}\n")
                accs2orgs[record.id] = sanitize_pat.sub(
                    "", record.annotations["organism"]
                )
    except Exception as e:
        logging.error(f"Error occured around GenBank ID: {record.id}")
        logging.error(f"Error: {e}")
        exit(1)

    return accs2orgs


def from_genbank_alt(gbk: os.PathLike) -> dict:
    """
    Method to parse GenBank file and return
    organism to latest accession mapping without
    using BioPython's GenBank Scanner
    """
    accs2orgs = dict()
    accs = dict()
    orgs = dict()
    acc = False
    acc_pat = re.compile(r"^VERSION\s+(.+)")
    org_pat = re.compile(r"^\s+ORGANISM\s+(.+)")
    sanitize_pat = re.compile(r"[\,\"]")

    if not (os.path.exists(gbk) or os.path.getsize(gbk) > 0):
        logging.info(
            f"The GenBank file [{os.path.basename(gbk)}] does not exist"
            + "\nor is of size 0."
        )
        exit(1)

    logging.info(
        f"Indexing {os.path.basename(gbk)} without using\nBioPython's GenBank Scanner..."
    )

    with open(gbk, "r") as gbk_fh:
        for line in gbk_fh:
            line = line.rstrip()
            if line.startswith("VERSION") and acc_pat.match(line):
                acc = acc_pat.match(line).group(1)
                accs[acc] = 1
            if org_pat.match(line):
                if acc and acc not in orgs.keys():
                    orgs[acc] = sanitize_pat.sub("", org_pat.match(line).group(1))
                elif acc and acc in orgs.keys():
                    logging.error(f"Duplicate VERSION line: {acc}")
                    exit(1)
        if len(accs.keys()) != len(orgs.keys()):
            logging.error(
                f"Got unequal number of organisms ({len(orgs.keys())})\n"
                + f"and accessions ({len(accs.keys())})"
            )
            exit(1)
        else:
            for acc in accs.keys():
                if acc not in orgs.keys():
                    logging.error(f"ORAGANISM not found for accession: {acc}")
                    exit(1)
                accs2orgs[acc] = orgs[acc]

    gbk_fh.close()
    return accs2orgs


def write_fasta(seq: str, id: str, basedir: os.PathLike, suffix: str) -> None:
    """
    Write sequence with no description to specified file.
    """
    SeqIO.write(
        SeqRecord(Seq(seq), id=id, description=str()),
        os.path.join(basedir, id + suffix),
        "fasta",
    )


# Main
def main() -> None:
    """
    This script takes:
        1. Downloads the RefSeq Mitochrondrial GenBank and FASTA format files.
        2. Takes as input and output .csv.gz or .csv file generated by `ncbitax2lin`.

    and then generates a folder containing individual FASTA sequence files
    per organelle, and a corresponding lineage file in CSV format.
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
        "-csv",
        dest="csv",
        default=False,
        required=True,
        help="Absolute UNIX path to .csv or .csv.gz file which is generated "
        + "\nby the `ncbitax2lin` tool.",
    )
    parser.add_argument(
        "-cols",
        dest="lineage_cols",
        default="tax_id,domain,phylum,class,order,family,genus,species,strain",
        required=False,
        help="Taxonomic lineage will be built using these columns from the output of"
        + "\n`ncbitax2lin` tool.",
    )
    parser.add_argument(
        "-url",
        dest="url",
        default="https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion",
        required=False,
        help="Base URL from where NCBI RefSeq mitochondrion files will be downloaded\nfrom.",
    )
    parser.add_argument(
        "-out",
        dest="out_folder",
        default=os.path.join(os.getcwd(), "organelles"),
        required=False,
        help="By default, the output is written to this folder.",
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
        "--fna-suffix",
        dest="fna_suffix",
        default=".fna",
        required=False,
        help="Suffix of the individual organelle FASTA files that will be saved.",
    )
    parser.add_argument(
        "-ml",
        dest="fa_min_len",
        default=200,
        required=False,
        help="Minimum length of the FASTA sequence for it to be considered for"
        + "\nfurther processing",
    )
    parser.add_argument(
        "--gen-per-fa",
        dest="gen_per_fa",
        default=False,
        required=False,
        action="store_true",
        help="Generate per sequence FASTA file.",
    )
    parser.add_argument(
        "--alt-gb-parser",
        dest="alt_gb_parser",
        default=False,
        required=False,
        action="store_true",
        help="Use alternate GenBank parser instead of BioPython's.",
    )

    # Parse defaults
    args = parser.parse_args()
    csv = args.csv
    out = args.out_folder
    overwrite = args.force_write_out
    fna_suffix = args.fna_suffix
    url = args.url
    tax_cols = args.lineage_cols
    gen_per_fa = args.gen_per_fa
    alt_gb_parser = args.alt_gb_parser
    min_len = int(args.fa_min_len)
    tcols_pat = re.compile(r"^[\w\,]+?\w$")
    mito_fna_suffix = re.compile(r".*?\.genomic\.fna")
    mito_gbff_suffix = re.compile(r".*?\.genomic\.gbff")
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

    if not tcols_pat.match(tax_cols):
        logging.error(
            f"Supplied columns' names {tax_cols} should only have words (alphanumeric) separated by a comma."
        )
        exit(1)
    else:
        tax_cols = re.sub("\n", "", tax_cols).split(",")

    # Get .fna and .gbk files
    fna = dl_mito_seqs_and_flat_files(url, mito_fna_suffix, out)
    gbk = dl_mito_seqs_and_flat_files(url, mito_gbff_suffix, out)

    # Get  taxonomy from ncbitax2lin
    lineages, raw_recs = get_lineages(csv, tax_cols)

    # Get parsed organisms and latest accession from GenBank file.
    if alt_gb_parser:
        accs2orgs = from_genbank_alt(gbk)
    else:
        accs2orgs = from_genbank(gbk, min_len)

    # # Finally, read FASTA and create individual FASTA if lineage exists.
    logging.info(f"Creating new sequences and lineages...")

    l_fh = open(final_lineages, "w")
    ln_fh = open(lineages_not_found, "w")
    l_fh.write(
        "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
    )
    ln_fh.write("fna_id,gbk_org\n")
    passed_lookup = 0
    failed_lookup = 0
    gbk_recs_missing = 0
    skipped_len_short = 0

    if gen_per_fa and not os.path.exists(base_fasta_dir):
        os.makedirs(base_fasta_dir)

    for record in SeqIO.parse(fna, "fasta"):
        if len(record.seq) < min_len:
            skipped_len_short += 1
            continue
        elif record.id in accs2orgs.keys():
            org_words = accs2orgs[record.id].split(" ")
        else:
            gbk_recs_missing += 1
            continue

        genus_species = (
            " ".join(org_words[0:2]) if len(org_words) > 2 else " ".join(org_words[0:])
        )

        if gen_per_fa:
            write_fasta(record.seq, record.id, base_fasta_dir, fna_suffix)

        if record.id in accs2orgs.keys() and accs2orgs[record.id] in lineages.keys():
            l_fh.write(",".join([record.id, lineages[accs2orgs[record.id]]]) + "\n")
            passed_lookup += 1
        elif record.id in accs2orgs.keys() and genus_species in lineages.keys():
            if len(org_words) > 2:
                l_fh.write(
                    ",".join(
                        [
                            record.id,
                            lineages[genus_species].rstrip(","),
                            accs2orgs[record.id],
                        ]
                    )
                    + "\n"
                )
            else:
                l_fh.write(",".join([record.id, lineages[genus_species]]) + "\n")
            passed_lookup += 1
        else:
            if len(org_words) > 2:
                l_fh.write(
                    ",".join(
                        [
                            record.id,
                            "",
                            "",
                            "",
                            "",
                            "",
                            org_words[0],
                            org_words[0] + " " + org_words[1],
                            accs2orgs[record.id],
                        ]
                    )
                    + "\n"
                )
            else:
                l_fh.write(
                    ",".join(
                        [
                            record.id,
                            "",
                            "",
                            "",
                            "",
                            "",
                            org_words[0],
                            accs2orgs[record.id],
                            "",
                        ]
                    )
                    + "\n"
                )
            ln_fh.write(",".join([record.id, accs2orgs[record.id]]) + "\n")
            failed_lookup += 1

    logging.info(
        f"No. of raw records present in `ncbitax2lin` [{os.path.basename(csv)}]: {raw_recs}"
        + f"\nNo. of valid records collected from `ncbitax2lin` [{os.path.basename(csv)}]: {len(lineages.keys())}"
        + f"\nNo. of sequences skipped (Sequence length < {min_len}): {skipped_len_short}"
        + f"\nNo. of records in FASTA [{os.path.basename(fna)}]: {passed_lookup + failed_lookup}"
        + f"\nNo. of records in GenBank [{os.path.basename(gbk)}]: {len(accs2orgs.keys())}"
        + f"\nNo. of FASTA records for which new lineages were created: {passed_lookup}"
        + f"\nNo. of FASTA records for which only genus, species and/or strain information were created: {failed_lookup}"
        + f"\nNo. of FASTA records for which no GenBank records exist: {gbk_recs_missing}"
    )

    if (passed_lookup + failed_lookup) != len(accs2orgs.keys()):
        logging.error(
            f"The number of FASTA records written [{len(accs2orgs.keys())}]"
            + f"\nis not equal to number of lineages created [{passed_lookup + failed_lookup}]!"
        )
        exit(1)
    else:
        logging.info("Succesfully created lineages and FASTA records! Done!!")

    l_fh.close()
    ln_fh.close()


if __name__ == "__main__":
    main()
