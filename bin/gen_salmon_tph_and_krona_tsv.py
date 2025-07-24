#!/usr/bin/env python3

# Kranti Konganti
# 03/06/2024

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
    The succesful execution of this script requires access to properly formatted
    lineages.csv file which has no more than 9 columns.

    It takes the lineages.csv file, the *_hits.csv results from `sourmash gather`
    mentioned with -smres option and and a root parent directory of the
    `salmon quant` results mentioned with -sal option and generates a final
    results table with the TPM values and a .krona.tsv file for each sample
    to be used by KronaTools.
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
        "-sal",
        dest="salmon_res_dir",
        default=False,
        required=True,
        help="Absolute UNIX path to the parent directory that contains the\n"
        + "`salmon quant` results. For example, if path to\n"
        + "`quant.sf` is in /hpc/john_doe/test/salmon_res/sampleA/quant.sf, then\n"
        + "use this command-line option as:\n"
        + "-sal /hpc/john_doe/test/salmon_res",
    )
    required.add_argument(
        "-lin",
        dest="lin",
        default=False,
        required=True,
        help="Absolute UNIX Path to the lineages CSV file.\n"
        + "This file should have only 9 columns.",
    )
    required.add_argument(
        "-smres",
        dest="sm_res_dir",
        default=False,
        required=True,
        help="Absolute UNIX path to the parent directory that contains the\n"
        + "filtered `sourmas gather` results. For example, if path to\n"
        + "`sampleA.csv` is in /hpc/john_doe/test/sourmash_gather/sampleA.csv,\n"
        + "then use this command-line option as:\n"
        + "-sal /hpc/john_doe/test",
    )
    parser.add_argument(
        "-op",
        dest="out_prefix",
        default="nowayout.tblsum",
        required=False,
        help="Set the output file(s) prefix for output(s) generated\n"
        + "by this program.",
    )
    parser.add_argument(
        "-sf",
        dest="scale_down_factor",
        default=float(10000),
        required=False,
        help="Set the scaling factor by which TPM values are scaled\ndown.",
    )
    parser.add_argument(
        "-smres-suffix",
        dest="sm_res_suffix",
        default="_hits.csv",
        required=False,
        help="Find the `sourmash gather` result files ending in this\nsuffix.",
    )
    parser.add_argument(
        "-failed-suffix",
        dest="failed_suffix",
        default="_FAILED.txt",
        required=False,
        help="Find the sample names which failed classification stored\n"
        + "inside the files ending in this suffix.",
    )
    parser.add_argument(
        "-num-lin-cols",
        dest="num_lin_cols",
        default=int(9),
        required=False,
        help="Number of columns expected in the lineages CSV file.",
    )
    parser.add_argument(
        "-lin-acc-regex",
        dest="lin_acc_regex",
        default=re.compile(r"\w+[\-\.]{1}[0-9]+"),
        required=False,
        help="The pattern of the lineage's accession.",
    )

    args = parser.parse_args()
    salmon_res_dir = args.salmon_res_dir
    sm_res_dir = args.sm_res_dir
    sm_res_suffix = args.sm_res_suffix
    failed_suffix = args.failed_suffix
    out_prefix = args.out_prefix
    lin = args.lin
    num_lin_cols = args.num_lin_cols
    acc_pat = args.lin_acc_regex
    scale_down = float(args.scale_down_factor)
    no_hit = "Unclassified"
    no_hit_reads = "reads mapped to the database"
    tpm_const = float(1000000.0000000000)
    round_to = 10
    all_samples = set()
    (
        lineage2sample,
        unclassified2sample,
        lineage2sm,
        sm2passed,
        reads_total,
        per_taxon_reads,
        lineages,
    ) = (
        defaultdict(defaultdict),
        defaultdict(defaultdict),
        defaultdict(defaultdict),
        defaultdict(defaultdict),
        defaultdict(defaultdict),
        defaultdict(defaultdict),
        defaultdict(int),
    )

    salmon_comb_res = os.path.join(os.getcwd(), out_prefix + ".txt")
    # salmon_comb_res_reads_mapped = os.path.join(
    #     os.getcwd(), re.sub(".tblsum", "_reads_mapped.tblsum", out_prefix) + ".txt"
    # )
    salmon_comb_res_indiv_reads_mapped = os.path.join(
        os.getcwd(),
        re.sub(".tblsum", "_indiv_reads_mapped.tblsum", out_prefix) + ".txt",
    )
    salmon_res_files = glob.glob(
        os.path.join(salmon_res_dir, "*", "quant.sf"), recursive=True
    )
    sample_res_files_failed = glob.glob(
        os.path.join(salmon_res_dir, "*" + failed_suffix), recursive=True
    )
    sm_res_files = glob.glob(
        os.path.join(sm_res_dir, "*" + sm_res_suffix), recursive=True
    )

    # Basic checks
    if lin and not (os.path.exists(lin) and os.path.getsize(lin) > 0):
        logging.error(
            "The lineages file,\n"
            + f"{os.path.basename(lin)} does not exist or is empty!"
        )
        exit(1)

    if salmon_res_dir:
        if not os.path.isdir(salmon_res_dir):
            logging.error("UNIX path\n" + f"{salmon_res_dir}\n" + "does not exist!")
            exit(1)
        if len(salmon_res_files) <= 0:
            with open(salmon_comb_res, "w") as salmon_comb_res_fh, open(
                salmon_comb_res_indiv_reads_mapped, "w"
            ) as salmon_comb_res_indiv_reads_mapped_fh:
                salmon_comb_res_fh.write(f"Sample\n{no_hit} reads in all samples\n")
                salmon_comb_res_indiv_reads_mapped_fh.write(
                    f"Sample\nNo {no_hit_reads} from all samples\n"
                )
            salmon_comb_res_fh.close()
            salmon_comb_res_indiv_reads_mapped_fh.close()
            exit(0)

    # Only proceed if lineages.csv exists.
    if lin and os.path.exists(lin) and os.path.getsize(lin) > 0:
        lin_fh = open(lin, "r")
        _ = lin_fh.readline()

        # Index lineages.csv
        for line in lin_fh:
            cols = line.strip().split(",")

            if len(cols) < num_lin_cols:
                logging.error(
                    f"The file {os.path.basename(lin)} seems to\n"
                    + "be malformed. It contains less than required 9 columns."
                )
                exit(1)

            if cols[0] in lineages.keys():
                continue
                # logging.info(
                #     f"There is a duplicate accession [{cols[0]}]"
                #     + f" in the lineages file {os.path.basename(lin)}!"
                # )
            elif acc_pat.match(cols[0]):
                lineages[cols[0]] = ",".join(cols[1:])

        lin_fh.close()

        # Index each samples' filtered sourmash results.
        for sm_res_file in sm_res_files:
            sample_name = re.sub(sm_res_suffix, "", os.path.basename(sm_res_file))

            with open(sm_res_file, "r") as sm_res_fh:
                _ = sm_res_fh.readline()
                for line in sm_res_fh:
                    acc = acc_pat.findall(line.strip().split(",")[9])

                    if len(acc) == 0:
                        logging.info(
                            f"Got empty lineage accession: {acc}"
                            + f"\nRow elements: {line.strip().split(',')}"
                        )
                        exit(1)
                    if len(acc) not in [1]:
                        logging.info(
                            f"Got more than one lineage accession: {acc}"
                            + f"\nRow elements: {line.strip().split(',')}"
                        )
                        logging.info(f"Considering first element: {acc[0]}")
                    if acc[0] not in lineages.keys():
                        logging.error(
                            f"The lineage accession {acc[0]} is not found in {os.path.basename(lin)}"
                        )
                        exit(1)
                    lineage2sm[lineages[acc[0]]].setdefault(sample_name, 1)
                    sm2passed["sourmash_passed"].setdefault(sample_name, 1)
            sm_res_fh.close()

        # Index each samples' salmon results.
        for salmon_res_file in salmon_res_files:
            sample_name = re.match(
                r"(^.+?)((\_salmon\_res)|(\.salmon))$",
                os.path.basename(os.path.dirname(salmon_res_file)),
            )[1]
            salmon_meta_json = os.path.join(
                os.path.dirname(salmon_res_file), "aux_info", "meta_info.json"
            )

            if (
                not os.path.exists(salmon_meta_json)
                or not os.path.getsize(salmon_meta_json) > 0
            ):
                logging.error(
                    "The file\n"
                    + f"{salmon_meta_json}\ndoes not exist or is empty!\n"
                    + "Did `salmon quant` fail?"
                )
                exit(1)

            if (
                not os.path.exists(salmon_res_file)
                or not os.path.getsize(salmon_res_file) > 0
            ):
                logging.error(
                    "The file\n"
                    + f"{salmon_res_file}\ndoes not exist or is empty!\n"
                    + "Did `salmon quant` fail?"
                )
                exit(1)

            # Initiate all_tpm, rem_tpm and reads_mapped
            # all_tpm
            reads_total[sample_name].setdefault("all_tpm", []).append(float(0.0))
            # rem_tpm
            reads_total[sample_name].setdefault("rem_tpm", []).append(float(0.0))
            # reads_mapped
            reads_total[sample_name].setdefault("reads_mapped", []).append(float(0.0))

            with open(salmon_res_file, "r") as salmon_res_fh:
                for line in salmon_res_fh.readlines():
                    if re.match(r"^Name.+", line):
                        continue
                    cols = line.strip().split("\t")
                    ref_acc = cols[0]
                    tpm = cols[3]
                    num_reads_mapped = cols[4]

                    (
                        reads_total[sample_name]
                        .setdefault("all_tpm", [])
                        .append(
                            round(float(tpm), round_to),
                        )
                    )

                    (
                        reads_total[sample_name]
                        .setdefault("reads_mapped", [])
                        .append(
                            round(float(num_reads_mapped), round_to),
                        )
                    )

                    if lineages[ref_acc] in lineage2sm.keys():
                        (
                            lineage2sample[lineages[ref_acc]]
                            .setdefault(sample_name, [])
                            .append(round(float(tpm), round_to))
                        )
                        (
                            per_taxon_reads[sample_name]
                            .setdefault(lineages[ref_acc], [])
                            .append(round(float(num_reads_mapped)))
                        )
                    else:
                        (
                            reads_total[sample_name]
                            .setdefault("rem_tpm", [])
                            .append(
                                round(float(tpm), round_to),
                            )
                        )

            salmon_res_fh.close()

        # Index each samples' complete failure results i.e., 100% unclassified.
        for sample_res_file_failed in sample_res_files_failed:
            sample_name = re.sub(
                failed_suffix, "", os.path.basename(sample_res_file_failed)
            )
            with open("".join(sample_res_file_failed), "r") as no_calls_fh:
                for line in no_calls_fh.readlines():
                    if line in ["\n", "\n\r", "\r"]:
                        continue
                    unclassified2sample[sample_name].setdefault(no_hit, tpm_const)
            no_calls_fh.close()

        # Finally, write all results.
        for sample in sorted(reads_total.keys()) + sorted(unclassified2sample.keys()):
            all_samples.add(sample)

        # Check if sourmash results exist but salmon `quant` failed
        # and if so, set the sample to 100% Unclassified as well.
        for sample in sm2passed["sourmash_passed"].keys():
            if sample not in all_samples:
                unclassified2sample[sample].setdefault(no_hit, tpm_const)
                all_samples.add(sample)

        # Write total number of reads mapped to nowayout database.
        # with open(salmon_comb_res_reads_mapped, "w") as nowo_reads_mapped_fh:
        #     nowo_reads_mapped_fh.write(
        #         "\t".join(
        #             [
        #                 "Sample",
        #                 "All reads",
        #                 "Classified reads",
        #                 "Unclassified reads (Reads failed thresholds )",
        #             ]
        #         )
        #     )

        #     for sample in all_samples:
        #         if sample in reads_total.keys():
        #             nowo_reads_mapped_fh.write(
        #                 "\n"
        #                 + "\t".join(
        #                     [
        #                         f"\n{sample}",
        #                         f"{int(sum(reads_total[sample]['reads_mapped']))}",
        #                         f"{int(reads_total[sample]['reads_mapped'])}",
        #                         f"{int(reads_total[sample]['rem_tpm'])}",
        #                     ],
        #                 )
        #             )
        #         else:
        #             nowo_reads_mapped_fh.write(f"\n{sample}\t{int(0.0)}")
        # nowo_reads_mapped_fh.close()

        # Write scaled down TPM values for each sample.
        with open(salmon_comb_res, "w") as salmon_comb_res_fh, open(
            salmon_comb_res_indiv_reads_mapped, "w"
        ) as salmon_comb_res_indiv_reads_mapped_fh:
            salmon_comb_res_fh.write("Lineage\t" + "\t".join(all_samples) + "\n")
            salmon_comb_res_indiv_reads_mapped_fh.write(
                "Lineage\t" + "\t".join(all_samples) + "\n"
            )

            # Write *.krona.tsv header for all samples.
            for sample in all_samples:
                krona_fh = open(
                    os.path.join(salmon_res_dir, sample + ".krona.tsv"), "w"
                )
                krona_fh.write(
                    "\t".join(
                        [
                            "fraction",
                            "superkingdom",
                            "phylum",
                            "class",
                            "order",
                            "family",
                            "genus",
                            "species",
                        ]
                    )
                )
                krona_fh.close()

            # Write the TPM values (TPM/scale_down) for valid lineages.
            for lineage in lineage2sm.keys():
                salmon_comb_res_fh.write(lineage)
                salmon_comb_res_indiv_reads_mapped_fh.write(lineage)

                for sample in all_samples:
                    krona_fh = open(
                        os.path.join(salmon_res_dir, sample + ".krona.tsv"), "a"
                    )

                    if sample in unclassified2sample.keys():
                        salmon_comb_res_fh.write(f"\t0.0")
                        salmon_comb_res_indiv_reads_mapped_fh.write(f"\t0")
                    elif sample in lineage2sample[lineage].keys():
                        reads = sum(per_taxon_reads[sample][lineage])
                        tpm = sum(lineage2sample[lineage][sample])
                        tph = round(tpm / scale_down, round_to)
                        lineage2sample[sample].setdefault("hits_tpm", []).append(
                            float(tpm)
                        )

                        salmon_comb_res_fh.write(f"\t{tph}")
                        salmon_comb_res_indiv_reads_mapped_fh.write(f"\t{reads}")
                        krona_lin_row = lineage.split(",")

                        if len(krona_lin_row) > num_lin_cols - 1:
                            logging.error(
                                "Taxonomy columns are more than 8 for the following lineage:"
                                + f"{krona_lin_row}"
                            )
                            exit(1)
                        else:
                            krona_fh.write(
                                "\n"
                                + str(round((tpm / tpm_const), round_to))
                                + "\t"
                                + "\t".join(krona_lin_row[:-1])
                            )
                    else:
                        salmon_comb_res_fh.write(f"\t0.0")
                        salmon_comb_res_indiv_reads_mapped_fh.write(f"\t0")
                    krona_fh.close()

                salmon_comb_res_fh.write("\n")
                salmon_comb_res_indiv_reads_mapped_fh.write(f"\n")

            # Finally write TPH (TPM/scale_down) for Unclassified
            # Row = Unclassified / No reads mapped to the database ...
            salmon_comb_res_fh.write(f"{no_hit}")
            salmon_comb_res_indiv_reads_mapped_fh.write(f"Total {no_hit_reads}")

            for sample in all_samples:
                krona_ufh = open(
                    os.path.join(salmon_res_dir, sample + ".krona.tsv"), "a"
                )
                # krona_ufh.write("\t")
                if sample in unclassified2sample.keys():
                    salmon_comb_res_fh.write(
                        f"\t{round((unclassified2sample[sample][no_hit] / scale_down), round_to)}"
                    )
                    salmon_comb_res_indiv_reads_mapped_fh.write(f"\t0")
                    krona_ufh.write(
                        f"\n{round((unclassified2sample[sample][no_hit] / tpm_const), round_to)}"
                    )
                else:
                    trace_tpm = tpm_const - sum(reads_total[sample]["all_tpm"])
                    trace_tpm = float(f"{trace_tpm:.{round_to}f}")
                    if trace_tpm <= 0:
                        trace_tpm = float(0.0)
                    tph_unclassified = float(
                        f"{(sum(reads_total[sample]['rem_tpm']) + trace_tpm) / scale_down:{round_to}f}"
                    )
                    krona_unclassified = float(
                        f"{(sum(reads_total[sample]['rem_tpm']) + trace_tpm) / tpm_const:{round_to}f}"
                    )
                    salmon_comb_res_fh.write(f"\t{tph_unclassified}")
                    salmon_comb_res_indiv_reads_mapped_fh.write(
                        f"\t{int(sum(sum(per_taxon_reads[sample].values(), [])))}"
                    )
                    krona_ufh.write(f"\n{krona_unclassified}")
                krona_ufh.write("\t" + "\t".join(["unclassified"] * (num_lin_cols - 2)))
                krona_ufh.close()

        salmon_comb_res_fh.close()
        salmon_comb_res_indiv_reads_mapped_fh.close()
        # ppp.pprint(lineage2sample)
        # ppp.pprint(lineage2sm)
        # ppp.pprint(reads_total)


if __name__ == "__main__":
    main()
