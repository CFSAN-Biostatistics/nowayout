#!/usr/bin/env python3

# Kranti Konganti

import argparse
import gzip
import inspect
import logging
import os
import pprint
import re

# Set logging.
logging.basicConfig(
    format="\n" + "=" * 55 + "\n%(asctime)s - %(levelname)s\n" + "=" * 55 + "\n%(message)s\n\n",
    level=logging.DEBUG,
)

# Debug print.
ppp = pprint.PrettyPrinter(width=50, indent=4)

# Multiple inheritence for pretty printing of help text.
class MultiArgFormatClasses(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def write_failures(prefix: str, file: os.PathLike) -> None:
    with open(file, "w") as outfile_failed_fh:
        outfile_failed_fh.write(f"{prefix}\n")
    outfile_failed_fh.close()


def main() -> None:
    """
    This script will take the CSV output of `sourmash search` and `sourmash gather`
    and will return a column's value filtered by requested column name and its value
    """

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
        help="Absolute UNIX path to CSV file containing output from\n"
        + "`sourmash gather` or `sourmash search`",
    )
    required.add_argument(
        "-extract",
        dest="extract",
        required=False,
        default="name",
        help="Extract this column's value which matches the filters.\n"
        + "Controlled by -fcn and -fcv.",
    )
    parser.add_argument(
        "-all",
        dest="alllines",
        required=False,
        default=False,
        action="store_true",
        help="Instead of just the column value, print entire row.",
    )
    parser.add_argument(
        "-fcn",
        dest="filter_col_name",
        default="f_match",
        required=False,
        help="Column name by which the filtering of rows\nshould be applied.",
    )
    parser.add_argument(
        "-fcv",
        dest="filter_col_val",
        default="0",
        required=False,
        help="Only rows where the column (defined by --fcn)\nsatisfies this value will be\n"
        + "will be considered. This can be numeric, regex\nor a string value.",
    )
    parser.add_argument(
        "-gt",
        dest="gt",
        default=True,
        required=False,
        action="store_true",
        help="Apply greater than or equal to condition on\nnumeric values of --fcn column.",
    )
    parser.add_argument(
        "-lt",
        dest="lt",
        default=False,
        required=False,
        action="store_true",
        help="Apply less than or equal to condition on\nnumeric values of --fcn column.",
    )

    args = parser.parse_args()
    csv = args.csv
    ex = args.extract
    all_lines = args.alllines
    fcn = args.filter_col_name
    fcv = args.filter_col_val
    gt = args.gt
    lt = args.lt
    hits = set()
    hit_lines = set()
    empty_lines = 0

    outfile_prefix = re.sub(r"(^.*?)\.csv\.gz", r"\1", os.path.basename(csv))
    outfile_failed = os.path.join(os.getcwd(), "_".join([outfile_prefix, "FAILED.txt"]))

    if csv and (not os.path.exists(csv) or not os.path.getsize(csv) > 0):
        logging.error(
            "The CSV file,\n" + f"{os.path.basename(csv)} does not exists or\nis of size zero."
        )
        write_failures(outfile_prefix, outfile_failed)
        exit(0)

    if all_lines:
        outfile = os.path.join(os.getcwd(), "_".join([outfile_prefix, "hits.csv"]))
    else:
        outfile = os.path.join(os.getcwd(), "_".join([outfile_prefix, "template_hits.txt"]))

    with gzip.open(csv, "rb") as csv_fh:
        header_cols = dict(
            [
                (col, ele)
                for ele, col in enumerate(csv_fh.readline().decode("utf-8").strip().split(","))
            ]
        )

        if fcn and ex not in header_cols.keys():
            logging.info(
                f"The header row in file\n{os.path.basename(csv)}\n"
                + "does not have a column whose names are:\n"
                + f"-fcn: {fcn} and -extract: {ex}"
            )
            exit(1)

        for line in csv_fh:
            line = line.decode("utf-8")

            if line in ["\n", "\n\r"]:
                empty_lines += 1
                continue

            cols = [x.strip() for x in line.strip().split(",")]
            investigate = float(format(float(cols[header_cols[fcn]]), '.10f'))
            fcv = float(fcv)

            if re.match(r"[\d\.]+", str(investigate)):
                if gt and investigate >= fcv:
                    hits.add(cols[header_cols[ex]])
                    hit_lines.add(line.strip())
                elif lt and investigate <= fcv:
                    hits.add(cols[header_cols[ex]])
                    hit_lines.add(line.strip())
            elif investigate == fcv:
                hits.add(cols[header_cols[ex]])
                hit_lines.add(line.strip())

        csv_fh.close()

        if len(hits) >= 1:
            with open(outfile, "w") as outfile_fh:
                outfile_fh.write(",".join(header_cols.keys()) + "\n")
                if all_lines:
                    outfile_fh.write("\n".join(hit_lines) + "\n")
                else:
                    outfile_fh.writelines("\n".join(hits) + "\n")
            outfile_fh.close()
        else:
            write_failures(outfile_prefix, outfile_failed)

        if empty_lines > 0:
            empty_lines_msg = f"Skipped {empty_lines} empty line(s).\n"

            logging.info(
                empty_lines_msg
                + f"File {os.path.basename(csv)}\n"
                + f"written in:\n{os.getcwd()}\nDone! Bye!"
            )
        exit(0)


if __name__ == "__main__":
    main()
