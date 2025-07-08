"""VCFv4.2 Validator.

This module implements a minimal validator for VCF files used at Congenica.  It
parses a ``.vcf`` or bgzipped ``.vcf.gz`` file and applies a set of structural
checks based on the VCF specification as well as some Congenica specific rules:

* Contig identifiers must not start with ``"chr"``.
* Records must contain a ``FORMAT`` column.
* ``INFO`` must include ``SVTYPE=CNV``.
* ``ALT`` must be ``<CNV>`` for copy number variants.
* ``ID`` values must contain ``LOSS`` or ``GAIN``.
* ``FORMAT`` column must contain ``CN``.

The script exits with ``sys.exit(1)`` on the first validation failure and prints
an error describing the offending line.
"""

# Change History
# 1.0.0 : WMF : Created.
# 1.0.1 : WMF : Updated with Congenica rules.
# 1.0.2 : WMF : Added support for bgzipped files. Updated error message for Alternate Alleles.
# 1.0.3 : WMF : Added header validations and duplicate sample detection.
#
# Header line syntax
# ------------------
# The header line names the 8 fixed, mandatory columns. These columns are as follows:
#   1. #CHROM
#   2. POS
#   3. ID
#   4. REF
#   5. ALT
#   6. QUAL
#   7. FILTER
#   8. INFO
# If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number
# of sample IDs. Duplicate sample IDs are not allowed. The header line is tab-delimited.
#
# Congenica Strict Rules
# contig = ID=1-22, ID=chr* - this fails
# Must have a FORMAT field (9 columns)
# INFO Must contain SVTYPE=CNV
# ALT must be <CNV> for CNV types
# ID must contain "LOSS" or "GAIN"
# FORMAT field must have "CN"

import sys
import re
import gzip
import argparse


def validate_chrom(value: str, line_number: int, line: str) -> None:
    """Validate the CHROM field."""
    if not re.match(r"^[0-9A-Za-z_]+$", value):
        print(f"Error: Invalid chromosome on line {line_number}: {line.strip()}")
        sys.exit(1)


def validate_pos(value: str, line_number: int, line: str) -> None:
    """Validate the POS field."""
    if not re.match(r"^[0-9]+$", value):
        print(f"Error: Invalid position on line {line_number}: {line.strip()}")
        sys.exit(1)


def validate_id(value: str, line_number: int, line: str) -> None:
    """Validate the ID field."""
    if not re.match(r"^([A-Za-z0-9:_.]+(;[A-Za-z0-9_.]+)*)?$", value):
        print(f"Error: Invalid ID on line {line_number}: {line.strip()}")
        sys.exit(1)
    if "LOSS" not in value and "GAIN" not in value:
        print(
            f"Error: ID field doesn't contain 'LOSS' or 'GAIN' on line {line_number}: {line.strip()}"
        )
        sys.exit(1)


def validate_ref(value: str, line_number: int, line: str) -> None:
    """Validate the REF field."""
    if not re.match(r"^[ACGTN]+$", value):
        print(f"Error: Invalid reference allele on line {line_number}: {line.strip()}")
        sys.exit(1)


def validate_alt(value: str, line_number: int, line: str) -> None:
    """Validate the ALT field."""
    if value != "<CNV>":
        print(f"Error: Invalid alternate allele on line {line_number}: {line.strip()}")
        print("ALT must be <CNV> for copy number variants.")
        sys.exit(1)


def validate_qual(value: str, line_number: int, line: str) -> None:
    """Validate the QUAL field."""
    if not re.match(r"^[0-9]+(\.[0-9]+)?$", value) and value != ".":
        print(f"Error: Invalid quality on line {line_number}: {line.strip()}")
        sys.exit(1)


def validate_filter(value: str, line_number: int, line: str) -> None:
    """Validate the FILTER field."""
    if not re.match(r"^([A-Za-z0-9_]+(;[A-Za-z0-9_]+)*)?$|^\.$", value):
        print(f"Error: Invalid filter on line {line_number}: {line.strip()}")
        sys.exit(1)


def validate_info(value: str, line_number: int, line: str) -> None:
    """Validate the INFO field."""
    if "SVTYPE=CNV" not in value:
        print(f"Error: Missing SVTYPE=CNV in INFO field on line {line_number}: {line.strip()}")
        sys.exit(1)


def validate_format(value: str, line_number: int, line: str) -> None:
    """Validate the FORMAT field."""
    if "CN" not in value:
        print(f"Error: Missing 'CN' in FORMAT field on line {line_number}: {line.strip()}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Validate a VCF file")
    parser.add_argument("vcf_file", help="Path to VCF or bgzipped file")
    parser.add_argument("--strict", action="store_true",
                        help="Enable Congenica strict rule checks")
    parser.add_argument("--report", action="store_true",
                        help="Print a summary report when validation completes")

    args = parser.parse_args()

    validate_vcf(args.vcf_file, strict=args.strict, report=args.report)

def validate_vcf(vcf_file, strict=False, report=False):
    if vcf_file.endswith('.vcf'):
        open_func = open
    elif vcf_file.endswith('.gz'):
        open_func = gzip.open
    else:
        print("File type not recognised.")
        print("Usage: python vcf_validation.py <*.vcf|*.gz>")
        sys.exit(1)
    
    with open_func(vcf_file, 'rt') as file:

        line_number = 0
        fileformat_found = False
        header_found = False
        for line in file:
            line_number += 1
            if line.startswith("##"):
                if line.startswith("##fileformat"):
                    fileformat_found = True
                if line.startswith("##contig"):
                    contig_info = line.split('<',1)[1].split('>')[0]
                    id_info = [x for x in contig_info.split(',') if x.startswith('ID=')]
                    if id_info:
                        contig_id = id_info[0].split('=')[1]
                        if contig_id.startswith("chr"):
                            print(f"Error: Contig ID starts with 'chr' on line {line_number}: {line.strip()}")
                            sys.exit(1)
                continue
            elif line.startswith("#CHROM"):
                header_found = True
                header_fields = line.strip().split('\t')
                if len(header_fields) > 8 and header_fields[8] != "FORMAT":
                    print("Error: FORMAT column missing from header line")
                    sys.exit(1)
                if len(header_fields) > 9:
                    sample_names = header_fields[9:]
                    if len(sample_names) != len(set(sample_names)):
                        print("Error: Duplicate sample names in header line")
                        sys.exit(1)
                continue
            else:
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    print(f"Error: Incorrect number of columns on line {line_number}: {line.strip()}")
                    sys.exit(1)

                validate_chrom(fields[0], line_number, line)
                validate_pos(fields[1], line_number, line)
                validate_id(fields[2], line_number, line)
                validate_ref(fields[3], line_number, line)
                validate_alt(fields[4], line_number, line)
                validate_qual(fields[5], line_number, line)
                validate_filter(fields[6], line_number, line)
                validate_info(fields[7], line_number, line)
                validate_format(fields[8], line_number, line)

        if not fileformat_found:
            print("Error: Missing ##fileformat header")
            sys.exit(1)
        if not header_found:
            print("Error: Missing #CHROM header line")
            sys.exit(1)

    if report:
        print("VCF file validation completed. No structural errors found.")

if __name__ == "__main__":
    main()
