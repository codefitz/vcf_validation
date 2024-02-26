# VCFv4.2 Validator
#
# Author: Wes Moskal-Fitzpatrick
#
# This script was built based on schema defined in https://samtools.github.io/hts-specs/VCFv4.2.pdf and assistance
# from ChatGPT. It is not meant for direct customer use.
# 
# Change History
# 1.0.0 : WMF : Created.
# 1.0.1 : WMF : Updated with Congenica rules.
# 1.0.2 : WMF : Added support for bgzipped files. Updated error message for Alternate Alleles.
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

def main():
    if len(sys.argv) != 2:
        print("Usage: python vcf_validation.py <*.vcf|*.gz>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    validate_vcf(vcf_file)

def validate_vcf(vcf_file):
    if vcf_file.endswith('.vcf'):
        open_func = open
    elif vcf_file.endswith('.gz'):
        open_func = gzip.open
    else:
        print("File type not recognised.")
        print("Usage: python vcf_validation.py <*.vcf|*.gz>")
        sys.exit(1)
    
    with open_func(vcf_file, 'rt') as file:
        ...

        line_number = 0
        header_found = False
        for line in file:
            line_number += 1
            if line.startswith("##"):
                if line.startswith("##contig"):
                    contig_info = line.split('<',1)[1].split('>')[0]
                    id_info = [x for x in contig_info.split(',') if x.startswith('ID=')]
                    if id_info:
                        contig_id = id_info[0].split('=')[1]
                        if contig_id.startswith("chr"):
                            print(f"Error: Contig ID starts with 'chr' on line {line_number}: {line.strip()}")
                            sys.exit(1)
                            
                if not header_found:
                    header_found = True
                continue
            elif line.startswith("#CHROM"):
                header_found = False
            else:
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    print(f"Error: Incorrect number of columns on line {line_number}: {line.strip()}")
                    sys.exit(1)
                
                # CHROM - Chromosome. An identifier from the reference genome or an angle-bracketed ID String (“<ID>”).
                # String, no whitespace permitted, Required.
                if not re.match("^[0-9A-Za-z_]+$", fields[0]):
                    print(f"Error: Invalid chromosome on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # POS - Position. The reference position, with the 1st base having position 1. Positions are sorted numerically, in increasing order,
                # within each reference sequence CHROM. Integer, Required.
                if not re.match("^[0-9]+$", fields[1]):
                    print(f"Error: Invalid position on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # ID - Identifier. Semicolon-separated list of unique identifiers where available.
                # String, no whitespace or semicolons permitted. Missing values denoted by ".".
                if not re.match("^([A-Za-z0-9:_.]+(;[A-Za-z0-9_.]+)*)?$", fields[2]):
                    print(f"Error: Invalid ID on line {line_number}: {line.strip()}")
                    sys.exit(1)
                # ID field must contain in the string somewhere "LOSS" or "GAIN"
                if "LOSS" not in fields[2] and "GAIN" not in fields[2]:
                    print(f"Error: ID field doesn't contain 'LOSS' or 'GAIN' on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # REF - Reference base(s). Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
                # String, Required.
                if not re.match("^[ACGTN]+$", fields[3]):
                    print(f"Error: Invalid reference allele on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # ALT - Altnerate base(s). Comma-separated list of alternate non-reference alleles.  Strings made up of the bases A,C,G,T,N,*, (case insensitive)
                # or an angle-bracketed ID String (“<ID>”) or a breakend replacement string. String; no whitespace, commas, or angle-brackets are permitted
                # in the ID String itself. Missing values denoted by ".".
                if fields[4] != "<CNV>":
                #if not re.match("^([ACGTN]+|<[^>]+>)(,([ACGTN]+|<[^>]+>))*$", fields[4]):
                    print(f"Error: Invalid alternate allele on line {line_number}: {line.strip()}")
                    print(f"ALT must be \<CNV\> for copy number variants.")
                    sys.exit(1)

                # QUAL - Quality. Phred-scaled quality score for the assertion made in ALT. Numeric, missing values are denoted by ".".
                if not re.match("^[0-9]+(\.[0-9]+)?$", fields[5]) and fields[5] != ".":
                    print(f"Error: Invalid quality on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # FILTER - Filter status. PASS if this position has passed all filters, i.e., a call is made at this position. Otherwise, if the site has not
                # passed all filters, a semicolon-separated list of codes for filters that fail. String, no whitespace or semicolons permitted. Missing
                # values denoted by ".".
                if not re.match("^([A-Za-z0-9_]+(;[A-Za-z0-9_]+)*)?$|^\\.$", fields[6]):
                    print(f"Error: Invalid filter on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # INFO - Additional information. encoded as a semicolon-separated series of short keys with optional values in the format: <key>=<data>[,data].
                # String, no whitespace, semicolons, or equals-signs permitted; commas are permitted only as delimiters for lists of values). Missing values
                # are denoted by ".".
                # INFO field (field 7) must contain SVTYPE=CNV
                if "SVTYPE=CNV" not in fields[7]:
                    print(f"Error: Missing SVTYPE=CNV in INFO field on line {line_number}: {line.strip()}")
                    sys.exit(1)

                # FORMAT field (field 8) must contain "CN"
                if "CN" not in fields[8]:
                    print(f"Error: Missing 'CN' in FORMAT field on line {line_number}: {line.strip()}")
                    sys.exit(1)

    print("VCF file validation completed. No structural errors found.")

if __name__ == "__main__":
    main()
