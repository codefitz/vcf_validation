# VCF File Validation

Version 1.0.1

Author: Wes Moskal-Fitzpatrick

This script was built based on schema defined in https://samtools.github.io/hts-specs/VCFv4.2.pdf and assistance from ChatGPT. It is not meant for direct customer use. This has not been unit tested at scale.
 
## Header line syntax

The header line names the 8 fixed, mandatory columns. These columns are as follows:
  1. #CHROM
  2. POS
  3. ID
  4. REF
  5. ALT
  6. QUAL
  7. FILTER
  8. INFO
  9. FORMAT

If genotype data is present in the file, these are followed by a FORMAT column header, then an arbitrary number of sample IDs. Duplicate sample IDs are not allowed. The header line is tab-delimited.

Test files obtained from NIH and are publicly available.

## Congenica Strict Rules

- contig = ID=1-22, ID=chr* - this fails
- Must have a FORMAT field (9 columns)
- INFO Must contain SVTYPE=CNV
- ALT must be <CNV> for CNV types
- ID must contain "LOSS" or "GAIN"
- FORMAT field must have "CN"

## Quickstart

Script accepts a VCF file or a compressed bgzipped file. You need python 3.6+ installed.

`python vcf_validation.py <*.vcf|*.gz>`