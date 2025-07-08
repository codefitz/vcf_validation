import tempfile
import os
import unittest
from vcf_validation import validate_vcf

class TestVCFValidation(unittest.TestCase):
    def _write_temp_vcf(self, content):
        tmp = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.vcf')
        tmp.write(content)
        tmp.flush()
        tmp.close()
        self.addCleanup(lambda: os.remove(tmp.name))
        return tmp.name

    def test_valid_file(self):
        content = """##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t1\tGAIN1\tA\t<CNV>\t30\tPASS\tSVTYPE=CNV\tCN:GT\t0/1\t1/1
"""
        path = self._write_temp_vcf(content)
        # Should not raise SystemExit
        validate_vcf(path)

    def test_missing_fileformat(self):
        content = """##contig=<ID=1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t1\tGAIN1\tA\t<CNV>\t30\tPASS\tSVTYPE=CNV\tCN:GT\t0/1
"""
        path = self._write_temp_vcf(content)
        with self.assertRaises(SystemExit):
            validate_vcf(path)

    def test_missing_chrom_header(self):
        content = """##fileformat=VCFv4.2
##contig=<ID=1>
1\t1\tGAIN1\tA\t<CNV>\t30\tPASS\tSVTYPE=CNV\tCN:GT\t0/1
"""
        path = self._write_temp_vcf(content)
        with self.assertRaises(SystemExit):
            validate_vcf(path)

    def test_missing_format_with_genotypes(self):
        content = """##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1
1\t1\tGAIN1\tA\t<CNV>\t30\tPASS\tSVTYPE=CNV\t0/1
"""
        path = self._write_temp_vcf(content)
        with self.assertRaises(SystemExit):
            validate_vcf(path)

    def test_duplicate_sample_names(self):
        content = """##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS1
1\t1\tGAIN1\tA\t<CNV>\t30\tPASS\tSVTYPE=CNV\tCN:GT\t0/1\t0/1
"""
        path = self._write_temp_vcf(content)
        with self.assertRaises(SystemExit):
            validate_vcf(path)

if __name__ == '__main__':
    unittest.main()
