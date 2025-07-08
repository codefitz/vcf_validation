import pathlib
import pytest
from vcf_validation import validate_vcf

test_dir = pathlib.Path(__file__).parent / "data"


def test_valid_vcf(capsys):
    valid_path = test_dir / "valid.vcf"
    validate_vcf(str(valid_path))
    captured = capsys.readouterr()
    assert "VCF file validation completed" in captured.out


def test_invalid_alt():
    invalid_path = test_dir / "invalid_alt.vcf"
    with pytest.raises(SystemExit):
        validate_vcf(str(invalid_path))
