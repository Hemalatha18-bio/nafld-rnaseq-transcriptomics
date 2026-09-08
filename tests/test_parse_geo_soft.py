import gzip
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from parse_geo_soft import derive_analysis_group, parse_soft


def test_derive_analysis_group():
    assert derive_analysis_group("Control", "0") == "Control"
    assert derive_analysis_group("NAFLD", "0") == "NAFLD_F0_F1"
    assert derive_analysis_group("NAFLD", "1") == "NAFLD_F0_F1"
    assert derive_analysis_group("NAFLD", "2") == "NAFLD_F2"
    assert derive_analysis_group("NAFLD", "3") == "NAFLD_F3_F4"
    assert derive_analysis_group("NAFLD", "4") == "NAFLD_F3_F4"
    assert derive_analysis_group("NAFLD", "unknown") == "Unknown"


def test_parse_soft_extracts_characteristics(tmp_path):
    text = """^SAMPLE = GSM0000001
!Sample_title = Liver control sample
!Sample_source_name_ch1 = Control_control_liver biopsy
!Sample_characteristics_ch1 = nas score: 0
!Sample_characteristics_ch1 = fibrosis stage: 0
!Sample_characteristics_ch1 = group in paper: control
!Sample_characteristics_ch1 = disease: Control
!Sample_characteristics_ch1 = Stage: control
!Sample_description = control-demo
^SAMPLE = GSM0000002
!Sample_title = Liver patient
!Sample_source_name_ch1 = NAFLD_advanced_liver biopsy
!Sample_characteristics_ch1 = nas score: 6
!Sample_characteristics_ch1 = fibrosis stage: 4
!Sample_characteristics_ch1 = group in paper: advanced
!Sample_characteristics_ch1 = disease: NAFLD
!Sample_characteristics_ch1 = Stage: advanced
!Sample_description = disease-demo
"""
    path = tmp_path / "family.soft.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(text)

    frame = parse_soft(path)
    assert frame["sample_id"].tolist() == ["GSM0000001", "GSM0000002"]
    assert frame["fibrosis_stage"].tolist() == ["0", "4"]
    assert frame["analysis_group"].tolist() == ["Control", "NAFLD_F3_F4"]
    assert frame.loc[1, "nas_score"] == "6"
