import sys
from pathlib import Path

project_root = Path("/Users/wuhhh/Desktop/徐老师/代码/metagenomic")
sys.path.insert(0, str(project_root))

from metagenomic_refactor.virus_analysis import prepare_influenza_reference_set

result = prepare_influenza_reference_set(
    pre="flu_demo",
    species="Influenza virus",
    requested_ref="noref",
    fq1="/Users/wuhhh/Desktop/徐老师/代码/Baiyiapp_example/inf/NGS_Inf/inf_H1N1_91_R1.fq",
    fq2="/Users/wuhhh/Desktop/徐老师/代码/Baiyiapp_example/inf/NGS_Inf/inf_H1N1_91_R2.fq",
    long_type="",
    threads=8,
)

print(result)
