#!/bin/bash
/usr/bin/Rscript /home/boyuanli/bashscript/bin/ruvb2/chip-seq/diff_peak/diff_peak_pipeline.r -d /WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/Pol2 -f sample.csv
/usr/bin/Rscript /home/boyuanli/bashscript/bin/ruvb2/chip-seq/diff_peak/diff_peak_pipeline.r -d /WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/SMC1 -f sample.csv
/usr/bin/Rscript /home/boyuanli/bashscript/bin/ruvb2/chip-seq/diff_peak/diff_peak_pipeline.r -d /WORK/lbyybl/WH/rvb/cor_uniq/sample_Ruvbl/CTCF -f sample.csv

# https://github.com/winston-lab/mnase-seq