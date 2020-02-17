#! /bin/bash
# aggregate cellranger snrnaseq counts
# TO RUN:
# bash cellrangerRnaAggr.sh path/to/aggregation.csv

cellranger aggr \
--id=cellranger_rna_aggr_control \
--csv=crRnaAggr_control.csv \
--normalize=none \
--nosecondary