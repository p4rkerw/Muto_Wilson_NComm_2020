#! /bin/bash
# aggregate cellranger snrnaseq counts
# TO RUN:
# bash cellrangerRnaCount.sh path/to/aggregation.csv

cellranger aggr \
--id=cellranger_rna_aggr_control \
--csv=crRnaAggr.csv \
--normalize=none \
--nosecondary