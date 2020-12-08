#! /bin/bash
# aggregate cellranger snrnaseq counts

cellranger aggr \
--id=cellranger_rna_aggr_control \
--csv=crRnaAggr_control.csv \
--normalize=none \
--nosecondary