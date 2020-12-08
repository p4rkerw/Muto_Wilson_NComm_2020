#!/bin/bash
# this script will aggregate snAtacseq output using cellranger-atac 
# TO RUN: bash cellrangerAtacAggr.sh path/to/aggregation.csv

# aggregate the output
cellranger-atac aggr \
--id=cellranger_atac_aggr_control \
--csv=crAtacAggr_control.csv \
--reference=$HOME/reference/refdata-cellranger-atac-GRCh38-1.2.0 \
--normalize="none" \
--nosecondary