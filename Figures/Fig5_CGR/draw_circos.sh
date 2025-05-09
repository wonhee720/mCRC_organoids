#!/usr/bin/bash

python /mCRC_organoids/FiguresFig5_CGR/circos_prep.py {sample_id} {input.snv} {input.indel} {input.cnv} {input.sv} {output.snv} {output.indel} {output.cnv} {output.sv}

Rscript /mCRC_organoids/FiguresFig5_CGR/circor.R {output.snv} {output.indel} {output.cnv} {output.sv} {final.circos_png} {sample_id}