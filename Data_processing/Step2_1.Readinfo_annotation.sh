#!/usr/bin/bash

/mCRC/inhouse_scripts/readinfoAnnot.wrapper.v10-lwh.sh \
-i /home/users/wonhee720/Projects/26_Sox9_Org_clonal/03_vcf/union/snv_merged/snvs.Proj20_Proj26_merged.v2.vcf.gz \
-o /home/users/wonhee720/Projects/26_Sox9_Org_clonal/03_vcf/union/snv_merged/snvs.Proj20_Proj26_merged.readinfo.v2.vcf.gz \
-b {bams_from_all_samples_seperated_by_comma}\
-f 19 \
-s {sample_names_from_all_samples_seperated_by_comma}