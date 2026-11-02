# GSE71378 real cancer cfDNA (Snyder et al. 2016) - authors' GRCh37 BAMs, no realignment, no spike-ins
Folder on the cluster: $P/00_metadata/gse71378/full
1. bash 1_submit_bams.sh      -> mirror BAM -> ~12x / 1x / 0.1x analysis BAMs -> task lists
2. bash R0_get_references.sh  -> (login node, meanwhile) TCGA / OncoKB / gnomAD-SV / segdup / Hallmark references
3. bash 3_check_p3.sh         -> when "stage 1: complete":
4. bash 2_submit_stage2.sh    -> callers (GRCh37 resources from Simulation 1)
5. bash 4_collect_p3.sh       -> ~/cfcnv_gse71378_for_claude.tar.gz (calls + QC + prepared references)
Old scripts 1_submit_stage1.sh, 5_resume_stage1.sh, 6_switch_to_mirror.sh, 10/11/12, M1 are obsolete.
