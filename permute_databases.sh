#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -l h_rt=24:00:00
#$ -t 1-536
scriptDir=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/PFOCRInPathwayAnalyses/
containerDir=/wynton/group/gladstone/biocore/containers
dataDir=/wynton/group/gladstone/biocore/projects/PFOCR/
dataDir1=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/permuted_geneset_databases_results_updated_pfocr_gene_set_sizes_may_2021/
export SINGULARITY_BINDPATH="$containerDir,$scriptDir,$dataDir,$dataDir1"


f=$1
singularity exec $containerDir/pathway_enrichment_pathway_databases_evaluation_latest.sif Rscript permuted_gene_sets.R $SGE_TASK_ID

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


