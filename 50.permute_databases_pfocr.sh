#!/bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l mem_free=5G
#$ -l scratch=50G

scriptDir=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/PFOCRInPathwayAnalyses/
containerDir=/wynton/group/gladstone/biocore/containers
dataDir=/wynton/group/gladstone/biocore/projects/PFOCR/
dataDir1=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/
export SINGULARITY_BINDPATH="$containerDir,$scriptDir,$dataDir,$dataDir1"


f=$1
singularity exec $containerDir/pathway_enrichment_pathway_databases_evaluation_latest.sif Rscript permuted_gene_sets_pfocr.R $f

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


