# EAC_brainMets


This repository contains the scripts used in the EAC brain metastasis manuscript


# Scripts in Figure 1 folder were used for WGS analysis and visualization.

Comutplot_with_allsamples_figure.ipynb

DriverGene_proportions_bytype_Figure.R

OncogeneCN_comparison_Figure.R

refphase_Figure.R developed from https://bitbucket.org/schwarzlab/refphase/src/master/
execution of refphase R script: 

```
SAMPLE_IDS="SAMPLEID"
PATIENT_ID=PATIENTID

ASCAT_PATH=/results/ascat/
REFPHASE_DIR=/results/refphase/$PATIENT_ID/
MEDICC_DIR=/results/medicc2/$PATIENT_ID/

. ~/.bashrc
module load R/4.2.1
RSCRIPT=/scripts/sv/refphase.R
mkdir -p $REFPHASE_DIR && cd $REFPHASE_DIR
Rscript --vanilla $RSCRIPT $ASCAT_PATH "$SAMPLE_IDS" $REFPHASE_DIR


conda activate /miniconda3/envs/medicc2
mkdir -p $MEDICC_DIR && cd $MEDICC_DIR
medicc2 --input-type tsv $REFPHASE_DIR/refphase-segmentation.tsv $MEDICC_DIR
```





# Scripts in Figure 2 folder were used for spatial transcriptomics analysis and visualization.

Xenium_Annotations_plotting_Script.ipynb, 

Xenium_ERBB2_EGFR_expression_boxplots_Figure.R,

TLS-Finder [script from GITHUB Page](https://github.com/AAKoksoy/TLS-Finder)




# Scripts in Figure 3 folder were used for ecDNA and expression analysis and visualization.

CycleViz [_script from GITHUB Page](https://github.com/AmpliconSuite/CycleViz)

execution of CycleViz script:

```
python3 CycleViz.py -g ../SampleID-WG01_amplicon1_graph.txt --cycles_file ../SampleID-WG01_amplicon1_cycles.txt 
    --cycle 3 --ref GRCh38 --gene_subset_file ../list.txt --annotate_structure genes --gene_fontsize 15 --tick_fontsize 7
```

Xenium_ERBB2_EGFR_expression_boxplots_Figure.R



# Scripts in Figure 4 folder were used for scWGS analysis and visualization. 

Heatmap scripts developed from the Navin Lab copy number pipeline https://github.com/navinlabcode/copykit

scWGS_heatmap_wholegenome_Figure.R

Pseudobulk_coverage_Figure.R

scWGS_Heatmap_bychromosome_Figure.R

