# EAC_brainMets


This repository contains the scripts used in the manuscript: Recurrent ERBB2 alterations are associated with esophageal adenocarcinoma brain metastases



# Scripts in Figure 1 folder were used for downstream analysis and visualization.

A: Comutplot_with_allsamples_figure.ipynb

B: DriverGene_proportions_bytype_Figure.R

C: OncogeneCN_comparison_Figure.R

D: jabba_visuzalize_v3_Figure.R developed from https://github.com/mskilab-org/JaBbA

E: refphase_Figure.R developed from https://bitbucket.org/schwarzlab/refphase/src/master/
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

F: MRI images



# Scripts in Figure 2 folder were used for downstream analysis and visualization.

A: Graphical summary 

B: Xenium_Annotations_plotting_Script.ipynb, Xenium_ERBB2_EGFR_expression_boxplots_Figure.R

C: Spatial Maps produced by Xenium Explorer, TLS-Finder [script from GITHUB Page](https://github.com/AAKoksoy/TLS-Finder)

D: Spatial Maps produced by Xenium Explorer,  MRI images, Xenium_comparing_celltype_primaryvsmet_barplot.R



# Scripts in Figure 3 folder were used for downstream analysis and visualization.

A: Xenium Explorer images and microscope images

B: ecDNAcountscomparison_Figure.R

C: CycleViz [_script from GITHUB Page](https://github.com/AmpliconSuite/CycleViz)

execution of CycleViz script:

```
python3 CycleViz.py -g ../SampleID-WG01_amplicon1_graph.txt --cycles_file ../SampleID-WG01_amplicon1_cycles.txt 
    --cycle 3 --ref GRCh38 --gene_subset_file ../list.txt --annotate_structure genes --gene_fontsize 15 --tick_fontsize 7
```

D: Xenium_ERBB2_EGFR_expression_boxplots_Figure.R



# Scripts in Figure 4 folder were used for downstream analysis and visualization. 

Heatmap scripts developed from the Navin Lab copy number pipeline https://github.com/navinlabcode/copykit

A: scWGS_heatmap_wholegenome_Figure.R

B: Pseudobulk_coverage_Figure.R

C: scWGS_Heatmap_bychromosome_Figure.R

D: scWGS_heatmap_wholegenome_Figure.R

E: Pseudobulk_coverage_Figure.R

F: scWGS_Heatmap_bychromosome_Figure.R
