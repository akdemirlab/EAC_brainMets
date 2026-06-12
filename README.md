# EAC_brainMets

## Plotting and Analysis Scripts

This repository contains the scripts used for analysis and plotting of figures in the associated paper, EAC brain metastasis manuscript.

Software versions and dependencies are described in the Methods section of the paper.


## Figure 1

Scripts in the `Figure1` folder were used for analysis and plotting of the data shown in this figure.

The [Signal website](https://signal.mutationalsignatures.com/) was used for COSMIC mutational signature analysis.

## Figure 2

Scripts in the `Figure2` folder were used for analysis and plotting of the data shown in this figure.

Additionally, scripts from [AMPLIFISH](https://github.com/akdemirlab/AMPLIFISH) were used for analysis and plotting of this data.

## Figure 3

Scripts in the `Figure3` folder were used for analysis and plotting of the data shown in this figure.

- Heatmap scripts were developed from the [Navin Lab copy number pipeline (CopyKit)](https://github.com/navinlabcode/copykit).
- Additional plotting and analysis used the scripts described below.

### CycleViz

Plots were generated using [CycleViz](https://github.com/AmpliconSuite/CycleViz).

Example execution:

```bash
python3 CycleViz.py -g ../SampleID-WG01_amplicon1_graph.txt --cycles_file ../SampleID-WG01_amplicon1_cycles.txt \
    --cycle 3 --ref GRCh38 --gene_subset_file ../list.txt --annotate_structure genes --gene_fontsize 15 --tick_fontsize 7
```

### refphase

`refphase_Figure.R` was developed from [refphase](https://bitbucket.org/schwarzlab/refphase/src/master/).

Example execution:

```bash
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

## Figure 4

Scripts in the `Figure4` folder were used for analysis and plotting of the data shown in this figure.

Spatial niche plotting was performed using `spatial_niche_plot.py` in order to zoom in on particular regions of the tumor section.

## Figure 5

Scripts in the `Figure5` folder were used for analysis and plotting of the data shown in this figure.

[CellCharter](https://github.com/CSOgroup/cellcharter) (v0.3.7) was used for neighborhood analysis prior to running `cellcharter-neighborhood_spatial_analysis.py`.

