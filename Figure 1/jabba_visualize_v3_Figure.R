library(tidyverse)
library(gTrack)
library(rtracklayer)
library(kableExtra)    
library(magrittr)
library(tidyr)
library(gGnome)
library(gUtils)

jab_v3_dir = "/results/jba_v3"
frag_dir = "/results/frag"
svaba_dir = "/results/svaba"
script_dir = "/scripts"

figures_dir = "/results/jba_figures_v3"

options(stringsAsFactors = FALSE)
setwd("/results");


##====================================================
## Read metadata and plot whole genome and every each chromosome
##
library(tidyverse)
data_meta = read_csv(file.path(script_dir, "meta_jabba_v3.csv"))
epgap = tibble(name=character(), epgap=numeric())

for(row in 1:nrow(data_meta))
{
  normal = data_meta[[row, "normal"]]
  tumor = data_meta[[row, "tumor"]]
  mrn = data_meta[[row, "patient"]]
  print(row)

  if(!file.exists(figures_dir)) 
  {
      dir.create(figures_dir)
  }

  sample_dir = file.path(figures_dir, paste(mrn, tumor, sep="_"))
  whole_genome_dir = file.path(figures_dir, "whole_genome")
  if(!file.exists(sample_dir)) 
  {
    dir.create(sample_dir)
  }
  if(!file.exists(whole_genome_dir)) 
  {
    dir.create(whole_genome_dir)
  }

  jab_v3_file = file.path(jab_v3_dir, tumor, "jabba.simple.rds")
  cov_normal_file = file.path(frag_dir, normal, "cov_50k.rds")
  cov_tumor_file = file.path(frag_dir, tumor, "cov_50k.rds")
  jab_v3 = gG(jabba = jab_v3_file)

  epgap = epgap %>% add_row(name=tumor, epgap=jab_v3$nodes$dt$epgap[1])

  gt_jab_v3 = jab_v3$gt
  gt_jab_v3$name = "jabba"

  cov_normal = readRDS(cov_normal_file)
  cov_tumor = readRDS(cov_tumor_file)

  cov_normal_tibble = cov_normal %>% 
    as_tibble() %>%
    group_by(seqnames) %>%
    mutate(reads.median = median(reads.corrected, na.rm=TRUE))

  gt_cov_normal_50k = gTrack(cov_normal, y.field = 'reads.corrected', name = 'cov_normal', y0=0, y1=2)
  gt_cov_tumor_50k = gTrack(cov_tumor, y.field = 'reads.corrected', name = 'cov_tumor', y0=0, y1=2)

  ###################
  ## generate pictures
  ##
  wg_file_name = paste(mrn, "_", tumor, ".jpg", sep = "")
  jpeg(file=file.path(sample_dir, wg_file_name), width = 4400, height = 2400, res=300)
  plot(c(gt_cov_normal_50k, gt_cov_tumor_50k, gt_jab_v3),
       c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
         '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
         '21', '22', 'X', 'Y'))
  title(main = paste(mrn, tumor, sep = "_"))
  dev.off()
  file.copy(file.path(sample_dir, wg_file_name), whole_genome_dir)
  
  for (chrom in c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                  '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                  '21', '22', 'X', 'Y'))
  {
    tmp_file_name = paste("chr", chrom, ".jpg", sep = "")
    jpeg(file=file.path(sample_dir, tmp_file_name), width = 3600, height = 2400, res=300)
    plot(c(gt_cov_normal_50k, gt_cov_tumor_50k, gt_jab_v3),
         c(chrom))
    title(main = paste(mrn, tumor, chrom, sep = "_"))
    dev.off()
  }
}


