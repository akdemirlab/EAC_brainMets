library(tidyverse)
library(refphase)

read_data_snps = function(files, sample_ids, prefix){
  
  data_snps = tibble()
  
  for(row in 1:length(files)){
    tsv = read_tsv(files[row])
    tmp = tsv %>% dplyr::rename(chrom="Chromosome", pos="Position", baf="BAF", logR=`Log R`) %>%
      dplyr::select(chrom, pos, baf, logR) %>%
      dplyr::filter(!is.na(logR)) %>%
      dplyr::filter(baf!=1 & baf!=0)
    tmp$sample_id = sample_ids[row]
    data_snps = rbind(data_snps, tmp)
  }
  data_snps$germline_zygosity = 'het' #v0.3.2
  data_snps = data_snps %>% dplyr::filter(chrom!="chrX" & chrom!="chrY")
  write_tsv(data_snps, paste0(prefix, "snps.tsv"))
}

read_segments = function(files, sample_ids, prefix){
  segments = tibble()
  
  for(row in 1:length(files)){
    colnames=c('seg_number', 'chrom', 'start', 'end', 'norm_cn_major', 'norm_cn_minor', 'cn_major', 'cn_minor')
    csv = read_csv(files[row], col_names=colnames)
    tmp = csv 
    tmp$sample_id = sample_ids[row]
    tmp = tmp %>% dplyr::select(sample_id, chrom, start, end, cn_major, cn_minor)
    segments = rbind(segments, tmp)
  }
  segments = segments %>% dplyr::filter(chrom!="chrX" & chrom!="chrY")
  write_tsv(segments, paste0(prefix, "segments.tsv"))
}

read_purity_ploidy = function(files, sample_ids, prefix){
  purity_ploidy = tibble()
  
  for(row in 1:length(files)){
    colnames=c('type', 'value')
    tsv = read_delim(files[row], col_names=colnames, " ")
    purity = tsv %>% dplyr::filter(type=="rho") %>% pull(value) %>% as.double()
    ploidy = tsv %>% dplyr::filter(type=="Ploidy") %>% pull(value) %>% as.double()
    #print(sample_ids)
    #print(sample_ids[row])
    tmp = data.frame(sample=sample_ids[row], purity=purity, ploidy=ploidy) 
    purity_ploidy = rbind(purity_ploidy, tmp)
  }
  write_tsv(purity_ploidy, paste0(prefix, "purity_ploidy.tsv"))
}

refphage_plot = function(sample_ids, prefix){
  
  refphase_data_tsv <- refphase_load(data_format = "tsv", samples = sample_ids,
                                     tsv_prefix = prefix)
  
  refphase_output = refphase(refphase_data_tsv)
  
  results = refphase_output
  
  seg = as.data.frame(results$phased_segs) %>% dplyr::rename(chrom=seqnames, sample_id=group_name) #v0.3.2
  write_tsv(seg, file.path(prefix, "refphase-segmentation.tsv"))
  seg2 = seg %>% dplyr::filter(chrom!="X" & chrom!="Y")
  write_tsv(seg2, file.path(prefix, "refphase-segmentation_nxy.tsv"))
  #write_segs(results$phased_segs, file = file.path(prefix, "refphase-segmentation.tsv"))
  
  # (optional) output the SNPs, including phasing information
  write_snps(results$phased_snps, file = file.path(prefix, "refphase-phased-snps.tsv.gz"))
  
  # Ploidy might have changed, if we have updated any copy numbers
  write.table(results$sample_data, file = file.path(prefix, "refphase-sample-data-updated.tsv"), sep = "\t", row.names = FALSE)

  png(file.path(prefix, "refphase-genome.png"), width = 40, height = 40, units = "cm", res = 300, family = "Sans")
  plot(results, what = "genome")
  dev.off()
  
  png(file.path(prefix, "refphase-summary.png"), width = 40, height = 40, units = "cm", res = 300, family = "Sans")
  plot(results, what = "summary")
  dev.off()
  
  png(file.path(prefix, "refphase-copy_numbers.png"), width = 40, height = 40, units = "cm", res = 300, family = "Sans")
  plot(results, what = "copy_numbers")
  dev.off()
  
  png(file.path(prefix, "refphase-BAF.png"), width = 40, height = 40, units = "cm", res = 300, family = "Sans")
  plot(results, what = "BAF")
  dev.off()
  
  # Note the %02d in the file name which will be replaced with the chromosome number
  # This function creates one plot per chromosome
  png(file.path(prefix, "refphase-chromosome%02d.png"), width = 40, height = 40, units = "cm", res = 300, family = "Sans")
  plot(results, what = "chromosome")
  dev.off()
}


refphage_from_ascat = function(ascat_path, sample_ids, result_prefix){
  #print(result_prefix)
  dir.create(result_prefix, showWarnings = FALSE)
  path = ascat_path
  #print(sample_ids)
  cpnb_files = paste0(path, sample_ids, '/', sample_ids, '.copynumber.txt.gz')
  seg_files = paste0(path, sample_ids, '/', sample_ids, '.copynumber.caveman.csv')
  ploidy_files = paste0(path, sample_ids, '/', sample_ids, '.samplestatistics.txt')
  read_data_snps(cpnb_files, sample_ids, result_prefix)
  read_segments(seg_files, sample_ids, result_prefix)
  read_purity_ploidy(ploidy_files, sample_ids, result_prefix)
  refphage_plot(sample_ids, result_prefix)
  
}

args = commandArgs(trailingOnly=TRUE)

#print(paste0('args[1] is:', args[1]))
#print(paste0('args[2] is:', args[2]))
#print(paste0('args[3] is:', args[3]))

ascat_path = args[1]
sample_ids = unlist(strsplit(args[2], "\\s+") )
result_dir = args[3]

refphage_from_ascat(ascat_path,sample_ids, result_dir)