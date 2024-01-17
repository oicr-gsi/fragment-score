##version 1.0

####packages####
library(optparse)

process_results <- function(results){
  results$label <- "CONTROLS"
  results$label[1] <- "THIS SAMPLE"
  this_sample_coverage <- results$reads_with_variants_checked[1]
  
  results$sites_detected_byreads <- results$sites_detected / results$reads_with_variants_checked
  results$sites_detected_normalized <- results$sites_detected_byreads * this_sample_coverage
  
  results$reads_detected_byreads <- results$fragments_detected / results$reads_with_variants_checked
  results$reads_detected_normalized <- results$reads_detected_byreads * this_sample_coverage
  
  results$ct_sites_detected_byreads <- results$ct_sites_detected / results$reads_with_variants_checked
  results$ct_sites_detected_normalized <- results$ct_sites_detected_byreads * this_sample_coverage
  
  results$ct_fragments_detected_byreads <- results$ct_fragments_detected / results$reads_with_variants_checked
  results$ct_fragments_detected_normalized <- results$ct_fragments_detected_byreads * this_sample_coverage
  
  return(results)
}

get_mrd_stats <- function(results, pval_cutoff){
  
  sites_zscore <- (results$sites_detected_normalized[results$label == "THIS SAMPLE"] - mean(results$sites_detected_normalized))/ sd(results$sites_detected_normalized)
  sites_pvalue <- pnorm(sites_zscore,lower.tail=F)
  sites_cutoff <- (qnorm(pval_cutoff,lower.tail = F) * sd(results$sites_detected_normalized)) +  mean(results$sites_detected_normalized[results$label == "CONTROLS"])
  
  if(sites_pvalue < pval_cutoff){significantly_more_sites = "TRUE"}else if(sites_pvalue >= pval_cutoff){significantly_more_sites = "FALSE"}
  
  reads_zscore <- (results$reads_detected_normalized[results$label == "THIS SAMPLE"] - mean(results$reads_detected_normalized))/ sd(results$reads_detected_normalized)
  reads_pvalue <- pnorm(reads_zscore,lower.tail=F)
  reads_cutoff <- (qnorm(pval_cutoff,lower.tail = F) * sd(results$reads_detected_normalized)) +  mean(results$reads_detected_normalized[results$label == "CONTROLS"])
  
  if(reads_cutoff < pval_cutoff){significantly_more_reads = "TRUE"}else if(reads_cutoff >= pval_cutoff){significantly_more_reads = "FALSE"}
  
  ##
  ct_sites_zscore <- (results$ct_sites_detected_normalized[results$label == "THIS SAMPLE"] - mean(results$ct_sites_detected_normalized))/ sd(results$ct_sites_detected_normalized)
  ct_sites_pvalue <- pnorm(ct_sites_zscore,lower.tail=F)
  ct_sites_cutoff <- (qnorm(pval_cutoff,lower.tail = F) * sd(results$ct_sites_detected_normalized)) +  mean(results$ct_sites_detected_normalized[results$label == "CONTROLS"])
  
  if(ct_sites_pvalue < pval_cutoff){significantly_more_ct_sites = "TRUE"}else if(ct_sites_pvalue >= pval_cutoff){significantly_more_ct_sites = "FALSE"}
  
  ct_fragments_zscore <- (results$ct_fragments_detected_normalized[results$label == "THIS SAMPLE"] - mean(results$ct_fragments_detected_normalized))/ sd(results$ct_fragments_detected_normalized)
  ct_fragments_pvalue <- pnorm(ct_fragments_zscore,lower.tail=F)
  ct_fragments_cutoff <- (qnorm(pval_cutoff,lower.tail = F) * sd(results$ct_fragments_detected_normalized)) +  mean(results$ct_fragments_detected_normalized[results$label == "CONTROLS"])
  
  if(ct_fragments_pvalue < pval_cutoff){significantly_more_ct_fragments = "TRUE"}else if(ct_fragments_pvalue >= pval_cutoff){significantly_more_ct_fragments = "FALSE"}
  
  mrd_stats <- list(
    "all_reads" = results$all_reads[results$label == "THIS SAMPLE"],
    
    "sites_checked" =  results$sites_checked[results$label == "THIS SAMPLE"],
    "reads_with_variants_checked" =  results$reads_with_variants_checked[results$label == "THIS SAMPLE"],
    "fragments_checked" =  results$fragments_checked[results$label == "THIS SAMPLE"],
    
    "sites_detected" =  results$sites_detected[results$label == "THIS SAMPLE"],
    "controls_mean_sites_detection" = mean(results$sites_detected[results$label == "CONTROLS"]),
    "sites_cutoff" = sites_cutoff,
    "sites_zscore" = sites_zscore,
    "sites_pvalue" = sites_pvalue,
    "significantly_more_sites" = significantly_more_sites,
    
    "ct_sites_detected" =  results$ct_sites_detected[results$label == "THIS SAMPLE"],
    "controls_mean_ct_sites_detection" = mean(results$ct_sites_detected[results$label == "CONTROLS"]),
    "ct_sites_cutoff" = ct_sites_cutoff,
    "ct_sites_zscore" = ct_sites_zscore,
    "ct_sites_pvalue" = ct_sites_pvalue,
    "significantly_more_ct_sites" = significantly_more_ct_sites,
    
    "fragments_detected" =  results$fragments_detected[results$label == "THIS SAMPLE"],
    "controls_mean_fragments_detection" = mean(results$fragments_detected[results$label == "CONTROLS"]),
    "fragments_cutoff" = reads_cutoff,
    "fragments_zscore" = reads_zscore,
    "fragments_pvalue" = reads_pvalue,
    "significantly_more_fragments" = significantly_more_reads,
    
    "ct_fragments_detected" =  results$ct_fragments_detected[results$label == "THIS SAMPLE"],
    "controls_mean_ct_fragments_detection" = mean(results$ct_fragments_detected[results$label == "CONTROLS"]),
    "ct_fragments_cutoff" = ct_fragments_cutoff,
    "ct_fragments_zscore" = ct_fragments_zscore,
    "ct_fragments_pvalue" = ct_fragments_pvalue,
    "significantly_more_ct_fragments" = significantly_more_ct_fragments,
    
    "WFS" =  results$WFS[results$label == "THIS SAMPLE"],
    "VFS" =  results$VFS[results$label == "THIS SAMPLE"]
    
  )
  return(mrd_stats)
}

options(scipen=999)
options(digits = 5)

# get options
option_list = list(
  make_option(c("-r", "--results"), type="character", default=NULL, help="results file path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./", help="output directory", metavar="character"),
  make_option(c("-p", "--pval"), type="numeric", default=0.001, help="p-value cutoff", metavar="numeric"),
  make_option(c("-j", "--print_json"), type="character", default="FALSE", help="p-value cutoff", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

results_path <- opt$results
pval_cutoff <- opt$pval
outdir <- opt$outdir
print_json <- opt$print_json

results_df <- read.table(results_path, header = T)
results <- process_results(results_df)
mrd_stats <- get_mrd_stats(results, pval_cutoff)

if(print_json == "TRUE"){
  
  library(jsonlite)
  mrd.json <- jsonlite::toJSON(mrd_stats, pretty=TRUE, auto_unbox=TRUE)
  write(mrd.json, file = paste0(outdir, "mrd.fragment.json"))
  
}
mrd_stats = as.data.frame(mrd_stats)
write.table(
  mrd_stats,
  file = paste0(outdir, "mrd.fragment.txt"),
  append = F, quote = FALSE, sep = "\t", 
  eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
  col.names = TRUE
)
