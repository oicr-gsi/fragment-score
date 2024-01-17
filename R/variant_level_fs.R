#! /usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(Rsamtools))
suppressMessages(library(vcfR))
suppressMessages(library(DescTools))

option_list <- list(
  make_option(c("--bam"), type = "character", help = "Path to BAM file."),
  make_option(c("--vcf"), type = "character", help = "Path to VCF file."),
  make_option(c("--sampleid"), type = "character", help = "Prefix for output files"),
  make_option(c("--ref"), default="ref/vessies_reference_set.txt", type = "character", help = "Path to reference set file."),
  make_option(c("--libdir"), default=NULL, type = "character", help = "Path to scripts directory."),
  make_option(c("--outdir"), default="./", type = "character", help = "Path to output directory.")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
options(scipen=0, stringsAsFactors=F)

sampleid <- opt$sampleid
bam_file <- opt$bam
vcf_file <- opt$vcf
ref_set_file <- opt$ref
outdir <- opt$outdir


if(is.null(opt$libdir)){
  libdir <- paste(Sys.getenv(c("BASE_DIR")), sep='/')
}else{
  libdir <- opt$libdir
}
source(paste0(libdir,"/functions.R"))

cat(format(Sys.time())," Getting Variants from VCF...\n")
variants <- GetVariantsFromVCF(vcf_file)
ref_set <- FormatRefSet(ref_set_file)

bam_analysis <- GetBamReadsAtVariantSites(bam_file, variants)

bam_ranges <- bam_analysis$bam_ranges
variant_scores <- GenerateVariantLevelFS(bam_ranges, variants, ref_set)
patient_scores <- GenerateTumourInformedFS(bam_ranges, variants, ref_set)

fragment_summary <- summarize_fragments(variant_scores, patient_scores$scores, bam_analysis$all_reads)
fragment_summary$sample_id <- sampleid

write.table(variant_scores, 
            file.path(outdir,  paste0("variant_scores.txt")),
            append = F, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
            col.names = TRUE
)

write.table(fragment_summary, 
            file.path(outdir,  paste0("fragment_summary.txt")),
            append = F, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA",dec = ".", row.names = FALSE, 
            col.names = TRUE
)

  
cat(format(Sys.time())," All analyses complete.")
