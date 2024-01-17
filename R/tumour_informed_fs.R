#! /usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(vcfR))
suppressMessages(library(dplyr))
suppressMessages(library(Rsamtools))
suppressMessages(library(DescTools))

option_list <- list(
  make_option(c("--bam"), type = "character", help = "Path to BAM file."),
  make_option(c("--vcf"), type = "character", help = "Path to VCF file."),
  make_option(c("--outfilePrefix"), type = "character", help = "Prefix for output files"),
  make_option(c("--ref"), default="ref/vessies_reference_set.txt", type = "character", help = "Path to reference set file."),
  make_option(c("--libdir"), default=NULL, type = "character", help = "Path to scripts directory."),
  make_option(c("--outdir"), default="./", type = "character", help = "Path to output directory.")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
options(scipen=0, stringsAsFactors=F)

outfilePrefix <- opt$outfilePrefix
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

output <- GenerateTumourInformedFS(bam_file, variants, ref_set)

write.table(output$scores, file.path(outdir, paste0(outfilePrefix, "_scores.txt")),
            sep = "\t", row.names = FALSE)
write.table(output$var_lengths, file.path(outdir, paste0(outfilePrefix, "_var_lengths.txt")),
            sep = "\t", row.names = FALSE)
write.table(output$wt_lengths, file.path(outdir, paste0(outfilePrefix, "_wt_lengths.txt")),
            sep = "\t", row.names = FALSE)
write.table(output$variant_info, file.path(outdir, paste0(outfilePrefix, "_variant_info.txt")),
            sep = "\t", row.names = FALSE)

cat(format(Sys.time())," All analyses complete.")