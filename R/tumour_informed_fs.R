library(optparse)

option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID."),
  make_option(c("--bam"), type = "character", help = "Path to BAM file."),
  make_option(c("--vcf"), type = "character", help = "Path to VCF file."),
  make_option(c("--ref"), type = "character", help = "Path to reference set file."),
  make_option(c("--libdir"), type = "character", help = "Path to scripts directory."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

id <- opt$id
bam_file <- opt$bam
vcf_file <- opt$vcf
ref_set_file <- opt$ref
libdir <- opt$libdir
outdir <- opt$outdir

source(paste0(libdir,"/functions.R"))

variants <- GetVariantsFromVCF(vcf_file)
ref_set <- FormatRefSet(ref_set_file)
output <- GenerateTumourInformedFS(bam_file, variants, ref_set)

dir.create(file.path(outdir, id), showWarnings = FALSE)
write.table(output$scores, file.path(outdir, id, paste0(id, "_scores.txt")),
            sep = "\t", row.names = FALSE)
write.table(output$var_lengths, file.path(outdir, id, paste0(id, "_var_lengths.txt")),
            sep = "\t", row.names = FALSE)
write.table(output$wt_lengths, file.path(outdir, id, paste0(id, "_wt_lengths.txt")),
            sep = "\t", row.names = FALSE)
write.table(output$variant_info, file.path(outdir, id, paste0(id, "_variant_info.txt")),
            sep = "\t", row.names = FALSE)

q('no')


