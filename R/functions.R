
CalculateScoreFromLengths <- function(readlengths, ref_set) {
  ## "readlengths" must be passed as a frequency table
  
  if (sum(readlengths) == 0) {
    return(NA)
  }
  score <- weighted.mean(x = ref_set[as.numeric(names(readlengths))], w = readlengths)
  return(score)
}

ConvertBaseQualSymbol <- function(base_qual_symbol) {
  ## converts base quality symbol (ASCII of base qual + 33) to its numeric value
  
  library(DescTools)
  
  base_qual_score <- c()
  for (symbol in base_qual_symbol) {
    if (symbol == "") {
      base_qual_score <- append(base_qual_score, 0)  # no quality given
    } else {
      base_qual_score <- append(base_qual_score, DescTools::CharToAsc(symbol)-33)
    }
  }
  return(base_qual_score)
}

CountReadsAboveCutoff <- function(readlengths, ref_set, cutoff=0){
  if (sum(readlengths) == 0) {
    return(NA)
  }
  names(ref_set) <- seq(1,length(ref_set))
  ct_ref_set <- ref_set[ref_set > 0 ]
  return(sum(readlengths[names(ct_ref_set)]))
}

FormatRefSet <- function(ref_set_file) {
  ref_set <- read.delim(ref_set_file, header = FALSE)
  ref_set <- as.vector(ref_set$V1)
  return(ref_set)
}

GeneratePatientLevelFS <- function(bam_file, ref_set) {
  ## generates a tumour-agnostic score from all reads in the BAM

  lengths <- GetAllFragmentLengths(bam_file)
  fs <- CalculateScoreFromLengths(lengths, ref_set)
  score <- data.frame(FS = fs)

  lengths_table <- data.frame(lengths)
  names(lengths_table) <- c("LENGTH", "FREQ")

  return(list("score" = score,
              "lengths" = lengths_table))
}

GenerateReferenceSet <- function(healthy_lengths, tumour_lengths) {
  ## generates tumour:healthy log-ratio reference score dataset
  ## inputs must be passed as frequency tables
  
  scores <- NULL
  
  # bootstrap to smooth out noisy areas and reduce impact of sampling errors
  for (i in 1:1000) {
    
    # randomly sample reads and calculate probability density
    tlen <- table(factor(sample(as.numeric(names(tumour_lengths)), size = 10000,
                                prob = tumour_lengths, replace = TRUE),
                         levels = 1:900))/10000
    hlen <- table(factor(sample(as.numeric(names(healthy_lengths)), size = 10000,
                                prob = healthy_lengths, replace = TRUE),
                         levels = 1:900))/10000
    
    # exclude lengths with few reads
    tlen[(tlen + hlen) < 20/10000] <- 0
    hlen[(tlen + hlen) < 20/10000] <- 0
    
    # log2 ratio, maxed at -5 and +5
    score <- log(tlen/hlen, 2)
    score[is.na(score)] <- 0
    score[score == Inf] <- 5
    score[score == -Inf] <- -5
    
    scores <- cbind(scores, score)
  }
  score <- rowMeans(scores)
  return(score)
}

GenerateTumourInformedFS <- function(bam_ranges, variants, ref_set, fragment_score_cutoff=0) {
  ## generates a variant (VFS) and wildtype (WFS) score from all variant and
  ## wildtype fragments in the sample; only SNVs are currently supported
  
  require(dplyr)

  cat(format(Sys.time())," Getting fragment lengths...\n")
  df <- mapply(
    GetReadsAtThisSite,
    chr = variants$CHR,
    pos = as.numeric(variants$POS),
    ref = variants$REF,
    alt = variants$ALT,
    bam_range = bam_ranges,
    SIMPLIFY = FALSE
  )
  df <- do.call(cbind, df)
  
  cat(format(Sys.time())," Making table...\n")
  variant_info <- data.frame(dplyr::bind_rows(df["variant_info",]))
  variant_info <- variant_info %>%
    dplyr::arrange(as.numeric(sub("chr", "", CHR)), POS)
  
  var_fragments <- dplyr::distinct(dplyr::bind_rows(df["var_reads",]))
  wt_fragments <- dplyr::distinct(dplyr::bind_rows(df["wt_reads",]))
  
  var_lengths <- table(factor(as.integer(var_fragments$isize), levels = 1:900))
  wt_lengths <- table(factor(as.integer(wt_fragments$isize), levels = 1:900))
  
  cat(format(Sys.time())," Calculating Scores...\n")
  scores <- data.frame(WT_FRAG_CNT = sum(wt_lengths),
                       VAR_FRAG_CNT = sum(var_lengths),
                       CT_FRAG_CNT = CountReadsAboveCutoff(var_lengths, ref_set, cutoff=fragment_score_cutoff),
                       WFS = CalculateScoreFromLengths(wt_lengths, ref_set),
                       VFS = CalculateScoreFromLengths(var_lengths, ref_set))
  
  var_lengths_table <- data.frame(var_lengths)
  wt_lengths_table <- data.frame(wt_lengths)
  names(var_lengths_table) <- c("LENGTH", "FREQ")
  names(wt_lengths_table) <- c("LENGTH", "FREQ")
  
  return(list("scores" = scores,
              "var_lengths" = var_lengths_table,
              "wt_lengths" = wt_lengths_table,
              "variant_info" = variant_info))
}

GenerateVariantLevelFS <- function(bam_ranges, variants, ref_set) {
  ## generates a variant (VFS) and wildtype (WFS) score for each variant
  ## SNVs and indels are supported

  require(dplyr)
  
  cat(format(Sys.time())," Getting fragment scores...\n")
  scores <- mapply(GetScoresAtThisSite,
                   chr = variants$CHR,
                   pos = variants$POS,
                   ref = variants$REF,
                   alt = variants$ALT,
                   bam_range = bam_ranges,
                   MoreArgs = list(ref_set = ref_set))
  scores <- data.frame(t(scores))
  rownames(scores) <- NULL

  scores <- scores %>%
    dplyr::arrange(as.numeric(sub("chr", "", CHR)), POS)
  scores <- data.frame(lapply(scores, as.character), stringsAsFactors = FALSE)

  scores$POS <- as.numeric(scores$POS)
  scores$WT_READ_CNT <- as.numeric(scores$WT_READ_CNT)
  scores$VAR_READ_CNT <- as.numeric(scores$VAR_READ_CNT)
  scores$MED_WT_LEN <- as.numeric(scores$MED_WT_LEN)
  scores$MED_VAR_LEN <- as.numeric(scores$MED_VAR_LEN)
  scores$WFS <- as.numeric(scores$WFS)
  scores$VFS <- as.numeric(scores$VFS)

  return(scores)
}

GetAllFragmentLengths <- function(bam_file, min_mapq = 30) {
  ## returns frequency table of fragment lengths in the BAM
  
  library(Rsamtools)
  
  bamfile <- BamFile(bam_file, yieldSize = 1e7)  # downsample to 10 million reads
  param <-  ScanBamParam(what = scanBamWhat(),
                         flag = scanBamFlag(isPaired = TRUE,
                                            isDuplicate = FALSE,
                                            isSecondaryAlignment = FALSE,
                                            isUnmappedQuery = FALSE),
                         mapqFilter = min_mapq)
  bam <- scanBam(bamfile, param = param)
  bam <- GetProperPairReads(bam[[1]])
  
  lengths <- as.integer(abs(bam$isize))
  lengths <- table(factor(lengths, levels = 1:900))
  
  return(lengths)
}

GetBamReadsAtVariantSites <- function(bam_file, variants, min_mapq = 30) {
  ## returns a list of length(variants)
  ## each element contains the reads that map to its corresponding variant site
  
  require(Rsamtools)
  variants$POS <- as.numeric(variants$POS)
  
  cat(format(Sys.time())," Getting reads from variant sites...\n")
  
  cat(format(Sys.time())," ...loading bam...\n")
  bamfile <- BamFile(bam_file)
  
  counted_bam <- countBam(bamfile)
  all_reads <- counted_bam$records
  ## counting how many reads pass filter and are within intervals is TOO slow

  cat(format(Sys.time())," ...scanning bam...\n")
  param <- ScanBamParam(what = scanBamWhat(),
                        which = GRanges(variants$CHR,
                                        IRanges(start = variants$POS,
                                                end = variants$POS)),
                        flag = scanBamFlag(isPaired = TRUE,
                                           isDuplicate = FALSE,
                                           isSecondaryAlignment = FALSE,
                                           isUnmappedQuery = FALSE),
                        mapqFilter = min_mapq)
  bam_ranges <- scanBam(bamfile, param = param)
  
  cat(format(Sys.time())," ...getting pairs...\n")
  bam_ranges <- lapply(bam_ranges, GetProperPairReads)
  bam_analysis <- 
    list(
      "all_reads" = all_reads,
      "bam_ranges" = bam_ranges
    )
  return(bam_analysis)
}

GetPosInRead <- function(pos, read_start, cigar) {
  ## returns base position relative to read (1-based)
  ## accounts for indels and soft-clippings ahead of the base
  ## adapted from MRDetect (Zviran et al. 2020)
  
  library(GenomicAlignments)
  
  # cigar operations: M = match, I = insertion, D = deletion, S = soft-clipping
  cigar_ops <- explodeCigarOps(cigar)
  cigar_lengths <- explodeCigarOpLengths(cigar)
  parsed_cigar <- data.frame(cigar_ops, cigar_lengths)
  names(parsed_cigar) <- c("op", "length")
  
  if (parsed_cigar[1,]$op == "S") {  # start of read is soft-clipped
    read_start <- read_start - parsed_cigar[1,]$length
  }
  
  raw_pos_in_read <- pos-read_start+1
  if (!any(c("I","D") %in% parsed_cigar$op)) {  # no indels
    return(raw_pos_in_read)
  }
  
  pos_in_read <- raw_pos_in_read
  cigar_pos <- 1
  for (i in 1:nrow(parsed_cigar)) {
    op <- parsed_cigar[i,]$op
    length <- parsed_cigar[i,]$length
    
    if (pos_in_read >= cigar_pos + length) {
      if (op == "I") {
        pos_in_read <- pos_in_read + length
      } else if (op == "D") {
        pos_in_read <- pos_in_read - length
      }
    } else if (pos_in_read < cigar_pos) {
      break
    } else {
      if (op == "I") {
        pos_in_read <- pos_in_read + length
      } else if (op != "M") {
        return(-1)  # base not in read
      }
    }
    if (op != "D") {
      cigar_pos <- cigar_pos+length
    }
  }
  
  return(pos_in_read)
}

GetProperPairReads <- function(bam_range) {
  ## "proper pairs" map to opposite strands of the same chromosome
  ## may be inward or outward oriented; distance between reads is not considered

  require(Rsamtools)

  flags <- data.frame(bamFlagAsBitMatrix(bam_range$flag))
  opposite_strands <- (!flags$isMinusStrand & flags$isMateMinusStrand) |
    (flags$isMinusStrand & !flags$isMateMinusStrand)

  pairs_on_same_chromosome <- bam_range$rname == bam_range$mrnm

  is_proper_pair <- opposite_strands & pairs_on_same_chromosome
  proper_pairs <- lapply(bam_range, `[`, is_proper_pair)
  return(proper_pairs)
}

GetReadsAtThisSite <- function(chr, pos, ref, alt, bam_range,
                               min_base_qual = 25, min_reads_per_allele = 1,
                               min_total_reads = 8) {
  ## returns the unique fragment ID (qname) and fragment size for each read
  
  require(dplyr)
  require(DescTools)
  
  var_reads <- data.frame(qname = character(), isize = numeric())
  wt_reads <- data.frame(qname = character(), isize = numeric())
  variant_info <- data.frame(CHR = character(), POS = numeric(),
                             REF = character(), ALT = character(),
                             WT_READ_CNT = numeric(), VAR_READ_CNT = numeric(),
                             MED_WT_LEN = numeric(), MED_VAR_LEN = numeric())
  read_start <- bam_range$pos
  
  if (length(read_start) > 0) {  # reads exist at this site
    cigar <- bam_range$cigar
    pos_in_read <- mapply(GetPosInRead, read_start = read_start, cigar = cigar,
                          MoreArgs = list(pos = pos))
    base <- substring(as.character(bam_range$seq), pos_in_read, pos_in_read)
    
    base_qual_symbol <- substring(as.character(bam_range$qual), pos_in_read, pos_in_read)
    base_qual_score <- ConvertBaseQualSymbol(base_qual_symbol)
    
    supports_alt <- base == alt & base_qual_score >= min_base_qual
    supports_ref <- base == ref & base_qual_score >= min_base_qual
    
    var_read_count <- sum(supports_alt)
    wt_read_count <- sum(supports_ref)
    
    if (var_read_count >= min_reads_per_allele &
        wt_read_count >= min_reads_per_allele &
        var_read_count + wt_read_count >= min_total_reads) {
      
      # paired reads have same qname and absolute isize
      var_qname <- bam_range$qname[supports_alt]
      wt_qname <- bam_range$qname[supports_ref]
      
      var_isize <- abs(bam_range$isize[supports_alt])
      wt_isize <- abs(bam_range$isize[supports_ref])
      
      # deduplicate overlapping paired reads with same base
      var_reads <- dplyr::distinct(data.frame(qname = var_qname, isize = var_isize))
      wt_reads <- dplyr::distinct(data.frame(qname = wt_qname, isize = wt_isize))
      variant_info <- data.frame(CHR = chr, POS = pos, REF = ref, ALT = alt,
                                 WT_READ_CNT = wt_read_count,
                                 VAR_READ_CNT = var_read_count,
                                 MED_WT_LEN = median(wt_reads$isize),
                                 MED_VAR_LEN = median(var_reads$isize))
    }
  }
  
  return(list("var_reads" = var_reads,
              "wt_reads" = wt_reads,
              "variant_info" = variant_info))
}

GetScoresAtThisSite <- function(chr, pos, ref, alt, bam_range, ref_set) {
  ## returns a VFS and WFS for fragments mapping to this site
  
  library(dplyr)
  
  var_reads <- data.frame(qname = character(), isize = numeric())
  wt_reads <- data.frame(qname = character(), isize = numeric())
  
  score <- data.frame(CHR = chr, POS = pos, REF = ref, ALT = alt,
                      WT_READ_CNT = 0, VAR_READ_CNT = 0,
                      MED_WT_LEN = NA, MED_VAR_LEN = NA,
                      WFS = NA, VFS = NA)
  
  read_start <- bam_range$pos
  
  if (length(read_start) > 0) {  # reads exist at this site
    cigar <- bam_range$cigar
    pos_in_read <- mapply(GetPosInRead, read_start = read_start, cigar = cigar,
                          MoreArgs = list(pos = pos))
    seq <- bam_range$seq
    
    sort_reads <- SortReads(seq, pos_in_read, ref, alt, cigar)
    supports_alt <- sort_reads$supports_alt
    supports_ref <- sort_reads$supports_ref
    
    # get number of var and wt reads
    var_read_count <- sum(supports_alt)
    wt_read_count <- sum(supports_ref)
    
    if (var_read_count > 0 | wt_read_count > 0) {
      
      # paired reads have same qname and absolute isize
      var_qname <- bam_range$qname[supports_alt]
      wt_qname <- bam_range$qname[supports_ref]
      
      var_isize <- abs(bam_range$isize[supports_alt])
      wt_isize <- abs(bam_range$isize[supports_ref])
      
      # deduplicate overlapping paired reads with same base
      var_reads <- dplyr::distinct(data.frame(qname = var_qname, isize = var_isize))
      wt_reads <- dplyr::distinct(data.frame(qname = wt_qname, isize = wt_isize))
      
      var_lengths <- table(factor(as.integer(var_reads$isize), levels = 1:900))
      wt_lengths <- table(factor(as.integer(wt_reads$isize), levels = 1:900))
      
      # use fragment lengths to generate score
      vfs <- CalculateScoreFromLengths(var_lengths, ref_set)
      wfs <- CalculateScoreFromLengths(wt_lengths, ref_set)
      
      score$WT_READ_CNT <- wt_read_count
      score$VAR_READ_CNT <- var_read_count
      score$MED_WT_LEN <- median(wt_reads$isize)
      score$MED_VAR_LEN <- median(var_reads$isize)
      score$VFS <- vfs
      score$WFS <- wfs
    }
  }
  return(score)
}

GetVariantsFromVCF <- function(vcf_file) {
  require(vcfR)
  require(dplyr)
  
  vcf <- read.vcfR(vcf_file)
  vcf <- cbind(as.data.frame(vcf@fix))
  vcf$POS <- as.numeric(vcf$POS)
  variants <- dplyr::distinct(data.frame(CHR = vcf$CHROM, POS = vcf$POS,
                                         REF = vcf$REF, ALT = vcf$ALT))
  return(variants)
}

IsTrueInsertion <- function(query_ranges, pos_in_read, alt_len) {
  return(nrow(query_ranges[query_ranges$op == "I" &
                           query_ranges$start == pos_in_read+1 &
                           query_ranges$width == alt_len-1,]) > 0)
}

IsTrueDeletion <- function(query_ranges, reference_ranges, pos_in_read, ref_len) {
  deletion_exists <- which(query_ranges$op == "D" & query_ranges$start == pos_in_read+1)
  if (identical(deletion_exists, integer(0))) {
    return(FALSE)
  }
  return(reference_ranges[deletion_exists,]$width == ref_len-1)
}

IsTrueMatch <- function(query_ranges, pos_in_read, ref_len) {
  return(nrow(query_ranges[query_ranges$op == "M" &
                           query_ranges$start <= pos_in_read &
                           query_ranges$end >= pos_in_read+ref_len-1,]) > 0)
}

SortReads <- function(seq, pos_in_read, ref, alt, cigar) {
  ## sorts reads into alternate v.s. reference, accounting for indels

  library(GenomicAlignments)

  alt_len <- nchar(alt)
  ref_len <- nchar(ref)
  num_bases <- if (alt_len > ref_len) alt_len else ref_len
  bases <- substr(as.character(seq), pos_in_read, pos_in_read+num_bases-1)

  if (alt_len == ref_len) {  # SNV
    supports_alt <- bases == alt
    supports_ref <- bases == ref

  } else {  # indel
    cigar_ops <- explodeCigarOps(cigar)
    query_ranges <- lapply(cigarRangesAlongQuerySpace(cigar), as.data.frame)
    reference_ranges <- lapply(cigarRangesAlongReferenceSpace(cigar), as.data.frame)

    for (i in 1:length(query_ranges)) {
      query_ranges[[i]]$op <- cigar_ops[[i]]
    }

    # search for reads supporting alt
    if (alt_len > ref_len) {  # insertion
      supports_alt <- bases == alt &
        mapply(IsTrueInsertion,
               query_ranges = query_ranges, pos_in_read = pos_in_read,
               MoreArgs = list(alt_len = alt_len))

    } else if (alt_len < ref_len) {  # deletion
      supports_alt <- startsWith(bases, alt) &
        mapply(IsTrueDeletion,
               query_ranges = query_ranges, reference_ranges = reference_ranges,
               pos_in_read = pos_in_read, MoreArgs = list(ref_len = ref_len))
    }

    # search for reads supporting ref
    supports_ref <- !supports_alt & startsWith(bases, ref) &
      mapply(IsTrueMatch,
             query_ranges = query_ranges, pos_in_read = pos_in_read,
             MoreArgs = list(ref_len = ref_len))

  }
  return(list("supports_alt" = supports_alt,
              "supports_ref" = supports_ref))
}

summarize_fragments <- function(variant_scores, patient_scores, all_reads){
  
  sample_candidate_SNVs = length(variant_scores$REF)
  reads_with_variants_checked = sum(variant_scores$WT_READ_CNT, na.rm = T) + sum(variant_scores$VAR_READ_CNT, na.rm = T)
  
  cf_scores_checked <- variant_scores %>% filter(WT_READ_CNT > 0)
  sites_checked = length(cf_scores_checked$REF)

  cf_scores_detected <- cf_scores_checked %>% filter(VAR_READ_CNT > 0)
  sites_detected = length(cf_scores_detected$REF)

  if(sites_detected == 0){

    ct_sites_detected = 0

  }else{

    ct_scores_detected <- cf_scores_detected %>% filter(VFS > 0)
    ct_sites_detected = length(ct_scores_detected$REF)

  }
  
  summary.table <- as.data.frame(cbind(
    "all_reads" = all_reads,
    "sample_candidate_SNVs" = sample_candidate_SNVs,
    "sites_checked"  = sites_checked ,
    "reads_with_variants_checked" = reads_with_variants_checked,
    "sites_detected" = sites_detected,
    "ct_sites_detected" = ct_sites_detected,
    "fragments_checked" = patient_scores$WT_FRAG_CNT + patient_scores$VAR_FRAG_CNT,
    "fragments_detected" = patient_scores$VAR_FRAG_CNT,
    "ct_fragments_detected" = patient_scores$CT_FRAG_CNT,
    "WFS" =  patient_scores$WFS,
    "VFS" =  patient_scores$VFS
  ))
  return(summary.table)
}
