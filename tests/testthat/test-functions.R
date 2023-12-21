library(testthat)

fragment_score_repo <- Sys.getenv("FRAGMENT_SCORE_REPO")
source(paste0(fragment_score_repo, "/R/functions.R"))


### PATIENT-LEVEL FRAGMENTATION SCORE ----

test_that("test GeneratePatientLevelFS", {
  bam_file <- paste0(fragment_score_repo, "/tests/test_data/test.bam")
  ref_set_file <- paste0(fragment_score_repo, "/ref/vessies_reference_set.txt")

  ref_set <- FormatRefSet(ref_set_file)
  output <- GeneratePatientLevelFS(bam_file, ref_set)

  lengths <- c(156, 156, 158, 158, 159, 166, 166, 167, 168, 168, 168, 185)

  expect_equal(output$lengths$FREQ, as.vector(table(factor(as.integer(lengths),
                                                           levels = 1:900))))
})


### VARIANT-LEVEL FRAGMENTATION SCORE ----

test_that("test GenerateVariantLevelFS", {
  bam_file <- paste0(fragment_score_repo, "/tests/test_data/test.bam")
  vcf_file <- paste0(fragment_score_repo, "/tests/test_data/test.vcf")
  ref_set_file <- paste0(fragment_score_repo, "/ref/vessies_reference_set.txt")

  variants <- GetVariantsFromVCF(vcf_file)
  ref_set <- FormatRefSet(ref_set_file)
  output <- GenerateVariantLevelFS(bam_file, variants, ref_set)

  wt_reads <- c(156, 156, 158, 158, 159, 167, 168, 168, 168)
  var_reads <- c(166, 166, 185)
  all_reads <- c(wt_reads, var_reads)

  wt_frags <- c(156, 158, 159, 167, 168, 168)
  var_frags <- c(166, 185)
  all_frags <- c(wt_frags, var_frags)

  var_1 <- unname(unlist(output[1,]))
  var_2 <- unname(unlist(output[2,]))
  var_3 <- unname(unlist(output[3,]))

  expect_equal(var_1[5:8], as.character(c(length(wt_reads), length(var_reads),
                                          median(wt_frags), median(var_frags))))
  expect_equal(var_2[5:8], as.character(c(length(all_reads), 0,
                                          median(all_frags), NA)))
  expect_equal(var_3[5:8], as.character(c(0, 0, NA, NA)))
})


### TUMOUR-INFORMED FRAGMENTATION SCORE ----

test_that("test GenerateTumourInformedFS", {
  bam_file <- paste0(fragment_score_repo, "/tests/test_data/test.bam")
  vcf_file <- paste0(fragment_score_repo, "/tests/test_data/test.vcf")
  ref_set_file <- paste0(fragment_score_repo, "/ref/vessies_reference_set.txt")

  variants <- GetVariantsFromVCF(vcf_file)
  ref_set <- FormatRefSet(ref_set_file)
  output <- GenerateTumourInformedFS(bam_file, variants, ref_set)

  wt_reads <- c(156, 156, 158, 158, 159, 167, 168, 168)
  var_reads <- c(166, 166, 185)

  wt_frags <- c(156, 158, 159, 167, 168, 168)
  var_frags <- c(166, 185)

  expect_equal(unname(unlist(output$scores)[1:2]), c(length(wt_frags), length(var_frags)))

  expect_equal(output$var_lengths$FREQ, as.vector(table(factor(as.integer(var_frags),
                                                          levels = 1:900))))
  expect_equal(output$wt_lengths$FREQ, as.vector(table(factor(as.integer(wt_frags),
                                                         levels = 1:900))))

  expect_equal(unname(unlist(output$variant_info)[5:8]),
               as.character(c(length(wt_reads), length(var_reads),
                              median(wt_frags), median(var_frags))))
})


### HELPER FUNCTIONS ----

test_that("test GetPosInRead", {
  pos <- 7
  read_start <- 5
  cigar_list <- c("3M", "2S4M", "2M1I1M", "1M1D1M", "2M1D1M")

  output <- lapply(cigar_list, GetPosInRead, pos = pos, read_start = read_start)
  expect_equal(unlist(output), c(3, 5, 4, 2, -1))
})

test_that("test SortReads, SNV", {
  ref <- "A"
  alt <- "G"
  seq <- c("A", "G", "G")
  pos_in_read <- c(1, 1, 1)
  cigar <- c("1M", "1M", "1M")

  output <- SortReads(seq, pos_in_read, ref, alt, cigar)

  expect_equal(output$supports_ref, c(TRUE, FALSE, FALSE))
  expect_equal(output$supports_alt, c(FALSE, TRUE, TRUE))
})

test_that("test SortReads, insertion", {
  ref <- "A"
  alt <- "AT"
  seq <- c("GACGG", "GATGG", "GATGG", "GATGG", "GATGG", "GGATG", "GATGG")
  pos_in_read <- c(2, 2, 2, 2, 2, 3, 2)
  cigar <- c("5M",      # ref
             "5M",      # matches alt but no insertions
             "2M1I2M",  # insertion (at pos 2)
             "2M2I1M",  # insertion different length
             "3M1I1M",  # insertion elsewhere
             "3M1I1M",  # insertion (at pos 3)
             "1M1I3M")  # neither ref nor alt

  output <- SortReads(seq, pos_in_read, ref, alt, cigar)

  expect_equal(output$supports_ref, c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(output$supports_alt, c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that("test SortReads, deletion", {
  ref <- "AT"
  alt <- "A"
  seq <- c("GATGG", "GATGG", "GATGG", "GATGG", "GACGG", "GACGG")
  pos_in_read <- c(2, 2, 2, 2, 2, 2)
  cigar <- c("5M",      # ref
             "2M1D3M",  # matches ref but deletion
             "2M2D3M",  # deletion but wrong length
             "3M1D2M",  # deletion but elsewhere
             "5M",      # does not match ref, no deletion
             "2M1D3M")  # deletion

  output <- SortReads(seq, pos_in_read, ref, alt, cigar)

  expect_equal(output$supports_ref, c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(output$supports_alt, c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE))
})

