# The "Fragmentation Score" Method

Several studies have shown that circulating tumour DNA (ctDNA) fragments differ in length from healthy cell-free DNA (cfDNA). This method, adapted from [Vessies et al. (2022)](https://github.com/DCLVessies/Fragmentomics) and [Wong et al. (2023)](https://github.com/pughlab/fragmentomics/tree/main/fragment_score), explores the use of a length-derived "fragmentation score" to detect ctDNA. The higher the score, the higher the likelihood that ctDNA is present.

The Fragmentation Score (FS) has three different uses:

Test | Description | Purpose | Requires
--- | --- | --- | ---
Patient-Level FS   | Generates a score for a sample using all reads. Also returns a fragment length frequencies. | Distinguish between ctDNA-positive and -negative samples. | BAM
Variant-Level FS   | Generates a score for each given variant. | Distinguish between tumour and nontumour variants. | BAM, list of variants in VCF format
Tumour-Informed FS | Generates a score for a sample using variant-containing reads. Also returns variant and wildtype fragment length frequencies. | Distinguish between ctDNA-positive and -negative samples. | BAM, VCF of primary tumour


## Parameters
Flag | Description
--- | ---
`--id` | sample id
`--bam` | BAM file
`--vcf` | VCF file (not used for Patient-Level FS)
`--ref` | reference set
`--libdir` | scripts directory
`--outdir` | output directory

A reference set by Vessies et al. is provided, but to generate your own, use the `GenerateReferenceSet` function in `/R/functions.R`.


## Usage
The `/R/` directory contains a script for each test. All functions are located in `/R/functions.R`. To run a script, pass the required options to the following command:

```
Rscript /R/tumour_informed_fs.R --id sample_id --bam /path/to/bam_file --vcf /path/to/vcf_file --ref /ref/vessies_reference_set.txt --libdir /R/ --outdir path/to/output_directory
```


## Example Outputs

Patient-Level FS:
```
FS
0.123456789
```

Variant-Level FS:
```
CHR	POS	REF	ALT	WT_READ_CNT	VAR_READ_CNT	MED_WT_LEN	MED_VAR_LEN	WFS		VFS
chr1	1000000	C	CA	2000		50		167		145		-0.123456789	0.987654321
```

Tumour-Informed FS:
```
WT_FRAG_CNT	VAR_FRAG_CNT	WFS		VFS
20000		1000		-0.123456789	-0.987654321
```

Fragment Length Frequency Table (for Patient-Level FS and Tumour-Informed FS):
```
LENGTH	FREQ
1	0
2	2
3	1
4	5
5	4
6	9
7	11
...
```
