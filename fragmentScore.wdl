version 1.0

workflow fragmentScore {
	input {
		File? plasmabam
		File? plasmabai
		String outputFileNamePrefix
		String tumorVcfSampleName
		File tumorvcf
		File tumorvcfindex
		Boolean full_analysis_mode = true
		String controlFileList = "/.mounts/labs/gsi/src/pwgs_hbc/1.0/HBC.bam.list"
	}

	parameter_meta {
		plasmabam: "plasma input .bam file"
		plasmabai: "plasma input .bai file"
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		outputFileNamePrefix: "Prefix for output file"
		tumorVcfSampleName: "ID for WGS tumor sample, must match .vcf header"
		controlFileList: "tab seperated list of bam and bai files for healthy blood controls"
	}

	call filterVCF {
		input:
		tumorvcf = tumorvcf,
		tumorvcfindex = tumorvcfindex,
		tumorSampleName = tumorVcfSampleName
	}

	if(full_analysis_mode) {
		
		call parseControls {
			input:
			controlFileList = controlFileList
		}

		scatter (control in parseControls.controlFiles) {
			call detectSNVs as detectControl {
				input:
				plasmabam = control[0],
				plasmabai = control[1],
				plasmaSampleName = basename(control[0], ".bam"),
				tumorvcf = filterVCF.filteredvcf
			}
		}

		call detectSNVs as detectSample {
			input:
			plasmabam = plasmabam,
			plasmabai = plasmabai,
			tumorvcf = filterVCF.filteredvcf,
			plasmaSampleName = outputFileNamePrefix
		}

		call snvDetectionSummary {
			input:
			controlCalls = select_all(detectControl.fragment_summary),
			sampleCalls = detectSample.fragment_summary
		}

	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Workflow for fragment-score"
		dependencies:
		[
			{
				name: "bcftools/1.9",
				url: "https://github.com/samtools/bcftools"
			}
		]
		output_meta: {
			snpcount: "number of SNPs in vcf after filtering",
			filteredvcf: "filtered vcf",
			final_call: "final_call",
			all_calls: "all_calls"
		}
	}
	output {
		File snpcount = filterVCF.snpcount
		File? filteredvcf = filterVCF.filteredvcf
		File? final_call = snvDetectionSummary.final_call
		File? all_calls = snvDetectionSummary.all_calls
	}
}

task filterVCF {
	input {
		File tumorvcf
		File tumorvcfindex
		String tumorSampleName
		String tumorVCFfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'slippage' | FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' | FILTER~'multiallelic' | FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"
		String tumorVAF = "0.1"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String? varType = "snps"
		String difficultRegions = "--regions-file /.mounts/labs/CGI/scratch/fbeaudry/wdl/fragment-score/ref/hg38-norepeats.bed"
		String modules = "bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 10
	}

	parameter_meta {
		tumorvcf: "tumor vcf file, bgzip"
		tumorvcfindex: "tumor vcf index file"
		tumorSampleName: "ID for WGS tumor sample"
		tumorVCFfilter: "set of filter calls to incl. in tumor VCF (any line with these flags will be included"
		tumorVAF: "Variant Allele Frequency for tumor VCF"
		genome: "Path to loaded genome .fa"
		difficultRegions: "Path to .bed excluding difficult regions, string must include the flag --regions-file "
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"

	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools view -s ~{tumorSampleName} ~{difficultRegions} ~{tumorvcf} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{tumorVCFfilter}" |\
		~{"$BCFTOOLS_ROOT/bin/bcftools filter -i \"TYPE='" + varType + "'\" |"} \
		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{tumorVAF}" > ~{tumorSampleName}.SNP.vcf

		awk '$1 !~ "#" {print}' ~{tumorSampleName}.SNP.vcf | wc -l >SNP.count.txt

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File filteredvcf = "~{tumorSampleName}.SNP.vcf"
		File snpcount = "SNP.count.txt"
	}

	meta {
		output_meta: {
			filteredvcf: "Filtered vcf",
			snpcount: "Number of SNPs after filtering"
		}
	}
}

task parseControls {
	input {
		String controlFileList 
		Int jobMemory = 4
		Int timeout = 12
	}

	parameter_meta {
		controlFileList: "file with list of control files"
		jobMemory: "Memory for this task in GB"
		timeout: "Timeout in hours, needed to override imposed limits"
	}

	command <<<
		python <<CODE
		import os, re

		with open("~{controlFileList}") as f:
			for line in f:
				line = line.rstrip()
				tmp = line.split("\t")
				r = tmp[0] + "\t" + tmp[1]
				print(r)
		f.close()
		CODE
	>>>

	runtime {
		memory:  "~{jobMemory} GB"
		timeout: "~{timeout}"
	}

	output {
		Array[Array[File]] controlFiles = read_tsv(stdout())
	}
}

task detectSNVs {
	input {
		File? plasmabam
		File? plasmabai
		String plasmaSampleName
		File tumorvcf
		String libdir = "/.mounts/labs/CGI/scratch/fbeaudry/wdl/fragment-score/R/"
		String fragmentScoreScript = "/.mounts/labs/CGI/scratch/fbeaudry/wdl/fragment-score/R/variant_level_fs.R"
		String vessies_reference_set = "/.mounts/labs/CGI/scratch/fbeaudry/wdl/fragment-score/ref/vessies_reference_set.txt"
		String modules = "fragmentomics/0.1"
		Int jobMemory = 64
		Int threads = 4
		Int timeout = 20
	}

	parameter_meta {

		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"

	}

	command <<<
		set -euo pipefail

		Rscript ~{fragmentScoreScript} \
			--libdir ~{libdir} \
			--sampleid ~{plasmaSampleName} \
			--bam ~{plasmabam} \
			--vcf ~{tumorvcf} \
			--ref ~{vessies_reference_set} \
			--outdir ./

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? snvDetectionReadsScored = "variant_scores.txt" 
		File? fragment_summary = "fragment_summary.txt" 
	}

	meta {
		output_meta: {
			snvDetectionReadsScored: "Reads with potential for tumor, with their scores"
		}
	}
}


task snvDetectionSummary {
	input {
		File? sampleCalls
		Array[File] controlCalls
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
		String modules = "fragmentomics/0.1"
		String collateScoreScript = "/.mounts/labs/CGI/scratch/fbeaudry/wdl/fragment-score/R/collate_results.R"
		String print_json = "FALSE"
		String pvalue = 0.05
	}

	parameter_meta {
		sampleCalls: "file of detection rate call for sample"
		controlCalls: "array of file of detection rate calls for HBCs"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		cat ~{sep=' ' controlCalls} | awk '$1 !~ "all_reads" {print}' > HBCs.txt

		cat ~{sampleCalls} HBCs.txt >cohort_summary.txt

		Rscript ~{collateScoreScript} \
			--results cohort_summary.txt \
			--pval ~{pvalue} \
			--print_json ~{print_json}

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File? all_calls = "cohort_summary.txt"
		File? final_call = "mrd.fragment.txt"
	}

	meta {
		output_meta: {
			all_calls : "HBC and sample mrdetect results",
			final_call: "Comparison of sample with controls"
		}
	}
}