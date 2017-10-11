## Bam-Merge-Preprocessing
a split from BFMC with Indel Realignment and Base Quality Recalibration (gatk3)

![bmpp](docs/bamFilterMergeSort.png)

## Workflow Options

This workflow has a lot of options, but user will be using mostly the ones listed below:

### Basic Options

Parameter | Value | Description
---|---|---
output_dir | seqware-results | A standard SeqWare parameter specifying the sub-directory where the output files will be provisioned
output_prefix	| ./ | A standard SeqWare parameter specifying the root directory where the output files will be provisioned
manual_output	| false | Whether or not to use manual output. When false, a random integer will be inserted into the path of the file in order to ensure uniqueness. When true, the output files will be moved to the location of output_prefix/output_dir
input_files\* | | The input files to be used by the workflow.
output_identifiers | | Identifier(s) used for producing the output file names of the result file, may be a donor name but may also include information about template type, tissue and other pieces of information.
java | | path to java binary
tmp_dir	| tmp | temporary directory
r_dir | | path to R home directory
do_sam_filter	| true | flag for enabling/disabling samtools filtering (quality and other thing defined by samtools_flags
do_remove_duplicates | true | flag for enabling/disabling duplicate removal with Picard tools
do_mark_duplicates | true | flag for enabling/disabling marking duplicates with Picard tools
aligner_name | | This information, although optional is useful since it allows to identify alignment software used to produce .bam files
sort_order | | sort order for sorting with
chr_sizes | | A comma separated list of intervals to parallelize by.  Note: intervals missing from this list will be filtered out
ref_fasta | hg19_random.fa | reference fasta file (path)
queue	| | queue is usually NOT set (but is still available as a parameter)
gatk_dbsnp_vcf | | dbsnp_142.CAF.chr.vcf.gz | path to gatk vcf file
do_split_trim_reassign_quality | | flag for enabling/disabling split and trim step for RNA Seq BAM files

\* If there are multiple input/output groups, provide the output_identifiers as an ordered list with ";" as the group delimiter and "," as the delimiter for multiple files. For example:
```
input_files=/path/to/group1/input1.bam,/path/to/group1/input2.bam,/path/to/group1/input3.bam;/path/to/group2/input1.bam
```
If there are multiple input/output groups, provide the output_identifiers as an ordered list with ";" as the delimiter. For example:
```
output_identifiers=GROUP1_OUTPUT_NAME;GROUP2_OUTPUT_NAME
```
Intervals can be combined together to run within the same job by using the "+" delimiter. For example:
```
chr_sizes=chr1,chr2,chr3+chr4+chr5,chr6+chr7+ch8
would result in 4 parallel jobs and all reads chr9+ being filtered out from the final bam.
```
### Picard

Parameter | Value | Description
---|---|---
picard_mark_duplicates | | path to Picard's Mark Duplicate.jar
picard_mark_duplicates_mem_mb | 10000 | memory allocated to MarkDuplicate operation in Mb
picard_merge_sam | | path to Picard's MergeSam.jar
picard_sort_sam | | path to Picard's SortSam.jar
picard_merge_sam_mem_mb | 10000 | memory allocated to MergeSam operation
samtools_filter_other_params | | additional parameters for siltering with samtools
picard_mark_duplicates_other_params | | additional parameters for MarkDuplicates (Picard tools)
picard_merge_other_params | | additional parameters for merging with Picard
picard_dir | | picard directory
picard_sort_other_params | | additional parameters for sorting with Picard tools
picard_merge_use_threading | | flag that determines if threading should be used when running Picard tools

### GATK

Parameter | Value | Description
---|---|---
gatk_realign_target_creator_xmx | 12 | Memory allocated to JVM when running RealignerTargetCreator in Gb
gatk_indel_realigner_xmx | 12 | Memory allocated to JVM when running IndelRealigner walker in Gb
gatk_print_reads_xmx | 12 | Memory allocated to JVM when running PrintReads walker in Gb
gatk_base_recalibrator_xmx | 8192 | Memory allocated to JVM when running BaseRecalibrator walker in Mb
gatk_base_recalibrator_mem | 9728 | Memory allocated to the SeqWare job when running BaseRecalibrator in Mb
gatk_base_recalibrator_nct | 24 | CpuThreadsPerDataThread setting for BaseRecalibrator
gatk_base_recalibrator_smp | 20 | Nodes requested for BaseRecalibrator
gatk_sched_overhead_mem | 4 | Memory overhead for GATK SeqWare jobs, Mb
downsampling_coverage | | downsampling coverage for IndelRealign job
preserve_qscores_less_than | | q-score filtering parameter
interval_padding | | interval padding parameter that is used by IndelRealign
gatk_realigner_target_creator_params | | additional parameters for RealignerTargetCreator walker
gatk_indel_realigner_params | | additional parameters for IndelRealigner walker
gatk_base_recalibrator_params | | additional parameters for BaseRecalibrator walker
gatk_print_reads_params | | additional parameters for PrintReads walker
downsampling_type | | downsampling used by indel realignment step
gatk_jar | | Path to GATK jar
gatk_key | | GATK key file
interval_files | | intervals to use with GATK's -L option to limit the search space to certain regions
bqsr_covariates | | Base quality recalibration covariates to consider when doing base quality recalibration
split_cigar_Xmxg | 16 | Memory allocated to JVM when running RealignerTargetCreator walker in Gb
split_cigar_RMQF | 255 | reassign mapping quality from
split_cigar_RMQT | 60 | reassign mapping quality to
split_cigar_mem | 16 | Memory allocated to the SplitNCigarReads job in Gb
reassign_One_Mapping_Quality | false | additional parameter for SplitNCigarReads program

## Decider Options

The Decider for BMPP allows automatic scheduling of workflow runs as more of new data becomes available

Parameter | Value | Description
---|---|---
group-by | | Opionally specify how input files should be separated into workflow runs groups.  For example, if "ROOT_SAMPLE_NAME" was provided as the argument, one workflow run per "ROOT_SAMPLE_NAME" would be scheduled for all input files related to that "ROOT_SAMPLE_NAME".
output-path | ./ | define output prefix
output-folder | seqware-results | define output folder
library-template-type | | Restrict the processing to samples of a particular template type, e.g. WG, EX, TS
tissue-type | | Restrict the processing to samples of particular tissue types, e.g. P, R, X, C. Multiple values can be comma-separated. Tissue types are processed individually (all R's together, all C's together)
use-tissue-prep | true | flags if we need to use tissue prep when grouping data
use-tissue-region | true | flags if we need to use tissue region when grouping data
sam-filter-flag | | bit flag to use for filtering with samtools
min-map-quality | | Minimal mapping quality to use for filtering with samtools
do-filter | true | flag to indicate if we want filtering with samtools
do-mark-duplicates | true | flag to indicate if we want marking of duplicates with Picard tools
do-remove-duplicates | true |	flag to indicate if we want to remove duplicates with Picard tools
group-by-aligner | true | flag to indicate if we want to use aligner info for grouping data
chr-sizes | | Ids of reference chromosomes, comma-separated
interval-padding | 100 | Padding for intervals used by GATK's HaplotypeCaller
stand-emit-conf | 1 | Standard emit confidence (higher means more confidence)
stand-call-conf | 30 | Standard call confidence (higher means more confidence)
downsampling | | Downsampling type, one of "NONE", "ALL_READS", "BY_SAMPLE"
dbsnp | | Path to dbSNP file
disable-bqsr | false | Disable Base quality calibaration
verbose | false | Set logging to verbose (aka debug mode)
do-split-and-trim | false | flag to indicate if we want to run splitNCigarReads for RNASeq BAM files
