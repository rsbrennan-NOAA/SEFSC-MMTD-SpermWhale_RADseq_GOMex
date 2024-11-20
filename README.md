# Sperm whale radseq population structure

This repository contains scripts to analyze RADseq data from 73 sperm whale samples across the Gulf of Mexico and Atlantic.

Raw data can be found at: TBD

All scripts are found in the `scripts` directory and their function is described below.

For questions contact Reid Brennan: reid.brennan@noaa.gov or reid.brennan@gmail.com



### Filtering

1. Fastqc raw files: `fastqc.nf`
2. rename:
	- `rename.sh SNP052_sperm_whale_plateB_75_indexes_HS100_437.txt`
3. Trim with fastp: 
	- `nextflow run fastp.nf --fastq "raw_data/*.R1.fq.gz"`
4. Fastqc trimmed files: 
	- `nextflow run fastqc.nf --input_dir ./analysis/trimmed_files --outdir ./analysis/fastqc_post_trim`
5. Align with bwa mem2: 
	- `nextflow run bwa.nf -c nextflow.config -profile bwamem2`
6. calculate alignment stats: 
	- `nextflow run bam_stats.nf -c nextflow.config -profile bamstats`, outputs: `count.aligned.txt`
7. Call variants with freebayes: 
	- `freebayes_parallel.sh`
8. Filter SNPs:
	- `filter.1.sh`: First pass depth, missingness, bialleleic and snps only, depth.
	- `filter.2.sh`: Remove `Pmac075`, truncated file. max missing 0.7, maf 0.05, min mean depth 10x, allelic balance
	- `filter.R`: requires `HDplot.R`. Max depth, HDplot for paralogs
	- `filter.3.sh`: remove sites found in `filter.R`

`filtered.final.vcf.gz`: 4441 snps


### Analysis

- PCA with all individuals:
	- `plink_pca.sh`
	- plot results: `pca.R`
		- note that this includes related individuals
- Relatedness
	- `ngsrelate.sh`
	- `relatedness.R`
	- `filter_related.sh`
- PCA unrelated individuals only
	- `plink_unrelated.R`
	- `pca_unrelated.R`
- Admixture
	- `admixture.sh`
	- `admixture.R`
- DAPC
	- `dapc.R`
- genetic diversity
	- Pi
	- heterozygosity
	- Fst
- Effective population size

Mitochondrial analysis

- add here 




---

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
