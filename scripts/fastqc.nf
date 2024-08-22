#!/usr/bin/env nextflow

// Default parameters
params.input_dir = "$projectDir/data"
params.outdir = "$projectDir/results"

log.info """\
         FASTQC AND MULTIQC WORKFLOW
         ===========================
         Input directory : ${params.input_dir}
         Output directory: ${params.outdir}
         """
         .stripIndent()

// Create a channel from input files
input_files = Channel.fromPath("${params.input_dir}/*.{fastq,fq,fastq.gz,fq.gz}")
                     .ifEmpty { error "Cannot find any fastq files in ${params.input_dir}" }

process fastqc {
    tag "$reads.simpleName"
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    path reads

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc -t 16 $reads
    """
}

process multiqc {
    publishDir params.outdir, mode: 'move'
    
    input:
    path '*', stageAs: 'fastqc_results/'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}

workflow {
    fastqc(input_files)
    multiqc(fastqc.out.fastqc_results.collect())
}
