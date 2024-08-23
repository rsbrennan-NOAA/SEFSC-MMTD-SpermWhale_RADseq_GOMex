#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.indir = "trimmed_reads"
params.outdir = "aligned"
params.refDir = "/home/rbrennan//reference_genomes/sperm_whale"
params.refGenome = "${params.refDir}/GCF_002837175.3_ASM283717v5_genomic.fna"

// Process to align reads
process bwamem2 {
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(reads)
    path(refGenome)
    path(ref_0123)
    path(ref_amb)
    path(ref_ann)
    path(ref_bwt)
    path(ref_pac)

    output:
    path("${sampleId}.bam"), emit: bam
    path("${sampleId}.bam.bai"), emit: bai

    script:
    """
    # Load modules
    module load aligners/bwa-mem2/2.2.1
    module load bio/samtools/1.19

    # Set read group
    RG="@RG\\tID:${sampleId}\\tSM:${sampleId}"

    # Align reads
    bwa-mem2 mem -t 4 -R "\$RG" $refGenome $reads | \
    samtools view -S -h -u - | \
    samtools sort -  > ${sampleId}.bam

    # Index BAM file
    samtools index ${sampleId}.bam
    """
}

// Workflow
workflow {
    // Input channel with sample ID extraction
    reads_ch = channel
        .fromPath("${params.indir}/*.fastq.gz")
        .map { file ->
            def sampleId = file.name.split('\\.')[0]
            return tuple(sampleId, file)
        }

    // Reference genome and associated files
    refGenome = file(params.refGenome)
    ref_0123 = file("${params.refDir}/GCF_002837175.3_ASM283717v5_genomic.fna.0123")
    ref_amb = file("${params.refDir}/GCF_002837175.3_ASM283717v5_genomic.fna.amb")
    ref_ann = file("${params.refDir}/GCF_002837175.3_ASM283717v5_genomic.fna.ann")
    ref_bwt = file("${params.refDir}/GCF_002837175.3_ASM283717v5_genomic.fna.bwt.2bit.64")
    ref_pac = file("${params.refDir}/GCF_002837175.3_ASM283717v5_genomic.fna.pac")

    // Run alignment
    bwamem2(reads_ch, refGenome, ref_0123, ref_amb, ref_ann, ref_bwt, ref_pac)
}
