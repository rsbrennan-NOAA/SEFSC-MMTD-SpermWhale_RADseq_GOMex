#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.indir = "trimmed_reads"
params.outdir = "aligned"
params.refGenome = "~/reference_genomes/sperm_whale/GCF_002837175.3_ASM283717v5_genomic.fna"

// Process to align reads
process alignReads {
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(reads)
    path(refGenome)

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
    samtools sort -T /scratch/${sampleId} -o ${sampleId}.bam -

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

    // Reference genome
    refGenome = file(params.refGenome)

    // Run alignment
    alignReads(reads_ch, refGenome)
}
