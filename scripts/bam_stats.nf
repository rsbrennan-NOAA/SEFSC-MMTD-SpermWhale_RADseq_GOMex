#!/usr/bin/env nextflow

// Define input parameters
params.input_dir = 'aligned'
params.output_file = 'count.aligned.txt'
params.output_dir = '.'

// Define the process to calculate alignment statistics
process BAMSTATS {
    input:
    file bam

    output:
    stdout

    script:
    """
    module load bio/samtools/1.19

    samtools flagstat ${bam} > tmp.${bam.baseName}.txt

    TOTAL=\$(grep "in total" tmp.${bam.baseName}.txt | cut -f 1 -d " ")
    MAPPED=\$(grep "mapped (" tmp.${bam.baseName}.txt | cut -f 1 -d " " | head -n 1)
    PERCENT=\$(grep "mapped (" tmp.${bam.baseName}.txt | cut -f 2 -d "(" | cut -f 1 -d "%" | head -n 1)
    MAPPED_q=\$(samtools view -F 4 -q 20 ${bam} | wc -l)
    PERCENT_q=\$(echo "scale=2 ; \$MAPPED_q / \$TOTAL" | bc)

    echo "${bam.baseName},\$TOTAL,\$MAPPED,\$PERCENT,\$MAPPED_q,\$PERCENT_q"
    rm tmp.${bam.baseName}.txt
    """
}

workflow {
    // Create a channel from input BAM files
    bam_files = Channel.fromPath("${params.input_dir}/*.bam")
    
    println "Input directory: ${params.input_dir}"
    bam_files.count().view { count ->
        println "Number of BAM files found: $count"
    }

    // Run BAMSTATS
    alignment_stats = BAMSTATS(bam_files)

    // Collect and process results
    alignment_stats.collectFile(
        name: "${params.output_dir}/${params.output_file}",
        seed: "sample,total_reads,total_mapped,map_percent,mapped_q20,map_percent_q20",
        newLine: false
    )

    println "Output file: ${params.output_file}"
}
