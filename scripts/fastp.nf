#!/usr/bin/env nextflow

params.fastq = false

// Validate inputs
if (!params.fastq) {
    error "specify the path to fq files with --fastq"
}

process fastp {
  publishDir 'analysis/trimmed_files', mode: 'move', pattern: '*.clip.fastq.gz'
  publishDir 'analysis/fastp_logs', mode: 'move', pattern: '*.html'
  tag "${sample_id}"
  cpus 8

  input:
    tuple val(sample_id), path(read)

  output:
    val(sample_id)
    path("${sample_id}.fastp.json"), emit: fastp_json
    path("${sample_id}.fastp.html")
    path("${sample_id}.clip.fastq.gz")
    
  script:
  def threads = task.cpus - 4 // -w just defines worker threads, there are others for reading etc
  """
  fastp -i ${read} \
        -o ${sample_id}.clip.fastq.gz \
        -w ${threads} --cut_right \
        -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json  
  """
}

process multiqc {
  publishDir 'analysis/multiqc_report_trimmed', mode: 'copy'

  input:
    path fastp_reports // note that fastp_reports doesn't correspond to anything. we pass the reports directly in the workflow

  output:
    path "multiqc_report.html"

  """
  multiqc .
  """
}

workflow {
    Channel
        .fromPath(params.fastq, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
        .set { reads_ch }
    ch_fastp = fastp(reads_ch)
    ch_multiqc = multiqc(ch_fastp.fastp_json.flatten().collect())
}


