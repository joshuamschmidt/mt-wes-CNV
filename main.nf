#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.inputFile = 'inputFile.txt'
params.batch = 'batch01'
params.reference_fasta = 'ref.fasta'
params.reference_fasta_index = 'ref.fasta.fai'
params.regions_bed = 'regions.bed'
batch=params.batch
reference_fasta=params.reference_fasta
reference_fasta_index=params.reference_fasta_index
regions_bed=params.regions_bed

Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample_id, file(row.input_cram), file(row.input_crai)) }
    .set { coverageInChannel }

process cramCoverage {
    publishDir "$params.outdir/CoverageSummary", pattern: "*.summary.txt"

    label 'bamTasks'

    input:
    set sample_id, file(input_cram), file(input_crai) from coverageInChannel

    output:
    file "*regions.bed.gz" into coverageOutChannel
    file "${sample_id}.mosdepth.summary.txt" into coverageSummaryChannel

    script:
    """
    cp $reference_fasta_index .
    mosdepth --fasta $reference_fasta \
    --by $regions_bed \
    --no-per-base \
    --mapq 25 \
    --threads $task.cpus \
    $sample_id \
    $input_cram
    """
}
