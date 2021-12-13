#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.inputFile = 'inputFile.txt'
params.batch = 'batch01'
params.reference_fasta = 'ref.fasta'
params.reference_fasta_index = 'ref.fai'
params.target_bed = 'targets.bed'
params.target_cov_txt = 'targets.txt'
params.target_picard = ''
params.bait_bed = 'baits.bed'
params.bait_picard = ''
batch=params.batch
reference_fasta=params.reference_fasta
reference_fasta_index=params.reference_fasta_index
target_bed=params.target_bed
target_cov_txt=params.target_cov_txt
target_picard_list=params.target_picard
bait_bed=params.bait_bed
bait_picard_list=params.bait_picard

/////////////////////////////////////////////////////////
/* --          VALIDATE INPUT FILES                 -- */
/////////////////////////////////////////////////////////

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
    --by $target_bed \
    --no-per-base \
    --mapq 25 \
    --threads $task.cpus \
    $sample_id \
    $input_cram
    """
}
