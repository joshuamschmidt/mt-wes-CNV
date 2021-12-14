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
    .set { cramCountsInChannel }

process cramCounts {
    publishDir "$params.outdir/CoverageCounts", pattern: "*.summary.txt"

    label 'bamTasks'

    input:
    set sample_id, path(input_cram), path(input_crai) from cramCountsInChannel

    output:
    file "*.cpt.bed.gz" into estimateMTcnInChannel

    script:
    """
    cp $reference_fasta_index .
    hts_nim_tools count-reads \
    --fasta $reference_fasta \
    --mapq 25 \
    --threads $task.cpus \
    $regions_bed $input_cram \
    | sort -k1,1 -k2,2n \
    | gzip > "$sample_id".cpt.bed.gz
    """
}

process estimateMTcn {
    publishDir "$params.outdir/CombinedCov/", pattern: "*MT-DNA-CN*"

    label 'combineTasks'

    input:
    file input_files from estimateMTcnInChannel.collect()

    output:
    file "${batch}.coverage.bed.gz"

    script:
    """
    relativeMT_cn.py $input_files \
    --suffix ".regions.bed.gz" \
    | gzip > "$batch".MT-DNA-CN.txt.gz
    """
}
