#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads      = "$baseDir/data/*/*_{1,2}.fastq.gz"
params.genomeDir  = "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
params.gtf        = "gencode.v38.annotation.gtf"
params.outdir     = "$baseDir/results"
params.threads    = 8

workflow.onComplete {
    println "✅ Workflow complete. Sorted BAMs in: ${params.outdir}"
}

// Paired-end FASTQ files grouped by sample ID (based on folder and filename)
Channel
    .fromFilePairs(params.reads, flat: false)
    .set { read_pairs }

// STAR alignment and BAM processing
process ALIGN_AND_SORT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}.sorted.bam", emit: sorted_bam

    publishDir params.outdir, mode: 'copy'

    script:
    """
    STAR \
        --genomeDir ${params.genomeDir} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --runThreadN ${params.threads} \
        --outFileNamePrefix ${sample_id}_ \
        --outSAMtype SAM \
        --quantMode TranscriptomeSAM GeneCounts \
        --sjdbGTFfile ${params.gtf} \
        --alignEndsType EndToEnd \
        --alignIntronMax 1000000 \
        --readFilesCommand zcat

    samtools view -@ ${params.threads} -bS ${sample_id}_Aligned.out.sam -o ${sample_id}.bam
    samtools sort -@ ${params.threads} ${sample_id}.bam -o ${sample_id}.sorted.bam
    """
}

workflow {
    read_pairs | ALIGN_AND_SORT
}


STAR --runThreadN 90 --runMode genomeGenerate --genomeDir Gen_Dir_hg38 --genomeFastaFiles GRCh38.p14.genome.fa --sjdbGTFfile gencode.v48.chr_patch_hapl_scaff.basic.annotation.gtf --sjdbOverhang 100
