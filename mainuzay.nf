#!/usr/bin/env nextflow
params.reads  = "/archive/binokayl/Uzay/nftest/data/*_{1,2}.fastq.gz"
params.fasta  = "/archive/binokayl/Gencode/GRCh38.primary_assembly.genome.fa"
params.gtf    = "/archive/binokayl/Gencode/gencode.v48.primary_assembly.basic.annotation.gtf"
params.outdir = "/archive/binokayl/Uzay/nftest/results_2"

process TrimGalore {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}_1_trim.fastq.gz"),
          path("${sample_id}_2_trim.fastq.gz"),
          path("fastqc_before/*_fastqc.*"),
          path("fastqc_after/*_fastqc.*")
    publishDir "${params.outdir}/TrimGalore", mode:'copy'
    script:
    """
    mkdir -p fastqc_before fastqc_after trim

    fastqc --outdir fastqc_before --threads 4 ${reads[0]} ${reads[1]}

    trim_galore --fastqc --paired --cores 4 --length 20 --quality 20 \
        --output_dir trim \
        ${reads[0]} ${reads[1]}

    mv trim/*val_1.fq.gz ${sample_id}_1_trim.fastq.gz
    mv trim/*val_2.fq.gz ${sample_id}_2_trim.fastq.gz

    mv trim/*_fastqc.zip fastqc_after/ || true
    mv trim/*_fastqc.html fastqc_after/ || true
    """
}

process IndexGenome {
  input:
    tuple path(fasta), path(gtf)

    output:
    path "star_index"
    publishDir "${params.outdir}/IndexGenome", mode:'copy'
    script:
    """
    mkdir -p star_index
    STAR --runMode genomeGenerate \
         --runThreadN 8 \
         --genomeDir star_index \
         --genomeFastaFiles $fasta \
         --sjdbGTFfile $gtf \
         --sjdbOverhang 100
    """
}

process AlignReads {
    tag "$sample_id"

    memory '72 GB'
    cpus 16

input:
  tuple path(star_index),
        val(sample_id),
        path(read1),
        path(read2),
        path("fastqc_before/*_fastqc.*"),
        path("fastqc_after/*_fastqc.*")


    output:
    path "${sample_id}.sorted.bam"
    publishDir "${params.outdir}/Alignment", mode:'copy'
    script:
    """
    STAR --runThreadN 4 \
         --genomeDir ${star_index} \
         --readFilesIn ${read1} ${read2} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${sample_id}.

    mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    """
}

process RNAseqAnalysis {
    publishDir "${params.outdir}/RNAseqAnalysis", mode: 'copy'
   
    output:
    path "*.tsv"
    path "*.png"
    path "*.pdf"
    path "*.csv"

    script:
    """
    Rscript /archive/binokayl/Uzay/nftest/counts_and_tests.R
    """
}

workflow {

    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    fasta_ch = Channel.value(params.fasta)
    gtf_ch   = Channel.value(params.gtf)

    fasta_gtf_ch = Channel.of([file(params.fasta), file(params.gtf)])

    trimmed_reads_ch = reads_ch | TrimGalore
    indexed_ch       = fasta_gtf_ch | IndexGenome

    indexed_ch.combine(trimmed_reads_ch) | AlignReads

    RNAseqAnalysis()
}

