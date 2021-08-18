#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

def helpMessage() {
    log.info"""
    # Tools | Analysis of gene expression using Tag-Seq (v2.0)
    A pipeline to do transcript counting with Tag-Seq (v2.0) `https://eli-meyer.github.io/TagSeq_utilities`

    ## Examples
    ```bash
    nextflow run -resume quantify_Tagseq.nf \
        --transcriptome "FASTQDIR/*.fastq" \
        --cores 14 \
        --adaptors adaptors.fa \
        --outdir results_TagSeq
    ```
    ## Parameters
    --transcriptome <glob>
        Required
        A glob of the fastq reads to be quantified.
        The basename of the file is used as the sample name.

    --cores <int>
        Optional
        Default: 14

    --outdir <path>
        Default: `results_TagSeq`
        The directory to store the results in.
    
    --adaptors <path>
        Default: `adaptors.fa`
        The path to adaptor sequences.

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}
params.help = false
params.transcriptome = false
params.reference = 'reference.fasta'
params.gtf = 'reference.gtf'
params.gff3 = 'reference.gff3'
params.outdir = "results_TagSeq"
params.adaptors = 'adapters.fa'
params.cores = 14

// All reads were renamed to remove redundant extension _R1_001 and \
// gunzipped in parallel using "parallel --gnu gunzip  ::: *gz" before analysis

println """\
=====================================================================
Analysis of gene expression using Tag-Seq (v2.0):
`https://eli-meyer.github.io/TagSeq_utilities`
=====================================================================
Workflow Information:
---------------------
  Project Directory:          ${workflow.projectDir}
  Work Directory:             ${workflow.launchDir}/work
  Launch Directory:           ${workflow.launchDir}
  Output Directory:           ${workflow.launchDir}/${params.outdir}
Samples:
--------
  Local sample glob:          ${params.transcriptome}
  Adaptor sequence:           ${params.adaptors}
  Reference genome:           ${params.reference}
  Reference GFF3:           ${params.gff3}
  Reference GTF:           ${params.gtf}
  
"""
if ( params.transcriptome) {
    fqreads = Channel
    .fromPath(params.transcriptome, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{readsForQualFilter}
} else {
    log.info "No sequencing reads supplied."
    exit 1
}


// Quality Filter
// First, we exclude low quality reads likely to contain sequencing errors.
process fastq_quality_filter {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fq'

    input:
    set val(sampleID), file('fqreads.fastq') from readsForQualFilter

    output:
    set val(sampleID), file("${sampleID}.qf.fq") into readsForHomopolymerRepeatFilter
    
    script:
    """
    #docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    QualFilterFastq.pl -i fqreads.fastq -m 20 -x 10 -o ${sampleID}.qf.fq
    """
}


//Homopolymer repeat filter

process fastq_homopolymer_repeat_filter {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fq'

    input:
    set val(sampleID), file("${sampleID}.qf.fq") from readsForHomopolymerRepeatFilter

    output:
    set val(sampleID), file("${sampleID}.hf.fq") into readsForAdaptorFilter

    script:
    """
    #docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    HRFilterFastq.pl -i "${sampleID}.qf.fq" -n 30 -o ${sampleID}.hf.fq
    """
}

//Adaptor filter
process fastq_adaptor_filter {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fq'
    //publishDir "${params.outdir}", mode: 'copy', pattern: '*.log'

    input:
    set val(sampleID), file("${sampleID}.hf.fq") from readsForAdaptorFilter

    output:
    set val(sampleID), file("${sampleID}.af.fq") into readsForPCRdereplication

    script:
    """
    #docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    bbduk.sh in="${sampleID}.hf.fq" ref=${params.adaptors} k=12 stats="${sampleID}adaptorTrim.log" out="${sampleID}.af.fq" overwrite=t
    """
}

//Remove PCR Duplicates process is optional.. 

process fastq_PCR_dereplication {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fq'

    input:
    set val(sampleID), file("${sampleID}.af.fq") from readsForPCRdereplication

    output:
    set val(sampleID), file("${sampleID}.pcr.fq") into readsForTrimTag

    script:
    """
    #docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    RemovePCRDups.pl -i "${sampleID}.af.fq" -o "${sampleID}.pcr.fq" -s 1 -e 4 -j 9
    """
}

//Trim tags
process fastq_trim_tags {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fq'

    input:
    set val(sampleID), file("${sampleID}.pcr.fq") from readsForTrimTag

    output:
    set val(sampleID), file("${sampleID}.tt.fq") into readsForAligner

    script:
    """
    #docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    TagTrimmer.pl -i "${sampleID}.pcr.fq" -b 1 -e 8 -o "${sampleID}.tt.fq"
    """
}

//Align reads to reference geome using BWA
process align_reads {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.bam'

    input:
    set val(sampleID), file("${sampleID}.tt.fq") from readsForAligner

    output:
    set val(sampleID), file("${sampleID}.bwa.bam")  into bamForCounting

    script:
    //##gmapper --qv-offset 33 -Q --strata -o 3 -N 1 "${sampleID}.tt.fq" ${params.reference} > '${sampleID}.bwa.bam'
    """#docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    bwa mem -t 8 '${params.reference}' '${sampleID}.tt.fq' | samtools sort -@8 -o '${sampleID}.bwa.bam' -
    samtools index '${sampleID}.bwa.bam'
    """
}

process count_alignments {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.tsv'

    input:
    set val(sampleID), file("${sampleID}.bwa.bam") from bamForCounting

    output:
    file "${sampleID}.ht_seq.tsv" into combineCounts
    // After the test run.. BWA->HTseqCount had the most consistent results

    script:
    //##SAMFilterByGene.pl -i "${sampleID}.bwa.bam" -m 40 -o "${sampleID}.filtered.sam"
    """
    #docker run --user $(id -u):$(id -g) -v $(pwd):/work_dir fmobegi/tagseq_utils:1.0.0-beta \
    htseq-count --format=bam \
        --stranded=no \
        --type=CDS \
        --order=pos \
        --idattr=Name \
        "${sampleID}.bwa.bam" '${params.gff3}' \
        > "${sampleID}.ht_seq.tsv"
    """
}

combineCounts.collect().set{all_counts_tab}
// combineCounts.collect().view()

/* 
Waiting for all samples to be processed before combining them into  single count matrix!
""" */
 
process combine_counts{
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    all_counts_tab

    output:
    file 'abundance.matrix.tsv'

    script:
    """
    CombineExpression.pl all_counts_tab > abundance.matrix.tsv
    """
}


