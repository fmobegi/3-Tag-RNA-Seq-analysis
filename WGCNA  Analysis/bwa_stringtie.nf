#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    # Tools | Mapping reads to reference using hisat2 and quantifying transcripts with stringtie

    ## Examples
    ```bash
    nextflow run -resume hisat_stringtie.nf \
        --transcriptome "FASTQDIR/*.fq" \
        --cores 14 \
        --outdir results_hisat2
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
        Default: `results_hisat2`
        The directory to store the results in.

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

//params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Ecoracana_560_v1.0.fa'
//params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Ecoracana_560_v1.1.gene_exons.gff3'
//..................................E coracana.................................................

//params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Osativa_204_v7.0.fa'
//params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Osativa_204_v7.0.gene_exons.gff3'
//..................................O sativa...................................................

//params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Striga_asiatica_UVA1.fasta'
//params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Striga_asiatica_UVA1.gff3'
//.................................S aciatica..................................................

//params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Sbicolor_454_v3.0.1.fa'
//params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Sbicolor_454_v3.1.1.gene_exons.gff3'
//...................................S bicolor.................................................

//params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Zmays_B73_v5.fa'
//params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Zmays_B73_v5.gff3'
//...................................Z mays.................................................

//params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Brachiaria_ruziziensis.fna'
//params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Brachiaria_ruziziensis.gff3'
//...................................Brachiaria..............................................

params.reference = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Cenchrus_americanus_23D2B1.fna'
params.gff3 = '/data/01_raw_data/FASTQ_Generation_2021-06-10_13_41_28Z-426511272/reference_genomes/Cenchrus_americanus_23D2B1.gff3'
//...................................C.americanus23D2B1.................................................

//params.outdir = "results_Osativa"
//params.outdir = "results_Ecorana"
//params.outdir = "results_Sasiatica"
//params.outdir = "results_Sbicolor"
//params.outdir = "results_Zmays"
//params.outdir = "results_Brachiaria"
params.outdir = "results_Camericanus"

params.cores = 14

// All reads were renamed to remove redundant extension _R1_001 and \
// gunzipped in parallel using "parallel --gnu gunzip  ::: *gz" before analysis

println """\
=====================================================================
Mapping cleaned Tag-Seq reads to new reference using BWA:
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
  Reference genome:           ${params.reference}
  Reference GFF3:           ${params.gff3}
  
"""
if ( params.transcriptome) {
    fqreads = Channel
    .fromPath(params.transcriptome, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap{readsForHisat2}
} else {
    log.info "No reads supplied."
    exit 1
}


// BWA align
process bwa_align {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.sam'
    //publishDir "${params.outdir}", mode: 'copy', pattern: '*.bam'

    input:
    set val(sampleID), file('fqreads.af.fq.gz') from readsForHisat2

    output:
    set val(sampleID), file("${sampleID}.sam") into bamForStringtie
    //set val(sampleID), file("${sampleID}.bam") into bamForStringtie

    script:
    """
    bwa mem '${params.reference}' fqreads.af.fq.gz > '${sampleID}.sam'
    #bwa mem -t 1 '${params.reference}' fqreads.af.fq.gz | \
    #samtools sort -@8 -o '${sampleID}.bam' -
    """
}


//Stringtie
process stringtie_count {
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*'

    input:
    set val(sampleID), file("${sampleID}.bam") from bamForStringtie

    output:
    set val(sampleID), file("${sampleID}.Hisat2.ga"), file("${sampleID}.Hisat2.gtf") into gtfForCounts

    script:
    """
    stringtie \
        -e \
        -v \
        -p ${params.cores} \
        -e \
        -o ${sampleID}.Hisat2.gtf \
        -G ${params.gff3} \
        -A ${sampleID}.Hisat2.ga \
        -l ${sampleID} \
        '${sampleID}.bam'
    """
}

process counts_tab {
    tag "${sampleID}"
    label "stringtie"

    publishDir "${params.outdir}", mode: 'copy', pattern: '*'

    input:
    set val(sampleID), file("${sampleID}.Hisat2.ga"), file("${sampleID}.Hisat2.gtf") from gtfForCounts

    output:
    set val(sampleID), file("${sampleID}.Hisat2.fpkm"), file("${sampleID}.Hisat2.tpm"),file("${sampleID}.Hisat2.raw") into combineCounts
    file "*.raw.pre"

    script:
    """
    awk -F"\\t" '{if (NR!=1) {print \$1, \$8}}' OFS='\\t' ${sampleID}.Hisat2.ga > ${sampleID}.Hisat2.fpkm
    awk -F"\\t" '{if (NR!=1) {print \$1, \$9}}' OFS='\\t' ${sampleID}.Hisat2.ga > ${sampleID}.Hisat2.tpm
    echo -e "${sampleID}\\t./${sampleID}.Hisat2.gtf" > ${sampleID}.gtf_files
    prepDE.py -i ${sampleID}.gtf_files -g ${sampleID}.raw.pre
    grep -v "gene_id" "${sampleID}.raw.pre" | sed 's/,/\t/g' > ${sampleID}.Hisat2.raw
    
    """
}

//combineCounts.basename().collect().set{all_counts_tab}
//combineCounts.collect().view()

/* process combine_counts{
    tag "${sampleID}"
    publishDir "${params.outdir}", mode: 'copy'

    when:
    counts_tab    
    
    input:
    "${params.outdir}/*.raw"
    "${params.outdir}/*.fpkm"
    "${params.outdir}/*.tpm"

    output:
    file '*.raw.txt'
    file '*.FPKM.txt'
    file '*.TPM.txt'

    script:
    """
    create-gem.py --sources ${params.outdir}/ --prefix Osativa_204_v7.countmatrix --type raw
    create-gem.py --sources ${params.outdir} --prefix Osativa_204_v7.countmatrix --type FPKM
    create-gem.py --sources ${params.outdir} --prefix Osativa_204_v7.countmatrix --type TPM
    """
}*/

// Do this process manually inside the ourput directory 

/* create-gem.py --sources . --prefix Ecoracana_560_v1.countmatrix --type raw
create-gem.py --sources . --prefix Ecoracana_560_v1.countmatrix --type FPKM
create-gem.py --sources . --prefix Ecoracana_560_v1.countmatrix --type TPM

create-gem.py --sources . --prefix Osativa_204_v7.countmatrix --type raw
create-gem.py --sources . --prefix Osativa_204_v7.countmatrix --type FPKM
create-gem.py --sources . --prefix Osativa_204_v7.countmatrix --type TPM  

create-gem.py --sources . --prefix SaciaticaUVA1.countmatrix --type raw
create-gem.py --sources . --prefix SaciaticaUVA1.countmatrix --type FPKM
create-gem.py --sources . --prefix SaciaticaUVA1.countmatrix --type TPM

create-gem.py --sources . --prefix Sbicolor.countmatrix --type raw
create-gem.py --sources . --prefix Sbicolor.countmatrix --type FPKM
create-gem.py --sources . --prefix Sbicolor.countmatrix --type TPM

create-gem.py --sources . --prefix Zmays.countmatrix --type raw
create-gem.py --sources . --prefix Zmays.countmatrix --type FPKM
create-gem.py --sources . --prefix Zmayscolor.countmatrix --type TPM

create-gem.py --sources . --prefix B.ruziziensis.countmatrix --type raw
create-gem.py --sources . --prefix B.ruziziensis.countmatrix --type FPKM
create-gem.py --sources . --prefix B.ruziziensis.countmatrix --type TPM

create-gem.py --sources . --prefix C.americanus23D2B1.countmatrix --type raw
create-gem.py --sources . --prefix C.americanus23D2B1.countmatrix --type FPKM
create-gem.py --sources . --prefix C.americanus23D2B1.countmatrix --type TPM

*/

