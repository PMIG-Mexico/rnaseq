#!/usr/bin/env nextflow
 
/*
 * Parámetros Default, pueden ser modificados en la línea de comando
 */
params.reads = "/mnt/a/Transdiferenciacion/RawData/120-1-NEG_{1,2}.fastq.gz"
params.genome_index = "/mnt/d/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/"
params.version = "hg38"
params.threads = 1

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }

/*
* Importar los índices
*/
Channel
    .fromPath( params.genome_index,type: 'dir' )
    .set { index_prefix } 

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {    
    input:
    file '*' from index_prefix
    set pair_id, file(reads) from read_pairs
  
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam_files
  
    """
    tophat2 --num-threads 32 Bowtie2Index/genome ${reads}
    """
}
  
/*
 * Step 3. Assembles the transcript by using the "cufflinks"
 * and publish the transcript output files into the `results` folder
 */
process makeTranscript {
    publishDir "results"
     
    input:
    set pair_id, bam_file from bam_files
      
    output:
    set pair_id, 'transcripts.gtf' into transcripts
  
    """
    cufflinks --num-threads 32 ${bam_file}
    """
}