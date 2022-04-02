#!/usr/bin/env nextflow
/*

========================================================================================
                         BWA Pipeline
========================================================================================
Simple bwa mapping
 #### Authors
 Anthony Underwood @bioinformant
----------------------------------------------------------------------------------------
*/

// Pipeline version
version = '1.0'

/***************** Setup inputs and channels ************************/
// Defaults for configurable variables
params.paired_read_dir = false
params.single_read_dir = false
params.bam_dir = false
params.output_dir = false
params.pattern_match = false
params.help = false


// print help if required
def helpMessage() {
    log.info"""
    =========================================
     Mykrobe Predictor Pipeline v${version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run bwa.nf -with-docker phebioinformatics/phenix --paired_read_dir /path/to/read_dir --pattern_match='*_{1,2}.fastq.gz' --output_dir /path/to/output_dir
    Mandatory arguments:
      --output_dir                       Path to output dir "must be surrounded by quotes"
      --pattern_match                    The regular expression that will match files e.g '*_{1,2}.fastq.gz' or '*.bam'
      --reference_file                   Path to a reference file which reads will be mapped to
    Options:
    One of these must be specified
      --paired_read_dir                  Path to directory containing paired fastq files
      --single_read_dir                  Path to directory containing non-paired fastq files
      --bam_dir                          Path to directory containing bam files
   """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

def check_parameter(params, parameter_name){
   if ( !params[parameter_name]){
      error "You must specifiy a " + parameter_name
   } else {
      variable = params[parameter_name]
      return variable
   }

}

def check_optional_parameters(params, parameter_names){
  if (parameter_names.collect{element -> params[element]}.every{element -> element == false}){
    error "You must specifiy at least one of these options: " + parameter_names.join(', ')
  }
}

// set up reference file
eference_file = file(check_parameter(params, "reference_file"))
reference_name = reference_file.baseName

// set up output directory
output_dir = file(check_parameter(params, "output_dir"))

// set up pattern_match
pattern_match = check_parameter(params, "pattern_match")

//check for read dir directory
check_optional_parameters(params, ['paired_read_dir', 'single_read_dir'])

// setup input channels
if ( params.paired_read_dir) {
    /*
     * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
     * three elements: the pair ID, the first read-pair file and the second read-pair file
     
          */
    fastqs = params.paired_read_dir + '/' + pattern_match
    Channel
      .fromFilePairs( fastqs )
      .ifEmpty { error "Cannot find any reads matching: ${fastqs}" }
      .set { read_pairs }
} else if ( params.single_read_dir ){
    fastqs = params.single_read_dir + '/' + pattern_match
    Channel
      .fromPath( fastqs )
      .ifEmpty { error "Cannot find any bam files matching: ${bams}" }
      .set {reads}
}


log.info "======================================================================"
log.info "                  BWA pipeline"
log.info "======================================================================"
log.info "Running version   : ${version}"
if ( params.paired_read_dir || params.single_read_dir ) {
    log.info "Fastq files             : ${fastqs}"
} else if ( params.bam_dir ) {
    log.info "Bam files             : ${bams}"
}
log.info "======================================================================"
log.info "Outputs written to path '${params.output_dir}'"
log.info "======================================================================"
log.info ""



// Building index for ref fasta
process bwa_index {
     publishDir "${output_dir}"
    input:
        file reference_file

    output:
        file "${reference_name}.*" into bwa_index

    script:
    """
    bwa index ${reference_file} -p ${reference_name}
    """
}


if ( params.paired_read_dir ) {
   process bwa_mem {
       publishDir "${output_dir}"

       tag { id }

       input:
       set id, file(reads) from read_pairs
       file bwa_index

       output:
       set val(id), file("${id}.sam") into sam_file

       script:
       """
       bwa mem -R '@RG\tID:${id}\tSM:null\tLB:null\tCN:null' ${reference_name} ${reads} > ${id}.sam
       """
   }
             }


process sam_to_sorted_bam_and_index {
  publishDir "${output_dir}", mode: 'copy'

  tag {id}

  input:
  set id, file(sam_file) from sam_file

  output:
  set val(id), file("${prefix}.sorted.bam"), file("${prefix}.sorted.bam.bai") into sorted_bam_and_index

  script:
  prefix = sam_file.baseName
  """
  samtools view -bS ${sam_file} | samtools sort - -o ${prefix}.sorted.bam
  samtools index ${prefix}.sorted.bam
  """
}

workflow.onComplete {
        log.info "Nextflow Version:  $workflow.nextflow.version"
  log.info "Command Line:      $workflow.commandLine"
        log.info "Container:         $workflow.container"
        log.info "Duration:         $workflow.duration"
        log.info "Output Directory:  $params.output_dir"
}

