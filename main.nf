#!/usr/bin/env nextflow 

log.info "\n"
log.info "========================================================================================="
log.info "Welcome to the WATERISK pipeline. For any questions or remarks, please contact the author \n"
log.info "=========================================================================================="
log.info "\n"

log.info """\
Pipeline parameters:\n
fastq input directory    : ${params.fastq_pass_dir}
c_size                   : ${params.c_size}
minimum read lenth       : ${params.read_min_length}
bases to trim at the end : ${params.trim_end_size}
"""

include { IDENTIFIED_SAMPLES } from './modules/waterisk_modules.nf'
include { GZIP_FASTQ } from './modules/waterisk_modules.nf'
include { REMOVE_BARCODES } from './modules/waterisk_modules.nf'
include { CLEAN_READS } from './modules/waterisk_modules.nf'
include { ASSEMBLE_GENOME } from './modules/waterisk_modules.nf'
include { IDENTIFY_AMR} from './modules/waterisk_modules.nf'

workflow {
    def fastq_pass_ch = Channel.fromPath(params.fastq_pass_dir)
    IDENTIFIED_SAMPLES(fastq_pass_ch)
    (id_fastq, fastq) = GZIP_FASTQ(IDENTIFIED_SAMPLES.out.flatten())
    (id_nobar, fastq_nobar) = REMOVE_BARCODES(id_fastq, fastq)
    CLEAN_READS(id_nobar, fastq_nobar)
    ASSEMBLE_GENOME(CLEAN_READS.out.barID,CLEAN_READS.out.trimmed_fastq)
    IDENTIFY_AMR(ASSEMBLE_GENOME.out.id,ASSEMBLE_GENOME.out.assembly_chr_fasta).view()
}
