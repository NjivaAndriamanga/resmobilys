
log.info "\n"
log.info "========================================================================================="
log.info "Welcome to the WATERISK pipeline. For any questions or remarks, please contact the author \n"
log.info "=========================================================================================="
log.info "\n"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IDENTIFIED_SAMPLES } from '../modules/waterisk_modules.nf'
include { GZIP_FASTQ } from '../modules/waterisk_modules.nf'
include { REMOVE_BARCODES } from '../modules/waterisk_modules.nf'
include { CLEAN_READS } from '../modules/waterisk_modules.nf'
include { ASSEMBLE_GENOME } from '../modules/waterisk_modules.nf'
include { IDENTIFY_AMR_PLASMID} from '../modules/waterisk_modules.nf'
include { IDENTIFY_AMR_CRM} from '../modules/waterisk_modules.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WATERISK {
    def fastq_pass_ch = Channel.fromPath(params.fastq_pass_dir)
    IDENTIFIED_SAMPLES(fastq_pass_ch)

    (id_fastq, fastq) = GZIP_FASTQ(IDENTIFIED_SAMPLES.out.flatten())
    
    (id_nobar, fastq_nobar) = REMOVE_BARCODES(id_fastq, fastq)
    
    CLEAN_READS(id_nobar, fastq_nobar)
    
    ASSEMBLE_GENOME(CLEAN_READS.out.barID,
                    CLEAN_READS.out.trimmed_fastq)
    
    IDENTIFY_AMR_PLASMID(ASSEMBLE_GENOME.out.id,
                        ASSEMBLE_GENOME.out.assembly_pls_fasta)
    
    IDENTIFY_AMR_CRM(ASSEMBLE_GENOME.out.id,
    ASSEMBLE_GENOME.out.assembly_chr_fasta)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
