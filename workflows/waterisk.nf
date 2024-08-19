
log.info "\n"
log.info "========================================================================================="
log.info "Welcome to the WATERISK pipeline. For any questions or remarks, please contact the author \n"
log.info "=========================================================================================="
log.info "\n"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS MANAGMENT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

//print help message, supply typical command line usage for the pipeline
if (params.help) {
    log.info paramsHelp("nextflow run waterisk --profile perso")
    exit 0
}

//
validateParameters()

//Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

//Check if the input dir exists

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IDENTIFIED_RAW_SAMPLES } from '../modules/waterisk_modules.nf'
include { IDENTIFIED_SAMPLES} from '../modules/waterisk_modules.nf'
include { MERGE_SEPARATE_FASTQ } from '../modules/waterisk_modules.nf'
include { REMOVE_BARCODES } from '../modules/waterisk_modules.nf'
include { CLEAN_READS } from '../modules/waterisk_modules.nf'
include {COUNT_BP} from '../modules/waterisk_modules.nf'
include {SAMPLE_FASTQ} from '../modules/waterisk_modules.nf'
include { ASSEMBLE_GENOME } from '../modules/waterisk_modules.nf'
include { IDENTIFY_AMR_PLASMID_COMPLETE} from '../modules/waterisk_modules.nf'
include { IDENTIFY_AMR_CHRM_COMPLETE} from '../modules/waterisk_modules.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WATERISK {
    
    //raw output from MINION
    if (params.raw == true){
        IDENTIFIED_RAW_SAMPLES(file(params.fastq_pass_dir), params.fastq_pass_dir)
        (id_fastq, fastq) = MERGE_SEPARATE_FASTQ(IDENTIFIED_SAMPLES.out.flatten())
    }
    else {
        def fastqs = Channel.fromPath(params.fastq_pass_dir + "*.fastq.gz")
        (id_fastq, fastq) = IDENTIFIED_SAMPLES(fastqs)
    }
    
    
    //if remove barcode is enabled
    if (params.remove_barcode == true){
        (id_nobar, fastq_nobar) = REMOVE_BARCODES(id_fastq, fastq)
        CLEAN_READS(id_nobar, fastq_nobar)
    }
    else {
        CLEAN_READS(id_fastq, fastq)
    }
        
    COUNT_BP(CLEAN_READS.out.barID,
                    CLEAN_READS.out.trimmed_fastq)

    
    SAMPLE_FASTQ(COUNT_BP.out.barID,
                    COUNT_BP.out.number_bp,
                        COUNT_BP.out.fastq)

    ASSEMBLE_GENOME(SAMPLE_FASTQ.out.barID,
                        SAMPLE_FASTQ.out.trimmed_fastq)

    
    IDENTIFY_AMR_PLASMID_COMPLETE(ASSEMBLE_GENOME.out.id,
                                    ASSEMBLE_GENOME.out.assembly_pls_fasta)
    
    IDENTIFY_AMR_CHRM_COMPLETE(ASSEMBLE_GENOME.out.id,
    ASSEMBLE_GENOME.out.assembly_chr_fasta)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
