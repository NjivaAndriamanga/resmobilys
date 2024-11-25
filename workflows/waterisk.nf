
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
include { validateParameters; paramsHelp; paramsSummaryLog;  samplesheetToList } from 'plugin/nf-schema'

//Validate input parameters
validateParameters()

//Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

/*
CHECKING PROFILE
*/
if ( (workflow.profile.contains('local') || workflow.profile.contains('slurm') || workflow.profile.contains('perso')) && (workflow.profile.contains('singularity') || workflow.profile.contains('conda')) ) 
    { "executer selected" }
else { exit 1, "No executer selected: executer must be suplied with -profile local/slurm,singularity/conda" }

/*
VALIDATE SAMPLE SHEET INDEX_FILE
*/
ch_input = Channel.fromList(samplesheetToList(params.index_file, "assets/schema_input.json"))
/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEF FUNCTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def check_nonCircularPlasmid(tsv_file) {
    def rows = tsv_file.readLines() // Read the file line by line
    def header = rows[0].split('\t') // Assuming the CSV is comma-separated
    // Iterate through each row (starting from the second, since the first is the header)
    for (def i = 1; i < rows.size(); i++) {
        def row = rows[i].split('\t')
        if (row[4] == 'False' && row[0].contains('plasmid')) {
            return true // Return true if both conditions are satisfied
        }
    }
    return false // Return false if no row matches the conditions
}

def check_plasmidAllCircular(tsv_file) {
    def rows = tsv_file.readLines() // Read the file line by line
    def header = rows[0].split('\t') // Assuming the CSV is comma-separated
    // Iterate through each row (starting from the second, since the first is the header)
    for (def i = 1; i < rows.size(); i++) {
        def row = rows[i].split('\t')
        if (row[4] == 'False' && row[0].contains('plasmid')) {
            return false // Return true if both conditions are satisfied
        }
    }
    return true // Return false if no row matches the conditions
}

def write_value = { value -> "test.txt" >> value + "\n" } 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { TEST_PROCESS }            from '../modules/waterisk_modules.nf'
include { DOWNLOAD_DATABASE }       from '../modules/waterisk_modules.nf'
include { IDENTIFIED_RAW_SAMPLES }  from '../modules/waterisk_modules.nf'
include { IDENTIFIED_SAMPLES}       from '../modules/waterisk_modules.nf'
include { MERGE_SEPARATE_FASTQ }    from '../modules/waterisk_modules.nf'
include { CLEAN_LONG_READS }        from '../modules/waterisk_modules.nf'
include { ASSEMBLE_GENOME }         from '../modules/waterisk_modules.nf'
include { FILTER_CIRCULAR_PLASMID } from '../modules/waterisk_modules.nf'
include { IDENTIFY_AMR_PLASMID }    from '../modules/waterisk_modules.nf'
include { IDENTIFY_AMR_CHRM }       from '../modules/waterisk_modules.nf'
include { PLASME_COMPLETE }         from '../modules/waterisk_modules.nf'
include { PLASME_INCOMPLETE }       from '../modules/waterisk_modules.nf'
include { ALIGN_READS_PLASMID }     from '../modules/waterisk_modules.nf'
include { ASSEMBLY_PLASMID }        from '../modules/waterisk_modules.nf'
include { ASSEMBLY_CHRM }           from '../modules/waterisk_modules.nf'
include { BUSCO }                   from '../modules/waterisk_modules.nf'
include { CHANGE_PLASMID_NAME }     from '../modules/waterisk_modules.nf'
include { MERGE_PLASMID }           from '../modules/waterisk_modules.nf'
include { MOB_TYPER }               from '../modules/waterisk_modules.nf'
include { MERGE_TYPE }              from '../modules/waterisk_modules.nf'
include { CREATE_TAXA }             from '../modules/waterisk_modules.nf'
include { MERGE_TAXA }              from '../modules/waterisk_modules.nf'
include { MOB_CLUSTER }             from '../modules/waterisk_modules.nf' 
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WATERISK {

    //download database
    DOWNLOAD_DATABASE().view()
    TEST_PROCESS().view()

    /* if (params.raw == true){
        IDENTIFIED_RAW_SAMPLES(file(params.long_reads_dir), params.long_reads_dir)
        (fastq) = MERGE_SEPARATE_FASTQ(IDENTIFIED_SAMPLES.out.flatten())
    }
    */

    /*
    IDENTIFIED_SAMPLES(ch_input)
    fastq_long_reads_ch = IDENTIFIED_SAMPLES.out.long_reads
    genome_size_ch = IDENTIFIED_SAMPLES.out.genome_size
    fastq_short_reads_ch = IDENTIFIED_SAMPLES.out.short_reads

    //Remove barcode then trim reads
    if (params.remove_barcode == true){
        (fastq_nobar) = REMOVE_BARCODES(fastq_long_reads_ch)
        CLEAN_LONG_READS(fastq_nobar)
    }
    else {
        CLEAN_LONG_READS(fastq_long_reads_ch)
    }

    //De novo assembly using Hybracter
    assembly_input_ch = CLEAN_LONG_READS.out.trimmed_fastq.join(genome_size_ch)
    assembly_input_ch = assembly_input_ch.join(fastq_short_reads_ch)
    ASSEMBLE_GENOME(assembly_input_ch)

    //Filtering complete where all plasmid is circular (1), complete but with non-circular plasmid (2) and incomplete (3)
    complete_assembly_ch = ASSEMBLE_GENOME.out.complete_assembly

    complete_circular_ch = ASSEMBLE_GENOME.out.complete_assembly //1 Will be directly analyzed
        .filter{ barID, fastq, contig_stats, plassembler, chromosome, plasmids -> 
            check_plasmidAllCircular(contig_stats)}
    complete_circular_chrm_ch = complete_circular_ch.map{ barID, fastq, contig, plassembler, chromosome, plasmid -> [barID, chromosome]}
    complete_circular_plasmid_ch = complete_circular_ch.map{ barID, fastq, contig, plassembler, chromosome, plasmid -> [barID, plasmid]}
    

    complete_non_circular_ch = complete_assembly_ch //2
        .filter{ barID, fastq, contig_stats, plassembler, chromosome, plasmids -> 
            check_nonCircularPlasmid(contig_stats)}
    
    incomplete_assembly_ch = ASSEMBLE_GENOME.out.incomplete_assembly //3
    
    //Complete assembly with non circular plasmid
    putitative_plasmid_ch = complete_non_circular_ch.map{ barID, fastq, contig_stats, plassembler, chromosome, plasmids -> 
                                                                                                                [ barID, contig_stats, chromosome,plasmids]}
    FILTER_CIRCULAR_PLASMID(putitative_plasmid_ch)
    PLASME_COMPLETE(FILTER_CIRCULAR_PLASMID.out, DOWNLOAD_DATABASE.out)
    complete_chrm_ch = PLASME_COMPLETE.out.map{ barID, chromosome, plasmid -> [barID, chromosome]}
    complete_plasmid_ch = PLASME_COMPLETE.out.map{ barID, chromosome, plasmid -> [barID, plasmid]}
    
    //Incomplete assembly, align reads and filter plasmid reads for incomplete assembly
    PLASME_INCOMPLETE(ASSEMBLE_GENOME.out.incomplete_assembly, DOWNLOAD_DATABASE.out)
    ALIGN_READS_PLASMID(PLASME_INCOMPLETE.out.inferred_plasmid)
    ASSEMBLY_PLASMID(ALIGN_READS_PLASMID.out.plasmid_reads)
    ASSEMBLY_CHRM(ALIGN_READS_PLASMID.out.chrm_reads)
    incomplete_plasmid_ch = ASSEMBLY_PLASMID.out
    incomplete_chrm_ch = ASSEMBLY_CHRM.out
    
    //AMR detection
    chrm_amr_ch = complete_circular_chrm_ch.concat(complete_chrm_ch).concat(incomplete_chrm_ch)
    plasmid_amr_ch = complete_circular_plasmid_ch.concat(complete_plasmid_ch).concat(incomplete_plasmid_ch) 

    IDENTIFY_AMR_PLASMID( plasmid_amr_ch )
    IDENTIFY_AMR_CHRM( chrm_amr_ch)

    //BUSCO
    BUSCO( chrm_amr_ch)

    // Plasmid typing and clustering
    //CREATE_MERGE_FILE()
    CHANGE_PLASMID_NAME( plasmid_amr_ch)
    MOB_TYPER(CHANGE_PLASMID_NAME.out)

    plasmid_to_merge = CHANGE_PLASMID_NAME.out.map{barID, plasmid -> plasmid }.collectFile()
    MERGE_PLASMID(plasmid_to_merge)

    type_to_merge = MOB_TYPER.out.map{barID, type -> type }.collectFile()
    MERGE_TYPE(type_to_merge)

    CREATE_TAXA(MOB_TYPER.out)
    taxa_to_merge = CREATE_TAXA.out.collectFile()
    MERGE_TAXA(taxa_to_merge)

    MOB_CLUSTER(MERGE_TAXA.out, MERGE_PLASMID.out, MERGE_TYPE.out)
    */
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
