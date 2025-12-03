
log.info "\n"
log.info "========================================================================================="
log.info "Welcome to the RESMOBILYS pipeline. For any questions or remarks, please contact the author \n"
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
if ( (workflow.profile.contains('local') || workflow.profile.contains('slurm') || workflow.profile.contains('test')) && (workflow.profile.contains('singularity') || workflow.profile.contains('apptainer')) ) 
    { "executer selected" }
else { exit 1, "No executer selected: executer must be suplied with -profile local/slurm/test,singularity/apptainer/docker" }

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
include { CHECK_PLASME_DATABASE }           from '../modules/waterisk_modules.nf'
include { DOWNLOAD_PLATON_DATABASE}         from '../modules/waterisk_modules.nf'
include { DOWNLOAD_KRAKEN_DATABASE }        from '../modules/waterisk_modules.nf'
include { DOWNLOAD_VF_DATABASE }            from '../modules/waterisk_modules.nf'
include { DOWNLOAD_CONJSCAN }               from '../modules/waterisk_modules.nf'
include { PREPROCESSING_DBSCANDB_CHMOD }    from '../modules/waterisk_modules.nf'
include { DOWNLOAD_RGI_DATABASE }           from '../modules/waterisk_modules.nf'
include { IDENTIFIED_SAMPLES}               from '../modules/waterisk_modules.nf'
include { CLEAN_LONG_READS }                from '../modules/waterisk_modules.nf'
include { ASSEMBLE_GENOME }                 from '../modules/waterisk_modules.nf'
include { FILTER_CIRCULAR_PLASMID }         from '../modules/waterisk_modules.nf'
include { ABRICATE }                        from '../modules/waterisk_modules.nf'
include { RGI}                              from '../modules/waterisk_modules.nf'
include { RGI2GFF}                          from '../modules/waterisk_modules.nf'
include { PLASME_COMPLETE }                 from '../modules/waterisk_modules.nf'
include { PLATON}                           from '../modules/waterisk_modules.nf'
include { PLASME }                          from '../modules/waterisk_modules.nf'
include { PLASME_INCOMPLETE }               from '../modules/waterisk_modules.nf'
include { PLATON_COMPLETE }                 from '../modules/waterisk_modules.nf'
include { BUSCO }                           from '../modules/waterisk_modules.nf'
include { CHANGE_PLASMID_NAME }             from '../modules/waterisk_modules.nf'
include { MOB_TYPER }                       from '../modules/waterisk_modules.nf'
include { MERGE_TYPE }                      from '../modules/waterisk_modules.nf'
include { CREATE_TAXA }                     from '../modules/waterisk_modules.nf'
include { MERGE_TAXA }                      from '../modules/waterisk_modules.nf'
include { MOB_CLUSTER }                     from '../modules/waterisk_modules.nf'
include { INTEGRON_FINDER }                 from '../modules/waterisk_modules.nf'
include { INTEGRON_FORMAT }                 from '../modules/waterisk_modules.nf'
include { KRAKEN }                          from '../modules/waterisk_modules.nf'
include { VF_BLAST }                        from '../modules/waterisk_modules.nf'
include { DBSCAN }                          from '../modules/waterisk_modules.nf'
include { DBSCAN2GFF }                      from '../modules/waterisk_modules.nf'
include { ISESCAN }                         from '../modules/waterisk_modules.nf'
include { PLASMID_CLUSTER_ARG }             from '../modules/waterisk_modules.nf'
include { VISUALIZATION_TABLE }             from '../modules/waterisk_modules.nf'
include { ARGS_MGES }                       from '../modules/waterisk_modules.nf'
include { ISESCAN2GFF }                     from '../modules/waterisk_modules.nf'
include { PROKKA }                          from '../modules/waterisk_modules.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RESMOBILYS {

    //download tools and databases
    CHECK_PLASME_DATABASE().view()
    DOWNLOAD_PLATON_DATABASE().view()
    DOWNLOAD_VF_DATABASE().view()
    PREPROCESSING_DBSCANDB_CHMOD().view()
    DOWNLOAD_RGI_DATABASE()
    DOWNLOAD_CONJSCAN().view()

    //Identified samples
    IDENTIFIED_SAMPLES(ch_input)
    fastq_long_reads_ch = IDENTIFIED_SAMPLES.out.long_reads
    genome_size_ch = IDENTIFIED_SAMPLES.out.genome_size
    fastq_short_reads_ch = IDENTIFIED_SAMPLES.out.short_reads

    //Long reads triming
    CLEAN_LONG_READS(fastq_long_reads_ch)

    //De novo assembly using Hybracter
    assembly_input_ch = CLEAN_LONG_READS.out.trimmed_fastq.join(genome_size_ch)
    assembly_input_ch = assembly_input_ch.join(fastq_short_reads_ch)
    ASSEMBLE_GENOME(assembly_input_ch)

    //Filtering complete where all plasmid is circular (1), complete but with non-circular plasmid (2) and incomplete (3)
    complete_assembly_ch = ASSEMBLE_GENOME.out.complete_assembly

    complete_circular_ch = ASSEMBLE_GENOME.out.complete_assembly //1 Will be directly analyzed if all contigs are circular
        .filter{ barID, contig_stats, plassembler, chromosome, plasmids -> 
            check_plasmidAllCircular(contig_stats)}
    complete_circular_chrm_ch = complete_circular_ch.map{ barID, contig, plassembler, chromosome, plasmid -> [barID, chromosome, "chromosome"]}
    complete_circular_plasmid_ch = complete_circular_ch.map{ barID, contig, plassembler, chromosome, plasmid -> [barID, plasmid, "plasmid"]}
    

    complete_non_circular_ch = complete_assembly_ch //2
        .filter{ barID, contig_stats, plassembler, chromosome, plasmids -> 
            check_nonCircularPlasmid(contig_stats)}
    
    incomplete_assembly_ch = ASSEMBLE_GENOME.out.incomplete_assembly //3
    
    //Complete assembly with non circular plasmid
    putitative_plasmid_ch = complete_non_circular_ch.map{ barID, contig_stats, plassembler, chromosome, plasmids -> 
                                                                                                                [ barID, contig_stats, chromosome,plasmids]}
    FILTER_CIRCULAR_PLASMID(putitative_plasmid_ch)

    if(params.plasmid == "plasme") {
        PLASME_COMPLETE(FILTER_CIRCULAR_PLASMID.out, CHECK_PLASME_DATABASE.out)
        plasme_complete_chrm_ch = PLASME_COMPLETE.out.map{ barID, chromosome, plasmid -> [barID, chromosome, "chromosome"]}
        plasme_complete_plasmid_ch = PLASME_COMPLETE.out.map{ barID, chromosome, plasmid -> [barID, plasmid, "plasmid"]}
        
        //For incomplete assembly, use plasme to infer chrm and plasmid contig
        PLASME(incomplete_assembly_ch, CHECK_PLASME_DATABASE.out)

        //AMR detection
        chrm_amr_ch = complete_circular_chrm_ch.concat(plasme_complete_chrm_ch).concat(PLASME.out.inferred_chrm)
        plasmid_amr_ch = complete_circular_plasmid_ch.concat(plasme_complete_plasmid_ch).concat(PLASME.out.inferred_plasmid)
    }
    if(params.plasmid == "platon") {
        PLATON_COMPLETE(FILTER_CIRCULAR_PLASMID.out, CHECK_PLATON_DATABASE.out)
        platon_complete_chrm_ch = PLATON_COMPLETE.out.map{ barID, chromosome, plasmid -> [barID, chromosome, "chromosome"]}
        platon_complete_plasmid_ch = PLATON_COMPLETE.out.map{ barID, chromosome, plasmid -> [barID, plasmid, "plasmid"]}
        
        //For incomplete assembly, use plasme to infer chrm and plasmid contig
        PLATON(incomplete_assembly_ch, DOWNLOAD_PLATON_DATABASE.out)
        
        //AMR detection
        chrm_amr_ch = complete_circular_chrm_ch.concat(platon_complete_chrm_ch).concat(PLATON.out.inferred_chrm)
        plasmid_amr_ch = complete_circular_plasmid_ch.concat(platon_complete_plasmid_ch).concat(PLATON.out.inferred_plasmid)
    }
    
    contig_ch = chrm_amr_ch.concat(plasmid_amr_ch)

    ABRICATE(contig_ch)
    RGI(DOWNLOAD_RGI_DATABASE.out, contig_ch)
    rgi_out = RGI.out
    RGI2GFF(rgi_out)
    
    //ISEScan
    ISESCAN( contig_ch )
    ISESCAN2GFF( ISESCAN.out )

    //Integron_finder
    INTEGRON_FINDER( contig_ch )
    INTEGRON_FORMAT( INTEGRON_FINDER.out)

    //BUSCO
    BUSCO( chrm_amr_ch.map{barID, contig, type -> [barID, contig]} )

    //DBSCAN
    DBSCAN( contig_ch , PREPROCESSING_DBSCANDB_CHMOD.out )

    //DBSCAN2GFF
    DBSCAN2GFF( DBSCAN.out )

    //ICE
    PROKKA( contig_ch)

    // Plasmid typing and clustering
    CHANGE_PLASMID_NAME( plasmid_amr_ch.map{barID, contig, type -> [barID, contig]} )
    MOB_TYPER(CHANGE_PLASMID_NAME.out)

    plasmid_merge = CHANGE_PLASMID_NAME.out.map{barID, plasmid -> plasmid }.collectFile(name:"plasmid_merge.fasta", storeDir:"${params.output_dir}plasmid_annotation/")

    type_to_merge = MOB_TYPER.out.map{barID, type -> type }.collectFile()
    MERGE_TYPE(type_to_merge)

    CREATE_TAXA(MOB_TYPER.out)
    taxa_to_merge = CREATE_TAXA.out.collectFile()
    MERGE_TAXA(taxa_to_merge)

    MOB_CLUSTER(MERGE_TAXA.out, plasmid_merge, MERGE_TYPE.out)

    //KRAKEN
    if (params.kraken_db != "null" && params.kraken_taxonomy == true) {
        DOWNLOAD_KRAKEN_DATABASE().view()
        KRAKEN(chrm_amr_ch.map{barID, contig, type -> [barID, contig]},DOWNLOAD_KRAKEN_DATABASE.out)
        kraken_ch = KRAKEN.out.map{ barID, kraken -> kraken}.collectFile(name:"kraken_summary.txt", storeDir:"${params.output_dir}kraken/")
    }
    
    //Virulence factor and heavy metals detection
    VF_BLAST(DOWNLOAD_VF_DATABASE.out, contig_ch)

    //formating output
    plasmid_args_ch = RGI.out.map { id, file, type -> file }.filter { file ->  file.name.contains('plasmid')}
    merge_plasmid_rgi = plasmid_args_ch.collectFile(name:'merged_plasmid_rgi.txt')
    PLASMID_CLUSTER_ARG(merge_plasmid_rgi, MOB_CLUSTER.out)

    rgi_amr = RGI2GFF.out.map{barID, rgi -> rgi}.collectFile(name:"merge_rgi.tsv")
    VISUALIZATION_TABLE(MOB_CLUSTER.out, rgi_amr)

    prohages_file = DBSCAN2GFF.out.map{ id, file -> file}.collectFile(name:"prophages.gff")
    integrons_file = INTEGRON_FORMAT.out.map{ id, file -> file }.collectFile(name:"integrons.gff")
    is_file = ISESCAN2GFF.out.map{ id, file -> file }.collectFile(name:"is.gff")
    ARGS_MGES(rgi_amr, integrons_file, prohages_file, is_file)
}

//TN3_FINDER( contig_ch, TNFINDER_CORRECTION.out )
//TNCOMP_FINDER( contig_ch , TNFINDER_CORRECTION.out )
//TN2GFF
//TNFINDER2GFF (TN3_FINDER.out)
//TNFINDERCOMP2GFF ( TNCOMP_FINDER.out )
// AMRFINDER( contig_ch )
//Transposan finder
//TNFINDER_CORRECTION()

//Incomplete assembly, align reads and filter plasmid reads for incomplete assembly
        // PLASME_INCOMPLETE(ASSEMBLE_GENOME.out.incomplete_assembly, DOWNLOAD_PLASME_DATABASE.out)
        // ALIGN_READS_PLASMID(PLASME_INCOMPLETE.out.inferred_plasmid)
        // ASSEMBLY_PLASMID(ALIGN_READS_PLASMID.out.plasmid_reads)
        // ASSEMBLY_CHRM(ALIGN_READS_PLASMID.out.chrm_reads)
        // incomplete_plasmid_ch = ASSEMBLY_PLASMID.out
        // incomplete_chrm_ch = ASSEMBLY_CHRM.out

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
