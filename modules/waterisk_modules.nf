
/*
List all barcodes contain in the input directory from miniON output
*/
process IDENTIFIED_RAW_SAMPLES {
    input:
    path fastq_dir
    val fastq_path

    output:
    path '*.txt'

    script:
    """
    ls ${fastq_dir}/ > path_list.txt
    while read -r line; do basename=\$(echo \$line | awk '{n=split(\$0,A,"/"); print A[n]}'); output_file="\${basename}.txt"; echo "${fastq_path}\$line" > \$output_file; done < path_list.txt
    rm path_list.txt
    """
}

process IDENTIFIED_SAMPLES {
    input:
    path fastq

    output:
    val barID
    path fastq

    script:
    barID = fastq.getSimpleName()
    """

    """
}

/*
Merge all seprates fastq.gz for each barcodes file into one file
*/
process MERGE_SEPARATE_FASTQ {
    input:
    path barcode_dir

    output:
    val barID
    path "${barID}.fastq.gz"

    script:
    barID = barcode_dir.getSimpleName()
    """
        while read -r line; do cat \${line}/*.fastq.gz > ${barID}.fastq.gz; done < $barcode_dir
    """
}

/*
Remove barcode identifiers (or not processed with guppy during demultiplexing)
*/
process REMOVE_BARCODES {
    label 'process_high'

    input:
    val barID
    path fastq

    output:
    val barID
    path "trimmed_${barID}.fastq.gz"

    script:
    """
    porechop --i $fastq -o trimmed_${barID}.fastq.gz -t ${task.cpus}
    """
}

/*
Reads trimming by length and quality score and filtering with fastp. FASTP 0.23.4 is not reproductible and some errors have not yet been fixed
NB: fastp can be ram greedy. 
*/
/* process CLEAN_READS {
    label "process_high"
    publishDir "${params.output_dir}trimmed_output/"
    maxRetries 5
    maxForks 2
    errorStrategy { if (task.attempt <= maxRetries) { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry'} else { return 'finish'} }
    
    input:
    val barID
    path query

    output:
    val barID, emit: barID
    path "${barID}_trimmed.fastq.gz", emit: trimmed_fastq
    path "${barID}_trimmming.html" //to save trimming report in publishDir

    script:
    """
    fastp -i $query -o ${barID}_trimmed.fastq.gz --thread ${task.cpus} --trim_front1 ${params.trim_end_size} --trim_tail1 ${params.trim_end_size} \
    -Q -A --cut_tail --cut_tail_window_size 5 --cut_tail_mean_quality 20 --cut_front --cut_front_window_size 5 \
    --cut_front_mean_quality 20 --length_required ${params.read_min_length} --html ${barID}_trimmming.html
    """
} */

/*
Reads trimming by length and quality score and filtering with cutadapt. Asses reads quality before and reads filtering with fastqc. The two reports are merged with multiqc
*/
process CLEAN_READS {
    label "process_high"
    publishDir "${params.output_dir}trimmed_output/"
    maxRetries 5
    
    input:
    val barID
    path query

    output:
    val barID, emit: barID
    path "${barID}Trimmed.fastq.gz", emit: trimmed_fastq
    path "${barID}_report.html" //to save trimming report in publishDir

    script:
    """
    fastqc --memory 2000 $query -t 1
    cutadapt --cut ${params.trim_end_size} --cut -${params.trim_end_size} -q 20,20 -o ${barID}Trimmed.fastq.gz $query -m ${params.read_min_length}
    fastqc --memory 2000 ${barID}Trimmed.fastq.gz
    multiqc .
    mv multiqc_report.html ${barID}_report.html
    """
}

/*
Count the number of base
*/
process COUNT_BP {
    
    input: 
    val barID
    path trimmed_fastq

    output:
    val barID, emit: barID
    env nbp, emit: number_bp
    path trimmed_fastq, emit: fastq

    script:
    """
    gzip -dc $trimmed_fastq > unzip.fastq
    nbp=\$(awk 'NR % 4 == 0' ORS="" unzip.fastq | wc -m)
    rm unzip.fastq
    """
}

/*
Remove the worst reads until only 500 Mbp remain (100x coverage), useful for very large read sets. If the input read set is less than 500 Mbp, this setting will have no effect.
*/
process SAMPLE_FASTQ {
    publishDir "${params.output_dir}trimmed_output/"

    input:
    val barID
    val nbp
    path query

    output:
    val barID, emit: barID
    path "${barID}sampleTrimmed.fastq.gz", emit: trimmed_fastq

    script:
    """
    if [ ${nbp} -gt 50000000]; then
        filtlong --min_length ${params.read_min_length} --target_bases 500000000 $query | gzip > "${barID}sampleTrimmed.fastq.gz"
    else
        cp $query "${barID}sampleTrimmed.fastq.gz"
    fi
    """
}

/*
Assembling genome and plasmid with hybracter
Hybracter also compare putative plasmid with PLSDB using MASH (see plassember_summary.tsv)
*/
process ASSEMBLE_GENOME { 
    label 'process_high'
    publishDir "${params.output_dir}genome_assembly/"

    input:
    val barID
    path fastq
    
    output:
    val barID, emit: id
    
    //For complete assembly
    path "${barID}_sample_per_contig_stats.tsv", emit: assembly_summary, optional: true
    path "${barID}_plassembler_summary.tsv", emit: plassember_summary, optional: true
    path "${barID}_sample_chromosome.fasta", emit: assembly_chr_fasta, optional: true
    path "${barID}_sample_plasmid.fasta", emit: assembly_pls_fasta, optional: true
    //For incomplete assembly. No contig size above the -c (minimal chrm size).
    path "${barID}_sample.fasta", emit: sample_fasta, optional: true

    script:
    if (params.medaka == true)
        """
        hybracter long-single -l $fastq -t ${task.cpus} --skip_qc --min_length ${params.read_min_length} --flyeModel --nano-hq -c ${params.c_size}

        [ ! -f hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ] || mv hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ${barID}_plassembler_summary.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_sample_plasmid.fasta

        [ ! -f hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ] || mv hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ${barID}_sample.fasta
        """
    if (params.medaka == false)
        """
        hybracter long-single -l $fastq -t ${task.cpus} --skip_qc --no_medaka --min_length ${params.read_min_length} --flyeModel --nano-hq -c ${params.c_size}
        
        [ ! -f hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ] || mv hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ${barID}_plassembler_summary.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_sample_plasmid.fasta

        [ ! -f hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ] || && mv hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ${barID}_sample.fasta
        """

}

/*
Identify AMR gene on plasmid and chromosome using abricate
*/
process IDENTIFY_AMR_PLASMID_COMPLETE {
    publishDir "${params.output_dir}plasmid_amr/"

    input:
    val barID
    path plasmid_fasta

    output:
    path "${barID}_plasmid_amr.txt", emit: plasmid_amr

    script:
    """
    abricate -db ${params.amr_db} ${plasmid_fasta} > ${barID}_plasmid_amr.txt
    """
}

process IDENTIFY_AMR_CHRM_COMPLETE {
    publishDir "${params.output_dir}chrm_amr/"

    input:
    val barID
    path chrm_fasta

    output:
    path "${barID}_chrm_amr.txt", emit: chrm_amr
    
    script:
    """
    abricate -db ${params.amr_db} ${chrm_fasta} > ${barID}_chrm_amr.txt
    """
}