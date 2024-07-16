
/*
List all barcodes from miniON output
*/
process IDENTIFIED_SAMPLES {
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

/*
Merge all seprates fastq.gz for each barcodes file into one file
*/
process GZIP_FASTQ {
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
Run only fastqc
*/
process FASTQC {
    input:
    val barID
    path fastq

    output: 
    val barID
    path fastq

    when:
    params.only_qc == true

    script:
    """
        //TODO
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
    if ( params.remove_barcode == true)
        """
        porechop --i $fastq -o trimmed_${barID}.fastq.gz -t ${task.cpus}
        """
    else
        """
        cp ${fastq} trimmed_${barID}.fastq.gz
        """
}

/*
Reads trimming and filtering with fastp: length < 50, headcrop and tailcrop score 20
NB: fastp can be ram greedy. 
*/
process CLEAN_READS {
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
}

/*
Assembling genome and identify plasmid with hybracter
Hybracter also compare putative plasmid with PLSDB using MASH (see plassember_summary.tsv)
*/
process ASSEMBLE_GENOME { 
    label 'process_high'
    publishDir "${params.output_dir}genome_assembly/"
    errorStrategy 'ignore' //ignore if the assembly failed (particularly with fly when depth and coverage are not enought and chromosome are not circularized)

    input:
    val barID
    path fastq
    
    output:
    val barID, emit: id
    path "${barID}_sample_per_contig_stats.tsv", emit: assembly_summary
    path "${barID}_sample_chromosome.fasta", emit: assembly_chr_fasta
    path "${barID}_sample_plasmid.fasta", emit: assembly_pls_fasta
    path "${barID}_plassembler_summary.tsv", emit: plassember_summary

    script:
    """
    hybracter long-single -l $fastq -t ${task.cpus} --skip_qc --min_length ${params.read_min_length} --flyeModel --nano-hq -c ${params.c_size}
    mv hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
    mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
    mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_sample_plasmid.fasta
    mv hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ${barID}_plassembler_summary.tsv
    """

    stub:
    """
    mkdir hybracter_out
    mkdir hybracter_out/FINAL_OUTPUT
    mkdir hybracter_out/FINAL_OUTPUT/complete
    touch hybracter_out/FINAL_OUTPUT/complete/sample_summary.tsv
    touch hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta
    touch hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta
    
    mv hybracter_out/FINAL_OUTPUT/complete/sample_summary.tsv ${barID}_sample_per_contig_stats.tsv
    mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
    mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_sample_plasmid.fasta
    """
}

/*
Identify AMR gene on plasmid and chromosome using abricate
*/
process IDENTIFY_AMR_PLASMID {
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

process IDENTIFY_AMR_CHRM {
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