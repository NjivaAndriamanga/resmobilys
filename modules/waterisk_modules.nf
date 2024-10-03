/*
This process will download the plasme database from github and unzip it in the same directory as the main script (main.nf)
But if the directory DB already exist, it will not be re-downloaded
*/
process DOWNLOAD_DATABASE {
    cpus 8
    label 'plasme'

    output:
    val 1

    script:
    """
    cd ${projectDir}
    if [ ! -d DB ]; then
        python PLASMe_db.py
        unzip DB.zip
    fi

    """
}


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
    tuple val(barID), path(fastq)

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
    tuple val(barID), path("${barID}.fastq.gz")

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
    tuple val(barID), path(fastq)

    output:
    tuple val(barID), path("trimmed_${barID}.fastq.gz")

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
    
    input:
    tuple val(barID), path(query)

    output:
    tuple val(barID), path("${barID}Trimmed.fastq.gz"), emit: trimmed_fastq
    path "${barID}_report.html" //to save trimming report in publishDir

    script:
    """
    fastqc --memory 2000 $query -t 1
    cutadapt --cut ${params.trim_end_size} --cut -${params.trim_end_size} -q ${params.quality_trim},${params.quality_trim} -o ${barID}Trimmed.fastq.gz $query -m ${params.read_min_length}
    fastqc --memory 2000 ${barID}Trimmed.fastq.gz
    multiqc .
    mv multiqc_report.html ${barID}_report.html
    """
}

/*
Remove the worst reads until only 500 Mbp remain (100x coverage for 5M genome size), useful for very large read sets. If the input read set is less than 500 Mbp, this setting will have no effect.
Alternative Rasusa
*/
process SAMPLING_FASTQ {
    debug true
    publishDir "${params.output_dir}trimmed_output/"

    input:
    tuple val(barID),path(query)

    output:
    tuple val(barID), path("${barID}sampleTrimmed.fastq.gz")

    script:
    """
    filtlong --min_length ${params.read_min_length} --target_bases ${params.target_bases} $query | gzip > "${barID}sampleTrimmed.fastq.gz"
    """
}

/*
Assembling genome and plasmid with hybracter
Hybracter also compare putative plasmid with PLSDB using MASH (see plassember_summary.tsv)
For incomplete assembly, contigs are written in sample_final.fasta
*/
process ASSEMBLE_GENOME { 
    //Error are mainly from medaka. Re run failed sample with no medaka
    label 'process_high'
    publishDir "${params.output_dir}genome_assembly/"

    input:
    tuple val(barID), path(fastq)
    
    output:
    tuple val(barID), path(fastq),path("${barID}_sample_per_contig_stats.tsv"), path("${barID}_plassembler_summary.tsv"), path("${barID}_sample_chromosome.fasta"), path("${barID}_sample_plasmid.fasta"), optional: true, emit: complete_assembly
    tuple val(barID), path(fastq),path("${barID}_sample_final.fasta"), optional: true, emit: incomplete_assembly

    script:
    if (params.medaka == true)
        """
        hybracter long-single -l $fastq -t ${task.cpus} --min_length ${params.read_min_length} --flyeModel --nano-hq -c ${params.c_size}
        
        [ ! -f hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ] || mv hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ${barID}_plassembler_summary.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_sample_plasmid.fasta

        [ ! -f hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ] || mv hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ${barID}_sample_final.fasta
        """
    else if (params.medaka == false)
        """
        hybracter long-single -l $fastq -t ${task.cpus} --no_medaka --min_length ${params.read_min_length} --flyeModel --nano-hq -c ${params.c_size}
        
        [ ! -f hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ] || mv hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ${barID}_plassembler_summary.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
        [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_sample_plasmid.fasta

        [ ! -f hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ] || mv hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ${barID}_sample_final.fasta

        """
}

/*
Identify AMR gene on plasmid and chromosome using abricate
*/
process IDENTIFY_AMR_PLASMID {
    //label 'amr_detection'
    publishDir "${params.output_dir}plasmid_amr/"

    input:
    tuple val(barID) ,path(plasmid_fasta)

    output:
    path "${barID}_plasmid_amr.txt", emit: plasmid_amr

    script:
    """
    abricate -db ${params.amr_db} ${plasmid_fasta} > ${barID}_plasmid_amr.txt
    """
}



process IDENTIFY_AMR_CHRM {
    //label 'amr_detection'
    publishDir "${params.output_dir}chrm_amr/"

    input:
    tuple val(barID), path(chrm_fasta)

    output:
    path "${barID}_chrm_amr.txt", emit: chrm_amr
    
    script:
    """
    abricate -db ${params.amr_db} ${chrm_fasta} > ${barID}_chrm_amr.txt
    """
}

//Filter circular in a fasta file from a tab file
process FILTER_NON_CIRCULAR {
    publishDir "${params.output_dir}genome_assembly/"

    input:
    tuple val(barID), path(fastq), path(tab_file), path(contig_fasta)
    
    output:
    tuple val(barID), path(fastq), path("non_circular_plasmid.fasta"),  emit: non_circular_plasmid
    path "${bar_id}_circular_plasmid.fasta"

    script: 
    """
    awk '\$5 == "True" { print \$1 }' barcode09_sample_per_contig_stats.tsv > hybracter_circular_plasmid.txt
    seqkit grep -f hybracter_circular_plasmid.txt barcode09_sample_plasmid.fasta -o ${bar_id}_circular_plasmid.fasta

    awk '\$5 == "False" { print \$1 }' barcode09_sample_per_contig_stats.tsv > hybracter_non_circular_plasmid.txt
    seqkit grep -f hybracter_non_circular_plasmid.txt barcode09_sample_plasmid.fasta -o ${bar_id}_non_circular_plasmid.fasta
    """    
}

//Infer contig from a fasta file
process PLASME {
    label 'plasme'
    //publishDir "${params.output_dir}plasme_output/"

    input:
    tuple val(barID), path(fastq), path(contig_fasta)
    val x

    output:
    tuple val(barID), path("${barID}_plasme.fasta"), path(fastq)

    """
    python plasme.py ${contig_fasta} ${barID}_plasme.fasta -d ${params.plasme_db}

    """
}

//Align and filtered reads on infered plasmid
process ALIGN_READS_PLASMID {

    input:
    tuple val(barID), path(inferred_plasmid_fasta), path(fastq)

    output:
    tuple val(barID), path("${barID}_mapped_reads.fastq") , emit: plasmid_reads
    tuple val(barID), path("${barID}_unmapped_reads.fastq"), emit: chrm_reads

    script:
    """
    minimap2 -ax map-ont ${inferred_plasmid_fasta} ${fastq} > aln.sam
    samtools view -Sb -o aln.bam aln.sam
    samtools sort aln.bam -o aln_sorted.bam
    samtools index aln_sorted.bam
    samtools view -b -F 4 aln_sorted.bam > mapped_reads.bam
    samtools view -b -f 4 aln_sorted.bam > unmapped_reads.bam
    samtools fastq mapped_reads.bam > ${barID}_mapped_reads.fastq
    samtools fastq unmapped_reads.bam > ${barID}_unmapped_reads.fastq
    """
}

//Plasmid assembly with unicycler
process ASSEMBLY_PLASMID {
    label 'process_high'
    publishDir "${params.output_dir}plasme_output/"
    errorStrategy "ignore" //When depth is low, assembly is not possible and there is no result

    input:
    tuple val(barID), path(mapped_reads)

    output:
    tuple val(barID), path("${barID}_plasme_plasmid.fasta")

    script:
    """
    if [ ! -s ${mapped_reads} ]
    then
        touch ${barID}_plasme_plasmid.fasta
    else
        unicycler -l ${mapped_reads} -o ${barID}_plasme_plasmid -t ${task.cpus}
        mv ${barID}_plasme_plasmid/assembly.fasta ${barID}_plasme_plasmid.fasta
    fi
    """
}

//chrm assembly with flye
process ASSEMBLY_CHRM {
    label 'process_high'
    publishDir "${params.output_dir}plasme_output/"

    input:
    tuple val(barID), path(mapped_reads)

    output:
    tuple val(barID), path("${barID}_plasme_chrm.fasta")

    script:
    """
    if [ ! -s ${mapped_reads} ]
    then
        touch ${barID}_plasme_chrm.fasta
    else
        flye --nano-hq ${mapped_reads} -t ${task.cpus} -o flye_output
        mv flye_output/assembly.fasta ${barID}_plasme_chrm.fasta
    fi
    """
} 