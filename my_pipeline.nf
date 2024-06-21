// Define the variable for the fastq_pass directory
params.fastq_pass_dir = 'fastq_pass/'

/*
List all barcodes from miniON output
*/
process identified_samples {
    input:
    path fastq_dir

    output:
    path '*.txt'

    script:
    """
    output_prefix="ID"
    find $PWD/$fastq_dir/ -mindepth 1 -type d -name bar* > path_list.txt
    while read -r line; do basename=\$(echo \$line | awk '{n=split(\$0,A,"/"); print A[n]}'); output_file="\${basename}.txt"; echo "\$line" > \$output_file; done < path_list.txt
    rm path_list.txt
    """
}

/*
Merge all seprates fastq.gz for each barcodes file into one file
*/
process gzip_fastq {
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
Remove barcodes if there is one left (or not processed with guppy during multiplexing)
*/
process remove_barcodes {
    cpus 12

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
Reads trimming and filtering with fastp: length < 50, headcrop and tailcrop score 20
*/
process clean_reads {

    cpus 2
    publishDir "trimmed_output/"

    input:
    val barID
    path query

    output:
    val barID
    path "${barID}_trimmed.fastq.gz"
    path "${barID}_trimmming.html" //to save trimming report in publishDir

    script:
    """
    fastp -i $query -o ${barID}_trimmed.fastq.gz --thread ${task.cpus}\
    -Q --cut_tail --cut_tail_window_size 5 --cut_tail_mean_quality 20 --cut_front --cut_front_window_size 5 \
    --cut_front_mean_quality 20 --length_required 50 --html ${barID}_trimmming.html
    """
}

//assembing genome and identify plasmid with hybracter
process assemble_genome { 
    cpus 10
    publishDir "genome_assembly/"

    input:
    val barID
    path fastq
    
    output:
    val barID, emit: id
    path "${barID}_sample_per_contig_stats.tsv", emit: assembly_summary
    path "${barID}_sample_chromosome.fasta", emit: assembly_chr_fasta
    path "${barID}_sample_plasmid.fasta", emit: assembly_pls_fasta

    script:
    """
    hybracter long-single -l $fastq -t ${task.cpus} --skip_qc --min_length 100 --flyeModel --nano-hq
    mv hybracter_out/FINAL_OUTPUT/complete ${barID}_assembly.fasta
    mv assembly_info.txt ${barID}_assembly_info.txt
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

process after_assemble{
    input:
    val barID
    path assembly

    output:
    stdout

    script:
    """
    echo $barID
    echo $assembly
    """
}

workflow {
    println "Welcome to the Waterisk pipeline. For any questions or remarks, please contact the author \n"

    def fastq_pass_ch = Channel.fromPath(params.fastq_pass_dir)
    identified_samples(fastq_pass_ch)
    (id_fastq, fastq) = gzip_fastq(identified_samples.out.flatten())
    (id_nobar, fastq_nobar) = remove_barcodes(id_fastq, fastq)
    (id_fastq_cleaned, fastq_cleaned) = clean_reads(id_nobar, fastq_nobar)
    assemble_genome(id_fastq_cleaned,fastq_cleaned)
    //after_assemble(assemble_genome.out.id, assemble_genome.out.assembly).view()
    
}