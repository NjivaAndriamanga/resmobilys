
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
Remove barcodes if there is one left (or not processed with guppy during demultiplexing)
*/
process remove_barcodes {
    label 'many_cpus'

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
    val barID, emit: barID
    path "${barID}_trimmed.fastq.gz", emit: trimmed_fastq
    path "${barID}_trimmming.html" //to save trimming report in publishDir

    script:
    """
    fastp -i $query -o ${barID}_trimmed.fastq.gz --thread ${task.cpus} --trim_front1 10 --trim_tail1 10 \
    -Q --cut_tail --cut_tail_window_size 5 --cut_tail_mean_quality 20 --cut_front --cut_front_window_size 5 \
    --cut_front_mean_quality 20 --length_required 50 --html ${barID}_trimmming.html
    """
}

/*
Assembling genome and identify plasmid with hybracter
Hybracter also compare putative plasmid with PLSDD using MASH (see plassember_summary.tsv)
Remarks: contig with size >= 500kb are considered as chromosome. If you want to change the lower-bound chrm length, modify -c parameters
*/
process assemble_genome { 
    label 'many_cpus'
    publishDir "genome_assembly/"
    errorStrategy 'ignore' //impossible assembly (impossible to construct big contigs)

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
    hybracter long-single -l $fastq -t ${task.cpus} --skip_qc --min_length 50 --flyeModel --nano-hq -c 500000
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
Identify AMR gene on plasmid (chromosome ?)
*/
process identify_AMR {
    input:
    val barID
    path fasta

    output:
    stdout

    script:
    """
    echo "$fasta"
    """

    stub:
    """
    echo "$barID"
    echo "$fasta"

    """
}

workflow {
    println "Welcome to the Waterisk pipeline. For any questions or remarks, please contact the author \n"

    def fastq_pass_ch = Channel.fromPath(params.fastq_pass_dir)
    identified_samples(fastq_pass_ch)
    (id_fastq, fastq) = gzip_fastq(identified_samples.out.flatten())
    (id_nobar, fastq_nobar) = remove_barcodes(id_fastq, fastq)
    clean_reads(id_nobar, fastq_nobar)
    assemble_genome(clean_reads.out.barID,clean_reads.out.trimmed_fastq)
    identify_AMR(assemble_genome.out.id,assemble_genome.out.assembly_chr_fasta).view()
}