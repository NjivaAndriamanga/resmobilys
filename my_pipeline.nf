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
Reads trimming and filtering with fastp: length < 50, headcrop and tailcrop score 20
*/
process clean_reads {

    cpus 5
    publishDir "trimming_output/"

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

//With flye
process assemble_genome { 
    cpus 5
    publishDir "genome_assembly/"
    errorStrategy 'ignore' //ignore flye error du to coverage. 

    input:
    val barID
    path fastq
    
    output:
    val barID, emit: id
    path "${barID}_assembly.fasta", emit: assembly
    path "${barID}_assembly_info.fasta"

    script:
    """
    flye --nano-hq $fastq --threads ${task.cpus} --genome-size 4.2m --threads ${task.cpus} -o .
    mv assembly.fasta ${barID}_assembly.fasta
    mv assembly_info.txt ${barID}_assembly_info.fasta
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
    (id1, fastq) = gzip_fastq(identified_samples.out.flatten())
    (id2, fastq_cleaned) = clean_reads(id1, fastq)
    assemble_genome(id2,fastq_cleaned)
    after_assemble(assemble_genome.out.id, assemble_genome.out.assembly).view()
    
}