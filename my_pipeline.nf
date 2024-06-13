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
    path "${barID}.fastq.gz"
    val barID

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

    publishDir "trimming_output/"

    input:
    path query
    val barID

    output:
    path "${barID}.txt"


    script:
    """
    fastp $query -o ${barID}_output_trimmed.fastq.gz \
    -Q --cut_tail --cut_tail_window_size 5 --cut_tail_mean_quality 20 --cut_front --cut_front_window_size 5 \
    --cut_front_mean_quality 20 --length_required 50 --html ${barID}_trimmming.html
    """

    stub:
    """
    echo "test" > ${barID}.txt

    """
}

workflow {
    println "Welcome to the Waterisk pipeline. For any questions or remarks, please contact the author \n"

    def fastq_pass_ch = Channel.fromPath(params.fastq_pass_dir)
    identified_samples(fastq_pass_ch)
    (pathtest, id) = gzip_fastq(identified_samples.out.flatten())
    clean_reads(pathtest, id)
}