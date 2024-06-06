// Define the variable for the fastq_pass directory
params.fastq_pass_dir = 'fastq_pass/'

/*
List all samples
*/
process identified_samples {
    input:
    path query

    output:
    path '*.txt'

    script:
    """
    output_prefix="ID"
    find $PWD/$query/ -mindepth 1 -type d -name bar* > path_list.txt
    sed -i 1d path_list.txt
    while read -r line; do basename=\$(echo \$line | awk '{n=split(\$0,A,"/"); print A[n]}'); output_file="\${basename}.txt"; echo "\$line" > \$output_file; done < path_list.txt
    rm path_list.txt
    """
}

/*
Extract fastq file from MiniION output for each sample (identified with a barcode)
*/
process gzip_fastq {
    input:
    path barcode_dir

    output:
    path 'output.fastq'
    val barID

    script:
    barID = barcode_dir.name
    """
        while read -r line; do gzip -dc \${line}/*.fastq.gz > output.fastq; done < $barcode_dir
    """
}

process clean_reads {
    input:
    path query

    output:
    stdout

    script:
    """
    echo "test"
    echo $query
    """
}

workflow {
    println "Welcome to the Waterisk pipeline. For any questions or remarks, please contact the author \n"

    def fastq_pass_ch = Channel.fromPath(params.fastq_pass_dir)
    identified_samples(fastq_pass_ch)
    (pathtest, id) = gzip_fastq(identified_samples.out.flatten())

    //clean_reads(gzip_fastq.out.flatten()).view()
}