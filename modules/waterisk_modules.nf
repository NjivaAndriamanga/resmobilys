
/*
This process will download the plasme database from github and unzip it in the same directory as the main script (main.nf)
Run PLasme.py script to unzip and install the dabatase (avoid conflict when accessing the database during plasme process)
Alternative: download the database from plasme and unzip it waterisk directory with the name DB

Database for kraken and for virulence factor detection
But if the directory DB already exist, it will not be re-downloaded

*/
process DOWNLOAD_PLASME_DATABASE {
    cache true
    label 'plasme'

    output:
    env output

    script:

    
    log.info "Downloading plasme database..."
    """
    cd ${projectDir}
    if [ ! -d DB ]; then 
        echo "Download plasme db in tar format"
    else
        output=" Plamse DB already exist"
    fi
    """

}

process DOWNLOAD_PLATON_DATABASE {
    cache true
    
    output:
    env output

    script:
    """
    cd ${projectDir}
    wget https://zenodo.org/record/4066768/files/db.tar.gz
    tar -xzf db.tar.gz
    rm db.tar.gz
    """
}

process DOWNLOAD_KRAKEN_DATABASE {
    cache true

    output:
    env output

    script:
    if (params.kraken_db.equals("${projectDir}/k2_standard_08gb_20240904") && params.kraken_taxonomy == true) {
        log.info "Downloading kraken database..."
        """
        cd ${projectDir}
        if [ ! -d k2_standard_08gb_20240904 ]; then
            mkdir k2_standard_08gb_20240904
            wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz
            tar -xvf k2_standard_08gb_20240904.tar.gz -C k2_standard_08gb_20240904/
            output="Kraken DB OK"
        else
            output=" Kraken DB already exist "
        fi
        """
    }
    else if (params.kraken_taxonomy == true && params.kraken_db != "${projectDir}/k2_standard_08gb_20240904" && params.kraken_db != null) {
        """
        output="Kraken DB are already provided"
        """
    }
}

process DOWNLOAD_VF_DATABASE {
    cache true

    output:
    env output

    script:
    log.info "Downloading VF database..."
    """
    cd ${projectDir}
    if [ ! -d VF_db ]; then
        mkdir VF_db
        wget https://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
        gzip -d VFDB_setB_nt.fas.gz
        mv VFDB_setB_nt.fas VF_db/
        cd VF_db
        makeblastdb -in VFDB_setB_nt.fas -dbtype nucl
        output="VF DB OK"
    else
        output="VF DB already exist"
    fi
    """
}

process PREPROCESSING_DBSCANDB_CHMOD {
    cache true

    output:
    env output

    script:
    log.info "Preparing DBSCAN tools and database..."
    """
    
    cd ${projectDir}
   
    chmod u+x -R bin/DBSCAN-SWA/bin
    chmod u+x -R bin/DBSCAN-SWA/software
    export PATH=$PATH:${projectDir}/bin/DBSCAN-SWA/software/blast+/bin
    export PATH=$PATH:${projectDir}/bin/DBSCAN-SWA/bin
    export PATH=$PATH:${projectDir}/bin/DBSCAN-SWA/software/diamond
    output="DBSCAN OK. "

    chmod u+x -R bin/PLASMe
    chmod u+x -R bin/GFF_parsing
    chmod u+x bin/tncomp_finder/TnComp_finder.py
    chmod u+x bin/tn3-ta_finder/Tn3+TA_finder.py

    cd ${projectDir}/bin/DBSCAN-SWA
    if [ ! -d db ]; then
        wget https://zenodo.org/records/10404224/files/db.tar.gz
        tar -xvf db.tar.gz
        output+="db for DBSCAN OK. "
    else
        output+="db for DBSCAN already exist. "
    fi
    
    """
}

process DOWNLOAD_RGI_DATABASE {
    cache true

    output:
    path "card.json"

    script:
    log.info "Downloading RGI database..."
    """
    wget https://card.mcmaster.ca/latest/data
    tar -xvf data ./card.json
    """
}

/*
DBSCAN output an exit status 1 when no prophage is found for some fasta. This message is ignored for now and an empty file is created
When fasta id doesn't start with chromosome, a "bac" column is added. This is removed by the awk script.
*/
process DBSCAN {
    tag "${barID}_${type}"
    label "process_high"

    input:
    tuple val(barID) ,path(fasta), val(type)
    val x

    output:
    tuple val(barID) ,path("${barID}_${type}_DBSCAN.txt")

    script:
    """
    touch ${barID}_${type}_DBSCAN.txt
    if [[ -s ${fasta} ]]; then
        set +e
        python ${projectDir}/bin/DBSCAN-SWA/bin/dbscan-swa.py --thread_num ${task.cpus} --input ${fasta} --output dbscan_output
        if [[ \$? -ne 0 ]]; then
            echo "DBSCAN script failed because 0 prophages were found, but continuing..."
        else
            mv dbscan_output/bac_DBSCAN-SWA_prophage_summary.txt DBSCAN.txt
            awk 'BEGIN {FS=OFS="\t"} NF==11 { \$1=""; sub(/^\t/, ""); print; next } { print }' DBSCAN.txt > ${barID}_${type}_DBSCAN.txt
        fi
        set -e
    fi
    """
}

process DBSCAN2GFF {
    tag "${barID}_${type}"
    label "process_single"
    
    input:
    tuple val(barID), path(dbscan)
    
    output:
    tuple val(barID), path("${id}.gff3")
    
    script:
    id = dbscan.getSimpleName()
    """
    awk -f ${projectDir}/bin/GFF_parsing/dbscan2gff.sh ${dbscan} > ${id}.gff3
    """
    
}

//TN3 script is outdated. Need to be updated (Temporary solution for the moment with sed)
process TNFINDER_CORRECTION {
    label "process_single"

    output:
    env output

    script:
    """
    sed -i '/Bio\\.Alphabet/d' ${projectDir}/bin/tn3-ta_finder/Tn3+TA_finder.py
    sed -i '/Bio\\.Alphabet/d' ${projectDir}/bin/tncomp_finder/TnComp_finder.py
    sed -i 's/gbk_record = SeqRecord(Seq(sequence, IUPAC\\.unambiguous_dna))/gbk_record = SeqRecord(Seq(sequence))/' ${projectDir}/bin/tncomp_finder/TnComp_finder.py
    sed -i '/whole_seq = get_sequence(/a\\
            extended_seq = ""' ${projectDir}/bin/tncomp_finder/TnComp_finder.py
    output="TNFINDER CORRECTION OK"
    """
}

process TN3_FINDER {
    tag "${barID}_${type}"
    cache true
    label "tnfinder","process_high"

    input:
    tuple val(barID) ,path(fasta), val(type)
    val x

    output:
    tuple val(barID) ,path("${barID}_${type}_tn3.txt")

    script:
    id = fasta.getSimpleName()
    """
    python3 ${projectDir}/bin/tn3-ta_finder/Tn3+TA_finder.py -f ${fasta} -o ${barID}_tn3 -t ${task.cpus}
    if [ -f ${barID}_tn3/${id}.txt ]; then
        mv ${barID}_tn3/${id}.txt ${barID}_${type}_tn3.txt
    else
        touch ${barID}_${type}_tn3.txt
    fi
    """
}

/*
Du to an output issue, multifasta file need to be split in single fasta files then the resuslt are merged
*/
process TNCOMP_FINDER {
    tag "${barID}_${type}"
    label "tnfinder","process_high"

    input:
    tuple val(barID) ,path(fasta), val(type)
    val x

    output:
    tuple val(barID) ,path("${barID}_${type}_tncomp_final.txt")

    script:
    id = fasta.getSimpleName()
    """
    seqkit split -i -f ${fasta}
    for f in ${fasta}.split; do
        python3 ${projectDir}/bin/tncomp_finder/TnComp_finder.py -f \${f} -o tncomp -p ${task.cpus}
        if find tncomp -name "*composite.txt" | grep -q .; then
            cat tncomp/*composite.txt > ${barID}_${type}_tncomp.txt
        else
            touch ${barID}_${type}_tncomp.txt
        fi
        cat ${barID}_${type}_tncomp.txt >> ${barID}_${type}_tncomp_final.txt
        rm ${barID}_${type}_tncomp.txt
        rm -r tncomp
    done
    """ 
}

process TNFINDER2GFF {
    tag "${barID}"
    label "tnfinder"

    input:
    tuple val(barID) ,path(tn_output)

    output:
    tuple val(barID), path("${id}_tnfinder.gff")

    script:
    id = tn_output.getSimpleName()
    """
    awk -f ${projectDir}/bin/GFF_parsing/tnfindershort.sh ${tn_output} > tnshort.txt
    awk -f ${projectDir}/bin/GFF_parsing/tn2gff.sh tnshort.txt > ${id}_tnfinder.gff
    """
}

process TNFINDERCOMP2GFF {
    tag "${barID}"
    label "tnfinder"

    input:
    tuple val(barID) ,path(tcomp_output)

    output:
    tuple val(barID), path("${id}_tnfindercomp.gff")

    script:
    id = tcomp_output.getSimpleName()
    """
    awk -f ${projectDir}/bin/GFF_parsing/tncomp_short.sh ${tcomp_output} > tcompshort.txt
    awk -f ${projectDir}/bin/GFF_parsing/tcomp2gff.sh tcompshort.txt > ${id}_tnfindercomp.gff
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

/*
    Identified samples from index_files and check the presence of short reads
*/
process IDENTIFIED_SAMPLES {
    cache true

    input:
    tuple path(fastq), val(genome_size), path(sr1), path(sr2)

    output:
    tuple val(barID), path(fastq), emit: long_reads
    tuple val(barID), val(genome_size), emit: genome_size
    tuple val(barID), path(sr1), path(sr2), emit: short_reads
    
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
Long reads trimming by length and quality score and filtering with cutadapt. Asses reads quality before and reads filtering with fastqc. The two reports are merged with multiqc
*/
process CLEAN_LONG_READS {
    tag "${barID}"
    cache true
    label "process_high"
    
    input:
    tuple val(barID), path(query)

    output:
    tuple val(barID), path("${barID}Trimmed.fastq.gz"), emit: trimmed_fastq

    script:
    """
    cutadapt --cut ${params.trim_end_size} --cut -${params.trim_end_size} -q ${params.quality_trim},${params.quality_trim} -o ${barID}Trimmed.fastq.gz $query -m ${params.read_min_length}
    """
    /* """
    fastqc --memory 2000 $query -t ${task.cpus}
    cutadapt --cut ${params.trim_end_size} --cut -${params.trim_end_size} -q ${params.quality_trim},${params.quality_trim} -o ${barID}Trimmed.fastq.gz $query -m ${params.read_min_length}
    fastqc --memory 2000 ${barID}Trimmed.fastq.gz -t ${task.cpus}
    multiqc .
    mv multiqc_report.html ${barID}_report.html
    """ */
}

/*
Assembling genome and plasmid with hybracter
Hybracter also compare putative plasmid with PLSDB using MASH (see plassember_summary.tsv)
For incomplete assembly, contigs are written in sample_final.fasta
*/
process ASSEMBLE_GENOME {
    tag "${barID}"
    cache true
    label 'process_high'
    cpus { task.attempt < 2 ? task.cpus : 1 } //If blastx in dnaapler doesn't found hit fot certain seq length, there is a segmentation fault (temporary fix: reduce cpus to 1)
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    tuple val(barID), path(fastq), val(genome_size), path(sr1), path(sr2)
    
    output:
    tuple val(barID), path(fastq),path("${barID}_sample_per_contig_stats.tsv"), path("${barID}_plassembler_summary.tsv"), path("${barID}_sample_chromosome.fasta"), path("${barID}_hybracter_plasmid.fasta"), optional: true, emit: complete_assembly
    tuple val(barID), path(fastq),path("${barID}_sample_final.fasta"), optional: true, emit: incomplete_assembly

    script:
    def args = " "
    if (sr1 == [] || sr2 == []){ //if no short reads
        args = "long-single -l $fastq -t ${task.cpus} --min_length ${params.read_min_length} --flyeModel ${params.flyeModel}"
    }
    else {
        args = "hybrid-single -l $fastq -1 $sr1 -2 $sr2 -t ${task.cpus} --min_length ${params.read_min_length} --flyeModel ${params.flyeModel}"        
    }

    if (params.medaka == false) {
        args = args + " --no_medaka"
    }
    if(genome_size == 0){
        args = args + " --auto"
    }
    if(genome_size > 0){
        args = args + " -c ${genome_size}"
    }

    """
    hybracter ${args}
    
    [ ! -f hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ] || mv hybracter_out/processing/plassembler/sample/plassembler_summary.tsv ${barID}_plassembler_summary.tsv
    [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
    [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_chromosome.fasta ${barID}_sample_chromosome.fasta 
    [ ! -f hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ] || mv hybracter_out/FINAL_OUTPUT/complete/sample_plasmid.fasta ${barID}_hybracter_plasmid.fasta

    [ ! -f hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ] || mv hybracter_out/FINAL_OUTPUT/incomplete/sample_final.fasta ${barID}_sample_final.fasta
    [ ! -f hybracter_out/FINAL_OUTPUT/incomplete/sample_per_contig_stats.tsv ] || mv hybracter_out/FINAL_OUTPUT/incomplete/sample_per_contig_stats.tsv ${barID}_sample_per_contig_stats.tsv
    """
}

/*
Busco assembly evaluation
*/
process BUSCO {
    tag "${barID}"
    label 'busco'

    input:
    tuple val(barID), path(fasta)

    output:
    path("${barID}_busco.txt")

    script:
    """
    busco -i ${fasta} -m genome -l ${params.lineage_db} -o ${barID}_busco
    cp ${barID}_busco/short*.txt ${barID}_busco.txt
    """
}

/*
Identify AMR gene on plasmid. Abricate with megares to identify heavy and efflux pump. AMR finder plus for ARG.
Resistance to drugs identified with abricate are remove with grep
*/
process ABRICATE {
    tag "${barID}_${type}"
    label 'abricate'

    input:
    tuple val(barID), path(fasta), val(type)

    output:
    tuple val(barID), path(fasta), path("${barID}_${type}_amr.txt"), val(type), emit: amr

    script:
    """
    abricate -db ${params.amr_db} ${fasta} | grep -v "Drugs" > ${barID}_${type}_amr.txt
    """
}

process AMRFINDER {
    tag "${barID}_${type}"
    label 'amrfinder'
    
    input:
    tuple val(barID) , path(fasta), val(type)
    
    output:
    tuple val(barID), path(fasta), path("${barID}_${type}_arg.txt"), val(type)

    script:
    """
    amrfinder --nucleotide ${fasta} > ${barID}_${type}_arg.txt
    """
}

process RGI {
    tag "${barID}_${type}"
    label 'rgi'
    
    input:
    path card_json
    tuple val(barID) , path(fasta), val(type)
    
    output:
    tuple val(barID), path("${barID}_${type}_rgi.txt")
    
    script:
    """
    if [ -s ${fasta} ]; then
        rgi load --card_json ${card_json} --local
        rgi main --input_sequence ${fasta} --output_file rgi --local --clean
        awk -F"\t" '{print "${barID}",\$2, \$3, \$4, \$5, \$9, \$15, \$16, \$17}' OFS="\t" rgi.txt > ${barID}_${type}_rgi.txt
    else
        touch ${barID}_${type}_rgi.txt
    fi
    """ 
}

process RGI2GFF {
    tag "${barID}_${type}"
    
    input:
    tuple val(barID), path(rgitxt)

    output:
    tuple val(barID), path("${id}.gff")

    script:
    id = rgitxt.getSimpleName()
    """
    touch test.txt
    awk -v pre=${barID} -f ${projectDir}/bin/GFF_parsing/rgi2gff.sh $rgitxt > ${id}.gff
    """
}

//Filter circular plasmid in a fasta file from a tab file
process FILTER_CIRCULAR_PLASMID {
    tag "${barID}"

    input:
    tuple val(barID), path(tab_file), path(chromosome),path(putative_plasmid)
    
    output:
    tuple val(barID), path(chromosome), path("${barID}_plasmid.fasta"), path("${barID}_putative_plasmid.fasta")

    script: 
    """
    awk '\$5 == "True" && \$2 == "plasmid" { print \$1 }' ${tab_file} > hybracter_circular_plasmid.txt
    seqkit grep -f hybracter_circular_plasmid.txt ${putative_plasmid} -o ${barID}_plasmid.fasta

    awk '\$5 == "False" && \$2 == "plasmid" { print \$1 }' ${tab_file} > hybracter_non_circular_plasmid.txt
    seqkit grep -f hybracter_non_circular_plasmid.txt ${putative_plasmid} -o ${barID}_putative_plasmid.fasta
    """    
}

process PLASME {
    tag "${barID}"
    label 'plasme'

    input:
    tuple val(barID), path(fastq), path(sample_fasta)
    val x

    output:
    tuple val(barID), path("${barID}_plasme_plasmid.fasta"), val("plasmid"), emit: inferred_plasmid
    tuple val(barID), path("${barID}_plasme_chrm.fasta"), val("chrm"), emit: inferred_chrm

    script:
    """
    python ${projectDir}/bin/PLASMe/PLASMe.py ${sample_fasta} ${barID}_plasme_plasmid.fasta -d ${params.plasme_db}
    awk ' { print \$1 }' ${barID}_plasme_plasmid.fasta_report.csv > chrm_contig.txt
    seqkit grep -v -f chrm_contig.txt ${sample_fasta} -o ${barID}_plasme_chrm.fasta
    """
}

//Infer contig from a fasta file
process PLASME_COMPLETE {
    tag "${barID}"
    label 'plasme'
    errorStrategy "ignore"

    input:
    tuple val(barID), path(chromosome), path(plasmid), path(putative_plasmid)
    val x

    output:
    tuple val(barID), path("${barID}_final_chrm.fasta"), path("${barID}_final_plasmid.fasta")
    
    """
    python ${projectDir}/bin/PLASMe/PLASMe.py ${putative_plasmid} ${barID}_plasme.fasta -d ${params.plasme_db}
    awk ' { print \$1 }' ${barID}_plasme.fasta_report.csv > chrm_contig.txt
    seqkit grep --invert-match -f chrm_contig.txt ${putative_plasmid} -o ${barID}_plasme_chrm.fasta
    cat ${barID}_plasme_chrm.fasta > ${barID}_final_chrm.fasta
    cat ${chromosome} >> ${barID}_final_chrm.fasta
    cat ${barID}_plasme.fasta > ${barID}_final_plasmid.fasta
    cat ${plasmid} >> ${barID}_final_plasmid.fasta
    touch test.txt
    """
}

process PLATON {
    tag "${barID}"
    label "platon","process_high"

    input:
    tuple val(barID), path(fastq), path(sample_fasta)
    val x

    output:
    tuple val(barID), path("${barID}_platon_plasmid.fasta"), val("plasmid"), emit: inferred_plasmid
    tuple val(barID), path("${barID}_platon_chrm.fasta"), val("chrm"), emit: inferred_chrm

    script:
    name = sample_fasta.getSimpleName()
    """
    platon --db ${params.platon_db} --output platon_results ${sample_fasta}
    mv platon_results/${name}.chromosome.fasta ${barID}_platon_chrm.fasta
    mv platon_results/${name}.plasmid.fasta "${barID}_platon_plasmid.fasta
    """
}

process PLATON_COMPLETE {
    tag "${barID}"
    label 'platon','process_high'
    errorStrategy "ignore"
    errorStrategy 'retry'
    maxRetries 1  // retry once if Platon fails

    input:
    tuple val(barID), path(chromosome), path(plasmid), path(putative_plasmid)
    val x

    output:
    tuple val(barID), path("${barID}_final_chrm.fasta"), path("${barID}_final_plasmid.fasta")

    script:
    name = putative_plasmid.getSimpleName()
    if ( task.attempt == 1)
        """
        platon --db ${params.platon_db} --output platon_results ${putative_plasmid}
        mv platon_results/${name}.chromosome.fasta ${barID}_final_chrm.fasta
        mv platon_results/${name}.plasmid.fasta ${barID}_final_plasmid.fasta
        cat ${chromosome} >> ${barID}_final_chrm.fasta
        cat ${plasmid} >> ${barID}_final_plasmid.fasta
        """
    else
        """
        # Retry attempt: skip Platon, just copy input
        cat ${chromosome} >> ${barID}_final_chrm.fasta
        cat ${plasmid} >> ${barID}_final_chrm.fasta
        touch ${barID}_final_plasmid.fasta
        """
}

process PLASME_INCOMPLETE {
    tag "${barID}"
    label 'plasme'

    input:
    tuple val(barID), path(fastq), path(sample_fasta)
    val x

    output:
    tuple val(barID), path("${barID}_plasme_plasmid.fasta"), emit: inferred_plasmid
    tuple val(barID), path("${barID}_plasme_chrm.fasta"), emit: inferred_chrms

    script:
    """
    ${projectDir}/bin/PLASMe/PLASMe.py ${sample_fasta} ${barID}_plasme_plasmid.fasta -d ${params.plasme_db}
    awk ' { print \$1 }' ${barID}_plasme_plasmid.fasta_report.csv > chrm_contig.txt
    seqkit grep -v -f chrm_contig.txt ${sample_fasta} -o ${barID}_plasme_chrm.fasta
    """
}   


//Align and filtered reads on infered plasmid.
process ALIGN_READS_PLASMID {
    tag "${barID}"
    label 'process_high'
    
    input:
    tuple val(barID), path(fastq),path(inferred_chrms_fasta),path(inferred_plasmid_fasta)

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
    samtools fastq mapped_reads.bam > ${barID}_mapped_reads.fastq

    samtools view -b -f 4 aln_sorted.bam > ${barID}_unmapped_reads.bam
    samtools fastq ${barID}_unmapped_reads.bam > ${barID}_unmapped_reads.fastq
    """
}

//Plasmid assembly with unicycler
process ASSEMBLY_PLASMID {
    tag "${barID}"
    label 'process_high'
    publishDir "${params.output_dir}plasme_assembly/"
    errorStrategy "ignore" //When depth is low, assembly is not possible and there is no result

    input:
    tuple val(barID), path(plasmid_reads)

    output:
    tuple val(barID), path("${barID}_plasme_plasmid.fasta")

    script:
    """
    if [ ! -s ${plasmid_reads} ]
    then
        touch ${barID}_plasme_plasmid.fasta
    else
        unicycler -l ${plasmid_reads} -o ${barID}_plasme_plasmid -t ${task.cpus}
        mv ${barID}_plasme_plasmid/assembly.fasta ${barID}.fasta
        awk '/^>/ {print \$0 "_plasmid"; next} {print \$0}' ${barID}.fasta > ${barID}_plasme_plasmid.fasta
    fi
    """
}

//chrm assembly with flye
process ASSEMBLY_CHRM {
    label 'process_high'
    publishDir "${params.output_dir}plasme_assembly/"

    input:
    tuple val(barID), path(chrm_reads)

    output:
    tuple val(barID), path("${barID}_plasme_chrm.fasta")

    script:
    """
    if [ ! -s ${chrm_reads} ]
    then
        touch ${barID}_plasme_chrm.fasta
    else
        flye --nano-hq ${chrm_reads} -t ${task.cpus} -o flye_output
        mv flye_output/assembly.fasta ${barID}.fasta
        awk '/^>/ {print  \$0 "_chromosome"; next} {print \$0}' ${barID}.fasta > ${barID}_plasme_chrm.fasta
    fi
    """
}

/*
Add BarID for each plasmids id
*/
process CHANGE_PLASMID_NAME {
    tag "${barID}"
    label 'process_single'

    input:
    tuple val(barID) ,path(plasmid_fasta)

    output:
    tuple val(barID) ,path("plasmids.fasta")

    script:
    """
    awk '/^>/ {\$0 = ">${barID}_" substr(\$0, 2); gsub(" ", "_", \$0)} 1' ${plasmid_fasta} > plasmids.fasta
    """
}

process MOB_TYPER {
    tag "${barID}"
    label 'mob'

    input:
    tuple val(barID) ,path(plasmid_fasta)

    output:
    tuple val(barID) ,path("plasmid_type.txt"), optional: true

    script:
    """
    mob_typer --multi --infile ${plasmid_fasta} --out_file plasmid_type.txt
    """
}

process MERGE_TYPE {
    label 'merge'

    input:
    path(all_plasmid_type)

    output:
    path("mergeAll_type.txt")

    script:
    """
    echo "sample_id	num_contigs	size	gc	md5	rep_type(s)	rep_type_accession(s)	relaxase_type(s)	relaxase_type_accession(s)	mpf_type	mpf_type_accession(s)	orit_type(s)	orit_accession(s)	predicted_mobility	mash_nearest_neighbor	mash_neighbor_distance	mash_neighbor_identification	primary_cluster_id	secondary_cluster_id	predicted_host_range_overall_rank	predicted_host_range_overall_name	observed_host_range_ncbi_rank	observed_host_range_ncbi_name	reported_host_range_lit_rank	reported_host_range_lit_name	associated_pmid(s)" > mergeAll_type.txt
    sed '/^sample_id/d' ${all_plasmid_type} >> "mergeAll_type.txt"
    """
}


process CREATE_TAXA {
    tag "${barID}"

    input:
    tuple val(barID) ,path(plasmid_type)

    output:
    path "plasmid_tax.txt"

    script:
    """
    tail -n +2 ${plasmid_type} | cut -f1 -d\$'\t' | awk  -F'\t' '{print \$1 "\t${barID}"}' >> plasmid_tax.txt
    """
}

process MERGE_TAXA {
    label 'merge'

    input:
    path(all_plasmid_tax)

    output:
    path("plasmidsAll_sample.txt")

    script:
    """
    echo "sample_id\ttaxonomy" > plasmidsAll_sample.txt
    cat ${all_plasmid_tax} >> plasmidsAll_sample.txt
    """
}

process MOB_CLUSTER {
    label 'mob'

    input:
    path(plasmids_tax)
    path(plasmids_fasta)
    path(plasmids_type)

    output:
    path("clusters.txt")

    script:
    """
    mob_cluster  --mode build -f ${plasmids_fasta} -t ${plasmids_tax} -p ${plasmids_type} -o output
    mv output/clusters.txt clusters.txt
    """
}

process INTEGRON_FINDER {
    tag "${barID}_${type}"
    label 'integron_finder'

    input:
    tuple val(barID) ,path(fasta), val(type)

    output:
    tuple val(barID) ,path("${barID}_${type}_integron.txt"), val(type)

    script:
    file_name = fasta.getSimpleName()
    """
    integron_finder --local-max ${fasta}
    mv Results_Integron_Finder_${file_name}/${file_name}.integrons ${barID}_${type}_integron.txt
    """
}

/*
Format integron_finder output to GFF3. 
*/
process INTEGRON_FORMAT {
    tag "${barID}_${type}"
    label 'integron_finder'

    input:
    tuple val(barID) ,path(integron), val(type)

    output:
    tuple val(barID) ,path("${file_name}_summary.gff"), val(type)

    script:
    file_name = integron.getSimpleName()
    """
    awk '
    BEGIN { 
        OFS="\t"; 
        print "##gff-version 3";  # Print GFF version header first
        }
    \$1 ~ /^integron_/ && \$11 ~ /(complete|CALIN)/ {
        if (!seen[\$1]++) {
            min_pos[\$1] = \$4
            max_pos[\$1] = \$5
            type[\$1] = \$11
            replicon[\$1] = \$2
        } else {
            if (\$4 < min_pos[\$1]) min_pos[\$1] = \$4
            if (\$5 > max_pos[\$1]) max_pos[\$1] = \$5
        }
    }
    END {
        for (id in seen) {
            print replicon[id], "integron_finder", id,  min_pos[id], max_pos[id], ".","+", "0",  type[id]
        }
    }' ${integron} > ${file_name}_summary.gff
    """
}

process KRAKEN {
    tag "${barID}"
    label 'process_high'

    input:
    tuple val(barID) ,path(chromosome_fasta)
    val x

    output:
    tuple val(barID) ,path("${barID}_kraken.txt")

    script:
    """
    kraken2 --db ${params.kraken_db} --threads ${task.cpus} --input ${chromosome_fasta} --use-names --report kraken.txt
    echo "${barID} \n" >> ${barID}_kraken.txt
    cat kraken.txt >> ${barID}_kraken.txt
    echo "\n" >> ${barID}_kraken.txt
    """
}

/*
Blast for virulence factors and keeps track of unique values (gene) in column 2
*/
process VF_BLAST {
    tag "${barID}_${type}"
    
    input:
    tuple val(barID) ,path(fasta), val(type)

    output:
    tuple val(barID) ,path("${sample}_vf_blast.txt")

    script:
    sample = fasta.getSimpleName()
    """
    blastn -db ${params.vf_db} -query ${fasta} -out vf_blast.txt -outfmt '6 stitle qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads ${task.cpus}
    awk '!seen[\$2]++ {print \$0}' vf_blast.txt > ${sample}_vf_blast.txt
    """
}
