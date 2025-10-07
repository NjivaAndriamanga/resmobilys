//DEVELOPMENT PROCESS

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
    if [ "\$(wc -l < ${fasta})" -eq 0 ]; then
        touch ${barID}_${type}_tncomp_final.txt
    else
        seqkit split -i -f ${fasta}
        for f in ${fasta}.split/*; do
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
    fi
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
