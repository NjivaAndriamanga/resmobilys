BEGIN {
    FS=OFS="\t";
    print "##gff-version 3"
}
NR > 1 {
    seqid = $1;
    start = $4;
    end = $5;
    id = "prophage_" NR-1;
    note = "key_proteins=" $6 ";best_hit=" $7 ";CDS_count=" $8;
    print pre"_"seqid, "DBSCAN-SWA", "prophage", start, end, ".", ".", ".", "ID=" id ";" note
}
