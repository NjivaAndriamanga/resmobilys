BEGIN {
    OFS = "\t"
    print "##gff-version 3"
}

NR > 1 {
    if (match($2, /[a-zA-Z]+[0-9]+/)) {
        seqid = substr($2, RSTART, RLENGTH)
    } else {
        seqid = $2
    }

    attr = "ID=" $6 ";Note=" $9 ";DRUG=" $7
    print seqid, "RGI", "CDS", $3, $4, ".", $5, "0", attr
}
