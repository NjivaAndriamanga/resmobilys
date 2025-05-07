BEGIN {
    OFS = "\t"
    FS = "\t"
    print "##gff-version 3"
}

NR > 1 {
    if (match($2, /[a-zA-Z]+[0-9]+/)) {
        seqid = substr($2, RSTART, RLENGTH)
    } else {
        seqid = $2
    }

    attr = "FAMILY=" $9 ",CLASS=" $7, ",MECANISM=" $8
    print pre"_"seqid, "RGI", $6, $3, $4, ".", $5, "0", attr
}
