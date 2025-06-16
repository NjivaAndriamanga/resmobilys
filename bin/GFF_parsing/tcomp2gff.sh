BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
    print "##gff-version 3"
    id = 1
}
/^QUERY:/ {
    split($2, parts, "=")
    seqid = parts[1]
}
/^[AB] transposon/ {
    split($3, coords, /\.\./)
    start = coords[1]
    end = coords[2]
    strand = $4
    name = $6
    print seqid, "tcomp", "transposable_element", start, end, ".", strand, ".", "ID=" id ";Name=" name
    id++
}
