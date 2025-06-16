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
/A transposon/ {
    split($3, coords, /\.\./)
    start = coords[1]
    strand = $4
    name = $6
}

/B transposon/ {
    split($3, coords, /\.\./)
    end = coords[2]
    print seqid, "tcomp", "transposable_element", start, end, ".", strand, ".", "ID=" id ";Name=" name
    id++
}
