BEGIN {
    print "##gff-version 3"
}

/^QUERY:/ {
    if (candidate != "") {
        # Print last candidate before switching QUERY
        print seqid "\tTNFINDER\ttransposon\t" start "\t" end "\t.\t.\t.\tID=" candidate
        candidate = ""
    }
    split($2, a, "=")
    seqid = a[1]
    c = 0
    next
}

/^\*Candidate/ {
    if (candidate != "") {
        # Print previous candidate before starting a new one
        print seqid "\tTNFINDER\ttransposon\t" start "\t" end "\t.\t.\t.\tID=" candidate
    }
    candidate = "Candidate_" ++c
    start = end = ""
    next
}

/^[0-9]+\.\.[0-9]+/ {
    match($1, /([0-9]+)\.\.([0-9]+)/, coords)
    s = coords[1] + 0
    e = coords[2] + 0
    # Adjust for reverse coordinates
    if (s > e) {
        tmp = s; s = e; e = tmp
    }
    if (start == "" || s < start) start = s
    if (end == "" || e > end) end = e
    next
}

/^$/ {
    if (candidate != "") {
        print seqid "\tTNFINDER\ttransposon\t" start "\t" end "\t.\t.\t.\tID=" candidate
        candidate = ""
    }
}

END {
    if (candidate != "") {
        print seqid "\tTNFINDER\ttransposon\t" start "\t" end "\t.\t.\t.\tID=" candidate
    }
}
