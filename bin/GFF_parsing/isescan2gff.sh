BEGIN { FS = ","; OFS = "\t" }
NR == 1 { next }   # skip header
{
    # Determine strand from column 18
    strand = $18

    # Build ID: ID=<family>_NR_<cluster>_<line>
    id = "ID=" $3 "_" NR

    # Print GFF3 fields:
    # 1: seqid  = pre "_" $1
    # 2: source = ISESCAN
    # 3: type   = IS
    # 4: start  = $4
    # 5: end    = $5
    # 6: score  = "."
    # 7: strand = $18
    # 8: phase  = 0
    # 9: attributes = id
    print pre"_"$1, "ISESCAN", "IS", $4, $5, ".", strand, "0", id

}