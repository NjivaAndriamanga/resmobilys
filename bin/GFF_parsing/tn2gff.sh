BEGIN {
    OFS="\t";
    print "##gff-version 3";
    id=""
}
NR==1{
    split($0,fields," ")
    seq=fields[2]
    min=10000000
    max=0
}

/Candidate/ {
   gsub(/\*/,"",$0)
   gsub(/ /,"_",$0)

   if (id != $0) {
       id = $0
   }
}

$0 ~ /\.\./ {
   split($0,orf,/\.\./)
   if(orf[1] < min) {
        min=orf[1]
   }
   if(orf[2] < min) {
        min=orf[2]
   }
   if(orf[1] > max) {
        max=orf[1]
   }
   if(orf[2] > max) {
       max=orf[2]
   }
}

NF == 0 {
    print seq, "TNFINDER", "transposon", min, max, ".", ".", ".", "ID=" id
    min=100000000
    max=0
}

END {
   print seq, "TNFINDER", "transposon", min, max, ".", ".", ".", "ID=" id
}
