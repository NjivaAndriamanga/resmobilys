BEGIN {FS=OFS="\t"}
{
  strand = ($8 < $9) ? "+" : "-";
  start = ($8 < $9) ? $8 : $9;
  end = ($8 < $9) ? $9 : 10;
  print $3, "BLAST", "match", start, end, $13, strand, ".", \
        "sample="pre";ID=" $2 ";Identity=" $4 ";Evalue=" $12 ";Length=" $5 ";Description=" $1;
}