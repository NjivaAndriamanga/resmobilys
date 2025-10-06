BEGIN {OFS="\t"}
{
  strand = ($9 < $10) ? "+" : "-";
  start = ($9 < $10) ? $9 : $10;
  end = ($9 < $10) ? $10 : $9;
  print $3, "BLAST", "match", start, end, $13, strand, ".", \
        "ID=" $2 ";Identity=" $4 ";Evalue=" $12 ";Length=" $5 ";Description=" $1;
}