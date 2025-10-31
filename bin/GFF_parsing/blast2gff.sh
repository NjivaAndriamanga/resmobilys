BEGIN {FS=OFS="\t"}
{
  strand = ($8 < $9) ? "+" : "-";
  start = ($8 < $9) ? $8 : $9;
  end = ($8 < $9) ? $9 : 10;
  print pre"_"$2, "BLAST", "MGE", start, end, $13, strand, ".", \
        "ID=" $1 ";Identity=" $4 ";Evalue=" $12 ";Length=" $5;
}
