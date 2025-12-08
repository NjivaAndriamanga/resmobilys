#!/usr/bin/env python3
import argparse

# ----------------------------
# 1. Parse command-line arguments
# ----------------------------
parser = argparse.ArgumentParser(
    description="Parse GFF3-like files and extract gene positions."
)
parser.add_argument("--gff", required=True, help="Path to gene annotation GFF3 file")
parser.add_argument("--system", required=True, help="Path to systems file")
parser.add_argument("--avg_size", required=True, help="ICE average size estimation")
parser.add_argument("--type", required=True, help="Chromosome or plasmid")
args = parser.parse_args()

if args.type == "plasmid":
    args.avg_size = 0

gene_positions = {}

with open(args.gff) as f:
    for line in f:
        if line.startswith("#"):
            continue  # skip header lines
        
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue  # skip malformed lines
        
        start = int(parts[3])
        end   = int(parts[4])
        
        attributes = parts[8]
        
        # extract ID=...
        gene_id = None
        for attr in attributes.split(";"):
            if attr.startswith("ID="):
                gene_id = attr.split("ID=")[1]
                break
        
        if gene_id:
            gene_positions[gene_id] = (start, end)

def parse_systems(path, gene_positions):
    systems = []
    current_genes = []

    with open(path) as f:
        
        for line in f:
            if line.startswith("#"):
                continue  # skip header lines

            line = line.strip()

            # empty line → end of a system
            if line == "":
                if current_genes:
                    systems.append(current_genes)
                    current_genes = []
                continue

            parts = line.split("\t")

            # col1 always present
            gene_id = parts[1]

            # col4 may override
            if len(parts) > 20 and parts[20].strip() != "":
                gene_id = parts[20].strip()

            current_genes.append(gene_id)

    # add last system if file does not end with newline
    if current_genes:
        systems.append(current_genes)

    # compute genomic coords
    results = []
    for sys_genes in systems:
        starts = []
        ends = []

        for gid in sys_genes:
            if gid not in gene_positions:
                #print(f"WARNING: {gid} not found in gene_positions")
                continue
            s, e = gene_positions[gid]
            starts.append(s)
            ends.append(e)

        if starts and ends:
            ice_start = min(starts)
            current_length = max(ends) - ice_start
            if (current_length) > int(args.avg_size):
                ice_end = max(ends)
            else:
                ice_end = max(ends) + (int(args.avg_size) - current_length )

            results.append({
                "start": ice_start,
                "end": ice_end,
                "n_genes": len(sys_genes)
            })

    return results

systems = parse_systems(args.system, gene_positions)

for i, sys in enumerate(systems, 1):
    if sys['end'] - sys['start'] > 150000 : #limit the size of the ice at 150kb
        continue
    print(f"CONJSCan\tICE\t{sys['start']}\t{sys['end']}\t.\t+\t0\tID=System_{i}")