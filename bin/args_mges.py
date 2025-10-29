#!/usr/bin/env python3
import csv
import re
import argparse

# ----------------------------
# 1. Parse command-line arguments
# ----------------------------
parser = argparse.ArgumentParser(
    description="Check if ARGs are located inside integrons, prophages, or ICEs."
)
parser.add_argument("--args", required=True, help="Path to ARGs input file (GFF-like TSV)")
parser.add_argument("--integrons", required=True, help="Path to integrons input file (GFF-like TSV)")
parser.add_argument("--prophages", required=True, help="Path to prophages input file (GFF-like TSV)")
parser.add_argument("--ices", required=True, help="Path to ICEs input file (GFF-like TSV)")
args = parser.parse_args()

# ----------------------------
# 2. Define helper function
# ----------------------------
def load_mges(file_path):
    """Parse GFF-like MGE file and return dict[(sample, loc_type)] = list of {id, start, end}."""
    mges = {}
    with open(file_path) as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            seq_id, source, feature, start, end, _, strand, _, attributes = cols

            # Extract sample and replicon
            match = re.match(r"(.+?)_(plasmid\d+|chromosome\d+)", seq_id)
            if not match:
                continue
            sample, loc_type = match.groups()

            # Extract element ID
            id_match = re.search(r"ID=([^;]+)", attributes)
            mge_id = id_match.group(1) if id_match else "NA"

            # Save
            mges.setdefault((sample, loc_type), []).append({
                "id": mge_id,
                "start": int(start),
                "end": int(end)
            })
    return mges

# ----------------------------
# 3. Load all MGEs
# ----------------------------
integrons = load_mges(args.integrons)
prophages = load_mges(args.prophages)
ices = load_mges(args.ices)

# ----------------------------
# 4. Parse ARGs and check overlaps
# ----------------------------
output_rows = []
with open(args.args) as f:
    for line in f:
        if line.strip().startswith("#") or not line.strip():
            continue

        cols = line.strip().split("\t")
        if len(cols) < 9:
            continue

        seq_id, source, gene, start, end, score, strand, phase, attributes = cols
        start, end = int(start), int(end)

        match = re.match(r"(.+?)_(plasmid\d+|chromosome\d+)", seq_id)
        if not match:
            continue
        sample, loc_type = match.groups()

        # Check inclusion in each MGE type
        mge_hit = []
        for mge_type, mge_dict in [
            ("Integron", integrons),
            ("Prophage", prophages),
            ("ICE", ices),
        ]:
            for region in mge_dict.get((sample, loc_type), []):
                if start >= region["start"] and end <= region["end"]:
                    mge_hit.append(f"{mge_type}:{region['id']}")

        # Add to output (new column)
        cols.append(",".join(mge_hit) if mge_hit else "NA")
        output_rows.append(cols)

# ----------------------------
# 5. Write output
# ----------------------------
with open("args_mges.tsv", "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["seq_id", "source", "gene", "start", "end", "score", "strand", "phase", "attributes", "MGE_hits"])
    writer.writerows(output_rows)
