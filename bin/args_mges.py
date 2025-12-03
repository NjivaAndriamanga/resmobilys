#!/usr/bin/env python3
import csv
import re
import argparse

# ----------------------------
# 1. Parse command-line arguments
# ----------------------------
parser = argparse.ArgumentParser(
    description="Annotate ARGs with Integrons, Prophages, ICEs, ISs, and composite transposons."
)
parser.add_argument("--args", required=True, help="Path to ARGs input file (GFF-like TSV)")
parser.add_argument("--integrons", required=True, help="Path to integrons input file (GFF-like TSV)")
parser.add_argument("--prophages", required=True, help="Path to prophages input file (GFF-like TSV)")
#parser.add_argument("--ices", required=True, help="Path to ICEs input file (GFF-like TSV)")
parser.add_argument("--ises", required=True, help="Path to IS input file (GFF-like TSV)")
args = parser.parse_args()

# ----------------------------
# 2. Helper function: load MGE files
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

            # Extract ID
            id_match = re.search(r"ID=([^;]+)", attributes)
            mge_id = id_match.group(1) if id_match else "NA"

            # Store region
            mges.setdefault((sample, loc_type), []).append({
                "id": mge_id,
                "start": int(start),
                "end": int(end)
            })
    return mges

# ----------------------------
# 3. Load MGEs
# ----------------------------
integrons = load_mges(args.integrons)
prophages = load_mges(args.prophages)
#ices = load_mges(args.ices)
ises = load_mges(args.ises)

# ----------------------------
# 4. Helper: extract IS family
# ----------------------------
def extract_is_family(is_id):
    # IS110_140_18 → "IS110_140"
    return "_".join(is_id.split("_")[:-1])

# ----------------------------
# 5. Process ARGs and add annotations
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

        # ----------------------------
        # Detect MGE hits (Integron, Prophage, ICE, IS)
        # ----------------------------
        mge_hit = []
        for mge_type, mge_dict in [
            ("Integron", integrons),
            ("Prophage", prophages),
            ("IS", ises),
        ]:
            for region in mge_dict.get((sample, loc_type), []):
                if start >= region["start"] and end <= region["end"]:
                    mge_hit.append(f"{mge_type}:{region['id']}")

        # ----------------------------
        # Detect composite transposons
        # ----------------------------
        composite_hit = "NA"

        is_regions = ises.get((sample, loc_type), [])
        if is_regions:
            # Group IS by family
            is_by_family = {}
            for region in is_regions:
                fam = extract_is_family(region["id"])
                is_by_family.setdefault(fam, []).append(region)

            # For each family, check adjacent IS pairs
            for fam, regions in is_by_family.items():
                if len(regions) < 2:
                    continue

                regions_sorted = sorted(regions, key=lambda r: r["start"])

                for i in range(len(regions_sorted) - 1):
                    left = regions_sorted[i]
                    right = regions_sorted[i + 1]

                    distance = abs(right["start"] - left["end"])
                    if distance > 20000:
                        continue

                    # Check if ARG lies strictly inside the IS pair
                    if start >= left["end"] and end <= right["start"]:
                        composite_hit = f"CompositeIS:{fam}:{left['id']},{right['id']}"
                        break

                if composite_hit != "NA":
                    break

        # Add composite IS to MGE hits
        if composite_hit != "NA":
            mge_hit.append(composite_hit)

        # Add final MGE_hits column
        cols.append(",".join(mge_hit) if mge_hit else "NA")
        output_rows.append(cols)

# ----------------------------
# 6. Write output file
# ----------------------------
with open("args_mges.tsv", "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow([
        "seq_id", "source", "gene", "start", "end",
        "score", "strand", "phase", "attributes",
        "MGE_hits"
    ])
    writer.writerows(output_rows)
