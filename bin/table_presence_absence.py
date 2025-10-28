from collections import defaultdict
import re
import argparse

parser = argparse.ArgumentParser(
    description="Build ARG presence/absence matrix with plasmid cluster info"
)
parser.add_argument("--amr", required=True, help="RGI AMR results file (TSV)")
parser.add_argument("--plasmids", required=True, help="Plasmid info file with clusters (TSV)")
args = parser.parse_args()

# ----- 1. Parse plasmid cluster ----
file_path = args.plasmids

plasmid_dict = defaultdict(dict)

with open(file_path, "r") as f:
    header = f.readline().strip().split("\t")  # read header
    col_index = {name: i for i, name in enumerate(header)}  # map column names to indices

    for line in f:
        line = line.strip()
        if not line:
            continue  # skip empty lines

        cols = line.split("\t")
        full_id = cols[col_index["sample_id"]]
        cluster = cols[col_index["primary_cluster_id"]]

        # Extract sample and plasmid ID using regex
        match = re.match(r"(.+?)_(plasmid\d+)_", full_id)
        if not match:
            continue
        sample, plasmid_id = match.groups()

        plasmid_dict[sample][plasmid_id] = cluster

# Convert defaultdict to normal dict
plasmid_dict = {k: dict(v) for k, v in plasmid_dict.items()}


# ---- 2. Parse RGI file ----
rgi_file = args.amr
output_file = "presence_absence_with_clusters.tsv"

sample_to_args = defaultdict(lambda: defaultdict(set))
all_args = set()

with open(rgi_file, "r") as f:
    for line in f:
        if not line.strip():
            continue
        cols = line.strip().split("\t")
        if len(cols) < 4:
            continue

        sample = cols[0]
        contig = cols[1].lower()
        arg = cols[3].strip()

        # Determine location (chromosome/plasmid + number)
        if "plasmid" in contig:
            match = re.search(r"plasmid0*([0-9]+)", contig)
            plasmid_num = match.group(1) if match else "?"
            plasmid_id = f"plasmid{plasmid_num.zfill(5)}"

            # Fetch cluster if available
            cluster = plasmid_dict.get(sample, {}).get(plasmid_id)
            if cluster:
                location = f"P{plasmid_num}({cluster})"
            else:
                location = f"P{plasmid_num}"
        elif "chromosome" in contig:
            match = re.search(r"chromosome0*([0-9]+)", contig)
            chrom_num = match.group(1) if match else "?"
            location = f"C{chrom_num}"
        else:
            location = "?"

        sample_to_args[sample][arg].add(location)
        all_args.add(arg)

# ---- 3. Write table ----
arg_list = sorted(all_args)

with open(output_file, "w") as out:
    out.write("sample\t" + "\t".join(arg_list) + "\n")

    for sample in sorted(sample_to_args.keys()):
        row = []
        for arg in arg_list:
            if arg in sample_to_args[sample]:
                row.append(",".join(sorted(sample_to_args[sample][arg])))
            else:
                row.append("0")
        out.write(sample + "\t" + "\t".join(row) + "\n")

print(f"✅ Table written to {output_file}")