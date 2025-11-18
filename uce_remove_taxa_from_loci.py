# 本脚本用于从 uce loci 中移除指定物种的序列

import os
from Bio import SeqIO

input_dir = "/home/salticidae/disk_computation/zy_chrysillini/5_3811loci.spruceup0.9.seqtool200/raw"   # 输入目录
output_dir = "/home/salticidae/disk_computation/zy_chrysillini/5_3811loci.spruceup0.9.seqtool200/no-1096-1370"	# 输出目录

# 要剔除的多个物种
remove_taxa = [
    "Tasa_davidi_JXZ1096",
    "Orienticius_chikunii_JXZ1370"
]

os.makedirs(output_dir, exist_ok=True)

remove_taxa_lower = [t.lower() for t in remove_taxa]

for fname in os.listdir(input_dir):
    if not fname.lower().endswith((".fasta", ".fa", ".fas", ".fna")):
        continue

    in_path = os.path.join(input_dir, fname)
    out_path = os.path.join(output_dir, fname)

    records = list(SeqIO.parse(in_path, "fasta"))

    filtered = []
    for r in records:
        header = r.description.lower()
        if not any(taxon in header for taxon in remove_taxa_lower):
            filtered.append(r)

    SeqIO.write(filtered, out_path, "fasta")

print("Done! Multiple taxa removed.")
