# This code was used to plot distribution of L1 insertions across B42 genome (supplementary fig. 10D)

import gzip
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

vcf_path = "B42tumor.vcf.gz"

positions = []

with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split('\t')
        chrom, pos, info = cols[0], int(cols[1]), cols[7]
        
        if "INS0" in cols[2] and "FAM_N=L1" in info:
            positions.append((chrom, pos))

positions.sort()

rank = []
distances = []

prev_pos = None
prev_chr = None
counter = 1

for chrom, pos in positions:
    if chrom == prev_chr and prev_pos is not None:
        distances.append(pos - prev_pos)
        rank.append(counter)
        counter += 1
    prev_chr = chrom
    prev_pos = pos

log_distances = np.log10([d if d > 0 else 1 for d in distances])

plt.figure(figsize=(5, 4))
plt.scatter(rank, log_distances, s=10, alpha=0.7)
plt.xlabel("L1 insertion rank (by genomic position)")
plt.ylabel("log10(distance to previous L1 insertion)")
plt.title("B42tumor L1 position Distribution")
plt.grid(False)
plt.tight_layout()
plt.show()