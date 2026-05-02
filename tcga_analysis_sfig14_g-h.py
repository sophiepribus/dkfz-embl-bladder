import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.stats import spearmanr

tcga_l1s_df = pd.read_csv("tcga_l1_paper_raw_data.csv", header=1)
tcga_ecdna_df = pd.read_csv("tcga_platinum_metadata.tsv", sep='\t')
tcga_l1s_bladder = tcga_l1s_df[tcga_l1s_df["subtype"] == "BLCA"]
tcga_ecdna_bladder = tcga_ecdna_df[tcga_ecdna_df["lineage"] == "Bladder"]
tcga_l1s_bladder["sample_id"] = tcga_l1s_bladder["sample_id"].str[:-1]
tcga_l1_ecdna = pd.merge(tcga_ecdna_bladder, tcga_l1s_bladder, left_on="sample_barcode", right_on="sample_id")
tcga_l1_ecdna = tcga_l1_ecdna[["sample_barcode", "sample_classification", "number_amp_region", "RNA", "WGS", 
                               "RT_burden", "RT_burden_log2", "RT_burden_adj", "RT_burden_log2_adj",
                               'all_loci_TPM', 'all_loci_TPM_log2', 'all_loci_TPM_adj', 'all_loci_TPM_log2_adj', 
                               'active_loci_TPM', 'active_loci_TPM_log2', 'active_loci_TPM_adj', 'active_loci_TPM_log2_adj']]
tcga_l1_ecdna["ecDNA+"] = (tcga_l1_ecdna["sample_classification"] == "Circular")

plot_df = tcga_l1_ecdna.dropna(subset=["active_loci_TPM_log2_adj"])

# split groups
group0 = plot_df.loc[plot_df["ecDNA+"] == False, "active_loci_TPM_log2_adj"]
group1 = plot_df.loc[plot_df["ecDNA+"] == True, "active_loci_TPM_log2_adj"]

# Mann–Whitney U test
stat, pval = mannwhitneyu(group0, group1, alternative="two-sided")

plt.figure(figsize=(4, 5))
ax = sns.boxplot(
    data=plot_df,
    x="ecDNA+",
    y="active_loci_TPM_log2_adj",
    showcaps=True,
    boxprops=dict(facecolor="none"),
    whiskerprops=dict(color="black"),
    medianprops=dict(color="black"),
    showfliers=False
)

sns.stripplot(
    data=plot_df,
    x="ecDNA+",
    y="active_loci_TPM_log2_adj",
    jitter=0.25, 
    size=4,   
    alpha=0.6,     
    color="black",  
    ax=ax
)

# counts per group
counts = plot_df["ecDNA+"].value_counts().sort_index()

# build new x-tick labels with n
new_labels = [
    f"{ec} (n={counts.loc[ec]})"
    for ec in counts.index
]

ax.set_xticklabels(new_labels)

# add p-value annotation
y_max = plot_df["active_loci_TPM_log2_adj"].max()
ax.plot([0, 1], [y_max * 1.1, y_max * 1.1], color="black", lw=1)
ax.text(
    0.5,
    y_max * 1.1,
    f"Mann–Whitney p = {pval:.3e}",
    ha="center",
    fontsize=10
)

ax.set_xlabel("ecDNA+")
ax.set_ylabel("Active L1 loci TPM (log2 adjusted)")
ax.set_title("TCGA: active L1 loci TPM by ecDNA status")

plt.tight_layout()
plt.show()

plot_df = tcga_l1_ecdna.dropna(subset=["number_amp_region", "RT_burden_log2_adj"])

# split groups
groups = {
    "ecDNA−": plot_df[plot_df["ecDNA+"] == False],
    "ecDNA+": plot_df[plot_df["ecDNA+"] == True]
}

# compute correlations
stats = {}
for label, g in groups.items():
    r, p = spearmanr(g["number_amp_region"], g["RT_burden_log2_adj"])
    stats[label] = (r, p)

plt.figure(figsize=(6, 5))
ax = sns.scatterplot(
    data=plot_df,
    x="number_amp_region",
    y="RT_burden_log2_adj",
    hue="ecDNA+",
    s=70
)

# add trend lines (visual aid only)
sns.regplot(
    data=groups["ecDNA−"],
    x="number_amp_region",
    y="RT_burden_log2_adj",
    scatter=False,
    ci=None,
    ax=ax
)
sns.regplot(
    data=groups["ecDNA+"],
    x="number_amp_region",
    y="RT_burden_log2_adj",
    scatter=False,
    ci=None,
    ax=ax
)

# annotate stats
y_max = plot_df["RT_burden_log2_adj"].max()
y_min = plot_df["RT_burden_log2_adj"].min()
y_range = y_max - y_min

for i, (label, (r, p)) in enumerate(stats.items()):
    ax.text(
        0.02,
        0.98 - i * 0.05,
        f"{label}: ρ = {r:.2f}, p = {p:.2e}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10
    )

ax.set_xlabel("Amplicon count")
ax.set_ylabel("L1 RT burden (log2 adjusted)")
ax.set_title("Correlation between amplicon count and L1 RT burden")

plt.tight_layout()
plt.show()
