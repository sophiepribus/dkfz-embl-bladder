# This code was used for L1 count plots (figure 3H, supplementary fig. 10F). 

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# load L1 information

l1_table_folder = "070425_L1_counts/"

all_l1s = []

for file in os.listdir(l1_table_folder):
    l1_file = os.path.join(l1_table_folder, file)
    if l1_file[-8:] != "L1.table":
        continue
    if os.stat(l1_file).st_size == 0:
        print(f"Skipping empty file: {l1_file}")
        continue
    l1_df = pd.read_csv(l1_file, 
                     sep='\t', 
                     header=None)
    sample_name = l1_file.split("tumor")[0].split("/")[1]
    if sample_name == "B060":
        sample_name = "B60"
    l1_df["sample"] = sample_name
    all_l1s.append(l1_df)
    
all_l1s_df = pd.concat(all_l1s, ignore_index=True)

# load amplicon information

amplicon_folder = "010825_ontbc_amplicons/"
amplicon_dict = {}

for filename in os.listdir(amplicon_folder):
    if filename.endswith(".tsv"):
        
        sample_name = filename.replace("tumor.amplicon.gene.tsv", "")
        
        filepath = os.path.join(amplicon_folder, filename)
        df = pd.read_csv(filepath, sep="\t")
        
        amplicon_dict[sample_name] = df

amplicon_counts = {sample: len(df) for sample, df in amplicon_dict.items()}

l1_counts = all_l1s_df["sample"].value_counts().to_dict()

all_samples = set(amplicon_counts) | set(l1_counts)

# load SV information

sv_folder = "031125_sv_tsvs"

data = []

for sv_file in os.listdir(sv_folder):
    if sv_file.endswith(".tsv"):
        sample = sv_file.replace("tumor.sv.tsv", "")
        sv_path = os.path.join(sv_folder, sv_file)

        sv_df = pd.read_csv(sv_path, sep="\t", header=None,
                            names=["chr1", "start1", "chr2", "end2", "svtype"])
        
        # filter out insertions
        sv_filtered = sv_df[sv_df["svtype"].str.upper() != "INS"]

        sv_count = len(sv_filtered)

        if sample not in l1_counts:
            l1_count = 0
        else:
            l1_count = l1_counts[sample]
        
        l1_count_log = np.log(l1_count + 1)
        data.append({"sample": sample, "sv_count": sv_count, "l1_count": l1_count, "l1_count_log": l1_count_log})
        
df = pd.DataFrame(data).dropna()

amp_high_samples = []

for sample, amp_count in amplicon_counts.items():
    if amp_count > 20:
        amp_high_samples.append(sample)
        
df["amplicon_high"] = (df["sample"].isin(amp_high_samples))

ecDNA_pos_samps = ["B42", "B123", "B125", "B4", "B5", "B39", "B74"]
df["ecDNA_pos"] = (df["sample"].isin(ecDNA_pos_samps))

### plot L1 counts (figure 3H)

L1ins_df_no_B42 = all_l1s_df[all_l1s_df["sample"] != "B42"]
sample_counts = L1ins_df_no_B42["sample"].value_counts().sort_values(ascending=False)  

# ordered by descending count for amplicon-high then amplicon-low samples (manually created)
custom_order_fig5H = ["B123", "B4", "B5", "B134", "B94", "B157", "B125", "B187", "B175", "B12", "B178", "B67", "B85", "B154", "B60", "B22", "B24", "B39", "B74", "B32", "B156", "B87"]

sample_counts_df = sample_counts.rename_axis("sample").reset_index(name="count")

sample_counts_df = sample_counts_df[sample_counts_df["sample"].isin(custom_order_fig5H)]
sample_counts_df["sample"] = pd.Categorical(sample_counts_df["sample"], categories=custom_order_fig5H, ordered=True)

# plot (no b42)
plt.figure(figsize=(15, 5))
sns.barplot(data=sample_counts_df, x="sample", y="count", color="steelblue", order=custom_order_fig5H)

plt.xticks(fontsize=16, rotation=0)
plt.xlabel("Sample", fontsize=16)
plt.ylabel("L1 insertions", fontsize=16)
plt.yticks(fontsize=16)
plt.ylim(0, 75, 10)
plt.title("L1 Insertion Counts per Sample", fontsize=16)
plt.tight_layout()
plt.show()

### plot L1 counts (figure 3H)

L1ins_df_no_B42 = all_l1s_df[all_l1s_df["sample"] != "B42"]
sample_counts = L1ins_df_no_B42["sample"].value_counts().sort_values(ascending=False)  

# ordered by descending count for ecDNA+ then ecDNA- samples (manually created)
custom_order_sfig10F = ["B123", "B4", "B5", "B125", "B39", "B74", "B134", "B94", "B157", "B187", "B175", "B12", "B178", "B67", "B85", "B154", "B60", "B22", "B24", "B32", "B156", "B87"]

sample_counts_df = sample_counts.rename_axis("sample").reset_index(name="count")

sample_counts_df = sample_counts_df[sample_counts_df["sample"].isin(custom_order_sfig10F)]
sample_counts_df["sample"] = pd.Categorical(sample_counts_df["sample"], categories=custom_order_sfig10F, ordered=True)

# plot (no b42)
plt.figure(figsize=(15, 5))
sns.barplot(data=sample_counts_df, x="sample", y="count", color="steelblue", order=custom_order_sfig10F)

plt.xticks(fontsize=16, rotation=0)
plt.xlabel("Sample", fontsize=16)
plt.ylabel("L1 insertions", fontsize=16)
plt.yticks(fontsize=16)
plt.ylim(0, 75, 10)
plt.title("L1 Insertion Counts per Sample", fontsize=16)
plt.tight_layout()
plt.show()

b42_count = all_l1s_df["sample"].value_counts().get("B42", 0)

# create b42 bar plot (separated for better scale visualization)
plt.figure(figsize=(.75, 5))
sns.barplot(x=["B42"], y=[b42_count], color="steelblue")

plt.ylabel("L1 Insertion Count", fontsize=16)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.title("L1 insertion Count for B42")
plt.ylim(0, max(b42_count + 100, 10))
plt.tight_layout()
plt.show()

### get stats for amplicon-high vs amplicon-low

group1 = df[df["amplicon_high"] == True]["l1_count"]
group2 = df[df["amplicon_high"] == False]["l1_count"]

stat, pval = mannwhitneyu(group1, group2, alternative="greater")
n1, n2 = len(group1), len(group2)
r_rb = 1 - (2 * stat) / (n1 * n2)

print(f'p-value = {pval}')
print(f'Rank-biserial effect size = {r_rb}')

### get stats for ecDNA+ vs ecDNA-

group1 = df[df["ecDNA_pos"] == True]["l1_count"]
group2 = df[df["ecDNA_pos"] == False]["l1_count"]

stat, pval = mannwhitneyu(group1, group2, alternative="greater")
n1, n2 = len(group1), len(group2)
r_rb = 1 - (2 * stat) / (n1 * n2)

print(f'p-value = {pval}')
print(f'Rank-biserial effect size = {r_rb}')