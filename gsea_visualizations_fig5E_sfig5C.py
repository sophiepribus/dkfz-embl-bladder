# This code was used for the GSEA analyses visualizations (L1s - figure 5E, LOY - supplementary fig. 5C). 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# load GSEA tsv results
pos_df = pd.read_csv("high_gsea.tsv", sep="\t")
neg_df = pd.read_csv("low_gsea.tsv", sep="\t")

pos_df["NOM p-val"] = pd.to_numeric(pos_df["NOM p-val"], errors="coerce")
neg_df["NOM p-val"] = pd.to_numeric(neg_df["NOM p-val"], errors="coerce")

# subset to gene sets with a nominal p-value < 0.05
df_pos_sig = pos_df[pos_df["NOM p-val"] < 0.05].copy()
df_neg_sig = neg_df[neg_df["NOM p-val"] < 0.05].copy()

# color gene sets by FDR q-value threshold at 0.05
df_pos_sig["Significant"] = df_pos_sig["FDR q-val"] < 0.05
df_neg_sig["Significant"] = df_neg_sig["FDR q-val"] < 0.05

df_all = pd.concat([df_pos_sig, df_neg_sig])

plt.figure(figsize=(4, 6))
ax = sns.barplot(
    data=df_all,
    y="NAME", 
    x="NES", 
    hue="Significant", 
    dodge=False,
    palette={True: "green", False: "orange"},
)

ax.set_xlabel("low group <-- NES --> high group", fontsize=10)
ax.set_ylabel("Pathway", fontsize=10)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.set_title("GSEA: high group versus low group", fontsize=10)
ax.axvline(0, color="black", linestyle="--", linewidth=1)

ax.legend(title="FDR q-value\n< 0.05", loc="lower right", fontsize=10, title_fontsize=10)

plt.show()