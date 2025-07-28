# This code was used for the structural variant clonality analysis (figure 5A calculations).

import os
import pandas as pd
import gzip
import shutil
import pysam

# load all somatic L1s (counts & ids)

l1_table_folder = "L1_counts/"

all_l1s_counts = []

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
    all_l1s_counts.append(l1_df)
    
all_l1s_counts_df = pd.concat(all_l1s_counts, ignore_index=True)

# load all variants

vcfs_path = "vcfs/"

all_l1s_raw = []

for file in os.listdir(vcfs_path):
    if file == ".ipynb_checkpoints":
        continue
    l1_file = os.path.join(vcfs_path, file)
    vcf_df = pd.read_csv(l1_file, 
                     compression='gzip', 
                     comment='#', 
                     sep='\t', 
                     header=None)
    all_l1s_raw.append(vcf_df)

all_l1s_raw_df = pd.concat(all_l1s_raw, ignore_index=True)
new_column_names = ['chr', 'pos', 'id', 'red', 'alt', 'qual', 'filter', 'info']
all_l1s_raw_df.columns = new_column_names
all_l1s_raw_df = all_l1s_raw_df[all_l1s_raw_df["filter"] == "PASS"]

# pull out L1, ALU and SVA type insertion IDs
all_L1_ids = list(all_l1s_df[3])
alu_ins = list(all_l1s_raw_df[all_l1s_raw_df['info'].str.contains('Alu', na=False)]["id"])
sva_ins = list(all_l1s_raw_df[all_l1s_raw_df['info'].str.contains('SVA', na=False)]["id"])

# load SV read data - for VAF calculation

bcfs_path = "sv_bcfs/"

bcf_dfs = []

for file in os.listdir(bcfs_path):
    
    if file.endswith('.somatic.final.bcf'):
    
        bcf_file = os.path.join(bcfs_path, file)
        sample_name = file.replace('tumor.somatic.final.bcf', '')

        bcf = pysam.VariantFile(bcf_file)

        bcf_results = []

        for record in bcf.fetch():
            variant_id = record.id
            chrom = record.chrom
            pos = record.pos

            # check if the FORMAT field contains 'RR' and 'RV' for the sample
            for sample in record.samples:
                # get RR (reference reads), default to 0 if missing
                rr = record.samples[sample].get('RR', 0)  
                # get RV (variant reads), default to 0 if missing
                rv = record.samples[sample].get('RV', 0)  

                # avoid division by 0 and calculate the fraction of variant reads if sum > 0
                if rr + rv > 0:
                    fraction = rv / (rr + rv)
                    bcf_results.append([sample_name, chrom, pos, variant_id, rr, rv, fraction])

        headers = ["sample", "chrom", "pos", "id", "RR", "RV", "VAF"]

        bcf_df = pd.DataFrame(bcf_results, columns=headers)
        bcf_dfs.append(bcf_df)
        
all_bcfs_df = pd.concat(bcf_dfs, ignore_index=True)

all_bcfs_df["SV type"] = all_bcfs_df["id"].str[:3]
all_bcfs_df.loc[all_bcfs_df['id'].isin(all_L1_ids), 'SV type'] = 'INS - L1'
all_bcfs_df.loc[all_bcfs_df['id'].isin(alu_ins), 'SV type'] = 'INS - Alu'
all_bcfs_df.loc[all_bcfs_df['id'].isin(sva_ins), 'SV type'] = 'INS - SVA'
all_bcfs_df.loc[all_bcfs_df['SV type'] == "INS", 'SV type'] = 'INS - other'

# assign clonality of each sample
threshold = 0.3
all_bcfs_df["Clonal"] = (all_bcfs_df["VAF"] > threshold)

all_bcfs_df['Clonal'] = pd.Categorical(all_bcfs_df['Clonal'], categories=[True, False], ordered=True)

# save for visualization in R
all_bcfs_df.to_csv("vafs.csv", index=False)