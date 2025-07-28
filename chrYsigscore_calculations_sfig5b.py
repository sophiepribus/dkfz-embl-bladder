import pandas as pd
import numpy as np
import glob
import gseapy as gp

# load all gene count data
gene_counts = pd.read_csv("salmon.merged.gene_counts.tsv", sep='\t')

# subset to male patients
male_gene_counts = gene_counts[["gene_id", "gene_name", "B123tumor", "B125tumor", "B134tumor", "B154tumor", "B157tumor", "B175tumor", "B178tumor", "B24tumor", "B39tumor", "B4tumor"]]

chrY_genes = ["RPS4Y1", "ZFY", "PCDH11Y", "TBL1Y", "USP9Y", "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "HSFY2", "KDM5D", "EIF1AY", "RBMY1A1", "PRY2", "DAZ1", "DAZ2", "DAZ3", "DAZ4"]

metadata = male_gene_counts.iloc[:, :2]
numeric_data = male_gene_counts.iloc[:, 2:]

# apply log2(x + 1)
numeric_data_log = numeric_data.applymap(lambda x: np.log2(x + 1))

gene_counts_log = pd.concat([metadata, numeric_data_log], axis=1)
gene_counts_log.index = gene_counts_log['gene_name'].str.upper()

expr_cols = gene_counts_log.columns.difference(['gene_id', 'gene_name'])

gene_counts_log = gene_counts_log[expr_cols]

gene_sets = {'Y_signature': chrY_genes}

# Run ssGSEA
ssgsea_result = gp.ssgsea(data=gene_counts_log,
                          gene_sets=gene_sets,
                          sample_norm_method='rank',
                          outdir=None,
                          permutation_num=0,
                          no_plot=True)

y_scores = ssgsea_result.res2d.T
y_scores = y_scores.T
y_scores = y_scores[['Name', 'ES']].copy()
y_scores.columns = ['Sample', 'Y_score']

# load all coverage data (from t2t)
cov_dir = 't2t_covs/'
cov_files = glob.glob(os.path.join(cov_dir, '*.cov'))

coverage_results = []

for filepath in cov_files:
    sample = os.path.basename(filepath).replace('.cov', '')
    
    df = pd.read_csv(filepath, sep='\t')
    
    # filter for chrY
    df_chrY = df[df['chr'].isin(['chrY', 'Y'])]
    
    cn_col = df_chrY.columns[-1]
    
    avg_cn = df_chrY[cn_col].mean()
    
    coverage_results.append({'Sample': sample, 'avg_chrY_CN': avg_cn})

avg_chrY_df = pd.DataFrame(coverage_results)

merged_df = pd.merge(y_scores, avg_chrY_df, on='Sample')
merged_df['avg_chrY_CN'] = pd.to_numeric(merged_df['avg_chrY_CN'], errors='coerce')
merged_df['Y_score'] = pd.to_numeric(merged_df['Y_score'], errors='coerce')
merged_df.to_csv("chrYCN_score.csv", index=False)