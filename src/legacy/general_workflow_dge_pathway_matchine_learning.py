import GEOparse
import pandas as pd

gse = GEOparse.get_GEO("GSE213001", destdir="./")

counts = pd.read_csv("GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz", sep=",", index_col=0)
print(counts.iloc[:3, :3])
pheno = gse.phenotype_data
print(pheno.iloc[:3, :3])
cols = [c for c in pheno.columns if ("age" in c.lower()) or ("sex" in c.lower()) or ("gender" in c.lower())]
print("Matching columns:", cols)
# Show those columns (first few rows)
print(pheno[cols].head())

rename_map = {}

pheno = pheno.rename(columns={"characteristics_ch1.11.age": "age"})
pheno = pheno.rename(columns={"characteristics_ch1.9.gender": "sex"})
pheno = pheno.rename(columns={"characteristics_ch1.8.diseasegroup": "condition"})

pheno['age'] = pd.to_numeric(pheno['age'], errors='coerce')
pheno['sex'] = pheno['sex'].astype('category')
pheno['condition'] = pheno['condition'].astype('category')

#print(pheno.columns)
#print(pheno['characteristics_ch1.8.diseasegroup'])
print("Renamed columns:", rename_map)
print(pheno[["age", "sex"]].head())

type(counts)
print("Omics Shape:", counts.shape)
print("Pheno Shape:", pheno.shape)
print(counts.iloc[:3, :3])
print(pheno.iloc[:3, :3])

pheno = pheno.set_index('title')
print(pheno.iloc[:3, :3])

common_samples = counts.columns.intersection(pheno.index)
counts = counts[common_samples]
pheno = pheno.loc[common_samples]
assert all(counts.columns == pheno.index)
print("✅ Counts and metadata are aligned!")

counts_t = counts.T
pheno = pheno.loc[counts_t.index]
print(pheno.shape)
print(counts_t.shape)
assert all(counts_t.index == pheno.index)

import pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Make sure counts are samples × genes

# Align metadata to counts
pheno = pheno.loc[counts_t.index]
assert all(counts_t.index == pheno.index), "Counts and metadata samples do not match!"

# Convert age to numeric and sex to categorical
pheno['age'] = pd.to_numeric(pheno['age'], errors='coerce')
pheno['sex'] = pheno['sex'].str.strip().str.upper().astype('category')
pheno = pheno.dropna(subset=['age', 'sex'])  # optional: drop missing

# Subset counts to matching samples again after dropping NAs
counts_t = counts_t.loc[pheno.index]

# Create DESeq2 dataset using formula (modern pyDESeq2)
dds = DeseqDataSet(
    counts=counts_t,
    metadata=pheno,
    design="~ age + sex + condition"  # main variable = condition
)

# Run DESeq2
dds.deseq2()

# Extract results for condition contrast
stat_res = DeseqStats(dds, contrast=["condition", "IPF", "NDC"])
stat_res.summary()
res_df = stat_res.results_df
print(res_df.head())

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Example res_df columns: 'gene', 'log2FoldChange', 'padj'
df = res_df.copy()

# -log10 adjusted p-value
df['neg_log10_padj'] = -np.log10(df['padj'])

# Optional: define significance thresholds
fc_thresh = 1      # log2 fold change threshold
pval_thresh = 0.05 # adjusted p-value threshold

# Define significance for coloring
df['significant'] = 'Not Sig'
df.loc[(df['padj'] < pval_thresh) & (df['log2FoldChange'] > fc_thresh), 'significant'] = 'Up'
df.loc[(df['padj'] < pval_thresh) & (df['log2FoldChange'] < -fc_thresh), 'significant'] = 'Down'

plt.figure(figsize=(8,6))
colors = {'Up':'red', 'Down':'blue', 'Not Sig':'grey'}

# Scatter plot
for sig in ['Up', 'Down', 'Not Sig']:
    subset = df[df['significant']==sig]
    plt.scatter(subset['log2FoldChange'], subset['neg_log10_padj'],
                c=colors[sig], label=sig, alpha=0.7, edgecolors='k', s=50)

# Labels and lines
plt.axhline(-np.log10(pval_thresh), color='black', linestyle='--', lw=1)
plt.axvline(fc_thresh, color='black', linestyle='--', lw=1)
plt.axvline(-fc_thresh, color='black', linestyle='--', lw=1)

plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(padj)')
plt.title('Volcano Plot')
plt.legend()
plt.show()

import seaborn as sns
sns.scatterplot(data=df, x='log2FoldChange', y='neg_log10_padj', hue='significant', palette={'Up':'red','Down':'blue','Not Sig':'grey'})

# gene set enrichment analyses

import pandas as pd
import numpy as np

res_df_output = res_df
res_df_output['gene'] = res_df_output.index

# Make sure you're using Series, not DataFrame slices
res_df_output['rank'] = (
    -np.log10(res_df_output['pvalue']) * 
    np.abs(res_df_output['log2FoldChange']) / 
    res_df_output['log2FoldChange']
)
# Make sure res_df has 'gene' and 'log2FoldChange' columns
ranked_genes = res_df_output[['gene', 'rank']].dropna()

# Sort by effect size (descending)
ranked_genes = ranked_genes.sort_values('rank', ascending=False)

# Convert to Series: index=gene, value=log2FC
ranked_series = ranked_genes.set_index('gene')['rank']

import mygene
import pandas as pd

# ranked_series: pd.Series, index = Ensembl IDs, values = signed log2FC or rank score
ensembl_ids = ranked_series.index.tolist()

# Initialize MyGeneInfo
mg = mygene.MyGeneInfo()

# Query Ensembl IDs for Entrez IDs
query_res = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

# Build mapping dictionary: Ensembl -> Entrez
ensembl_to_symbol = {r['query']: str(r.get('symbol')) for r in query_res if 'symbol' in r}

# Filter ranked_series to only IDs that successfully mapped
ranked_series = ranked_series[ranked_series.index.isin(ensembl_to_symbol.keys())]

# Replace index with Entrez IDs
ranked_series.index = ranked_series.index.map(ensembl_to_symbol)

# Optional: remove duplicates if multiple Ensembl IDs map to same Entrez
ranked_series = ranked_series.groupby(ranked_series.index).mean()
print(ranked_series.head(10))

import gseapy as gp

# Example using GO Biological Process 2021
pre_res = gp.prerank(
    rnk=ranked_series,
    gene_sets='KEGG_2019_Human',  # or use KEGG, Reactome, MSigDB collections
    processes=4,                            # number of cores
    permutation_num=100,                     # reduce for testing, default=1000
    outdir='gsea_results',                  # folder to save results
    format='png',                            # save figures as PNG
    seed=42
)

# Table of top enriched gene sets
gsea_res_df = pre_res.res2d
print(gsea_res_df.head())

# Plot top gene set enrichment
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Assuming you already ran prerank GSEA
# pre_res = gp.prerank(...)

# Extract results table
gsea_res_df = pre_res.res2d

# Sort by FDR or p-value and select top 10 pathways
gsea_res_df['abs_NES'] = gsea_res_df['NES'].abs()

# Sort by absolute NES and select top 20
top_terms = gsea_res_df.sort_values('abs_NES', ascending=False).head(20)
colors = ['red' if x > 0 else 'blue' for x in top_terms['NES']]
# Optional: reverse order so top pathway is on top in barplot
#top_terms = top_terms[::-1]

# Plot barplot of NES (Normalized Enrichment Score)
plt.figure(figsize=(8,6))
sns.barplot(x='NES', y='Term', data=top_terms, palette=colors)
plt.xlabel('Normalized Enrichment Score (NES)')
plt.ylabel('KEGG Pathway')
plt.title('Top 10 Enriched KEGG Pathways')
plt.tight_layout()
plt.show()

# residual preparation for machine learning
import pandas as pd
import numpy as np
import statsmodels.api as sm

# -----------------------------
# 1️⃣ Ensure numeric data
# -----------------------------
# counts: genes x samples or samples x genes
# pheno: metadata with 'age', 'sex', 'condition'

# Make counts numeric
counts_cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
counts = np.log2(counts_cpm + 1)
counts = counts.apply(pd.to_numeric, errors='coerce')

# Ensure phenotype covariates are numeric
pheno['age'] = pd.to_numeric(pheno['age'], errors='coerce')
print(pheno['age'])
# Convert sex to numeric (example: Male=1, Female=0)
pheno['sex'] = pheno['sex'].map({'M': 1, 'F': 0})
print(pheno['sex'])
# -----------------------------
# 2️⃣ Align samples

# -----------------------------
# Ensure counts and pheno have the same order of samples
if counts.shape[0] == pheno.shape[0]:
    X = counts.copy()
else:
    X = counts.T  # samples x genes
X = X.loc[pheno.index]  # align

# -----------------------------
# 3️⃣ Drop any samples with missing values
# -----------------------------
data = pd.concat([X, pheno[['age', 'sex']]], axis=1)
data = data.dropna()
X = data[X.columns]
covariates = data[['age', 'sex']]


# -----------------------------
# 4️⃣ Calculate residuals
# -----------------------------
residuals = pd.DataFrame(index=X.index, columns=X.columns, dtype=float)

for gene in X.columns:
    model = sm.OLS(X[gene], sm.add_constant(covariates))
    results = model.fit()
    residuals[gene] = results.resid

# -----------------------------
# 5️⃣ residuals are numeric and ready
# -----------------------------
print(residuals.head())
print(residuals.dtypes.unique())  # should all be float

# sklearn train-test split

mask = pheno['condition'].isin(['IPF', 'NDC'])
pheno_filtered = pheno.loc[mask]
X_filtered = residuals.loc[mask]  # align rows with filtered pheno

y = pheno_filtered['condition'].apply(lambda x: 1 if x == 'IPF' else 0)


from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X_filtered, y, test_size=0.2, random_state=42, stratify=y
)

# train model

from sklearn.ensemble import RandomForestClassifier

rf = RandomForestClassifier(n_estimators=500, max_depth=None, random_state=42)
rf.fit(X_train, y_train)

# ACC/AUC

from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

y_pred = rf.predict(X_test)
y_pred_proba = rf.predict_proba(X_test)[:, 1]
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(6,6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve - Random Forest Disease Prediction')
plt.legend(loc='lower right')
plt.show()

from sklearn.metrics import confusion_matrix

# y_test = true labels, y_pred = predicted labels
cm = confusion_matrix(y_test, y_pred)
print(cm)

# top 20 features

import seaborn as sns
import numpy as np

# Extract feature importance
feat_imp = pd.DataFrame({
    'gene': X_train.columns,
    'importance': rf.feature_importances_
})

# Sort by importance and take top 20
top_features = feat_imp.sort_values('importance', ascending=False).head(20)
top_features['gene'] = top_features.gene.map(ensembl_to_symbol)

# Plot
plt.figure(figsize=(8,6))
sns.barplot(x='importance', y='gene', data=top_features, palette='viridis')
plt.title('Top 20 Important Genes for Disease Prediction')
plt.tight_layout()
plt.show()

