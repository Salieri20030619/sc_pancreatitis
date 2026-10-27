import scanpy as sc, numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
sc.settings.figdir="results"
import os; os.chdir("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

adata = sc.read_h5ad("results/GSE181276_step1.h5ad")
print(">> 细胞:", adata.n_obs, "| 聚类:", adata.obs['leiden'].nunique())

# 胰腺 + 免疫 经典 marker(小鼠)
markers = {
  "Acinar":      ["Prss1","Cpa1","Cela1","Amy2a5","Try4","Ctrb1"],
  "Ductal":      ["Krt19","Sox9","Spp1","Krt8","Cftr"],
  "Endocrine":   ["Ins1","Ins2","Gcg","Sst","Chga"],
  "Stellate/Fib":["Col1a1","Col3a1","Pdgfrb","Acta2","Dcn"],
  "Endothelial": ["Pecam1","Cdh5","Cldn5","Flt1"],
  "Macrophage":  ["Adgre1","Csf1r","Lyz2","C1qa","Itgam"],
  "Neutrophil":  ["S100a8","S100a9","Ly6g","Retnlg"],
  "DC":          ["Flt3","Cd209a","Batf3","Itgax"],
  "T/NK":        ["Cd3e","Cd8a","Nkg7","Ccl5"],
  "B":           ["Cd79a","Cd79b","Ms4a1"],
}
# 只保留数据中存在的基因
markers = {k:[g for g in v if g in adata.var_names] for k,v in markers.items()}

# 各 cluster 的 top marker(rank_genes_groups)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=15)
top = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(10)
top.to_csv("results/cluster_top_markers.csv", index=False)
print(">> 各 cluster top10 marker 已存 results/cluster_top_markers.csv")

# marker 点图(cluster × 已知marker)
sc.pl.dotplot(adata, markers, groupby="leiden", save="_step2_markers.png",
              show=False, standard_scale="var", dendrogram=False)
print("\n== 各 cluster top10 marker(用于注释判断)==")
for cl in sorted(top.columns, key=int):
    print(f"cluster {cl}: {', '.join(top[cl].tolist())}")
