# ============================================================
# GSE181276 小鼠 cerulein 急性胰腺炎单细胞时间轴
# Step1: 10x稀疏读入 + QC + 归一化 + 聚类 + UMAP
# ============================================================
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
sc.settings.figdir = "results"
import os; os.chdir("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

# ---- 循环读入 4 个时间点(10x 稀疏, 秒级)----
tps = ["CTRL","Day1","Day7","Day28"]
adatas = []
for tp in tps:
    a = sc.read_10x_mtx(f"raw10x/{tp}", var_names="gene_symbols", cache=True)
    a.var_names_make_unique()
    a.obs_names = [f"{tp}_{bc}" for bc in a.obs_names]
    a.obs["timepoint"] = tp
    print(f"   {tp}: {a.n_obs} 细胞 (raw)")
    adatas.append(a)
adata = ad.concat(adatas, join="outer")
adata.obs["timepoint"] = pd.Categorical(adata.obs["timepoint"], categories=tps, ordered=True)
adata.obs["condition"] = np.where(adata.obs["timepoint"]=="CTRL","Control","AP")
print(">> 合并 raw:", adata.n_obs, "细胞 x", adata.n_vars, "基因")

# ---- QC ----
adata.var["mt"] = adata.var_names.str.startswith(("mt-","Mt-"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < 6000) & (adata.obs.pct_counts_mt < 20)].copy()
print(">> QC 后:", adata.n_obs, "细胞")
print(adata.obs["timepoint"].value_counts().sort_index())

# ---- 归一化 + log ----
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ---- HVG + PCA + (Harmony) + 聚类 + UMAP ----
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="timepoint")
ah = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(ah, max_value=10)
sc.tl.pca(ah, n_comps=30)
try:
    sc.external.pp.harmony_integrate(ah, "timepoint"); rep="X_pca_harmony"; print(">> Harmony 完成")
except Exception as e:
    rep="X_pca"; print(">> 跳过 Harmony(", str(e)[:40], ")")
sc.pp.neighbors(ah, n_neighbors=15, use_rep=rep)
sc.tl.leiden(ah, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(ah)
adata.obs["leiden"] = ah.obs["leiden"].values
adata.obsm["X_umap"] = ah.obsm["X_umap"]

sc.pl.umap(adata, color=["leiden","timepoint"], save="_step1_cluster_timepoint.png",
           show=False, ncols=2, wspace=0.3)
adata.write("results/GSE181276_step1.h5ad")
print("\n==== Step1 完成 ====")
print("聚类数:", adata.obs["leiden"].nunique(), "| 细胞:", adata.n_obs, "| 基因:", adata.n_vars)
