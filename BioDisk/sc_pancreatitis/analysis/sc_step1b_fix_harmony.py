import scanpy as sc, numpy as np
import harmonypy
import warnings; warnings.filterwarnings("ignore")
sc.settings.figdir="results"
import os; os.chdir("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

adata = sc.read_h5ad("results/GSE181276_step1.h5ad")
print(">> 读回:", adata.n_obs, "细胞")

ah = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(ah, max_value=10)
sc.tl.pca(ah, n_comps=30)
n_cells = ah.n_obs

# ---- 手动 Harmony,自动判断 Z_corr 方向 ----
ho = harmonypy.run_harmony(ah.obsm["X_pca"], ah.obs, ["timepoint"])
Z = np.asarray(ho.Z_corr)
print(">> Z_corr 原始形状:", Z.shape, "(应含", n_cells, "细胞)")
# 让行数 = 细胞数
if Z.shape[0] != n_cells:
    Z = Z.T
print(">> 校正后用于 obsm 的形状:", Z.shape)
assert Z.shape[0] == n_cells, "维度仍不匹配"
ah.obsm["X_pca_harmony"] = Z.copy()

sc.pp.neighbors(ah, n_neighbors=15, use_rep="X_pca_harmony")
sc.tl.leiden(ah, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(ah)

adata.obs["leiden"] = ah.obs["leiden"].values
adata.obsm["X_umap"] = ah.obsm["X_umap"]
adata.obsm["X_pca_harmony"] = ah.obsm["X_pca_harmony"]

sc.pl.umap(adata, color=["leiden","timepoint"], save="_step1_harmony_fixed.png",
           show=False, ncols=2, wspace=0.3)
adata.write("results/GSE181276_step1.h5ad")
print("\n==== 修复完成 ====")
print("聚类数:", adata.obs["leiden"].nunique(), "| 细胞:", adata.n_obs)
