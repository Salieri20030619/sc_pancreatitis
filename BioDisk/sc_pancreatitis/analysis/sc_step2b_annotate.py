import scanpy as sc, numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
sc.settings.figdir="results"; sc.settings.verbosity=1
import os; os.chdir("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

adata = sc.read_h5ad("results/GSE181276_step1.h5ad")

# ---- cluster -> 主要细胞类型映射(基于 marker 判定)----
mapping = {
 "0":"Fibroblast","1":"Fibroblast","2":"Stellate","3":"Fibroblast","4":"Acinar",
 "5":"Endothelial","6":"Endocrine","7":"Ductal","8":"Endothelial","9":"Proliferating",
 "10":"Ductal","11":"Endocrine","12":"Acinar","13":"Acinar","14":"Acinar",
 "15":"SmoothMuscle","16":"Immune","17":"Pericyte","18":"Glia","19":"Erythrocyte",
 "20":"Mesothelial","21":"Acinar",
}
adata.obs["celltype"] = adata.obs["leiden"].map(mapping).astype("category")
print(">> 主要细胞类型分布:\n", adata.obs["celltype"].value_counts())

# ---- UMAP 按细胞类型 ----
sc.pl.umap(adata, color=["celltype","timepoint"], save="_step2_celltype.png",
           show=False, ncols=2, wspace=0.35)

# ---- 免疫细胞单独细分 ----
imm = adata[adata.obs["celltype"]=="Immune"].copy()
print(f"\n>> 免疫细胞 {imm.n_obs} 个,单独再聚类")
sc.pp.highly_variable_genes(imm, n_top_genes=1500)
imm2 = imm[:, imm.var.highly_variable].copy()
sc.pp.scale(imm2, max_value=10); sc.tl.pca(imm2, n_comps=20)
sc.pp.neighbors(imm2, n_neighbors=15, use_rep="X_pca")
sc.tl.leiden(imm2, resolution=0.6, flavor="igraph", n_iterations=2, directed=False, key_added="imm_leiden")
sc.tl.umap(imm2)
imm.obs["imm_leiden"]=imm2.obs["imm_leiden"].values; imm.obsm["X_umap"]=imm2.obsm["X_umap"]

imm_markers={
 "Macrophage":["Adgre1","Csf1r","C1qa","C1qb","Lyz2"],
 "Neutrophil":["S100a8","S100a9","Retnlg","Ly6g"],
 "DC":["Flt3","Itgax","Cd209a","Batf3"],
 "T/NK":["Cd3e","Cd8a","Nkg7","Ccl5","Cd3d"],
 "B":["Cd79a","Cd79b","Ms4a1"],
 "Mast":["Cpa3","Mcpt4","Kit"],
}
imm_markers={k:[g for g in v if g in imm.var_names] for k,v in imm_markers.items()}
sc.tl.rank_genes_groups(imm,"imm_leiden",method="wilcoxon",n_genes=15)
itop=pd.DataFrame(imm.uns["rank_genes_groups"]["names"]).head(10)
sc.pl.dotplot(imm, imm_markers, groupby="imm_leiden", save="_step2_immune_markers.png",
              show=False, standard_scale="var")
sc.pl.umap(imm, color="imm_leiden", save="_step2_immune_umap.png", show=False)

print("\n== 免疫亚群 top10 marker ==")
for cl in sorted(itop.columns,key=int):
    print(f"imm cluster {cl}: {', '.join(itop[cl].tolist())}")

imm.write("results/GSE181276_immune.h5ad")
adata.write("results/GSE181276_step2.h5ad")
print("\n==== Step2 完成,免疫亚群数:", imm.obs['imm_leiden'].nunique(),"====")
