import scanpy as sc, numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
import os; os.chdir("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

adata = sc.read_h5ad("results/GSE181276_step2.h5ad")
tps=["CTRL","Day1","Day7","Day28"]

# 1. UMAP 坐标 + 注释(全细胞)
umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1","UMAP2"], index=adata.obs_names)
umap["celltype"]=adata.obs["celltype"].values
umap["timepoint"]=adata.obs["timepoint"].values
# 3.2万点太多,降采样到1.5万加速R渲染(分层抽样保比例)
umap_s = umap.groupby("celltype", group_keys=False).apply(lambda d: d.sample(min(len(d), int(len(d)*15000/len(umap))+1), random_state=1))
umap_s.to_csv("results/export_umap.csv")
print(">> UMAP 导出:", umap_s.shape[0], "点")

# 2. 细胞类型比例
prop = pd.crosstab(adata.obs["timepoint"], adata.obs["celltype"], normalize="index")*100
prop.loc[tps].to_csv("results/export_proportion.csv")

# 3. 抗菌肽在各细胞类型的平均表达(用于点图)
amp=[g for g in ["Reg3b","Reg3g","Reg3a","Reg2","Reg1","Lyz2","Dmbt1","Lcn2"] if g in adata.var_names]
expr = sc.get.obs_df(adata, keys=amp+["celltype"])
amp_mean = expr.groupby("celltype")[amp].mean()
# 表达比例(每类里表达该基因的细胞占比)
amp_pct = expr.groupby("celltype")[amp].apply(lambda d:(d>0).mean()*100)
amp_mean.to_csv("results/export_amp_mean.csv")
amp_pct.to_csv("results/export_amp_pct.csv")

# 4. 抗菌肽在腺泡内随时间
ac = adata[adata.obs["celltype"]=="Acinar"]
expr_t = sc.get.obs_df(ac, keys=amp+["timepoint"])
expr_t.groupby("timepoint")[amp].mean().loc[tps].to_csv("results/export_amp_acinar_time.csv")
print(">> 抗菌肽基因:", amp)
print("==== 导出完成 ====")
