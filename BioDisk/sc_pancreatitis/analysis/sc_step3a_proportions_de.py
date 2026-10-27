import scanpy as sc, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings("ignore")
sc.settings.figdir="results"; sc.settings.verbosity=1
import os; os.chdir("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

adata = sc.read_h5ad("results/GSE181276_step2.h5ad")
tps=["CTRL","Day1","Day7","Day28"]

# ================= 1. 细胞类型比例随时间 =================
prop = pd.crosstab(adata.obs["timepoint"], adata.obs["celltype"], normalize="index")*100
prop = prop.loc[tps]
prop.to_csv("results/celltype_proportion.csv")
print(">> 细胞类型比例(%)已存")

fig,ax=plt.subplots(figsize=(9,5))
prop.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.75)
ax.set_ylabel("Proportion (%)",fontweight="bold"); ax.set_xlabel("")
ax.set_title("Cell type composition over AP time course",fontweight="bold")
ax.legend(bbox_to_anchor=(1.01,1),loc="upper left",fontsize=8,frameon=False)
plt.xticks(rotation=0); plt.tight_layout()
plt.savefig("results/sc_celltype_proportion.png",dpi=200,bbox_inches="tight"); plt.close()

# 免疫比例单独折线(趋势更清晰)
imm_prop = prop["Immune"]
fig,ax=plt.subplots(figsize=(6,4))
ax.plot(range(len(tps)), imm_prop.values, marker="o", linewidth=2, color="#E16CFC")
ax.set_xticks(range(len(tps))); ax.set_xticklabels(tps,fontweight="bold")
ax.set_ylabel("Immune cell %",fontweight="bold")
ax.set_title("Immune infiltration over time",fontweight="bold")
plt.tight_layout(); plt.savefig("results/sc_immune_proportion.png",dpi=200); plt.close()
print(">> 免疫比例趋势:", dict(zip(tps, imm_prop.round(2).values)))

# ================= 2. 抗菌肽基因在各细胞类型/时间点表达 =================
amp_genes=[g for g in ["Reg3b","Reg3g","Reg3a","Reg2","Reg1","Lyz2","Dmbt1","Lcn2"] if g in adata.var_names]
# 各细胞类型点图
sc.pl.dotplot(adata, amp_genes, groupby="celltype", save="_step3_AMP_celltype.png",
              show=False, standard_scale="var", title="Antimicrobial peptides by cell type")
# 腺泡细胞内,抗菌肽随时间
ac = adata[adata.obs["celltype"]=="Acinar"].copy()
sc.pl.dotplot(ac, amp_genes, groupby="timepoint", save="_step3_AMP_acinar_time.png",
              show=False, standard_scale="var", title="AMP in Acinar over time")
print(">> 抗菌肽定位图完成; 腺泡细胞数:", ac.n_obs)

# ================= 3. 腺泡 AP vs CTRL 差异,导出给 R 富集 =================
ac.obs["grp"]=np.where(ac.obs["timepoint"]=="CTRL","CTRL","AP")
sc.tl.rank_genes_groups(ac, "grp", groups=["AP"], reference="CTRL", method="wilcoxon")
de = sc.get.rank_genes_groups_df(ac, group="AP")
de.columns=["gene","score","logFC","pval","padj"]
de.to_csv("results/acinar_AP_vs_CTRL_DE.csv", index=False)
sig = de[(de.padj<0.05)&(de.logFC.abs()>0.5)]
print(f">> 腺泡 AP vs CTRL: DEG={len(sig)} (up={sum(sig.logFC>0)}, down={sum(sig.logFC<0)})")

# 全体细胞 AP vs CTRL(备用,给富集更全的背景)
adata.obs["grp"]=np.where(adata.obs["timepoint"]=="CTRL","CTRL","AP")
sc.tl.rank_genes_groups(adata,"grp",groups=["AP"],reference="CTRL",method="wilcoxon")
de_all=sc.get.rank_genes_groups_df(adata,group="AP")
de_all.columns=["gene","score","logFC","pval","padj"]
de_all.to_csv("results/allcells_AP_vs_CTRL_DE.csv",index=False)
sig_all=de_all[(de_all.padj<0.05)&(de_all.logFC.abs()>0.5)]
print(f">> 全体 AP vs CTRL: DEG={len(sig_all)}")

# 导出显著上调基因列表(给 R 富集用)
sig[sig.logFC>0].sort_values("padj")["gene"].to_csv("results/acinar_AP_up_genes.txt",index=False,header=False)
print("\n==== Step3a 完成,DE 已导出供 R 富集 ====")
print("腺泡上调 top15:", sig[sig.logFC>0].sort_values("padj").head(15)["gene"].tolist())
