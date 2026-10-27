# ============================================================
# GSE181276 单细胞 6-panel 大图(WTA 同款风格)
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(ggplot2); library(tidyr); library(tibble)
  library(cowplot); library(RColorBrewer)
})
setwd("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")

# ---- WTA 同款配色 ----
PAL_BLUE<-"#00C5FD"; PAL_PURPLE<-"#E16CFC"; PAL_GREEN<-"#51E6AD"; PAL_MAGENTA<-"#EE3BA5"
PAL_RED<-"#D83010"; PAL_ORANGE<-"#FFC565"; PAL_PINK<-"#FF9EC4"; PAL_NAVY<-"#1F4E9C"
PAL8 <- c(PAL_MAGENTA,PAL_RED,PAL_ORANGE,PAL_PURPLE,PAL_BLUE,PAL_GREEN,PAL_PINK,PAL_NAVY)
# 13 种细胞类型配色(扩展色板)
ct_cols <- c(colorRampPalette(PAL8)(13))
BASE<-13; TITLE<-14
bold_theme <- theme(text=element_text(face="bold"), axis.text=element_text(face="bold",color="black"),
  axis.title=element_text(face="bold"), legend.text=element_text(face="bold"),
  legend.title=element_text(face="bold"), plot.title=element_text(face="bold",size=TITLE),
  plot.margin=margin(5,5,5,2))
tps <- c("CTRL","Day1","Day7","Day28")

# ============ A/B: UMAP ============
um <- read.csv("results/export_umap.csv", row.names=1, stringsAsFactors=FALSE)
um$timepoint <- factor(um$timepoint, levels=tps)
cts <- sort(unique(um$celltype)); names(ct_cols) <- cts

pA <- ggplot(um, aes(UMAP1,UMAP2,color=celltype)) +
  geom_point(size=0.35, alpha=0.7) +
  scale_color_manual(values=ct_cols, name="") +
  guides(color=guide_legend(override.aes=list(size=3), ncol=2)) +
  labs(title="A. Cell types (UMAP)") +
  theme_bw(base_size=BASE)+bold_theme+
  theme(legend.text=element_text(size=BASE-5,face="bold"), legend.key.height=unit(10,"pt"),
        legend.position="right", axis.text=element_blank(), axis.ticks=element_blank())

pB <- ggplot(um, aes(UMAP1,UMAP2,color=timepoint)) +
  geom_point(size=0.35, alpha=0.7) +
  scale_color_manual(values=c(PAL_BLUE,PAL_ORANGE,PAL_GREEN,PAL_RED), name="") +
  guides(color=guide_legend(override.aes=list(size=3))) +
  labs(title="B. Timepoints (UMAP)") +
  theme_bw(base_size=BASE)+bold_theme+
  theme(legend.position="right", axis.text=element_blank(), axis.ticks=element_blank())

# ============ C: 细胞比例堆叠 ============
prop <- read.csv("results/export_proportion.csv", row.names=1, check.names=FALSE)
pl <- prop %>% rownames_to_column("timepoint") %>%
  pivot_longer(-timepoint, names_to="celltype", values_to="pct") %>%
  mutate(timepoint=factor(timepoint,levels=tps), celltype=factor(celltype,levels=cts))
pC <- ggplot(pl, aes(timepoint,pct,fill=celltype)) +
  geom_col(width=0.72) + scale_fill_manual(values=ct_cols,name="") +
  guides(fill=guide_legend(ncol=2)) +
  labs(x=NULL,y="Proportion (%)",title="C. Cell composition over time") +
  theme_bw(base_size=BASE)+bold_theme+
  theme(legend.text=element_text(size=BASE-5,face="bold"), legend.key.size=unit(10,"pt"))

# ============ D: 免疫比例趋势 ============
imm <- data.frame(timepoint=factor(tps,levels=tps), pct=prop[tps,"Immune"])
pD <- ggplot(imm, aes(as.integer(timepoint),pct)) +
  geom_line(linewidth=1.3,color=PAL_PURPLE) + geom_point(size=3,color=PAL_PURPLE) +
  scale_x_continuous(breaks=1:4,labels=tps) +
  labs(x=NULL,y="Immune cell %",title="D. Immune proportion over time") +
  theme_bw(base_size=BASE)+bold_theme

# ============ E: 抗菌肽在各细胞类型(点图: 大小=表达比例, 颜色=平均表达) ============
am_mean <- read.csv("results/export_amp_mean.csv", row.names=1, check.names=FALSE)
am_pct  <- read.csv("results/export_amp_pct.csv", row.names=1, check.names=FALSE)
amp_genes <- colnames(am_mean)
dm <- am_mean %>% rownames_to_column("celltype") %>% pivot_longer(-celltype,names_to="gene",values_to="mean")
dp <- am_pct  %>% rownames_to_column("celltype") %>% pivot_longer(-celltype,names_to="gene",values_to="pct")
de <- left_join(dm,dp,by=c("celltype","gene")) %>%
  mutate(gene=factor(gene,levels=amp_genes), celltype=factor(celltype,levels=rev(cts)))
pE <- ggplot(de, aes(gene,celltype,size=pct,color=mean)) + geom_point() +
  scale_color_gradient(low=PAL_BLUE, high=PAL_PURPLE, name="mean\nexpr") +
  scale_size(range=c(0,7), name="% expr") +
  labs(x=NULL,y=NULL,title="E. Antimicrobial peptides by cell type") +
  theme_bw(base_size=BASE)+bold_theme+
  theme(axis.text.x=element_text(angle=45,hjust=1,face="bold"), axis.text.y=element_text(size=BASE-3,face="bold"))

# ============ F: 免疫/防御 GO 富集(复用 step3c 逻辑) ============
IMM <- "bacter|antimicrob|defense|immun|inflamm|microb|pathogen|wound|leukocyte|cytokine|neutrophil|chemotax|phagocy|complement|interleukin|interferon|tumor necrosis"
go <- read.csv("results/sc_acinar_GO_BP.csv", stringsAsFactors=FALSE)
gimm <- go %>% filter(grepl(IMM,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(15) %>%
  mutate(Description=str_wrap(Description,30), Description=factor(Description,levels=rev(Description)),
         neglogp=-log10(p.adjust))
pF <- ggplot(gimm, aes(neglogp,Description,fill=neglogp)) + geom_col() +
  scale_fill_gradient(low=PAL_BLUE, high=PAL_PURPLE, name="-log10\n(padj)") +
  labs(x="-log10(p.adjust)",y=NULL,title="F. Acinar AP-up: immune/defense GO") +
  theme_bw(base_size=BASE)+bold_theme+theme(axis.text.y=element_text(size=BASE-5,face="bold"))

# ============ 拼图 3x2 ============
dir.create("results/tiff_panels", showWarnings=FALSE)
.sp <- list(A=pA,B=pB,C=pC,D=pD,E=pE,F=pF)
for(.nm in names(.sp)){
  tiff(sprintf("results/tiff_panels/SC_panel_%s.tiff",.nm), width=7,height=6,units="in",res=600,compression="lzw")
  print(.sp[[.nm]]); dev.off()
}
cat(">> SC 6 panel TIFF 完成\n")
big <- plot_grid(pA,pB,pC,pD,pE,pF, ncol=2, align="none", labels=NULL)
title_gg <- ggplot()+labs(title="GSE181276 scRNA-seq - Cerulein acute pancreatitis: cellular landscape & antimicrobial/immune response")+
  theme_void()+theme(plot.title=element_text(face="bold",size=17,hjust=0.01))
big <- plot_grid(title_gg, big, ncol=1, rel_heights=c(0.03,1))
big <- ggdraw()+draw_grob(grid::rectGrob(gp=grid::gpar(fill="white",col=NA)))+draw_plot(big)
ggsave("results/Figure_GSE181276_6panel.png", big, width=17, height=18, dpi=600, bg="white")
ggsave("results/Figure_GSE181276_6panel.pdf", big, width=17, height=18)
cat("==== 单细胞 6-panel 大图完成 ====\n")
