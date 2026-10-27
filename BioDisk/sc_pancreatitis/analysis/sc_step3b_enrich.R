# ============================================================
# GSE181276 腺泡细胞 AP vs CTRL 上调基因 → 菌+免疫 GO/KEGG 富集
# 复用 WTA 配色与风格
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(ggplot2); library(tidyr); library(tibble)
  library(clusterProfiler); library(org.Mm.eg.db); library(patchwork); library(cowplot)
})
setwd("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")
dir.create("results", showWarnings=FALSE)

# ---- WTA 同款配色 ----
PAL_BLUE<-"#00C5FD"; PAL_PURPLE<-"#E16CFC"; PAL_ORANGE<-"#FFC565"; PAL_RED<-"#D83010"
BASE<-13; TITLE<-14
bold_theme <- theme(text=element_text(face="bold"), axis.text=element_text(face="bold",color="black"),
  axis.title=element_text(face="bold"), legend.text=element_text(face="bold"),
  legend.title=element_text(face="bold"), plot.title=element_text(face="bold",size=TITLE),
  plot.margin=margin(5,5,5,2))
retry <- function(expr,times=5,wait=6,what="op"){for(i in seq_len(times)){
  r<-tryCatch(force(expr),error=function(e){message(sprintf("[%s]%d失败:%s",what,i,conditionMessage(e)));NULL})
  if(!is.null(r))return(r); if(i<times){message(sprintf("[%s]%ds重试",what,wait));Sys.sleep(wait)}}; NULL}

# ---- 读腺泡上调基因 ----
genes <- readLines("results/acinar_AP_up_genes.txt")
genes <- genes[genes!=""]
cat(">> 腺泡上调基因:", length(genes), "\n")
eg <- bitr(genes, "SYMBOL","ENTREZID", org.Mm.eg.db)
cat(">> 成功映射:", nrow(eg), "\n")

# ================= GO:BP 富集 =================
ego <- enrichGO(eg$ENTREZID, org.Mm.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
write.csv(as.data.frame(ego), "results/sc_acinar_GO_BP.csv", row.names=FALSE)

god <- as.data.frame(ego) %>% arrange(p.adjust) %>% head(15) %>%
  mutate(Description=str_wrap(Description,30), Description=factor(Description,levels=rev(Description)),
         neglogp=-log10(p.adjust), bact=grepl("bact|antimicr|defense|immun|inflamm|wound|leukocyte|cytokine",Description,ignore.case=TRUE))
pC <- ggplot(god, aes(neglogp,Description,fill=bact)) + geom_col() +
  scale_fill_manual(values=c(`TRUE`=PAL_PURPLE,`FALSE`=PAL_BLUE),labels=c("other","microbe/immune"),name="") +
  labs(x="-log10(p.adjust)",y=NULL,title="Acinar AP-up: GO:BP") +
  theme_bw(base_size=BASE)+bold_theme+theme(axis.text.y=element_text(size=BASE-4,face="bold"),
    legend.position="top",legend.direction="horizontal")
bac_go <- as.data.frame(ego) %>% filter(grepl("bact|antimicr|defense",Description,ignore.case=TRUE))
cat("== 菌相关 GO 命中 ==\n"); print(head(bac_go[,c("Description","p.adjust","Count")],10))

# ================= KEGG 富集(带重试)=================
invisible(retry(clusterProfiler::download_KEGG("mmu"), times=5, wait=6, what="KEGG预取"))
ekegg <- retry(setReadable(enrichKEGG(eg$ENTREZID,"mmu",pvalueCutoff=0.1), org.Mm.eg.db,"ENTREZID"),
               times=5, wait=8, what="enrichKEGG")
if(!is.null(ekegg) && nrow(as.data.frame(ekegg))>0){
  kdf<-as.data.frame(ekegg); write.csv(kdf,"results/sc_acinar_KEGG.csv",row.names=FALSE)
  noise<-"Ribosome|Coronavirus|Estrogen|Alcoholism|Oxytocin|Proteoglycans in cancer|Pancreatic secretion"
  ikey<-"bact|infection|Salmonella|trap|Toll|NOD|phagos|Leukocyte|Antigen|Neutrophil|Complement|inflamm|cytokine|chemokine|immun|IL-17|TNF|NF-kappa|apopto|necro"
  kk<-kdf %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>%
    filter(grepl(ikey,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(14)
  if(nrow(kk)<6) kk<-kdf %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(14)
  kk<-kk %>% mutate(M=as.numeric(sub("/.*","",BgRatio)),RichFactor=Count/M,neglogp=-log10(p.adjust),
    patho=grepl("bact|Salmonella|infection|trap|phagos|Toll|NOD",Description,ignore.case=TRUE),
    Description=str_wrap(Description,30),Description=factor(Description,levels=rev(Description)))
  pD <- ggplot(kk, aes(RichFactor,Description,size=Count,color=neglogp)) + geom_point() +
    geom_text(aes(label=ifelse(patho,"*","")),color="red",size=6,fontface="bold",nudge_x=0.004) +
    scale_color_gradient(low=PAL_ORANGE,high=PAL_PURPLE,name="-log10\n(padj)") +
    scale_size(range=c(2.5,8),name="Count") +
    labs(x="Rich Factor",y=NULL,title="Acinar AP-up: KEGG (* pathogen)") +
    theme_bw(base_size=BASE)+bold_theme+theme(axis.text.y=element_text(size=BASE-4,face="bold"))
  cat(">> KEGG 成功, 通路:", nrow(kdf), "\n")
} else { pD<-ggplot()+annotate("text",1,1,label="KEGG unavailable")+theme_void(); cat(">> KEGG 失败\n") }

# ---- 拼图导出 ----
big <- plot_grid(pC, pD, ncol=2, align="none")
big <- ggdraw()+draw_grob(grid::rectGrob(gp=grid::gpar(fill="white",col=NA)))+draw_plot(big)
ggsave("results/Figure_sc_acinar_enrich.png", big, width=15, height=7, dpi=300, bg="white")
ggsave("results/Figure_sc_acinar_enrich.pdf", big, width=15, height=7)
cat("\n==== Step3b 完成 ====\n菌相关GO:",nrow(bac_go)," KEGG:",ifelse(is.null(ekegg),"失败","成功"),"\n")
