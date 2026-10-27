suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(ggplot2); library(cowplot)
})
setwd("/mnt/BioDisk/wangmingyuan/sc_pancreatitis")
PAL_BLUE<-"#00C5FD"; PAL_PURPLE<-"#E16CFC"; PAL_ORANGE<-"#FFC565"
BASE<-13; TITLE<-14
bold_theme <- theme(text=element_text(face="bold"), axis.text=element_text(face="bold",color="black"),
  axis.title=element_text(face="bold"), legend.text=element_text(face="bold"),
  legend.title=element_text(face="bold"), plot.title=element_text(face="bold",size=TITLE),
  plot.margin=margin(5,5,5,2))

# 免疫/防御主题关键词(放宽)
IMM <- "bacter|antimicrob|defense|immun|inflamm|microb|pathogen|wound|leukocyte|cytokine|neutrophil|chemotax|phagocy|complement|interleukin|interferon|tumor necrosis"

# ---- GO: 只展示免疫/防御通路 top15 ----
go <- read.csv("results/sc_acinar_GO_BP.csv", stringsAsFactors=FALSE)
gimm <- go %>% filter(grepl(IMM, Description, ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(15) %>%
  mutate(Description=str_wrap(Description,32), Description=factor(Description,levels=rev(Description)),
         neglogp=-log10(p.adjust))
cat(">> GO 免疫/防御通路命中:", nrow(go %>% filter(grepl(IMM,Description,ignore.case=TRUE))), "条\n")
pC <- ggplot(gimm, aes(neglogp,Description,fill=neglogp)) + geom_col() +
  scale_fill_gradient(low=PAL_BLUE, high=PAL_PURPLE, name="-log10\n(padj)") +
  labs(x="-log10(p.adjust)",y=NULL,title="Acinar AP-up: immune/defense GO:BP") +
  theme_bw(base_size=BASE)+bold_theme+theme(axis.text.y=element_text(size=BASE-4,face="bold"))

# ---- KEGG: 免疫/防御通路(若文件在)----
kf <- "results/sc_acinar_KEGG.csv"
if(file.exists(kf)){
  kg <- read.csv(kf, stringsAsFactors=FALSE)
  noise <- "Ribosome|Coronavirus|Estrogen|Alcoholism|Oxytocin|Proteoglycans in cancer|Pancreatic secretion|Protein digestion"
  kk <- kg %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>%
    filter(grepl(IMM, Description, ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(12)
  if(nrow(kk)<4) kk <- kg %>% filter(!grepl(noise,Description,ignore.case=TRUE)) %>% arrange(p.adjust) %>% head(12)
  kk <- kk %>% mutate(M=as.numeric(sub("/.*","",BgRatio)),RichFactor=Count/M,neglogp=-log10(p.adjust),
    patho=grepl("bacter|Salmonella|infection|trap|phagos|Toll|NOD",Description,ignore.case=TRUE),
    Description=str_wrap(Description,30),Description=factor(Description,levels=rev(Description)))
  pD <- ggplot(kk, aes(RichFactor,Description,size=Count,color=neglogp)) + geom_point() +
    geom_text(aes(label=ifelse(patho,"*","")),color="red",size=6,fontface="bold",nudge_x=0.004) +
    scale_color_gradient(low=PAL_ORANGE,high=PAL_PURPLE,name="-log10\n(padj)") +
    scale_size(range=c(2.5,8),name="Count") +
    labs(x="Rich Factor",y=NULL,title="Acinar AP-up: KEGG") +
    theme_bw(base_size=BASE)+bold_theme+theme(axis.text.y=element_text(size=BASE-4,face="bold"))
  cat(">> KEGG 展示通路:", nrow(kk), "条\n")
} else pD <- ggplot()+annotate("text",1,1,label="KEGG file missing")+theme_void()

big <- plot_grid(pC, pD, ncol=2, align="none", rel_widths=c(1,1))
big <- ggdraw()+draw_grob(grid::rectGrob(gp=grid::gpar(fill="white",col=NA)))+draw_plot(big)
ggsave("results/Figure_sc_acinar_enrich.png", big, width=16, height=7, dpi=300, bg="white")
ggsave("results/Figure_sc_acinar_enrich.pdf", big, width=16, height=7)
cat("\n== 展示的免疫/防御 GO ==\n"); print(as.character(gimm$Description))
cat("\n==== 重画完成 ====\n")
