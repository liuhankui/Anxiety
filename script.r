#require R-packages
library(EWCE)
library(ggplot2)

#download the full Karolinska dataset 
#reference: Skene N G, Bryois J, Bakken T E, et al. Genetic identification of brain cell types underlying schizophrenia[J]. Nature genetics, 2018, 50(6): 825-833.
download.file("http://www.hjerling-leffler-lab.org/data/scz_singlecell/ctdFiles.zip",destfile="ctdFiles.zip") 
unzip("ctdFiles.zip")
path <- "ctd_allKI.rda"
load(path)

#load the human orthologs of MGI genes
data("mouse_to_human_homologs")
m2h <- unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
mouse.bg <- unique(m2h$MGI.symbol)

#prepare gene lists
ADs.genes <- c('AC109466','ADCY2','AP003464','ASXL3','BTN1A1','C15orf62','CAMKMT','CNGB3','DNAJC17','EDIL3','ESR1','GCHFR','GIGYF2','IQCE','KNL1','MAD1L1','MYH15','NCAM1','NTRK2','NXPE3','OPRL1','PDE4B','PFND1','PLCG1','PREPL','PTPRT','RFWD2','RGS19','RMDN3','RN7SKP120','RRIQ3','SATB1','SDK1','TCEA2','TMEM106B','TRPV6','ZFYVE19')
ADs.specific.genes <- c('AC109466','AP003464','BTN1A1','C15orf62','CNGB3','DNAJC17','EDIL3','ESR1','GCHFR','GIGYF2','IQCE','KNL1','MYH15','NTRK2','NXPE3','OPRL1','PFND1','PREPL','PTPRT','RFWD2','RGS19','RMDN3','RN7SKP120','RRIQ3','SATB1','SDK1','TCEA2','TRPV6','ZFYVE19')
CBT.genes <- c('ARL5A','BDNF','C4orf51','C8orf86','CACNB4','COMT','FBLN2','FGFR1','GRIK4','GRIN2B','HTR2A','LETM2','MMAA','NEB','NGF','NTNG1','SLC6A4','SMAD1','TPH2','VAV3','WHSC1L1','WNT7A','ZNF827')

#exclude genes without orthologs
ADs.genes.mouse <- unique(m2h[m2h$HGNC.symbol %in% ADs.genes,"MGI.symbol"])
ADs.specific.genes.mouse <- unique(m2h[m2h$HGNC.symbol %in% ADs.specific.genes,"MGI.symbol"])
CBT.genes.mouse <- unique(m2h[m2h$HGNC.symbol %in% CBT.genes,"MGI.symbol"])

#calculate cell-type specificity and corresponding p-value for each gene list
set.seed(202104)
ADs.results <- bootstrap.enrichment.test(hits=ADs.genes.mouse,sct_data=ctd,bg=mouse.bg,reps=10000,annotLevel=1)$results[,3:5]
set.seed(202104)
ADs.specific.results <- bootstrap.enrichment.test(hits=ADs.specific.genes.mouse,sct_data=ctd,bg=mouse.bg,reps=10000,annotLevel=1)$results[,3:5]
set.seed(202104)
CBT.results <- bootstrap.enrichment.test(hits=CBT.genes.mouse,sct_data=ctd,bg=mouse.bg,reps=10000,annotLevel=1)$results[,3:5]

#save the cell types for x-axis of plot
ADs.results$cellTypes <- row.names(ADs.results)
ADs.specific.results$cellTypes <- row.names(ADs.specific.results)
CBT.results$cellTypes <- row.names(CBT.results)

#adjust p-values via fdr method 
ADs.results$p_adjust <- p.adjust(ADs.results$p,method='fdr')
ADs.specific.results$p_adjust <- p.adjust(ADs.specific.results$p,method='fdr')
CBT.results$p_adjust <- p.adjust(CBT.results$p,method='fdr')

#mark the names of the three gene lists
ADs.results$GeneSet <- 'ADs-related genes'
ADs.specific.results$GeneSet  <- 'ADs-specific genes'
CBT.results$GeneSet <- 'CBT-related genes'

#merge the results of the three gene lists
All.results <- rbind(ADs.results,ADs.specific.results,CBT.results)

#visualization
ggplot(All.results,aes(cellTypes,-log10(p_adjust)))+
    geom_histogram(stat = 'identity',position = "dodge",aes(fill = GeneSet))+
    geom_hline(yintercept = -log10(0.05),colour='red')+
    xlab('')+ylab(expression(-log10(paste("FDR ","P-Values"))))+ylim(0,2.2)+theme_bw()+
    theme(axis.text.x = element_text(angle=60,hjust=1),
          legend.title = element_blank(),
          legend.position = c(0.15,0.8),
          legend.background  = element_rect(colour = 'black'))
