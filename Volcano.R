res<-read.csv("~/Desktop/gene_DE.csv",sep=',',header=T, row.names=1)
genetype2<-c('CXCL1','CCL7','FCN1','C3','CD300E','MMP14','VCAN','VSIG4','CCL2','S100A8','S100A9','PPBP','IL6','CA12','CEBPB','IRAK3','IL1A','CD163','LAIR1','S100A12','CXCL3','VNN2','LILRB1','LILRB2','CXCL6','MMP8','IL10')
genetype1<-c('CD1A','CD1C','CCL17','MMP12','IFIT1','CD1B','RSAD2','ADAM19','HLA-DRB1','USP18','CD1E','CLEC10A','HLA-DPB1','LIPA','HLA-DMB','HLA-DQB1','HLA-DRA','STAT1','HLA-DQA1','HLA-DPA1','HLA-DMA','HLA-DQB2','ECM1')



genetype2<-c('CXCL1','CCL7','FCN1','C3','CD300E','MMP14','VCAN','VSIG4','CCL2','S100A8','S100A9','PPBP','IL6','CA12','IL1A','CD163','LAIR1','S100A12','CXCL3','VNN2','CXCL6','MMP8','IL10')
genetype1<-c('CD1A','CD1C','CCL17','MMP12','IFIT1','CD1B','RSAD2','ADAM19','HLA-DRB1','USP18','CD1E','CLEC10A','LIPA','HLA-DRA','STAT1','ECM1')

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'qValue',
                selectLab = c(genetype1,genetype2),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-6,
                FCcutoff = 1.0,
                pointSize = 4.0,
                labSize = 3.0,
                colAlpha = 0.7,
                drawConnectors=TRUE,
                legendLabSize = 12,
                legendIconSize = 4.0)
