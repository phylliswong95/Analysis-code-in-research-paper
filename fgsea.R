
library(data.table)
library(fgsea)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)

#读取差异分析好的数据
ranks1 <- read.csv("ranks.csv",head=T)
#读取基因名时的处理，跳过ID转换
ranks <- ranks1$log2FC
names(ranks) <- ranks1$SYMBOL

#ID转换，先制作ID对照表，再合并表格
keytypes(org.Hs.eg.db)
k=keys(org.Hs.eg.db, keytype="SYMBOL")
list=AnnotationDbi::select(org.Hs.eg.db, keys=k, columns = c("ENTREZID","SYMBOL"), keytype = "SYMBOL")
ID_list=list[match(ranks1$SYMBOL, list[,"SYMBOL"]),]
ranks2 <- inner_join(ranks1, ID_list, by="SYMBOL")
#去除重复值
table(duplicated(ranks2$ENTREZID))
#提取所需数据并保存
ranks3 <- ranks2[!duplicated(ranks2$ENTREZID), ]
ranks <- ranks3$log2FC
names(ranks) <- ranks3$ENTREZID
write.csv(ranks, file="ranks.csv")


#读取基因列表gmt文件-方法1
gmtfile <- "c7.all.v2023.1.Hs.symbols.gmt"
hallmark <- read.gmt(gmtfile)
hallmark$term <- gsub ('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply("[[",2)

#读取gmt文件-方法2
x <- readLines("c7.all.v2023.1.Hs.entrez.gmt")
c7 <- strsplit(x, "\t")
names(c7) <- vapply(c7, function(y) y[1], character(1))
c7 <- lapply(c7, "[", -c(1:2))

#用fgsea进行富集分析
fgseaRes <- fgsea(pathways=hallmark.list, stats = ranks, minSize=0, maxSize=10000)
sig <- fgseaRes[fgseaRes$padj<0.05]
sig <- sig[order(sig$NES, decrease=T),]

#数次尝试分析p值都很大，考虑用GSEA软件分析的GeneCounts结果直接对应
#读取GSEA分析结果
report1<-read.table("gsea_report.tsv",header = T, comment.char = "", sep = "\t")
report2 <- inner_join(report1, fgseaRes, by="pathway")
report3 <- report2[,-(7:11)]
colnames(report3) <- c("pathway","ES","NES","pval","padj", "log2err","size", "leadingEdge")


pic <- fgseaRes[c(3370,3375,1327,1404,475,3537,224),]

#GeneSet列表
plotGseaTable(hallmark.list[pic$pathway], ranks, report3)

#用clusterProfiler
ego3 <- gseGO(geneList     = ranks,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

GO_kk_entrez <- gseGO(geneList     = ranks,
                      ont          = "BP",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = "org.Hs.eg.db",#人类org.Hs.eg.db 鼠org.Mm.eg.db
                      keyType      = "SYMBOL",
                      pvalueCutoff = 0.25)   #实际为padj阈值可调整
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb="org.Hs.eg.db",
                           keyType='ENTREZID')#转化id 
save(GO_kk_entrez, file = "C5BP_result.RData")

#clusterprofiler本地基因集
hallmark <- read.gmt("c2.all.v2023.1.Hs.symbols.gmt")

c2res <- GSEA(ranks, TERM2GENE = hallmark, pvalueCutoff = 2)



#载入做图包
devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
#圆圈图
dotplotGsea(data = c5BPres,topn = 20,
            order.by = 'NES',
            add.seg = T)
#GSEA富集图
gseaNb(object = c2res,
       geneSetID = "BUCKANOVICH_T_LYMPHOCYTE_HOMING_ON_TUMOR_DN" )


