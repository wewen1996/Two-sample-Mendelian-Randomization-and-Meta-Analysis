library(tidyverse)

ID <- read.csv("id.csv",row.names = 1)
CeD <- read.csv("res_mendelian_CeD.csv")
DM <- read.csv("res_mendelian_DM.csv")
IBD <- read.csv("res_mendelian_IBD.csv")
RA <- read.csv("res_mendelian_RA.csv")
SLE <- read.csv("res_mendelian_SLE.csv")
T1D <- read.csv("res_mendelian_T1D.csv")

CeD$Metabolite.ID <- unlist(strsplit(CeD$exposure,".metal.pos.txt.gz"))
DM$Metabolite.ID <- unlist(strsplit(DM$exposure,".metal.pos.txt.gz"))
IBD$Metabolite.ID <- unlist(strsplit(IBD$exposure,".metal.pos.txt.gz"))
RA$Metabolite.ID <- unlist(strsplit(RA$exposure,".metal.pos.txt.gz"))
SLE$Metabolite.ID <- unlist(strsplit(SLE$exposure,".metal.pos.txt.gz"))
T1D$Metabolite.ID <- unlist(strsplit(T1D$exposure,".metal.pos.txt.gz"))

CeD <- merge(CeD, ID, by = "Metabolite.ID")
DM <- merge(DM, ID, by = "Metabolite.ID")
IBD <- merge(IBD, ID, by = "Metabolite.ID")
RA <- merge(RA, ID, by = "Metabolite.ID")
SLE <- merge(SLE, ID, by = "Metabolite.ID")
T1D <- merge(T1D, ID, by = "Metabolite.ID")

Total <- merge(CeD, DM, by = "Metabolite.ID")
Total <- merge(Total, IBD, by = "Metabolite.ID")
Total <- merge(Total, RA, by = "Metabolite.ID")
Total <- merge(Total, SLE, by = "Metabolite.ID")
Total <- merge(Total, T1D, by = "Metabolite.ID")
write.csv(Total, file = "total.csv")

P_heatmap <- read.csv("P_heatmap.csv")
P_heatmap <-  P_heatmap[!duplicated(P_heatmap$Metabolite.x),]
rownames(P_heatmap) <- P_heatmap[,1]
P_heatmap <- P_heatmap[,-1]
P_heatmap$CeD <- ifelse(P_heatmap$CeD > 0.05, 0, P_heatmap$CeD)
P_heatmap$DM <- ifelse(P_heatmap$DM > 0.05, 0, P_heatmap$DM)
P_heatmap$IBD <- ifelse(P_heatmap$IBD > 0.05, 0, P_heatmap$IBD)
P_heatmap$RA <- ifelse(P_heatmap$RA > 0.05, 0, P_heatmap$RA)
P_heatmap$SLE <- ifelse(P_heatmap$SLE > 0.05, 0, P_heatmap$SLE)
P_heatmap$T1D <- ifelse(P_heatmap$T1D > 0.05, 0, P_heatmap$T1D)

P_heatmap <- P_heatmap [which(rowSums(P_heatmap) > 0),]
P_heatmap <- as.data.frame(P_heatmap)
P_heatmap$Metabolite.x <- rownames(P_heatmap)
write.csv(P_heatmap, file = "P_heatmap_new.csv")

P_heatmap$CeD <- ifelse(P_heatmap$CeD > 0.05, " ",
                        ifelse(P_heatmap$CeD > 0.01, "*",
                               ifelse(P_heatmap$CeD > 0.001, "**", "***")))
P_heatmap$DM <- ifelse(P_heatmap$DM > 0.05, " ",
                        ifelse(P_heatmap$DM > 0.01, "*",
                               ifelse(P_heatmap$DM > 0.001, "**", "***")))
P_heatmap$IBD <- ifelse(P_heatmap$IBD > 0.05, " ",
                       ifelse(P_heatmap$IBD > 0.01, "*",
                              ifelse(P_heatmap$IBD > 0.001, "**", "***")))
P_heatmap$RA <- ifelse(P_heatmap$RA > 0.05, " ",
                        ifelse(P_heatmap$RA > 0.01, "*",
                               ifelse(P_heatmap$RA > 0.001, "**", "***")))
P_heatmap$SLE <- ifelse(P_heatmap$SLE > 0.05, " ",
                       ifelse(P_heatmap$SLE > 0.01, "*",
                              ifelse(P_heatmap$SLE > 0.001, "**", "***")))
P_heatmap$T1D <- ifelse(P_heatmap$T1D > 0.05, " ",
                        ifelse(P_heatmap$T1D > 0.01, "*",
                               ifelse(P_heatmap$T1D > 0.001, "**", "***")))
P_heatmap <- P_heatmap[,-7]
P_heatmap <- as.matrix(P_heatmap)

P_heatmap1 <- read.csv("P_heatmap_new.csv")
rownames(P_heatmap1) <- P_heatmap1[,1]
P_heatmap1 <- P_heatmap1[,-1]
# P_heatmap1 <- P_heatmap1[,-7]
# P_heatmap1$CeD <- ifelse(P_heatmap1$CeD > 0.05, " ",
#                         ifelse(P_heatmap1$CeD > 0.01, "*",
#                                ifelse(P_heatmap1$CeD > 0.001, "**", "***")))
# P_heatmap1$DM <- ifelse(P_heatmap1$DM > 0.05, " ",
#                        ifelse(P_heatmap1$DM > 0.01, "*",
#                               ifelse(P_heatmap1$DM > 0.001, "**", "***")))
# P_heatmap1$IBD <- ifelse(P_heatmap1$IBD > 0.05, " ",
#                         ifelse(P_heatmap1$IBD > 0.01, "*",
#                                ifelse(P_heatmap1$IBD > 0.001, "**", "***")))
# P_heatmap1$RA <- ifelse(P_heatmap1$RA > 0.05, " ",
#                        ifelse(P_heatmap1$RA > 0.01, "*",
#                               ifelse(P_heatmap1$RA > 0.001, "**", "***")))
# P_heatmap1$SLE <- ifelse(P_heatmap1$SLE > 0.05, " ",
#                         ifelse(P_heatmap1$SLE > 0.01, "*",
#                                ifelse(P_heatmap1$SLE > 0.001, "**", "***")))
# P_heatmap1$T1D <- ifelse(P_heatmap1$T1D > 0.05, " ",
#                         ifelse(P_heatmap1$T1D > 0.01, "*",
#                                ifelse(P_heatmap1$T1D > 0.001, "**", "***")))
# P_heatmap1 <- as.matrix(P_heatmap1)

# annotation <- read.csv("annotation_heatmap.csv")
# annotation <- merge(annotation, P_heatmap, by = "Metabolite.x")
# annotation <- annotation[,-c(3:8)]
# write.csv(annotation,file = "annotation_new.csv")

annotation <- read.csv("annotation_heatmap.csv")
annotation <- merge(annotation, P_heatmap1, by = "Metabolite.x")
annotation <- annotation[,-c(3:8)]
annotation <-  annotation[!duplicated(annotation$Metabolite.x),]
rownames(annotation) <- annotation$Metabolite.x
write.csv(annotation,file = "annotation_new.csv")
annotation2 <- read.csv("annotation_new.csv")
rownames(annotation2) <- annotation2[,1]
annotation2 <- annotation2[,-c(1), drop=FALSE]

beta_heatmap <- read.csv("beta_heatmap.csv")
beta_heatmap <-  beta_heatmap[!duplicated(beta_heatmap$Metabolite.x),]
rownames(beta_heatmap) <- beta_heatmap[,1]
beta_heatmap <- beta_heatmap[P_heatmap1$Metabolite.x,]
beta_heatmap <- beta_heatmap[,-1]
beta_heatmap$Metabolite.x <- rownames(beta_heatmap)
beta_heatmap <- beta_heatmap[,-7]

beta_heatmap <- beta_heatmap[match(annotation2$Metabolite.x, rownames(beta_heatmap)),] 
#annotation <- annotation[match(rownames(beta_heatmap),annotation$Metabolite.x),] 
a <- rownames(annotation2)
b <- rownames(beta_heatmap)
identical(a,b)
annotation <- annotation2[,-1, drop=FALSE]
colnames(annotation) <- c("Category")
annotation$Category <- as.factor(annotation$Category)




library(pheatmap)
# 绘制热图
p4 <- pheatmap(beta_heatmap,  
               color = colorRampPalette(c("navy","white","firebrick3"))(100), 
               border_color = "black",  
               scale = "column", 
               cluster_rows = FALSE, 
               cluster_cols = FALSE, 
               legend = TRUE, 
               legend_breaks = c(-2, 0, 2), 
               legend_labels = c("Low","","Heigh"), 
               show_rownames = TRUE, 
               show_colnames = TRUE, 
               fontsize = 8,
               gaps_row = c(22, 26, 27, 64, 65, 67),
               annotation_row  = annotation,   #添加列注释信息
               display_numbers = P_heatmap,  #将预先定义的矩阵传递给display_numbers参数
               number_color = "black", 
               fontsize_number = 6  
)
pdf("热图新.pdf",width = 12,height = 8)
p4
dev.off()

