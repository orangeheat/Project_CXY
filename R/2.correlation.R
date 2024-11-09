###setup
library(here)
library(dplyr)
rm(list = ls())
#######Pearson分析#########
matrix <- read.csv("Results/RawTPM.csv", row.names=1)
#提取lncRNA矩阵
lncRNA_matrix <- matrix[matrix$gene_type=='lncRNA',] |>
  select(-gene_type)

#提取m6A/m5C相关RNA
m6A_related_gene <-c('METTL3','METTL14','METTL16','WTAP','VIRMA','RBMX','RBM3','RBM10','ZC3H13','FTO','ALKBH5','YTHDC1','YTHDF1','YTHDF2','YTHDF3','HNRNPA2B1','IGFBP1', 'IGFBP3', 'IGFBP2','HNRNPC')
m5C_related_gene <- c('NOP2',"NSUN2","TRDMT1",'TET1',"TET2",'TET3',"ALYREF","YBX1")
MG <- unique(c(m6A_related_gene, m5C_related_gene))
MG_matrix <- matrix[rownames(matrix) %in% MG,] |>
  select(-gene_type)

#写下lncRNA和MG矩阵
write.csv(lncRNA_matrix,'Results/lncRNA_matrix.csv',row.names = TRUE)
write.csv(MG_matrix,'Results/MG_matrix.csv',row.names = TRUE)


# 改造矩阵，因为相关性分析需要逐行计算，所以转置
lncRNA_matrix <- t(lncRNA_matrix)
MG_matrix <- t(MG_matrix)

cor_matrix <- cor(lncRNA_matrix, MG_matrix, method = "pearson")
p_values_matrix <- matrix(NA, ncol = ncol(cor_matrix), nrow = nrow(cor_matrix))
colnames(p_values_matrix) <- colnames(MG_matrix)
rownames(p_values_matrix) <- colnames(lncRNA_matrix)

for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    correlation_test <- cor.test(lncRNA_matrix[, i], MG_matrix[, j])
    p_values_matrix[i, j] <- correlation_test$p.value
  }
}


#######筛选######
threshold <- 0.5
p_threshold <- 0.001
filtered_matrix <- cor_matrix * (abs(cor_matrix) > threshold & p_values_matrix < p_threshold)
# 提取满足相关性系数大于0.5，p小于0.001的lncRNA
MG_related_lncRNA <- rownames(filtered_matrix)[which(filtered_matrix != 0, arr.ind = TRUE)[, 1]]
MG_related_lncRNA <- unique(MG_related_lncRNA)
MGR_lncRNA_matrix <- lncRNA_matrix[,colnames(lncRNA_matrix) %in% MG_related_lncRNA]
write.csv(MGR_lncRNA_matrix,"Results/MGR_lncRNA.csv", row.names = T)

#### 做了一张表，最后用网站画的桑基图
lncRNAs <- rownames(filtered_matrix_1)
MGs <- colnames(filtered_matrix_1)
cor_coefficient <- as.vector(filtered_matrix_1)
# 创建数据框，过滤掉值为0的元素
plot_data <- data.frame(
  From = rep(lncRNAs, each = length(MGs)),
  To = rep(MGs, times = length(lncRNAs)),
  Value = cor_coefficient
) 
plot_data <- dplyr::filter(plot_data, Value!= 0)
write.csv(plot_data, "./RawData/初筛相关系数数据.csv", row.names = F)
