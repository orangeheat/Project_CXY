#setup
library(here)
library(jsonlite)
library(data.table)
library(rjson)
library(dplyr)
library(tibble)

json <- jsonlite::fromJSON("Data/metadata.cart.2024-01-03.json")
#提取TCGA ID aka submitter ID
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  

#读文件
count_file <- list.files('Data/13af6c0698b34f0a9d8a49a505ecc3f2',pattern = "\\.rna_seq\\.augmented_star_gene_counts\\.tsv$",recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

path = paste0('Data/13af6c0698b34f0a9d8a49a505ecc3f2//',count_file[1])
data0<- read.delim(path,fill = TRUE,header = F, sep = "\t", fileEncoding = "UTF-8", skipNul = TRUE, row.names = 1)
data0<-data0[-c(1:6),] #前六行是统计信息，删除
matrix <- data0[c(1,2)] # 前两列是symbol和genetype

#提取每个文件tpm列
for (i in 1:length(count_file_name)){
  path = paste0('Data/13af6c0698b34f0a9d8a49a505ecc3f2//',count_file[i])
  data<- read.delim(path,fill = TRUE,header = F, sep = "\t", fileEncoding = "UTF-8", skipNul = TRUE, row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[6] #取出tpm列
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#整理，取重复基因的最大值
data <- as.matrix(read.delim(path,fill = TRUE, header = F, row.names = 1))
gene_name <- data[-c(1:6), 1]
matrix0 <- matrix
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max) |>
  select(-V2) |>
  column_to_rownames(var = "gene_name") |>
  rename(gene_type = V3)

#去除80%样本中未表达的基因
matrix0[, 2:ncol(matrix0)] <- lapply(matrix0[, 2:ncol(matrix0)], as.numeric)
noexpression_80 <- rowMeans(matrix0[, 2:ncol(matrix0)] == 0) >= 0.8
matrix_filtered <- matrix0[!noexpression_80,]

#这里写好经处理的原始数据
write.csv(matrix_filtered,'Results/RawTPM.csv',row.names = TRUE)
