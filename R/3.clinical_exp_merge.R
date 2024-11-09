#setup
library(here)
library(readr)
library(data.table)
library(dplyr)
library(tibble)
rm(list = ls())

#临床数据是重复的，且一位患者可以取多个sample，后续将根据TCGA barcode筛选肿瘤样本
clinical <- read_tsv('Data/clinical.tsv')
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),]) |>
  select(case_id, case_submitter_id,age_at_index,gender,race,
         vital_status,days_to_death,days_to_last_follow_up,
         ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,
         ajcc_pathologic_n,treatment_type)
colnames(clinical) <- c("caseID","TCGA_ID","Age","Gender","Race",
               "Status","days_to_death","days_to_last_follow_up",
               "Stage","T","M","N","Treatment")

#删去没有OS的数据
cli_matrix = clinical[clinical$Status %in% c('Alive','Dead'),] 

columns_to_numeric <- c("days_to_last_follow_up", "days_to_death", "Age")
cli_matrix[columns_to_numeric] <- lapply(cli_matrix[columns_to_numeric], function(x) {
  x[is.na(x)] <- 0
  as.numeric(x)
})

cli_matrix$days <- ifelse(cli_matrix$Status=='Alive',cli_matrix$days_to_last_follow_up,cli_matrix$days_to_death)
cli_matrix <- cli_matrix[cli_matrix$days > 0, ]
cli_matrix$months <- ceiling(cli_matrix$days / 30)
# 不足一个月的记为1个月
cli_matrix$months[cli_matrix$months == 0] <- 1 
cli_matrix <- cli_matrix |>
  mutate(OS = ifelse(Status == "Alive", 0, 1),
         OS.time = round(months / 12, 1))

write.csv(cli_matrix, "./Results/Clinical_data.csv", row.names = TRUE)





#加载筛选后的lncRNA数据
exp <- read.csv("./Results/MGR_lncRNA.csv", row.names = 1) |>
  mutate(ID = gsub("\\.", "-", rownames(exp))) |>
  mutate(TCGA_ID = substr(ID, 1, 12),
    type_barcode = substr(ID, 14, 16),
    portion = as.numeric(substr(ID, 23, 24))) |>
  select(ID, TCGA_ID, type_barcode, portion, everything())
#ID中TCGA_D8_A1XO_01A_11R_A14M_07，第四部分代表样本类型，01是肿瘤样本，06是转移样本，11是正常样本
#我这里是有"01A" "11A" "01B" "11B" "06A" "01C"，只保留01的样本；如果还有重复，优先保留肿瘤组织（01A）和portion大的样本
exp01 <- filter(exp, grepl("^01", type_barcode))
exp_filtered <- exp01 |>
  group_by(TCGA_ID) |>
  arrange(desc(type_barcode == "01A"), desc(portion)) |>
  filter(row_number() == 1) |>
  ungroup() |>
  select(-c(ID, type_barcode, portion))
  
#连接临床数据和表达数据
cli_exp <- cli_matrix |>
  inner_join(exp_filtered, by = "TCGA_ID") |>
  column_to_rownames(var = "TCGA_ID")
write.csv(cli_exp, "./Results/clinical_expression.csv", row.names = T)


