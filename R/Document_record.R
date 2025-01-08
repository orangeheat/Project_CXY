#文件太多了，写着写着都忘了
#####1.tidyup########
#"Data/metadata.cart.2024-01-03.json"从TCGA网站上选取数据是得到的包含文件名的json
#"Data/clinical.tsv"从TCGA网站上下载得到的临床数据
#'Results/RawTPM.csv'经下载整合清洗过后的TPM数据, rather "raw",如果要复现可以从这个数据开始

#######2.correlation##########
#'Results/lncRNA_matrix.csv' “RawTPM”里所有的lncRNA表达数据
#'Results/MG_matrix.csv' “RawTPM”里所有的m5C/m6A相关基因表达数据
#'Results/MGR_lncRNA.csv' 与MG(m5C/m6A genes)相关性系数大于0.5，p小于0.001的lncRNA表达数据

######3.clincal_exp_merge##########
#'Results/Clinical_data.csv'清洗过的临床数据
#'Results/clinical_expression.csv'连接临床数据和表达数据，通过TCGAid筛选了肿瘤样本

######4.modelling&KM######
#'/Results/train.csv'通过caret随机用OS分的训练集，包含临床数据和MGR_lncRNA表达数据
#'/Results/test.csv'测试集
#'/Results/model_train.csv'多因素cox建模的矩阵，包含建模基因表达量，OS和风险评分
#'/Results/model_test.csv'包含预测测试集风险评分
#/Results/model_whole.csv'整集，包含风险评分
#'/Results/multiCox_results.csv'模型用到的基因相关系数
#'/Results/km_train.jpeg'训练集的KM曲线
#'/Results/km_test.jpeg'测试集的KM曲线
#'/Results/km_whole.jpeg'整集的KM曲线

######5.KM_by_groups######
# km_stage_I_II, km_stage_III_IV, Age_under65, Age_over65, T1_2, T3_4, M0, M1, N0, N123这些jpeg都是相关分组的KM曲线
#'./Results/independence_verification.csv' 在整集中对独立预后预测因子的鉴定








