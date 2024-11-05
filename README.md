# TCGA_BRCA

RawTPM.csv：只包含处理后的数据（提取了标准化后的TPM数据，提取重复基因的最大值，去除80%样本中未表达的基因）
lncRNA_matrix.csv：从RawTPM中提取了gene_type为lncRNA的基因表达矩阵
MG_matrix.csv：从RawTPM中提取了m6A/m5C相关的基因表达矩阵
脚本运行顺序
