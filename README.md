# TCGA_BRCA

复现数据，需要先做非常重要的setup。

在空文件夹中，创建一个新的项目，在这个文件夹中建三个名字分别是"Data", "R", "Results"的子文件夹（严格大小写，否则在运行时路径可能出错），然后按照下面指示把文件放到对应文件夹中
create a new project in a new folder -> create three folders named "Data", "R", "Results",respectively -> move documents to their folders
<img width="360" alt="Screen Shot 2024-11-09 at 16 54 43" src="https://github.com/user-attachments/assets/f6952074-e992-4765-a3e0-e1d7af8e165d">

后缀为.R的脚本放入R文件夹里，
RawTPM.csv：只包含处理后的数据（提取了标准化后的TPM数据，提取重复基因的最大值，去除80%样本中未表达的基因）
lncRNA_matrix.csv：从RawTPM中提取了gene_type为lncRNA的基因表达矩阵
MG_matrix.csv：从RawTPM中提取了m6A/m5C相关的基因表达矩阵
脚本运行顺序
