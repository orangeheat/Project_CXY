# TCGA_BRCA

<h3>setup环境</h3>

在空文件夹中，创建一个proj，在这个文件夹中建"Data", "R", "Results"三个子文件夹，然后把相关文件放到对应文件夹中
create a new project in a new folder -> create three folders named "Data", "R", "Results",respectively -> move documents to their folders
<img width="360" alt="Screen Shot 2024-11-09 at 16 54 43" src="https://github.com/user-attachments/assets/f6952074-e992-4765-a3e0-e1d7af8e165d">

<h3>文件内容</h3>
R文件夹里都是后缀为.R的脚本，Data文件夹是初始数据，Results文件夹里是结果
<p>
RawTPM.csv：只包含处理后的数据（提取了标准化后的TPM数据，提取重复基因的最大值，去除80%样本中未表达的基因）
RawCOUNT.csv：用来差异分析的COUNT数据
lncRNA_matrix.csv：从RawTPM中提取了gene_type为lncRNA的基因表达矩阵
MG_matrix.csv：从RawTPM中提取了m6A/m5C相关的基因表达矩阵</p>

