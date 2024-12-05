#Setup
library(here)
library(survival)
library(dplyr)
library(caret)
library(glmnet)
library(survminer)
rm(list = ls())
#read previously made matrix
meta <- read.csv("Results/clinical_expression.csv", row.names = 1)
MGR_lncRNA <- meta |>
  select(-c(1:20)) #去除临床数据，仅提取表达数据
gene <- colnames(MGR_lncRNA)

#分组
set.seed(18679)
index <- createDataPartition(meta$OS, p = 0.6, list = FALSE, times = 1)
train <- meta[index, ]
test <- meta[-index, ]

#连续变量的t检验
continuous_variables <- c("Age", "OS.time")
t_test_results <- lapply(continuous_variables, function(var) {
  t_test_result <- t.test(train[[var]], test[[var]])
  return(data.frame(Variable = var, Test = "T-Test", 
                    Statistic = t_test_result$statistic, 
                    PValue = t_test_result$p.value))
})

print(t_test_results)
write.csv(train, "./Results/train.csv", row.names = T)
write.csv(test, "./Results/test.csv", row.names = T)

#uni-cox in train
# 构建生存数据表
mysurv <- Surv(train$OS.time, train$OS)
# 计算P值的函数
Unicox_PValue <- function(exp, genes) {
  pvalues <- sapply(genes, function(gene) {
    fml <- as.formula(paste0('mysurv ~ ', gene))
    cox <- summary(coxph(fml, exp))
    round(cox$coefficients[1, 5], 4) #P值
  })
  data.frame(Gene = genes, PValue = pvalues)
}

# 提取信息的函数
Unicox_Detail <- function(exp, genes) {
  results <- lapply(genes, function(gene) {
    fml <- as.formula(paste0('mysurv ~ ', gene))
    cox <- summary(coxph(fml, exp))
    data.frame(
      Gene = gene,
      HR = round(cox$coefficients[1, 2], 2),
      lower95 = round(cox$conf.int[1, 3], 2),
      upper95 = round(cox$conf.int[1, 4], 2),
      PValue = round(cox$coefficients[1, 5], 4)
    )
  })
  do.call(rbind, results)
}

# 筛选P值小于0.05的基因
gene_pvalues <- Unicox_PValue(train, gene)
significant_genes <- gene_pvalues |>
  filter(PValue < 0.05) |>
  pull(Gene)
Univar <- Unicox_Detail(train, significant_genes)

#lasso after unicox
x <- as.matrix(meta[ ,colnames(meta) %in% gene])
y <- as.matrix(Surv(meta$OS.time,meta$OS))
set.seed(89435)
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1)
fit <-  glmnet(x, y, family = "cox", alpha = 1, label = TRUE)

# 使用最优 lambda 的系数
coef_lasso <- coef(fit, s = cv_fit$lambda.min)
lasso_genes <- rownames(coef_lasso)[as.numeric(coef_lasso) != 0]
lasso_coefficients <- as.numeric(coef_lasso)[as.numeric(coef_lasso) != 0]
print(data.frame(Gene = lasso_genes, Coefficient = lasso_coefficients))

jpeg("./Results/lasso_1.jpg", width = 2000, height = 1500, res = 300)
plot(cv_fit)
dev.off()
jpeg("./Results/lasso_2.jpg", width = 2000, height = 1500, res = 300)
plot(fit)
dev.off()


#multi-cox
lasso_exp <- train |>
  select(OS, OS.time, lasso_genes)
multiCox <- coxph(Surv(OS.time, OS) ~ ., data =  lasso_exp)
cox_summary <- summary(multiCox)
#绘制森林图
multicox_ggforest <- ggforest(multiCox,  #coxph得到的Cox回归结果
                              data = lasso_exp,  #数据集
                              main = 'Hazard ratio of multi cox',  #标题
                              cpositions = c(0.05, 0.15, 0.35),  #前三列距离
                              fontsize = 0.8, #字体大小
                              refLabel = 'reference', #相对变量的数值标签，也可改为1
                              noDigits = 3 #保留HR值以及95%CI的小数位数
)
jpeg("./Results/multicox.jpg", width = 3000, height = 2500, res = 300)
multicox_ggforest
dev.off()

#筛选基因
model_gene <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]
model_exp <- lasso_exp[, c("OS", "OS.time", rownames(model_gene))]
#建模
multiCox1 <- coxph(Surv(OS.time, OS) ~ ., data =  model_exp)
# 提取 summary 结果
multiCox_summary <- summary(multiCox1)

# 创建包含系数、置信区间和 p 值的表格
multiCox_results <- data.frame(
  Variable = rownames(multiCox_summary$coefficients),
  Coefficient = multiCox_summary$coefficients[, "coef"],
  HR = exp(multiCox_summary$coefficients[, "coef"]),
  Lower95 = exp(multiCox_summary$conf.int[, "lower .95"]),
  Upper95 = exp(multiCox_summary$conf.int[, "upper .95"]),
  PValue = multiCox_summary$coefficients[, "Pr(>|z|)"]
)

# 查看结果表格
print(multiCox_results)

# 保存到文件
write.csv(multiCox_results, "./Results/multiCox_results.csv", row.names = FALSE)

#predict函数计算训练集风险评分
model_train <- model_exp |>
  mutate(riskScore=predict(multiCox1,type="risk",newdata=model_exp),
         risk=ifelse(riskScore > median(riskScore), "High", "Low")) |>
  select(OS, OS.time, riskScore, risk, everything())

#计算测试集风险评分
model_test <- select(test, c(OS, OS.time, rownames(model_gene)))|>
  mutate(riskScore=predict(multiCox1,type="risk",newdata=test),
         risk=ifelse(riskScore > median(riskScore), "High", "Low")) |>
  select(OS, OS.time, riskScore, risk, everything())
#整集风险评分
model_whole <- select(meta, c(OS, OS.time, rownames(model_gene)))|>
  mutate(riskScore=predict(multiCox1,type="risk",newdata=meta),
         risk=ifelse(riskScore > median(riskScore), "High", "Low")) |>
  select(OS, OS.time, riskScore, risk, everything())

write.csv(model_train,'Results/model_train.csv',row.names = TRUE)
write.csv(model_test,'Results/model_test.csv',row.names = TRUE)
write.csv(model_whole,'Results/model_whole.csv',row.names = TRUE)

#KM survival
#规定函数，需输入绘图数据，图上标题，文件名和路径以及画面截止的范围
km_plot <- function(data, plot_title, file_name, x_limit = c(0, 15), file_path = "./Results/") {
  surv_fit <- survfit(Surv(OS.time, OS) ~ risk, data = data)
  plot <- ggsurvplot(
    surv_fit, 
    data = data,
    pval = TRUE,
    risk.table = TRUE, 
    surv.median.line = "hv",
    legend.labs = c("High risk", "Low risk"),
    legend.title = "Risk",
    title = plot_title,
    ylab = "Cumulative Survival (Percentage)", 
    xlab = "Time (Years)",
    censor.shape = 124, censor.size = 2, 
    conf.int = FALSE,
    break.x.by = 2,
    xlim = x_limit
  )
  jpeg(filename = paste0(file_path, file_name, ".jpeg"), width = 1500, height = 1500, res = 300)
  print(plot) 
  dev.off() 
}


km_plot(data = model_train, 
        plot_title = "Overall Survival in Training Set", 
        file_name = "km_train")
km_plot(data = model_test,
        plot_title = "Overall Survival in Test Set",
        file_name = "km_test")
km_plot(data = model_whole,
        plot_title = "Overall Survival in the Entire Set",
        file_name = "km_whole")

