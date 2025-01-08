# setup
library(here)
library(dplyr)
library(ggplot2)
library(survminer)
library(timeROC)
library(survival)
rm(list = ls())
# 读文件，合并matrix
model_whole <- read.csv('Results/model_whole.csv', row.names = 1)
meta <- read.csv("Results/clinical_expression.csv", row.names = 1) 
train <- read.csv('Results/model_train.csv', row.names = 1)
meta_1 <- cbind(model_whole, meta) 
meta_1 <- meta_1[, !duplicated(colnames(meta_1))]|>
  select(OS, OS.time, riskScore, risk, Age, Stage_num, T_num, M_num, N_num)
train <- meta_1[rownames(train),]
test <- meta_1[setdiff(rownames(meta_1), rownames(train)), ]

# 还是规定函数，输入绘图数据，标题，文件名和路径以及画面截止的范围
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



# 分组stage
stage_I_II <- subset(meta_1, Stage_num <=2)
km_plot(data = stage_I_II, 
        plot_title = "Stage I-II", 
        file_name = "km_stage_I_II")
stage_III_IV <- subset(meta_1, Stage_num >2)
km_plot(data = stage_III_IV, 
        plot_title = "Stage III-IV", 
        file_name = "km_stage_III_IV")

# Age
Age_under65 <- subset(meta_1, Age <=65)
km_plot(data = Age_under65, 
        plot_title = "Age under 65", 
        file_name = "Age_under65")
Age_over65 <- subset(meta_1, Age >65)
km_plot(data = Age_over65, 
        plot_title = "Age over 65", 
        file_name = "Age_over65")
# T
T1_2 <- subset(meta_1, T_num <=2)
km_plot(data = T1_2, 
        plot_title = "T1-2", 
        file_name = "T1_2")
T3_4 <- subset(meta_1, T_num >2)
km_plot(data = T3_4, 
        plot_title = "T3-4", 
        file_name = "T3_4")
# M
M0 <- subset(meta_1, M_num =0)
km_plot(data = M0, 
        plot_title = "No distant metastasis (M0)", 
        file_name = "M0")
M1 <- subset(meta_1, M_num =1)
km_plot(data = M1, 
        plot_title = "Confirmed distant metastasis (M1)", 
        file_name = "M1")

# N
N0 <- subset(meta_1, N_num =0)
km_plot(data = N0, 
        plot_title = "No regional lymph node involvement (N0)", 
        file_name = "N0")
N123 <- subset(meta_1, N_num >0)
km_plot(data = N123, 
        plot_title = "Lymph node involvement (N1-3)", 
        file_name = "N123")


#画ROC曲线验证风险评分和其他临床因素的预测性能write a function here, present every factor's ROC
time_ROC <- function(data, file_name, plot_title) {
  
  ROC_riskscore <- timeROC(
    T = data$OS.time,
    delta = data$OS,
    marker = data$riskScore,
    cause = 1,
    weighting = "marginal",
    times = c(1, 3, 5),
    ROC = TRUE,
    iid = TRUE
  )
  jpeg(paste0(file_name, ".jpeg"), width = 2000, height = 1500, res = 300)
  
  plot(ROC_riskscore, time = 1, col = "#FFBBFF", lwd = 2, add = FALSE, title = FALSE)
  plot(ROC_riskscore, time = 3, col = "#A4D3EE", lwd = 2, add = TRUE)
  plot(ROC_riskscore, time = 5, col = "#8B668B", lwd = 2, add = TRUE)
  title(main = plot_title)
  legend("bottomright", c("1-Year", "3-Year", "5-Year"), col = c("#FFBBFF", "#A4D3EE", "#8B668B"), lty = 1, lwd = 2)
  text(0.5, 0.2, paste("1-Year AUC =", round(ROC_riskscore$AUC[1], 3)))
  text(0.5, 0.15, paste("3-Year AUC =", round(ROC_riskscore$AUC[2], 3)))
  text(0.5, 0.1, paste("5-Year AUC =", round(ROC_riskscore$AUC[3], 3)))
  
  dev.off()
}

time_ROC(train, "./Results/train_time_ROC", "TimeROC in train set")
time_ROC(test, "./Results/test_time_ROC", "TimaROC in test set")
time_ROC(meta_1, "./Results/wholeset_time_ROC", "TimaROC in whole dataset")


# timeROC表现不好，从临床因素鉴中定其他独立预后因子
multiCox_train <- coxph(Surv(OS.time, OS) ~ riskScore + Age + Stage_num + T_num + M_num + N_num, data =  train)
multiCox_test <- coxph(Surv(OS.time, OS) ~ riskScore + Age + Stage_num + T_num + M_num + N_num, data =  test)
multiCox_all <- coxph(Surv(OS.time,OS) ~ riskScore + Age + Stage_num + T_num + M_num + N_num, data = meta_1)

summary(multiCox_train)
summary(multiCox_test)
summary(multiCox_all)
multiCox_summary <- summary(multiCox_all)
multiCox_results <- data.frame(
  Variable = rownames(multiCox_summary$coefficients),
  Coefficient = multiCox_summary$coefficients[, "coef"],
  HR = exp(multiCox_summary$coefficients[, "coef"]),
  Lower95 = exp(multiCox_summary$conf.int[, "lower .95"]),
  Upper95 = exp(multiCox_summary$conf.int[, "upper .95"]),
  PValue = multiCox_summary$coefficients[, "Pr(>|z|)"]
)

write.csv(multiCox_results, "./Results/independence_verification.csv", row.names = FALSE)
# 发现风险评分和年龄在训练集测试集和整集中都是独立预后预测因子


####### 把年龄也纳入multicox模型 ######

# 构建 Cox 模型
train_new <- read.csv('Results/model_train.csv', row.names = 1) |>
  mutate(Age= train$Age) |>
  select(-c(risk, riskScore))
cox_model <- cph(Surv(OS.time, OS) ~ ., 
                 data = train_new, x = TRUE, y = TRUE, surv = TRUE)
train_new <- train_new|>
  mutate(riskScore=predict(cox_model,type="risk",newdata=train_new))
test_new <- test |>
  mutate(riskScore=prediect(cox_model, type="risk", newdata= test))
time_ROC(train_new, "./Results/train_time_ROC_new", "TimeROC in train set")
time_ROC(test_new, "./Results/test_time_ROC_new", "TimaROC in test set")
time_ROC(meta_1, "./Results/wholeset_time_ROC", "TimaROC in whole dataset")

# 绘制 Nomogram
regplot(cox_model,
        observation = train_new[24,],
        title = "Nomogram",
        clickable = TRUE,  # 允许点击查看详情
        failtime = c(1, 3, 5))  # 显示风险比
regplot(cox_model, 
        observation = meta_1[500,],
        interval ="confidence", 
        title="Nomogram",
        plots=c("violin", "boxes"), 
        clickable = T,
        failtime = c(1,3,5)) #设置随访时间1年、3年和5年

