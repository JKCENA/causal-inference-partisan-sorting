library(boot)        
library(dplyr)
library(ggplot2)
library(ResourceSelection)
library(pROC)
library(car)

df <- read.csv("C:/Users/jihun/OneDrive/바탕 화면/데이터 프로그램/vscode/작업장/df_ANES.csv")
df <- df[!is.na(df$affective_polarization_2016), ]
df$age_group_4    <- as.factor(df$age_group_4)
df$race_collapsed <- as.factor(df$race_collapsed)
df$region         <- as.factor(df$region)


#Propensity Score Model
ps.formula <- sorted_2016 ~ pol_interest_index + participation + 
  social_political_days + scaled_channel_count + 
  age_group_4 + race_collapsed + region + affective_polarization_2016

ps.model   <- glm(ps.formula, family = binomial(), data = df)
df$ps    <- predict(ps.model, type = "response")
n <- nobs(ps.model)
cat("총 관측치 수:", n, "\n")
summary(ps.model)

#Wald Chi2
lh <- linearHypothesis(ps.model,
                       names(coef(ps.model))[-1],
                       test = "Chisq")

print(lh)
pearson_chi2 <- sum(residuals(ps.model, type = "pearson")^2)
df_resid <- df.residual(ps.model)
p_value_pearson <- 1 - pchisq(pearson_chi2, df_resid)
data.frame(
  Pearson_Chi2 = pearson_chi2,
  df           = df_resid,
  p_value      = p_value_pearson
)


# Hosmer-Lemeshow test
hl_test <- hoslem.test(ps.model$y, fitted(ps.model), g = 10)
print(hl_test)

# ROC curve 및 AUC
roc_obj <- roc(ps.model$y, fitted(ps.model))
plot(roc_obj, main = "ROC Curve for PS Model")
auc_value <- auc(roc_obj)
print(auc_value)

ggplot(df, aes(x = ps, fill = as.factor(sorted_2016))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(fill = "Sorted_2016", x = "Propensity Score", y = "Count") +
  ggtitle("Propensity Score Distribution by Treatment Group")

library(tidyverse)
df_trim <- df %>% 
  filter(ps >= 0.4, ps <= 0.8)
cat("Trimmed sample size:", nrow(df_trim), "\n")

ps.model   <- glm(ps.formula, family = binomial(), data = df_trim)
df_trim$ps    <- predict(ps.model, type = "response")
summary(ps.model)
n <- nobs(ps.model)
cat("총 관측치 수:", n, "\n")

ggplot(df_trim, aes(x = ps, fill = as.factor(sorted_2016))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(fill = "Sorted_2016", x = "Propensity Score", y = "Count") +
  ggtitle("Propensity Score Distribution by Treatment Group")

#VIF
vif_results <- vif(ps.model)
print(vif_results)


#Stab
pA <- mean(df_trim$sorted_2016)
df_trim$w_stab <- with(df_trim, ifelse(sorted_2016==1, pA/ps, (1-pA)/(1-ps)))
summary(df_trim$w_stab)
hist(df_trim$w_stab, breaks = 50, main = "Distribution of stab_IPTW", xlab = "IPTW")

cutoff <- quantile(df_trim$w_stab, 0.95)
df_trim$w_stab_trimmed <- ifelse(df_trim$w_stab > cutoff, cutoff, df_trim$w_stab)
df_trim$w_stab_trimmed <- with(df_trim, ifelse(sorted_2016==1, pA/ps, (1-pA)/(1-ps)))
summary(df_trim$w_stab_trimmed) 

#out.formula
out.formula <- affective_polarization_2020 ~ sorted_2016 + 
  pol_interest_index + participation +
  social_political_days + scaled_channel_count +
  age_group_4 + race_collapsed + region + affective_polarization_2016

#IPW 
ipw.model_stab <- lm(out.formula, data = df_trim, weights = w_stab)
coef(ipw.model_stab)["sorted_2016"]
library(sandwich)
library(lmtest)
confint(ipw.model_stab, "sorted_2016", level = 0.95)
ipw.model_stab <- lm(out.formula, data = df_trim, weights = w_stab)

robust_vcov <- vcovHC(ipw.model_stab, type = "HC3") 

se_robust <- sqrt(robust_vcov["sorted_2016", "sorted_2016"])

beta_hat <- coef(ipw.model_stab)["sorted_2016"]

df_resid <- df.residual(ipw.model_stab)

t_crit <- qt(0.975, df = df_resid)

lower_robust <- beta_hat - t_crit * se_robust
upper_robust <- beta_hat + t_crit * se_robust

cbind(
  Estimate = beta_hat,
  `Lower 95% (robust)` = lower_robust,
  `Upper 95% (robust)` = upper_robust
)

#SMD
library(cobalt)
bal <- bal.tab(ps.model,
               weights = df_trim$w_stab,
               estimand = "ATE",
               un = TRUE)
print(bal)        # SMD 표
love.plot(bal)    # 균형 시각화


out.model   <- lm(out.formula, data = df_trim)
summary(out.model)
nobs(out.model)
s <- summary(out.model)
cat("Root MSE:", s$sigma, "\n")
cat("R-squared:", s$r.squared, "\n")
cat("Adjusted R-squared:", s$adj.r.squared, "\n")
cat("F-statistic:", s$fstatistic[1], "on", s$fstatistic[2], "and", s$fstatistic[3], "DF\n")

f_val <- s$fstatistic[1]
df1   <- s$fstatistic[2]  # numerator df (설명변수 수)
df2   <- s$fstatistic[3]  # denominator df (잔차 자유도)


p_val <- pf(f_val, df1, df2, lower.tail = FALSE)


cat("F-statistic:", f_val, "on", df1, "and", df2, "DF\n")
cat("F-statistic p-value:", p_val, "\n")

#Parametic G-Formula
df_trim$mu1  <- predict(out.model, newdata = transform(df_trim, sorted_2016=1),  type="response")
df_trim$mu0  <- predict(out.model, newdata = transform(df_trim, sorted_2016=0),  type="response")
mean(df_trim$mu1 - df_trim$mu0)

#AIPW
df_trim$aipw <- with(df_trim,
                mu1 - mu0 +
                  ifelse(sorted_2016==1, (affective_polarization_2020 - mu1)/ps, (mu0 - affective_polarization_2020)/(1-ps))
)

# ATE
ate_aipw <- mean(df_trim$aipw)
ate_aipw

n    <- nrow(df_trim)
ate  <- ate_aipw                     
sdev <- sd(df_trim$aipw)           
se   <- sdev / sqrt(n)             


lower_norm <- ate - qnorm(0.975) * se
upper_norm <- ate + qnorm(0.975) * se

c(
  Estimate        = ate,
  `SE (analytic)` = se,
  `Lower 95%`     = lower_norm,
  `Upper 95%`     = upper_norm
)

library(EValue)
library(lmtest)
library(sandwich)
sdev_aipw <- sd(df_trim$aipw, na.rm = TRUE)
se_aipw   <- sdev_aipw / sqrt(n)
sd_y <- sd(df_trim$affective_polarization_2020, na.rm = TRUE)


e_val_results <- evalues.OLS(
  est = ate,
  se  = se_aipw,
  sd  = sd_y
)
print(e_val_results)

#Boot
library(sandwich)
dr_ate_stab <- function(data, indices) {
  d <- data[indices, ]  
  
  
  ps_mod <- glm(sorted_2016 ~ 
                  pol_interest_index +
                  participation +
                  social_political_days +
                  scaled_channel_count +
                  age_group_4 +
                  race_collapsed +
                  affective_polarization_2016 +
                  region,
                family = binomial(), data = d)
  d$ps <- predict(ps_mod, type = "response")
  
 
  pA <- mean(d$sorted_2016)  
  d$stab_ipw <- ifelse(d$sorted_2016 == 1,
                       pA / d$ps,
                       (1 - pA) / (1 - d$ps))
  
  
  out_mod <- lm(affective_polarization_2020 ~ sorted_2016 +
                   pol_interest_index +
                   participation +
                   social_political_days +
                   scaled_channel_count +
                   age_group_4 +
                   race_collapsed +
                   affective_polarization_2016 +
                   region,
                   data = d)
  
  
  mu0 <- predict(out_mod, newdata = transform(d, sorted_2016 = 0), type = "response")
  mu1 <- predict(out_mod, newdata = transform(d, sorted_2016 = 1), type = "response")
  
  
  dr_terms <- (mu1 - mu0) +
    ifelse(
      d$sorted_2016 == 1,
      pA * (d$affective_polarization_2020 - mu1) / d$ps,
      (1 - pA) * (d$affective_polarization_2020 - mu0) / (1 - d$ps)
    )
  
  return(mean(dr_terms))
}


set.seed(2025)
boot_res_stab <- boot(data = df_trim,
                      statistic = dr_ate_stab,
                      R = 5000,
                      parallel = "multicore",  # 병렬처리
                      ncpus = 4)

print(boot_res_stab)
# BCa CI
boot.ci(boot_res_stab, type = "bca")


boot_estimates <- boot_res_stab$t


point_estimate <- mean(boot_estimates)

ci_bca <- boot.ci(boot_res_stab, type = "bca")$bca[4:5]
ggplot(data = data.frame(ATE = boot_estimates), aes(x = ATE)) +
  geom_histogram(aes(y = ..density..), bins = 30, 
                 fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "blue", size = 1.2) +
  geom_vline(xintercept = point_estimate, color = "red", 
             linetype = "solid", size = 1.2) +
  geom_vline(xintercept = ci_bca, color = "darkgreen", 
             linetype = "dashed", size = 1.2) +
  annotate("text", x = point_estimate, y = 0, label = paste0("Mean: ", round(point_estimate, 2)),
           vjust = -1, color = "red", size = 4) +
  annotate("text", x = ci_bca[1], y = 0, label = paste0("Lower BCa: ", round(ci_bca[1], 2)),
           vjust = 1.5, hjust = 1.1, color = "darkgreen", size = 3.5) +
  annotate("text", x = ci_bca[2], y = 0, label = paste0("Upper BCa: ", round(ci_bca[2], 2)),
           vjust = 1.5, hjust = -0.1, color = "darkgreen", size = 3.5) +
  labs(title = "Bootstrap Distribution of ATE (Stabilized Weight, BCa CI)",
       x = "ATE Estimate",
       y = "Density") +
  theme_minimal(base_size = 14)

