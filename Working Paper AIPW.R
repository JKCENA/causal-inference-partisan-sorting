library(boot)        # 부트스트랩
library(dplyr)
library(ggplot2)
library(ResourceSelection)
library(pROC)
library(car)
library(tidyverse)

df <- read.csv("C:/Users/jihun/OneDrive/바탕 화면/데이터 프로그램/vscode/작업장/df_ANES.csv")
df <- df[!is.na(df$affective_polarization_2016), ]
df$age_group_4    <- as.factor(df$age_group_4)
df$race_collapsed <- as.factor(df$race_collapsed)
df$region         <- as.factor(df$region)
#Propensity Score 모델 만들기
ps.formula <- sorted_2016 ~ pol_interest_index + participation + 
  social_political_days + scaled_channel_count + 
  age_group_4 + race_collapsed + region 

ps.model   <- glm(ps.formula, family = binomial(), data = df)
df$ps    <- predict(ps.model, type = "response")

#Wald Chi2
lh <- linearHypothesis(ps.model,
                       names(coef(ps.model))[-1],
                       test = "Chisq")

print(lh)

#Gof 계산을 위한 pearsos chi2

pearson_chi2 <- sum(residuals(ps.model, type = "pearson")^2)

# 2) 자유도: 관측치 수 − 추정 파라미터 수
df_resid <- df.residual(ps.model)

# 3) p‑value 계산
p_value_pearson <- 1 - pchisq(pearson_chi2, df_resid)

# 4) 결과 출력
data.frame(
  Pearson_Chi2 = pearson_chi2,
  df           = df_resid,
  p_value      = p_value_pearson
)


# propensity score model 진단
# Hosmer-Lemeshow test
hl_test <- hoslem.test(ps.model$y, fitted(ps.model), g = 10)
print(hl_test)

# ROC curve 및 AUC
roc_obj <- roc(ps.model$y, fitted(ps.model))
plot(roc_obj, main = "ROC Curve for PS Model")
auc_value <- auc(roc_obj)
print(auc_value)

#  다중공선성 (VIF) 점검
vif_results <- vif(ps.model)
print(vif_results)


# Propensity Score 분포를 히스토그램이나 밀도함수로 시각화
ggplot(df, aes(x = ps, fill = as.factor(sorted_2016))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(fill = "Sorted_2016", x = "Propensity Score", y = "Count") +
  ggtitle("Propensity Score Distribution by Treatment Group")

# 1) PS 컷오프(트리밍)
df_trim <- df %>% 
  filter(ps >= 0.3, ps <= 0.7)
cat("Trimmed sample size:", nrow(df_trim), "\n")

#Trim 이후 분포도 재점검
ggplot(df_trim, aes(x = ps, fill = as.factor(sorted_2016))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(fill = "Sorted_2016", x = "Propensity Score", y = "Count") +
  ggtitle("Propensity Score Distribution by Treatment Group")

pA <- mean(df_trim$sorted_2016)
df_trim$w_stab <- with(df_trim, ifelse(sorted_2016==1, pA/ps, (1-pA)/(1-ps)))
summary(df_trim$w_stab)
hist(df_trim$w_stab, breaks = 50, main = "Distribution of stab_IPTW", xlab = "IPTW")

cutoff <- quantile(df_trim$w_stab, 0.95)
# FIX: use pmin so truncation isn't immediately overwritten
df_trim$w_stab_trimmed <- pmin(df_trim$w_stab, cutoff)
summary(df_trim$w_stab_trimmed) #결과상 차이 없음

#균형 확인 
library(cobalt)
bal <- bal.tab(
  x = df_trim[, c("pol_interest_index", "participation", "social_political_days", 
                  "scaled_channel_count", "age_group_4", "race_collapsed", "region")],
  treat = df_trim$sorted_2016,
  weights = df_trim$w_stab,
  estimand = "ATE",
  un = TRUE
)

print(bal)        # SMD 표
love.plot(bal)    # 균형 시각화






#AIPW_short
out.short <- lm(affective_polarization_2016 ~ sorted_2016 +
                  pol_interest_index + participation +
                  social_political_days + scaled_channel_count +
                  age_group_4 + race_collapsed + region,
                data = df_trim, weights = w_stab)
summary(out.short)
ate_short <- coef(out.short)["sorted_2016"]
ate_short

##Parametic G-Formula_Short
df_trim$mu1_short  <- predict(out.short, newdata = transform(df_trim, sorted_2016=1),  type="response")
df_trim$mu0_short  <- predict(out.short, newdata = transform(df_trim, sorted_2016=0),  type="response")
mean(df_trim$mu1_short - df_trim$mu0_short)

##AIPW_Stab_short
# FIX: control-term must have a negative sign so augmentation correctly subtracts control residuals
df_trim$aipw_stab_short <- with(df_trim,
                                mu1_short - mu0_short +
                                  ifelse(sorted_2016 == 1,
                                         pA * (affective_polarization_2016 - mu1_short) / ps,
                                         - (1 - pA) * (affective_polarization_2016 - mu0_short) / (1 - ps)
                                  )
)

## 전체 평균 ATE_short
ate_aipw_short <- mean(df_trim$aipw_stab_short)
ate_aipw_short


#AIPW_Total
out.total <- lm(affective_polarization_2020 ~ sorted_2016 +
                  pol_interest_index + participation +
                  social_political_days + scaled_channel_count +
                  age_group_4 + race_collapsed + region,
                data = df_trim, weights = w_stab)
ate_total <- coef(out.total)["sorted_2016"]
ate_total

##Parametic G-Formula_Total
df_trim$mu1_total  <- predict(out.total, newdata = transform(df_trim, sorted_2016=1),  type="response")
df_trim$mu0_total  <- predict(out.total, newdata = transform(df_trim, sorted_2016=0),  type="response")
mean(df_trim$mu1_total - df_trim$mu0_total)

##AIPW_Stab_Total
# FIX: control-term must have a negative sign here as well
df_trim$aipw_stab_total <- with(df_trim,
                                mu1_total - mu0_total +
                                  ifelse(sorted_2016 == 1,
                                         pA * (affective_polarization_2020 - mu1_total) / ps,
                                         - (1 - pA) * (affective_polarization_2020 - mu0_total) / (1 - ps)
                                  )
)

## 전체 평균 ATE_Total
ate_aipw_total <- mean(df_trim$aipw_stab_total)
ate_aipw_total





#Bootstrap
library(boot)
library(ggplot2)

#-----------------------------
# 1) Short (2016 outcome)
#-----------------------------
dr_ate_stab_short <- function(data, indices) {
  d <- data[indices, ]
  
  out_mod <- lm(affective_polarization_2016 ~ sorted_2016 +
                  pol_interest_index + participation +
                  social_political_days + scaled_channel_count +
                  age_group_4 + race_collapsed + region,
                data = d, weights = d$w_stab)
  
  mu0 <- predict(out_mod, newdata = transform(d, sorted_2016 = 0), type = "response")
  mu1 <- predict(out_mod, newdata = transform(d, sorted_2016 = 1), type = "response")
  
  pA <- mean(d$sorted_2016)
  dr_terms <- (mu1 - mu0) +
    ifelse(d$sorted_2016 == 1,
           pA * (d$affective_polarization_2016 - mu1) / d$ps,
           - (1 - pA) * (d$affective_polarization_2016 - mu0) / (1 - d$ps))  # FIX: added minus
  
  mean(dr_terms)
}

#-----------------------------
# 2) Total (2020 outcome)
#-----------------------------
dr_ate_stab_total <- function(data, indices) {
  d <- data[indices, ]
  
  out_mod <- lm(affective_polarization_2020 ~ sorted_2016 +
                  pol_interest_index + participation +
                  social_political_days + scaled_channel_count +
                  age_group_4 + race_collapsed + region,
                data = d, weights = d$w_stab)
  
  mu0 <- predict(out_mod, newdata = transform(d, sorted_2016 = 0), type = "response")
  mu1 <- predict(out_mod, newdata = transform(d, sorted_2016 = 1), type = "response")
  
  pA <- mean(d$sorted_2016)
  dr_terms <- (mu1 - mu0) +
    ifelse(d$sorted_2016 == 1,
           pA * (d$affective_polarization_2020 - mu1) / d$ps,
           - (1 - pA) * (d$affective_polarization_2020 - mu0) / (1 - d$ps))  # FIX: added minus
  
  mean(dr_terms)
}



#-----------------------------
# 4) Bootstrap run
#   - Windows: snow / Others: multicore
#-----------------------------
par_type <- ifelse(.Platform$OS.type == "windows", "snow", "multicore")
set.seed(2025)

boot_short <- boot(data = df_trim, statistic = dr_ate_stab_short,
                   R = 2000, parallel = par_type, ncpus = 4)
boot_total <- boot(data = df_trim, statistic = dr_ate_stab_total,
                   R = 2000, parallel = par_type, ncpus = 4)


print(boot_short); print(boot_total)

#-----------------------------
# 5) BCa 95% CI
#-----------------------------
ci_short_bca <- boot.ci(boot_short, type = "bca")$bca[4:5]
ci_total_bca <- boot.ci(boot_total, type = "bca")$bca[4:5]


cat("\n[BCa 95% CI]\n")
cat(sprintf("Short  : (%.3f, %.3f)\n", ci_short_bca[1], ci_short_bca[2]))
cat(sprintf("Total  : (%.3f, %.3f)\n", ci_total_bca[1], ci_total_bca[2]))


#-----------------------------
# 6) Quick plots (optional)
#-----------------------------
plot_boot <- function(boot_obj, ttl) {
  est <- boot_obj$t
  pe  <- mean(est)
  ci  <- boot.ci(boot_obj, type = "bca")$bca[4:5]
  
  # 분포 기준으로 y위치 자동 설정
  dens <- density(est)
  y_top <- max(dens$y) * 0.25  # Mean 숫자 위치
  y_ci  <- y_top * 0.6         # CI 숫자 위치
  
  ggplot(data = data.frame(ATE = est), aes(x = ATE)) +
    geom_histogram(aes(y = ..density..), bins = 30,
                   fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(size = 1.2, color = "blue") +
    geom_vline(xintercept = pe, linetype = "solid", size = 1.2, color = "red") +
    geom_vline(xintercept = ci, linetype = "dashed", size = 1.1, color = "darkgreen") +
    
    # 숫자 직접 표시
    annotate("text", x = pe, y = y_top,
             label = paste0("Mean: ", round(pe, 3)),
             color = "red", size = 5.5, fontface = "bold") +
    annotate("text", x = ci[1], y = y_ci,
             label = paste0("Lower: ", round(ci[1], 3)),
             color = "darkgreen", size = 4.8, hjust = 1.2) +
    annotate("text", x = ci[2], y = y_ci,
             label = paste0("Upper: ", round(ci[2], 3)),
             color = "darkgreen", size = 4.8, hjust = -0.1) +
    
    labs(title = paste0(ttl, " — Bootstrap AIPW Distribution"),
         x = "ATE", y = "Density") +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# 예시 실행
p_short <- plot_boot(boot_short, "Short")
p_total <- plot_boot(boot_total, "Total")


# 출력
print(p_short)
print(p_total)




#기술통계 부분적으로 확인해보기 
mean(df_trim$affective_polarization_2016[df_trim$sorted_2016 == 1], na.rm = TRUE)
mean(df_trim$affective_polarization_2016[df_trim$sorted_2016 == 0], na.rm = TRUE)
mean(df_trim$affective_polarization_2020[df_trim$sorted_2016 == 1], na.rm = TRUE)
mean(df_trim$affective_polarization_2020[df_trim$sorted_2016 == 0], na.rm = TRUE)

df_trim <- df_trim %>% mutate(delta_aff_16_20 = affective_polarization_2016 - affective_polarization_2020)
df_trim %>%
  filter(sorted_2016 == 1) %>%
  summarise(
    total = sum(!is.na(delta_aff_16_20)),
    positive = sum(delta_aff_16_20 > 0, na.rm = TRUE),
    zero     = sum(delta_aff_16_20 == 0, na.rm = TRUE),
    negative = sum(delta_aff_16_20 < 0, na.rm = TRUE)
  ) %>%
  mutate(
    pct_pos = round(positive / total * 100, 1),
    pct_zero = round(zero / total * 100, 1),
    pct_neg = round(negative / total * 100, 1)
  )

df_trim %>%
  filter(sorted_2016 == 0) %>%
  summarise(
    total = sum(!is.na(delta_aff_16_20)),
    positive = sum(delta_aff_16_20 > 0, na.rm = TRUE),
    zero     = sum(delta_aff_16_20 == 0, na.rm = TRUE),
    negative = sum(delta_aff_16_20 < 0, na.rm = TRUE)
  ) %>%
  mutate(
    pct_pos = round(positive / total * 100, 1),
    pct_zero = round(zero / total * 100, 1),
    pct_neg = round(negative / total * 100, 1)
  )



library(sandwich)
library(lmtest)
#HC3를 사용한 강건성 check
fit_short <- lm(affective_polarization_2016 ~ sorted_2016 +
            pol_interest_index +
            participation +
            social_political_days +
            scaled_channel_count +
            age_group_4 +
            race_collapsed +
            region,
          data = df_trim,
          weights = w_stab)
coeftest(fit_short, vcov = vcovHC(fit_short, type = "HC3"))

fit_total <- lm(affective_polarization_2020 ~ sorted_2016 +
                  pol_interest_index +
                  participation +
                  social_political_days +
                  scaled_channel_count +
                  age_group_4 +
                  race_collapsed +
                  region,
                data = df_trim,
                weights = w_stab)
coeftest(fit_total, vcov = vcovHC(fit_total, type = "HC3"))



# --- E-value: HC3 robust SE 사용해서 계산 ---
library(EValue)
library(sandwich)  # 이미 로드되어 있으면 중복 로드되어도 괜찮음

# 종속변수 표준편차 (이미 계산해두셨다면 재계산 불필요)
sd_y <- sd(df_trim$affective_polarization_2020, na.rm = TRUE)
sd_y_short <- sd(df_trim$affective_polarization_2016, na.rm = TRUE)

# HC3 견고 공분산행렬에서 표준오차 추출
vcov_total_hc3 <- vcovHC(fit_total, type = "HC3")
vcov_short_hc3 <- vcovHC(fit_short, type = "HC3")

se_robust_total <- sqrt(diag(vcov_total_hc3))["sorted_2016"]
se_robust_short <- sqrt(diag(vcov_short_hc3))["sorted_2016"]

# 만약 이름이 정확히 맞지 않아 NA가 나오면 위치로 뽑아도 됩니다:
if (is.na(se_robust_total)) {
  se_robust_total <- sqrt(diag(vcov_total_hc3))[which(names(coef(fit_total)) == "sorted_2016")]
}
if (is.na(se_robust_short)) {
  se_robust_short <- sqrt(diag(vcov_short_hc3))[which(names(coef(fit_short)) == "sorted_2016")]
}

# E-value 계산 (HC3 표준오차 사용)
ev_total <- evalues.OLS(
  est = coef(fit_total)["sorted_2016"],
  se  = se_robust_total,
  sd  = sd_y
)
ev_short <- evalues.OLS(
  est = coef(fit_short)["sorted_2016"],
  se  = se_robust_short,
  sd  = sd_y_short
)

print(ev_total)
print(ev_short)




##Sequential G-Fromula
# 공변량 벡터 (네가 쓰는 X들)
xvars <- c("pol_interest_index","participation","social_political_days",
           "scaled_channel_count","age_group_4","race_collapsed","region")

# Mediator: AP_2016 ~ A + X (가중 OLS)
fit_M <- lm(affective_polarization_2016 ~ sorted_2016 + 
              pol_interest_index + participation + social_political_days + 
              scaled_channel_count + age_group_4 + race_collapsed + region,
            data = df_trim, weights = w_stab)

# Outcome: AP_2020 ~ A + M + X (가중 OLS)
fit_Y <- lm(affective_polarization_2020 ~ sorted_2016 + affective_polarization_2016 +
              pol_interest_index + participation + social_political_days + 
              scaled_channel_count + age_group_4 + race_collapsed + region,
            data = df_trim, weights = w_stab)


set.seed(2025)
K <- 1000  # 시뮬레이션 반복 수 (정밀도 필요하면 1000~2000)
n <- nrow(df_trim)

# mediator 분포 파라미터
mu_M_A1 <- predict(fit_M, newdata = transform(df_trim, sorted_2016 = 1), type = "response")
mu_M_A0 <- predict(fit_M, newdata = transform(df_trim, sorted_2016 = 0), type = "response")
sigma_M <- sqrt(sum(residuals(fit_M)^2 / (n - length(coef(fit_M)))))  # OLS residual SD

# 헬퍼: outcome 예측 함수
predY <- function(A, Mvec, data) {
  newd <- data
  newd$sorted_2016 <- A
  newd$affective_polarization_2016 <- Mvec
  as.numeric(predict(fit_Y, newdata = newd, type = "response"))
}

# 저장 벡터
EY1      <- numeric(K)  # E[Y_{A=1}]
EY0      <- numeric(K)  # E[Y_{A=0}]
EY1_M1   <- numeric(K)  # E[Y_{1, M(1)}]
EY1_M0   <- numeric(K)  # E[Y_{1, M(0)}]
EY0_M0   <- numeric(K)  # E[Y_{0, M(0)}]

for (k in 1:K) {
  # M(1), M(0) 시뮬레이션 (정규잔차 가정)
  M1 <- rnorm(n, mean = mu_M_A1, sd = sigma_M)
  M0 <- rnorm(n, mean = mu_M_A0, sd = sigma_M)
  
  # E[Y_a] : M의 실제 분포는 A에 의존 → A=a일 때 M(a)로 적분
  EY1[k]    <- mean(predY(1, M1, df_trim))
  EY0[k]    <- mean(predY(0, M0, df_trim))
  
  # NDE/NIE 구성요소
  EY1_M1[k] <- mean(predY(1, M1, df_trim))  # Y_{1, M(1)}
  EY1_M0[k] <- mean(predY(1, M0, df_trim))  # Y_{1, M(0)}
  EY0_M0[k] <- mean(predY(0, M0, df_trim))  # Y_{0, M(0)}
}

# 평균화
TE  <- mean(EY1 - EY0)                 # Total Effect
NDE <- mean(EY1_M0 - EY0_M0)           # Natural Direct Effect
NIE <- mean(EY1_M1 - EY1_M0)           # Natural Indirect Effect
c(TE = TE, NDE = NDE, NIE = NIE, PropMediated = NIE/TE)


library(boot)
library(parallel)

seq_g_stat_param_resid <- function(data, indices){
  d <- data[indices, ]
  
  # mediator 모델 (가중 회귀)
  fit_M_b <- lm(affective_polarization_2016 ~ sorted_2016 +
                  pol_interest_index + participation + social_political_days +
                  scaled_channel_count + age_group_4 + race_collapsed + region,
                data = d, weights = d$w_stab)
  
  # outcome 모델 (가중 회귀)
  fit_Y_b <- lm(affective_polarization_2020 ~ sorted_2016 + affective_polarization_2016 +
                  pol_interest_index + participation + social_political_days +
                  scaled_channel_count + age_group_4 + race_collapsed + region,
                data = d, weights = d$w_stab)
  
  n_b <- nrow(d)
  
  # 예측된 M 평균들 (각 관측치별)
  mu_M_A1 <- predict(fit_M_b, newdata = transform(d, sorted_2016 = 1))
  mu_M_A0 <- predict(fit_M_b, newdata = transform(d, sorted_2016 = 0))
  
  # 가중 잔차 표준편차 추정 (WLS residual variance)
  resid_b <- residuals(fit_M_b)
  w_b <- d$w_stab
  # 안전장치: 결측치 처리
  if(any(is.na(resid_b))) resid_b[is.na(resid_b)] <- 0
  if(any(is.na(w_b))) w_b[is.na(w_b)] <- 0
  
  denom <- sum(w_b, na.rm = TRUE) - length(coef(fit_M_b))
  if(denom <= 0) {
    sigma_b <- 0
  } else {
    sigma_b <- sqrt(sum(w_b * resid_b^2, na.rm = TRUE) / denom)
  }
  if(is.na(sigma_b) || sigma_b < 0) sigma_b <- 0
  
  # parametric residual resampling: 한 번의 시뮬레이션으로 근사(속도/정확성 균형)
  M1_sim <- rnorm(n_b, mean = mu_M_A1, sd = sigma_b)
  M0_sim <- rnorm(n_b, mean = mu_M_A0, sd = sigma_b)
  
  predY <- function(A, Mvec){
    nd <- d
    nd$sorted_2016 <- A
    nd$affective_polarization_2016 <- Mvec
    as.numeric(predict(fit_Y_b, newdata = nd))
  }
  
  EY1    <- mean(predY(1, M1_sim))
  EY0    <- mean(predY(0, M0_sim))
  EY1M1  <- mean(predY(1, M1_sim))
  EY1M0  <- mean(predY(1, M0_sim))
  EY0M0  <- mean(predY(0, M0_sim))
  
  TE  <- EY1 - EY0
  NDE <- EY1M0 - EY0M0
  NIE <- EY1M1 - EY1M0
  c(TE=TE, NDE=NDE, NIE=NIE, PropMediated=ifelse(abs(TE) < .Machine$double.eps, NA, NIE/TE))
}

# -----------------------------
# Bootstrap 실행 (병렬 + RNG 재현성 처리)
# -----------------------------
set.seed(2027)
par_type <- ifelse(.Platform$OS.type == "windows", "snow", "multicore")
Rboots <- 2000
ncpus  <- 4  # 필요에 따라 변경

if(par_type == "snow"){
  cl <- makeCluster(ncpus)
  # RNG 스트림 설정해서 재현성 확보
  clusterSetRNGStream(cl, 2027)
  # 필요 패키지/변수 노출 (boot will export statistic env generally, 추가로 필요하면 clusterExport)
  boot_seq_param <- boot(data = df_trim, statistic = seq_g_stat_param_resid,
                         R = Rboots, parallel = "snow", ncpus = ncpus, cl = cl)
  stopCluster(cl)
} else {
  # multicore: set.seed before boot is typically enough on Unix-like systems
  boot_seq_param <- boot(data = df_trim, statistic = seq_g_stat_param_resid,
                         R = Rboots, parallel = "multicore", ncpus = ncpus)
}

# 결과 확인
print(boot_seq_param$t0)         # 원점 추정
apply(boot_seq_param$t, 2, sd)   # bootstrap 표준편차

# BCa 신뢰구간 (예: TE, NDE, NIE)
boot.ci(boot_seq_param, type = "bca", index = 1)  # TE
boot.ci(boot_seq_param, type = "bca", index = 2)  # NDE
boot.ci(boot_seq_param, type = "bca", index = 3)  # NIE

library(purrr)


param_names <- c("TE","NDE","NIE")

# BCa 신뢰구간 안전 추출 헬퍼 (기존과 동일)
get_bca <- function(boot_obj, idx, conf=.95){
  ci <- try(boot.ci(boot_obj, type="bca", index=idx, conf=conf), silent=TRUE)
  if (inherits(ci, "try-error") || is.null(ci$bca)) return(c(lower=NA_real_, upper=NA_real_))
  bc <- ci$bca
  k <- ncol(bc)                          # 다양한 R 버전 호환: 마지막 두 컬럼이 하/상한
  c(lower = as.numeric(bc[1, k-1]), upper = as.numeric(bc[1, k]))
}

# 점추정 + BCa CI 테이블 (boot_seq_param 사용)
ci_tbl <- map_dfr(seq_along(param_names), function(i){
  b <- get_bca(boot_seq_param, idx = i, conf = .95)
  tibble(param = param_names[i], lower = b["lower"], upper = b["upper"])
})

sum_tbl <- tibble(param = param_names,
                  t0 = as.numeric(boot_seq_param$t0[1:3])) %>%
  left_join(ci_tbl, by = "param")

print(sum_tbl)
# => param, t0(점추정), lower/upper(BCa 95% CI) 확인

# ---- (그림 1) Forest plot: boot.ci(BCa) 결과 (parametric residual resampling) ----
p_forest <- ggplot(sum_tbl, aes(y = param, x = t0)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  labs(x = "Effect size (BCa 95% CI)", y = NULL,
       title = "BCa intervals from boot.ci",
       subtitle = "TE, NDE, NIE (parametric residual resampling sequential g-formula)") +
  theme_minimal(base_size = 13)
print(p_forest)

# ---- (그림 2) 부트스트랩 분포 + 점추정/CI 선 오버레이 (parametric residual resampling) ----
# boot_seq_param$t: columns are (TE, NDE, NIE, PropMediated) — 선택해서 사용
draws_df <- as_tibble(boot_seq_param$t[, 1:3, drop = FALSE]) %>%
  setNames(param_names) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value")

p_dist <- ggplot(draws_df, aes(x = value)) +
  geom_histogram(aes(y = ..density..), bins = 40, alpha = 0.6, fill = "lightblue", color = "black") +
  geom_density(linewidth = 1) +
  facet_wrap(~ param, scales = "free", ncol = 1) +
  geom_vline(data = sum_tbl, aes(xintercept = t0), linetype = 1, color = "red") +
  geom_vline(data = sum_tbl, aes(xintercept = lower), linetype = 2, color = "darkgreen") +
  geom_vline(data = sum_tbl, aes(xintercept = upper), linetype = 2, color = "darkgreen") +
  geom_vline(xintercept = 0, linetype = 3) +
  labs(x = "Bootstrap draws", y = "Density",
       title = "Bootstrap distributions with BCa 95% intervals",
       subtitle = "Parametric residual resampling sequential g-formula") +
  theme_minimal(base_size = 13)
print(p_dist)


#CDE
fit_Y <- lm(affective_polarization_2020 ~ sorted_2016 + affective_polarization_2016 +
              pol_interest_index + participation + social_political_days +
              scaled_channel_count + age_group_4 + race_collapsed + region,
            data = df_trim, weights = w_stab)

# 1️⃣ mediator의 고정값 선택
m_fixed <- mean(df_trim$affective_polarization_2016, na.rm = TRUE)
# 또는 분위수 등: quantile(df_trim$affective_polarization_2016, 0.5)

# 2️⃣ Y 예측 (M 고정)
pred_Y_CDE <- function(A_value, m_value){
  newd <- df_trim
  newd$sorted_2016 <- A_value
  newd$affective_polarization_2016 <- m_value
  predict(fit_Y, newdata = newd, type = "response")
}

# 3️⃣ CDE 계산
E_Y1_m <- mean(pred_Y_CDE(1, m_fixed))
E_Y0_m <- mean(pred_Y_CDE(0, m_fixed))
CDE_est <- E_Y1_m - E_Y0_m
CDE_est

library(boot)

cde_stat <- function(data, indices){
  d <- data[indices, ]
  fit_Y_b <- lm(affective_polarization_2020 ~ sorted_2016 + affective_polarization_2016 +
                  pol_interest_index + participation + social_political_days +
                  scaled_channel_count + age_group_4 + race_collapsed + region,
                data = d, weights = d$w_stab)
  m_fixed <- mean(d$affective_polarization_2016, na.rm = TRUE)
  
  newd1 <- d; newd1$sorted_2016 <- 1; newd1$affective_polarization_2016 <- m_fixed
  newd0 <- d; newd0$sorted_2016 <- 0; newd0$affective_polarization_2016 <- m_fixed
  
  Y1 <- mean(predict(fit_Y_b, newdata = newd1))
  Y0 <- mean(predict(fit_Y_b, newdata = newd0))
  return(Y1 - Y0)
}

set.seed(2030)
par_type <- ifelse(.Platform$OS.type=="windows","snow","multicore")
boot_cde <- boot(data = df_trim, statistic = cde_stat, R = 1000, parallel = par_type, ncpus = 4)

boot_cde$t0  # 점추정
boot.ci(boot_cde, type = "bca")  # BCa 신뢰구간



# ---- Forest plot from BCa CIs for Short (2016) & Total (2020) ----
library(ggplot2)

# 1) 요약 테이블 구성 (원점추정: boot_*$t0, 구간: ci_*_bca)
sum_tbl <- data.frame(
  param = factor(c("Short (2016 outcome)", "Total (2020 outcome)"),
                 levels = c("Total (2020 outcome)", "Short (2016 outcome)")), # 위에 Total 오도록
  t0    = c(boot_short$t0, boot_total$t0),
  lower = c(ci_short_bca[1], ci_total_bca[1]),
  upper = c(ci_short_bca[2], ci_total_bca[2])
)

# 2) x축 범위 여백 설정
xmin <- min(sum_tbl$lower, 0)
xmax <- max(sum_tbl$upper)
pad  <- 0.2 * (xmax - xmin)
xmin <- xmin - pad
xmax <- xmax + pad

# 3) Forest plot
p_forest <- ggplot(sum_tbl, aes(y = param, x = t0)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  coord_cartesian(xlim = c(xmin, xmax)) +
  labs(
    x = "Effect size (BCa 95% CI)",
    y = NULL,
    title = "BCa intervals from boot.ci",
    subtitle = "Short (2016) and Total (2020) ATE via DR + stabilized weights"
  ) +
  theme_minimal(base_size = 13)

print(p_forest)

# 4) (옵션) 수치 라벨을 오른쪽에 붙이고 싶다면:
lab_df <- transform(sum_tbl,
                    label = sprintf("%.2f  [%.2f, %.2f]", t0, lower, upper))
p_forest +
  geom_text(data = lab_df, aes(x = upper, y = param, label = label),
            hjust = -0.1, size = 3.7) +
  coord_cartesian(xlim = c(xmin, xmax + pad))

# 5) (옵션) 파일 저장
# ggsave("forest_bca_short_total.png", p_forest, width = 7, height = 3, dpi = 300)



# ---- Forest plot: Bootstrap mean & 95% CI (Short / Total) ----
library(ggplot2)

# 1) 요약 테이블: 사용자가 제공한 점추정과 CI
sum_tbl <- data.frame(
  param = factor(c("Short (2016 outcome)", "Total (2020 outcome)"),
                 levels = c("Total (2020 outcome)", "Short (2016 outcome)")), # Total을 위에
  mean  = c(14.18, 8.39),
  lower = c(11.85, 5.71),
  upper = c(16.77, 11.14)
)

# 2) 라벨 문자열 (표시 자리수는 필요에 맞게 조절)
lab_df <- transform(
  sum_tbl,
  label = sprintf("%.2f  [%.2f, %.2f]", mean, lower, upper)
)

# 3) x축 패딩(라벨이 잘리지 않도록 오른쪽 여백 크게)
xmin <- min(sum_tbl$lower, 0)
xmax <- max(sum_tbl$upper)
pad  <- 0.5 * (xmax - xmin)

# 4) 그리기
p_forest <- ggplot(sum_tbl, aes(y = param, x = mean)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_text(data = lab_df, aes(x = upper, y = param, label = label),
            hjust = -0.1, size = 3.7) +
  coord_cartesian(xlim = c(xmin, xmax + pad)) +
  labs(
    x = "Effect size (Bootstrap 95% CI)",
    y = NULL,
    title = "Bootstrap intervals",
    subtitle = "Short (2016) and Total (2020) mean estimates with 95% CI"
  ) +
  theme_minimal(base_size = 13)

print(p_forest)


