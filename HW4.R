install.packages("ggplot2")
install.package("broom")
install.packages("gridExtra")
install.packages("multcomp")

library(ggplot2)
library(multcomp)
library(broom)
library(gridExtra)

#Problem 1 data
techniques <- factor(rep(1:4, each = 4))
y <- c(
  # Technique 1
  3129, 3000, 2865, 2890,
  # Technique 2
  3200, 3300, 2975, 3150,
  # Technique 3
  2800, 2900, 2985, 3050,
  # Technique 4
  2600, 2700, 2600, 2765
)
dat <- data.frame(tech = techniques, y = y)

group_summary <- aggregate(y ~ tech, dat, function(x) c(n = length(x), mean = mean(x), sd = sd(x)))
group_summary <- do.call(data.frame, group_summary)
names(group_summary) <- c("tech", "n", "mean", "sd")
print(group_summary)

grand_mean <- mean(dat$y)
grand_mean

#SSB = 489737.6875 in the solutions pdf
#SSE = 53911.625 in the solutions pdf

fit <- lm(y ~ tech, data = dat)
anova_tab <- anova(fit)
anova_tab$`Mean Sq` <- anova_tab$`Sum Sq` / anova_tab$Df
print(anova_tab)

# Manual calculation of SST, SSB, SSE
N <- nrow(dat)
k <- length(levels(dat$tech))
ni <- as.numeric(table(dat$tech))

# SST (total)
SST <- sum( (dat$y - grand_mean)^2 )

# SSB (between / treatment)
group_means <- tapply(dat$y, dat$tech, mean)
SSB <- sum( ni * (group_means - grand_mean)^2 )

# SSE (within)
SSE <- sum( (dat$y - rep(group_means, each = 4))^2 )

SST
SSB
SSE


MSB <- SSB / (k - 1)
MSE <- SSE / (N - k)
Fval <- MSB / MSE
pval <- pf(Fval, df1 = k - 1, df2 = N - k, lower.tail = FALSE)

cat("MSB =", MSB, "\n")
cat("MSE =", MSE, "\n")
cat("F   =", Fval, " with p-value =", pval, "\n\n")

anova_report <- data.frame(
  Source = c("Between", "Within", "Total"),
  SS = c(SSB, SSE, SST),
  df = c(k - 1, N - k, N - 1),
  MS = c(MSB, MSE, NA)
)
print(anova_report)

resid_vals <- residuals(fit)
fitted_vals <- fitted(fit)

qq_dat <- data.frame(resid = resid_vals)
qq_plot <- ggplot(qq_dat, aes(sample = resid)) +
  stat_qq() + stat_qq_line() +
  ggtitle("QQ-plot of residuals (Problem 1)") +
  theme_minimal()
qq_plot

rvf_plot <- ggplot(data.frame(fitted = fitted_vals, resid = resid_vals), aes(x = fitted, y = resid)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("Residuals vs Fitted (Problem 1)") +
  xlab("Fitted values") + ylab("Residuals") +
  theme_minimal()
rvf_plot

# compute pairwise differences manually and run t-tests using pooled MSE
pairwise <- combn(levels(dat$tech), 2, simplify = FALSE)
pair_results <- lapply(pairwise, function(p) {
  a <- as.numeric(p[1]); b <- as.numeric(p[2])
  idx_a <- dat$tech == p[1]
  idx_b <- dat$tech == p[2]
  mean_diff <- mean(dat$y[idx_a]) - mean(dat$y[idx_b])
  # standard error using common MSE
  se_diff <- sqrt(MSE * (1/ni[1] + 1/ni[2]))
  t_stat <- mean_diff / se_diff
  p_two <- 2 * pt(-abs(t_stat), df = N - k)
  list(pair = paste(p[1], "vs", p[2]), diff = mean_diff, t = t_stat, p = p_two)
})
pair_results

pair_df <- do.call(rbind, lapply(pair_results, as.data.frame))
rownames(pair_df) <- NULL
print("\nPairwise (unadjusted) t-test results:")
print(pair_df)

#Bonferroni adjustment
pair_df$padj_bonf <- pmin(1, pair_df$p * nrow(pair_df))
pair_df$signif_bonf <- pair_df$padj_bonf < 0.05
cat("\nPairwise results with Bonferroni adjustment:\n")
print(pair_df[, c("pair", "diff", "t", "p", "padj_bonf", "signif_bonf")])


#Problem 2
k2 <- 6; n2 <- 5; N2 <- 30
SSTotal <- 900.25; SSTreat <- 750.5
SSE2 <- SSTotal - SSTreat
MSE2 <- SSE2 / (N2 - k2)
MS_treat2 <- SSTreat / (k2 - 1)
F2 <- MS_treat2 / MSE2
R2 <- SSTreat / SSTotal
cat("SSE =", SSE2, "\n")
cat("MSE (estimate sigma^2) =", MSE2, "\n")
cat("Proportion explained (R^2) = ", R2)
cat("F statistic =", F2, " with df =", (k2 - 1), "and", (N2 - k2), "\n")

#Problem 3
estimate <- function(coefs, group_means, group_ns, MSE, df_resid, conf.level = 0.95) {
  Lhat <- sum(coefs * group_means)
  est_var <- MSE * sum( (coefs^2) / group_ns )
  se <- sqrt(est_var)
  tval <- Lhat / se
  alpha <- 1 - conf.level
  tcrit <- qt(1 - alpha/2, df_resid)
  CI <- c(Lhat - tcrit * se, Lhat + tcrit * se)
  pval <- 2 * pt(-abs(tval), df_resid)
  list(Lhat = Lhat, se = se, t = tval, p = pval, CI = CI)
}

group_means_p1 <- as.numeric(group_means)
group_ns_p1 <- rep(4, length(group_means_p1))
c_example <- c(0.5, 0.5, -0.5, -0.5)
res <- estimate(c_example, group_means_p1, group_ns_p1, MSE, N - k)
print(res)

print_<- function(cvec, name = "contrast") {
  res <- estimate(cvec, group_means_p1, group_ns_p1, MSE, N - k)
  cat("\nContrast:", name, "\n")
  cat("Coefficients:", paste(cvec, collapse = ", "), "\n")
  cat("Estimate (Lhat) =", res$Lhat, "\n")
  cat("Std. Error =", res$se, "\n")
  cat("t =", res$t, " df =", N - k, " p =", res$p, "\n")
  cat("95% CI = (", res$CI[1], ",", res$CI[2], ")\n")
}
print_(c_example, "Avg(1,2) - Avg(3,4)")

anova_report
pair_df
data.frame(Tukey = I(list(tukey$tech)))

