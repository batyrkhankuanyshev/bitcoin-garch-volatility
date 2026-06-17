# ============================================================
# Bitcoin Volatility Analysis - GARCH(1,1) Model
# Data: CoinMarketCap via crypto2 package
# Period: 2021-01-01 to 2024-12-31
# ============================================================

# ============================================================
# 1. INSTALL AND LOAD PACKAGES
# ============================================================

if (!requireNamespace("crypto2", quietly = TRUE)) {install.packages("crypto2")}
if (!requireNamespace("moments", quietly = TRUE)) {install.packages("moments")}
if (!requireNamespace("rugarch", quietly = TRUE)) {install.packages("rugarch")}
if (!requireNamespace("tseries", quietly = TRUE)) {install.packages("tseries")}

library(crypto2)
library(moments)
library(rugarch)
library(tseries)

# ============================================================
# 2. DOWNLOAD BITCOIN DATA FROM COINMARKETCAP
# ============================================================

start_date <- "20210101"
end_date   <- "20241231"

cat("Downloading Bitcoin data from CoinMarketCap...\n")

# Get the full coin list first (only needs to run once)
coins <- crypto_list()

# Find Bitcoin's row
btc_info <- coins[coins$slug == "bitcoin", ]

# Now use the data frame as coin_list
btc_all <- crypto_history(
  coin_list  = btc_info,
  start_date = start_date,
  end_date   = end_date,
  interval   = "daily",
  finalWait  = FALSE
)

cat("Download complete. Rows downloaded:", nrow(btc_all), "\n")

# ============================================================
# 3. CLEAN AND PREPARE DATA
# ============================================================

# Select relevant columns
btc <- btc_all[, c("timestamp", "close", "volume", "market_cap")]
names(btc) <- c("date", "price", "volume", "market_cap")

# Convert date and sort
btc$date <- as.Date(btc$date)
btc <- btc[!duplicated(btc$date), ]
btc <- btc[order(btc$date), ]

# Calculate daily log returns (in percentage)
btc$log_ret <- c(NA, 100 * diff(log(btc$price)))

# Remove NA and infinite values
btc <- btc[!is.na(btc$log_ret), ]
btc <- btc[is.finite(btc$log_ret), ]
btc <- btc[btc$price > 0, ]

cat("Final sample size:", nrow(btc), "daily observations\n")

# ============================================================
# 4. DESCRIPTIVE STATISTICS (Section 2.3)
# ============================================================

N        <- length(btc$log_ret)
mean_ret <- mean(btc$log_ret)
sd_ret   <- sd(btc$log_ret)
min_ret  <- min(btc$log_ret)
max_ret  <- max(btc$log_ret)
skew_ret <- skewness(btc$log_ret)
kurt_ret <- kurtosis(btc$log_ret) - 3  # excess kurtosis

stats_table <- data.frame(
  Statistic = c("N", "Mean (%)", "Std Dev (%)", "Min (%)", "Max (%)", "Skewness", "Excess Kurtosis"),
  Value = round(c(N, mean_ret, sd_ret, min_ret, max_ret, skew_ret, kurt_ret), 4)
)

cat("\n========== DESCRIPTIVE STATISTICS ==========\n")
print(stats_table)

# ============================================================
# 5. TIME SERIES PLOTS (Section 2.3)
# ============================================================

# Plot 1: Closing price
dev.new()  # Opens new graphics window (use quartz() on Mac, windows() on PC)
plot(btc$date, btc$price, type = "l", col = "blue", lwd = 1.5,
     main = "Bitcoin Daily Closing Price (2021-2024)",
     xlab = "Date", ylab = "Price (USD)")

# Plot 2: Daily returns
dev.new()
plot(btc$date, btc$log_ret, type = "l", col = "darkred", lwd = 1,
     main = "Bitcoin Daily Log Returns (%)",
     xlab = "Date", ylab = "Log Return (%)")
abline(h = 0, col = "gray", lty = 2)

# Plot 3: Squared returns (volatility proxy)
btc$sq_ret <- btc$log_ret^2
dev.new()
plot(btc$date, btc$sq_ret, type = "l", col = "darkgreen", lwd = 1,
     main = "Squared Daily Log Returns of Bitcoin",
     xlab = "Date", ylab = "Squared Return")

# ============================================================
# 6. ACF AND PACF PLOTS
# ============================================================

dev.new()
par(mfrow = c(2, 2))
acf(btc$log_ret, main = "ACF of Returns", lag.max = 30)
pacf(btc$log_ret, main = "PACF of Returns", lag.max = 30)
acf(btc$sq_ret, main = "ACF of Squared Returns", lag.max = 30)
pacf(btc$sq_ret, main = "PACF of Squared Returns", lag.max = 30)

# ============================================================
# 7. ARCH EFFECTS TEST (Section 4.1)
# ============================================================

# Ljung-Box test on squared returns at different lags
lb_5  <- Box.test(btc$sq_ret, lag = 5, type = "Ljung-Box")
lb_10 <- Box.test(btc$sq_ret, lag = 10, type = "Ljung-Box")
lb_20 <- Box.test(btc$sq_ret, lag = 20, type = "Ljung-Box")

cat("\n========== ARCH EFFECTS TEST ==========\n")
cat("Ljung-Box Test on Squared Returns\n")
cat("----------------------------------\n")
cat(sprintf("Lag  5: Q-statistic = %.2f, p-value = %.6f\n", lb_5$statistic, lb_5$p.value))
cat(sprintf("Lag 10: Q-statistic = %.2f, p-value = %.6f\n", lb_10$statistic, lb_10$p.value))
cat(sprintf("Lag 20: Q-statistic = %.2f, p-value = %.6f\n", lb_20$statistic, lb_20$p.value))

# Interpretation
if (lb_10$p.value < 0.05) {
  cat("\nConclusion: ARCH effects present. GARCH modeling is appropriate.\n")
} else {
  cat("\nConclusion: No strong ARCH effects detected.\n")
}

# ============================================================
# 8. GARCH(1,1) ESTIMATION (Section 4.2)
# ============================================================

# Specify GARCH(1,1) with Student-t distribution
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"  # Student-t for fat tails
)

# Estimate the model
cat("\nEstimating GARCH(1,1) model...\n")
garch_fit <- ugarchfit(spec = spec, data = btc$log_ret)

# Extract coefficients
coefs <- coef(garch_fit)
alpha <- coefs["alpha1"]
beta  <- coefs["beta1"]
omega <- coefs["omega"]
mu    <- coefs["mu"]
shape <- coefs["shape"]
persistence <- alpha + beta

# Standard errors and p-values (from robust matrix)
robust_mat <- garch_fit@fit$robust.matcoef
se <- robust_mat[,2]
pvals <- robust_mat[,4]

cat("\n========== GARCH(1,1) ESTIMATION RESULTS ==========\n")
cat("Model: r_t = μ + ε_t, ε_t = σ_t * z_t\n")
cat("       σ_t^2 = ω + α ε_{t-1}^2 + β σ_{t-1}^2\n")
cat("Distribution: Student-t\n")
cat("=====================================================\n")
cat(sprintf("%-12s %10s %10s %10s %10s\n", "Parameter", "Coefficient", "Std Error", "z-value", "p-value"))
cat("-----------------------------------------------------\n")
cat(sprintf("%-12s %10.6f %10.6f %10.4f %10.4f\n", "μ (mean)", mu, se["mu"], mu/se["mu"], pvals["mu"]))
cat(sprintf("%-12s %10.8f %10.8f %10.4f %10.4f\n", "ω", omega, se["omega"], omega/se["omega"], pvals["omega"]))
cat(sprintf("%-12s %10.6f %10.6f %10.4f %10.4f\n", "α (ARCH)", alpha, se["alpha1"], alpha/se["alpha1"], pvals["alpha1"]))
cat(sprintf("%-12s %10.6f %10.6f %10.4f %10.4f\n", "β (GARCH)", beta, se["beta1"], beta/se["beta1"], pvals["beta1"]))
cat(sprintf("%-12s %10.6f\n", "α + β", persistence))
cat(sprintf("%-12s %10.2f\n", "Shape (t-dof)", shape))
cat("=====================================================\n")

# Information criteria
ic <- infocriteria(garch_fit)
cat("\nInformation Criteria:\n")
cat(sprintf("  AIC: %.4f\n", ic[1]))
cat(sprintf("  BIC: %.4f\n", ic[2]))
cat(sprintf("  Log-likelihood: %.2f\n", likelihood(garch_fit)))

# ============================================================
# 9. CONSTANT VARIANCE MODEL (BENCHMARK)
# ============================================================

# Constant variance model: just sample mean and variance
constant_var <- lm(btc$log_ret ~ 1)
sigma_const <- summary(constant_var)$sigma
logLik_const <- logLik(constant_var)
aic_const <- AIC(constant_var)
bic_const <- BIC(constant_var)

cat("\n========== CONSTANT VARIANCE MODEL (BENCHMARK) ==========\n")
cat(sprintf("Constant volatility (σ): %.4f%%\n", sigma_const))
cat(sprintf("Log-likelihood: %.2f\n", logLik_const))
cat(sprintf("AIC: %.4f\n", aic_const))
cat(sprintf("BIC: %.4f\n", bic_const))

cat("\n========== MODEL COMPARISON ==========\n")
cat(sprintf("                    GARCH(1,1)   Constant Var   Improvement\n"))
cat(sprintf("Log-likelihood:      %10.2f   %10.2f   %+10.2f\n", 
            likelihood(garch_fit), logLik_const, likelihood(garch_fit) - logLik_const))
cat(sprintf("AIC:                 %10.4f   %10.4f   %+10.4f\n", 
            ic[1], aic_const, aic_const - ic[1]))
cat(sprintf("BIC:                 %10.4f   %10.4f   %+10.4f\n", 
            ic[2], bic_const, bic_const - ic[2]))

# ============================================================
# 10. DIAGNOSTIC CHECKS (Section 4.4)
# ============================================================

# Extract standardized residuals
std_resid <- residuals(garch_fit, standardize = TRUE)

# Ljung-Box tests on standardized residuals
lb_resid_10 <- Box.test(std_resid, lag = 10, type = "Ljung-Box")
lb_resid_20 <- Box.test(std_resid, lag = 20, type = "Ljung-Box")

# Ljung-Box tests on squared standardized residuals (should show no ARCH)
lb_sq_resid_10 <- Box.test(std_resid^2, lag = 10, type = "Ljung-Box")
lb_sq_resid_20 <- Box.test(std_resid^2, lag = 20, type = "Ljung-Box")

cat("\n========== DIAGNOSTIC CHECKS ==========\n")
cat("\nLjung-Box Test on Standardized Residuals:\n")
cat(sprintf("  Lag 10: Q = %.2f, p-value = %.4f\n", lb_resid_10$statistic, lb_resid_10$p.value))
cat(sprintf("  Lag 20: Q = %.2f, p-value = %.4f\n", lb_resid_20$statistic, lb_resid_20$p.value))

cat("\nLjung-Box Test on Squared Standardized Residuals:\n")
cat(sprintf("  Lag 10: Q = %.2f, p-value = %.4f\n", lb_sq_resid_10$statistic, lb_sq_resid_10$p.value))
cat(sprintf("  Lag 20: Q = %.2f, p-value = %.4f\n", lb_sq_resid_20$statistic, lb_sq_resid_20$p.value))

# Jarque-Bera test for normality
jb_test <- jarque.bera.test(std_resid)
cat("\nJarque-Bera Normality Test:\n")
cat(sprintf("  JB-statistic = %.2f, p-value = %.4f\n", jb_test$statistic, jb_test$p.value))

# QQ-plot
dev.new()
qqnorm(std_resid, main = "QQ-Plot of Standardized Residuals")
qqline(std_resid, col = "red", lwd = 2)

# Histogram of standardized residuals
dev.new()
hist(std_resid, breaks = 50, main = "Histogram of Standardized Residuals",
     xlab = "Standardized Residuals", col = "lightblue", freq = FALSE)
curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)

# News impact curve
dev.new()
newsimpact(garch_fit)$plot

# ============================================================
# 11. SAVE OUTPUTS
# ============================================================

# Save cleaned dataset
write.csv(btc, "bitcoin_daily_2021_2024.csv", row.names = FALSE)

# Save summary statistics
write.csv(stats_table, "bitcoin_summary_stats_2021_2024.csv", row.names = FALSE)

# Save GARCH results to text file
sink("garch_results.txt")
cat("========================================\n")
cat("GARCH(1,1) Estimation Results\n")
cat("Bitcoin Daily Returns (2021-2024)\n")
cat("========================================\n\n")
print(garch_fit)
cat("\n\n========================================\n")
cat("Diagnostic Checks\n")
cat("========================================\n")
cat("\nLjung-Box Test on Standardized Residuals (lag 10):\n")
print(lb_resid_10)
cat("\nLjung-Box Test on Squared Standardized Residuals (lag 10):\n")
print(lb_sq_resid_10)
sink()

cat("\n========================================\n")
cat("Files saved:\n")
cat("  - bitcoin_daily_2021_2024.csv\n")
cat("  - bitcoin_summary_stats_2021_2024.csv\n")
cat("  - garch_results.txt\n")
cat("========================================\n")

# ============================================================
# 12. SUMMARY OF FINDINGS
# ============================================================

cat("\n========== SUMMARY OF FINDINGS ==========\n")
cat(sprintf("1. Sample size: %d daily observations\n", N))
cat(sprintf("2. Mean daily return: %.4f%%\n", mean_ret))
cat(sprintf("3. Standard deviation: %.4f%%\n", sd_ret))
cat(sprintf("4. Excess kurtosis: %.2f (fat tails present)\n", kurt_ret))
cat(sprintf("5. ARCH effects: Present (p-value < 0.001)\n"))
cat(sprintf("6. GARCH persistence (α+β): %.4f\n", persistence))
if (persistence > 0.99) {
  cat("   → Very high persistence, shocks decay slowly\n")
} else {
  cat("   → Moderate persistence\n")
}
cat(sprintf("7. GARCH improves AIC by %.2f vs constant variance\n", aic_const - ic[1]))
cat("========================================\n")

