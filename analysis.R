N <- 120   # population size: 15-min intervals, 24 per day x 5 weekdays

# expansion estimator of a total under SRS with fpc and a t-based 95% CI
srs_total <- function(y, N) {
  n   <- length(y)
  tau <- N * mean(y)                             
  se  <- N * sqrt((1 - n / N) * var(y) / n) 
  t   <- qt(0.975, df = n - 1)            # t used since n < 50
  c(estimate = tau, se = se, lower = tau - t * se, upper = tau + t * se)
}

# SRS sample size giving Var(ybar) = V0:  n = N / (1 + N V0 / S^2)
srs_n <- function(N, S2, V0) N / (1 + N * V0 / S2)

# pilot survey
pilot <- read.csv("data/pilot_counts.csv")   # people = walkers + cyclists

# combine 5-min bins into 10 and 15-min counts within each session
counts5  <- pilot$people
counts10 <- as.vector(tapply(counts5, list(ceiling(pilot$bin / 2), pilot$date), sum))
counts15 <- as.vector(tapply(counts5, list(ceiling(pilot$bin / 3), pilot$date), sum))

# table 1 (left): longer intervals smooth the bursts, so the CV falls
interval_cv <- t(sapply(list(`5 min` = counts5, `10 min` = counts10, `15 min` = counts15),
                        function(y) c(n = length(y), mean = mean(y), sd = sd(y),
                                      cv = sd(y) / mean(y))))
round(interval_cv, 2)

# pilot estimate of the weekly total, using 15-min units to match the main survey
pilot_est <- srs_total(counts15, N)
round(pilot_est)                                  # 3240, SE 944

# sample size to halve the pilot SE (s^2 from the pilot in place of S^2)
V0 <- unname(pilot_est["se"] / 2 / N)^2 
n_req <- srs_n(N, var(counts15), V0)
n_req                                             # 20.9 -> n = 21

# cross-check: 30-min sessions as units (N = 60), then double
counts30 <- as.vector(tapply(counts5, pilot$date, sum))
V0_30 <- unname(srs_total(counts30, 60)["se"] / 2 / 60)^2
2 * srs_n(60, var(counts30), V0_30)               # 20.9 


# main survey
main <- read.csv("data/main_survey_counts.csv")   # people = walkers + cyclists
n <- nrow(main)                                   # 21

# tables 1 (right) and 2: weekly totals with SEs and 95% CIs
est <- rbind(people = srs_total(main$people, N),
             cars   = srs_total(main$cars, N))
round(est)                                        # people 2314 (264), cars 354 (40)
round(est / 5)                                    # per weekday
round(100 * est[, "se"] / est[, "estimate"])      # relative SE: 11%

# precision achieved vs expected from the pilot
N * sqrt((1 - n / N) * var(counts15) / n)         # expected SE, about 470
unname(est["people", "se"] / pilot_est["se"])     # achieved / pilot SE = 0.28

# difference between each day
days <- c("Mon", "Tue", "Wed", "Thu", "Fri")
day  <- factor(main$day, levels = days)
n_j  <- as.vector(table(day))                   
ybar_j <- tapply(main$people, day, mean)
s_j    <- tapply(main$people, day, sd)
round(cbind(n = n_j, mean = ybar_j, sd = s_j), 1)


# planning a better survey
# SRS size for +-10% of the total with 95% confidence: n = N / (1 + N theta^2 mu^2 / (4 S^2))
theta <- 0.1
N / (1 + N * theta^2 * mean(main$people)^2 / (4 * var(main$people)))   # about 63

# stratify by day: N_j = 24, W_j = 0.2, S_j estimated by the daily SDs
N_j <- rep(24, 5)
W_j <- N_j / N

# Neyman allocation (Theorem 5.5): n_j proportional to N_j S_j
round(n * N_j * s_j / sum(N_j * s_j), 1)          # 3.7 3.3 2.1 8.7 3.3 -> 4 3 2 9 3
var_neyman <- sum(W_j * s_j)^2 / n - sum(W_j * s_j^2) / N
N * sqrt(var_neyman)                              # SE of total about 165 (vs 264)

# proportional allocation for comparison: n_j = n W_j
n_prop <- n * W_j
N * sqrt(sum(W_j^2 * (1 - n_prop / N_j) * s_j^2 / n_prop))   # about 195


# figures

col_day <- c("blue", "orange", "green", "red", "purple")   # one colour per day
font <- "serif"

# Figure 1: pilot counts in each 5-min bin, grouped by session
fig1 <- function() {
  par(family = font, mar = c(3, 4.2, 1, 1))
  barplot(matrix(counts5, nrow = 6), beside = TRUE, col = "grey", border = "white",
          names.arg = c("Mon 12 Oct (pm)", "Wed 14 Oct (am)", "Fri 16 Oct (am)"),
          ylim = c(0, 40), ylab = "People per 5 minutes")
  abline(h = mean(counts5), lty = 2)                # overall mean
}

# Figure 2: (a) counts against time of day, (b) mean count for each day
fig2 <- function() {
  par(family = font, mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  hour <- as.numeric(substr(main$start_time, 1, 2)) + as.numeric(substr(main$start_time, 4, 5)) / 60
  plot(hour, main$people, pch = 19, col = col_day[day], ylim = c(0, 60),
       xlab = "Time of day (hours)", ylab = "People per interval", main = "(a) Counts by time of day")
  legend("top", legend = days, col = col_day, pch = 19, bty = "n", ncol = 5, cex = 0.7)
  barplot(ybar_j, col = col_day, border = NA, ylim = c(0, 60),
          ylab = "Mean people per interval", main = "(b) Mean count by day")
  abline(h = mean(main$people), lty = 2)            # overall mean
}

dir.create("figures", showWarnings = FALSE)
png("figures/fig1_pilot_bursts.png", width = 2100, height = 800, res = 300); fig1(); dev.off()
png("figures/fig2_main_survey.png",  width = 2100, height = 900, res = 300); fig2(); dev.off()
