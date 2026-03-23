
          ###### Gull attack joint model (GAPC, GAPM and GAF) ######


# Article title: Calf mortality in a whale population varies with seabird attacks
# Authors: María Piotto, Iván Barberá, Mariano Sironi, Victoria J. Rowntree,
# Marcela M. Uhart, Macarena Agrelo, Alejandro A. Fernández Ajó, Jon Seger,
#Carina F. Marón.

# Scrip author: Iván Barberá
# Questions can be addressed to ivanbarbera93@gmail.com

# Library -----------------------------------------------------------------

library(tidyverse)
library(plyr)
library(rstan)
library(logitnorm)
library(bayestestR) # hdi

# Functions ---------------------------------------------------------------

qmean <- function(x) {
  quants <- quantile(x, probs = c(0.025, 0.975), method = 8)
  res <- c(mean(x), quants)
  names(res) <- c("mean", "lower", "upper")
  return(res)
}

hdi_lower <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_low"]])
}

hdi_upper <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_high"]])
}

hdint <- function(x, ci = 0.95, name = "y") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, ci$CI_high)
  names(result) <- paste(rep(name, 2), c("lower", "upper"), sep = "_")
  result
}

# function to compute probability of a quotient being >1
qprob <- function(x) {
  length(x[x > 1]) / length(x)
}


# Data -----------------------------------------------------------------------

# Gull attack data belong to the Southern Right Whale Research Program (Instituto
# de Conservación de Ballenas, Argentina). Requests to use the data should be
# addressed to Mariano Sironi mariano.sironi@icb.org.ar

fd <- read.csv("GAF data.csv", sep = ",")

# About this database:
# year: year when observations were made
# golfo = gulf: where the data was collected (GN: Golfo Nuevo, GJS: Golfo San José)
# intt:total number of intervals of the day
# intw: number of daily 5-min observation intervals with at least one attack on
# either mother or calf
# Obs: day of observations

pd <- read.csv("GAP data.csv", sep = ",")

# year: year when observations were made
# golfo = gulf: where the data was collected (GN: Golfo Nuevo, GJS: Golfo San José)
# mc: mother-calf factor. M: mother, C: calf
# pa: daily number of attacks on mothers or calves
# interv: total number of 5-min intervals of the day
# h: number of observation hours of the day
# HighLow: class of year (High or low mortality years)

pd <- pd[pd$interv >0, ] # remove rows with 0 hours of observation

# set equal levels for factors

fd$golfo <- revalue(fd$golfo, replace = c("GSJ" = "SJ", "GN" = "N"))

years <- unique(fd$year)
#years == unique(pd$year) # OK

fd$year <- factor(fd$year, levels = as.character(years))
pd$year <- factor(pd$year, levels = as.character(years))

fd$golfo <- factor(fd$golfo, levels = c("N", "SJ"))
pd$golfo <- factor(pd$golfo, levels = c("N", "SJ"))

pd$mc <- factor(pd$mc, levels = c("M", "C"))

# order

fd <- fd[with(fd, order(year, golfo)), ]
pd <- pd[with(pd, order(year, golfo, mc)), ]

# create interaction factors
# gulf-year
pd$gy <- paste(pd$golfo, pd$year, sep = "_")
fd$gy <- paste(fd$golfo, fd$year, sep = "_")

pd$gy <- factor(pd$gy, levels = unique(pd$gy))
fd$gy <- factor(fd$gy, levels = unique(fd$gy))

# gulf-mc (for gap)
pd$gm <- paste(pd$golfo, pd$mc, sep = "_")
pd$gm <- factor(pd$gm, levels = unique(pd$gm))

# gulf-mc-year (for gap)
pd$gmy <- paste(pd$golfo, pd$mc, pd$year, sep = "_")
pd$gmy <- factor(pd$gmy, levels = unique(pd$gmy))


# create design matrices
# start with X, followed by the factor (g = gulf, y = year, m = mother-calf) and
# p or f for pressure or frequency

# gulf matrix
X_g_p <- model.matrix(pa ~ golfo - 1, pd)
X_g_f <- model.matrix(intw ~ golfo - 1, fd)

# gulf - year matrix
X_gy_p <- model.matrix(pa ~ gy - 1, pd)
X_gy_f <- model.matrix(intw ~ gy - 1, fd)

# gulf-mc matrix (for gap)
X_gm_p <- model.matrix(pa ~ gm - 1, pd)
# gulf-mc-year matrix (for gap)
X_gmy_p <- model.matrix(pa ~ gmy - 1, pd)


# Prediction matrices at year level (start with Z)
# gap
pd$rate <- pd$pa / pd$h
pd_agg_gmy <- aggregate(rate ~ golfo + mc + year + gm + gy + gmy, pd, mean)
str(pd_agg_gmy)

Z_gm_p <- model.matrix(rate ~ gm - 1, pd_agg_gmy)
Z_gy_p <- model.matrix(rate ~ gy - 1, pd_agg_gmy)
Z_gmy_p <- model.matrix(rate ~ gmy - 1, pd_agg_gmy)

# gaf
fd$prop <- fd$intw / fd$intt
fd_agg_gy <- aggregate(prop ~ golfo + year + gy, fd, mean)
str(fd_agg_gy)

Z_g_f <- model.matrix(prop ~ golfo - 1, fd_agg_gy)
Z_gy_f <- model.matrix(prop ~ gy - 1, fd_agg_gy)

# Prediction matrices at gulf level (start with V)
# gap
pd_agg_gm <- aggregate(rate ~ golfo + mc + gm, pd_agg_gmy, mean)
V_gm_p <- model.matrix(rate ~ gm - 1, pd_agg_gm)
V_g_p <- model.matrix(rate ~ golfo - 1, pd_agg_gm)

# gaf
fd_agg_g <- aggregate(prop ~ golfo, fd_agg_gy, mean)
V_g_f <- model.matrix(prop ~ golfo - 1, fd_agg_g)


sdata <- list(
  N_f = nrow(fd),
  N_pred_f = nrow(fd_agg_gy),
  N_p = nrow(pd),
  N_pred_p = nrow(pd_agg_gmy),

  with = fd$intw,
  total = fd$intt,
  gan = pd$pa,
  hours = pd$h,

  N_g = length(levels(pd$golfo)),
  N_gm = length(levels(pd$gm)),
  N_gy = length(levels(pd$gy)),
  N_gmy = length(levels(pd$gmy)),
  gulf_id_gy = as.numeric(fd_agg_gy$golfo),

  X_g_p = X_g_p,
  X_g_f = X_g_f,
  X_gy_p = X_gy_p,
  X_gy_f = X_gy_f,
  X_gm_p = X_gm_p,
  X_gmy_p = X_gmy_p,

  Z_gm_p = Z_gm_p,
  Z_gy_p = Z_gy_p,
  Z_gmy_p = Z_gmy_p,

  Z_g_f = Z_g_f,
  Z_gy_f = Z_gy_f,

  V_g_p = V_g_p,
  V_gm_p = V_gm_p,

  V_g_f = V_g_f,

  prior_mu_f_sd = 2,
  prior_mu_p_sd = 10,
  prior_mu_p_mean = log(7),
  prior_sd_gy_f_sd = 1.5,
  prior_sd_gy_p_sd = 2,
  prior_sd_gmy_sd = 2,
  prior_sd_obs_sd = 1.5,
  prior_scale_sd = 5
)


# Model fit ---------------------------------------------------------------

smodel <- stan_model("GAP-GAF joint model.stan", verbose = TRUE)
m1 <- sampling(
  smodel, sdata, seed = 15498,
  cores = 10, chains = 10, iter = 2000, refresh = 10,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)

sm1 <- summary(m1)[[1]]

# saveRDS(m1, "GAP-GAF joint model_samples.rds")
# saveRDS(sm1, "GAP-GAF joint model_summary.rds")
# m1 <- readRDS("GAP-GAF joint model_samples.rds")
# sm1 <- readRDS("GAP-GAF joint model_summary.rds")

summary(sm1[, "n_eff"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 1260    4407    8924    8417   11866   16605       4
summary(sm1[, "Rhat"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.9992  0.9998  1.0002  1.0008  1.0013  1.0087       4

# View(sm1)
# plot(density(as.matrix(m1, "rho"), from = -1, to = 1), xlim = c(0, 1))


# GAP predictions ---------------------------------------------------------

# function to extract variables from Stan summary
takeout <- function(table, name) {
  string <- paste("\\b", name, "\\b", sep = "")
  x <- table[grep(string, rownames(table)), c("mean", "2.5%", "97.5%")] %>% as.data.frame()
  names(x) <- c("mean", "lower", "upper")
  rownames(x) <- NULL
  return(x)
}

# predictions by year
gap_pred_y <- cbind(pd_agg_gmy[, c("golfo", "mc", "year", "rate")], takeout(sm1, "gap_pred"))
gap_pred_y$predtype <- "years"

# predicted averages
gap_pred_avg <- cbind(pd_agg_gm[, c("golfo", "mc", "rate")], takeout(sm1, "gap_means"))
gap_pred_avg$year <- "avg"
gap_pred_avg$predtype <- "avg"
gap_pred_avg <- gap_pred_avg[, names(gap_pred_y)]

gap_pred <- rbind(gap_pred_y, gap_pred_avg)

# plot
ggplot(gap_pred, aes(x = year, y = mean, ymin = lower, ymax = upper,
                      colour = predtype)) +
  geom_point() +
  geom_errorbar() +
  facet_grid(rows = vars(golfo), cols = vars(mc)) +
  geom_point(gap_pred, mapping = aes(x = year, y = rate),
             inherit.aes = FALSE, colour = "black")

# GAF predictions ---------------------------------------------------------

# predictions by year
# Get mu and sd by year
gaf_pred_mu <- as.matrix(m1, "gaf_pred") %>% t
sd_obs <- as.matrix(m1, "sd_obs")
#dim(gaf_pred_mu)

gaf_pred_y_samples <- matrix(NA, nrow(gaf_pred_mu), ncol(gaf_pred_mu))
for(i in 1:nrow(gaf_pred_mu)) {
  print(i)
  for(j in 1:ncol(gaf_pred_mu)) {
    gaf_pred_y_samples[i, j] <- momentsLogitnorm(gaf_pred_mu[i, j], sd_obs[j])[1]
  }
}

## takes long, so I save it
#saveRDS(gaf_pred_y_samples, "gaf_pred_y_samples.rds")
# gaf_pred_y_samples <- readRDS("gaf_pred_y_samples.rds")
gaf_pred_temp <- apply(gaf_pred_y_samples, 1, qmean) %>% t %>% as.data.frame

gaf_pred_y <- cbind(fd_agg_gy[, c("golfo", "year", "prop")], gaf_pred_temp)
gaf_pred_y$predtype <- "years"

# means of each gulf

gaf_means_logit <- as.matrix(m1, "gaf_means_logit") %>% t
gaf_sds_logit <- as.matrix(m1, "gaf_sds_logit") %>% t

gaf_means_samples <- matrix(NA, nrow(gaf_means_logit), ncol(gaf_means_logit))
#gaf_sds_samples <- matrix(NA, nrow(gaf_means_logit), ncol(gaf_means_logit))

for(i in 1:nrow(gaf_means_samples)) {
  print(i)
  for(j in 1:ncol(gaf_means_samples)) {
    moments <- momentsLogitnorm(gaf_means_logit[i, j], gaf_sds_logit[i, j])
    gaf_means_samples[i, j] <- moments["mean"]
    #gaf_sds_samples[i, j] <- moments["var"] %>% sqrt
  }
}

gaf_means_temp <- apply(gaf_means_samples, 1, qmean) %>% t %>% as.data.frame

# predicted averages
gaf_pred_avg <- cbind(fd_agg_g[, c("golfo", "prop")], gaf_means_temp)
gaf_pred_avg$year <- "avg"
gaf_pred_avg$predtype <- "avg"
gaf_pred_avg <- gaf_pred_avg[, names(gaf_pred_y)]

gaf_pred <- rbind(gaf_pred_y, gaf_pred_avg)

# plot
ggplot(gaf_pred, aes(x = year, y = mean, ymin = lower, ymax = upper,
                     colour = predtype)) +
  geom_point() +
  geom_errorbar() +
  facet_wrap(vars(golfo)) +
  geom_point(gaf_pred, mapping = aes(x = year, y = prop),
             inherit.aes = FALSE, colour = "black")




# GAP means and sds comparison (relative diffs) ------------------------------

gap_pred_means <- gap_pred_avg[, c("golfo", "mc", "mean", "lower", "upper")]
names(gap_pred_means)[3:5] <- paste("mu", names(gap_pred_means)[3:5], sep = "_")

gap_means <- as.matrix(m1, "gap_means") %>% t

# M-C within gulf

# GN
qmean((gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] -
         gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ]) /
        gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ]) * 100
# At GN, the mean GAP on calves was 193.6939 % [108.4381, 280.4063] larger than on mothers.

# GSJ
qmean((gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ] -
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) /
        gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) * 100
# At GSJ, GAP on calves was 175.37423 % [71.80143, 329.75241] larger than on mothers.


# M-C averaging gulfs
tab_temp <- cbind(
 gsj = ((gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ] -
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) /
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) * 100,
 gn = ((gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] -
        gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ]) /
        gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ]) * 100
)
qmean(rowMeans(tab_temp))
# mean     lower    upper
# 184.5341 116.9412 270.6478

# GULF comparisons
# Calves
qmean((gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] -
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ]) /
        gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ]) * 100
# At GN, GAP on calves was 54.98471 % [-10.28851, 145.65526] larger than at GSJ

# Mothers
qmean((gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ] -
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) /
        gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) * 100
# At GN, GAP on mothers was 43.73623 % [-21.00909, 138.03287] larger than at GSJ

# Gulfs averaging M-C
tab_temp2 <- cbind(
  m = ((gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ] -
          gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) /
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) * 100,
  c = ((gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] -
          gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ]) /
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ]) * 100
)
qmean(rowMeans(tab_temp2))
#   mean       lower      upper
# 49.360467  -5.763266 127.517173


# Compute sd among years for observed estimated years
gap_pred_samples <- as.matrix(m1, "gap_pred") %>% t
gap_sd_samples <- aggregate(gap_pred_samples ~ golfo + mc, pd_agg_gmy, sd)
gap_sd_samples <- gap_sd_samples[with(gap_sd_samples, order(golfo, mc)), ]
gap_sds <- apply(as.matrix(gap_sd_samples[, -c(1, 2)]), 1, qmean) %>% t %>% as.data.frame
names(gap_sds) <- c("sd_mean", "sd_lower", "sd_upper")
gap_pred_mean_sd <- cbind(gap_pred_means, gap_sds)
aggregate(rate ~ golfo + mc, pd_agg_gmy, sd); gap_pred_mean_sd
#write.csv(gap_pred_mean_sd, "GAP predicted means and sds among years.csv")
# This matches the observed sds

# golfo mc  mu_mean mu_lower mu_upper   sd_mean  sd_lower sd_upper
# 1     N  M 1.977424 1.404243 2.974459 0.8599600 0.5993101 1.217442
# 3     N  C 5.709069 4.320433 8.126048 2.1105634 1.4226339 2.980445
# 2    SJ  M 1.429932 1.014456 2.196425 0.6319039 0.3857415 1.040849
# 4    SJ  C 3.832377 2.642629 5.837276 1.5887576 1.1839292 2.091190


# gap_pred_mean_sd shows averages and sds among years for each level at the
# response scale (attacks/h). (This table goes in the paper)

# The variability among years in the mean GAP was higher for calves in both gulfs.
# Comparing between gulfs, both calves and mothers showed higher variability among
# years in GN.

# Data ranges
range(pd$rate) # 0.00000 24.03101
# Estimated means ranges
range(gap_pred_y$mean) # 0.5315009 9.8004194


# GAP means comparison (quotients) -------------------------------------------

gap_pred_means <- gap_pred_avg[, c("golfo", "mc", "mean", "lower", "upper")]
names(gap_pred_means)[3:5] <- paste("mu", names(gap_pred_means)[3:5], sep = "_")

gap_means <- as.matrix(m1, "gap_means") %>% t

# M-C within gulf

# GN
qmean(gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] /
      gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ]) %>% round(digits = 2)
# At GN, GAP on calves was larger than on mothers by a factor of
# mean lower upper
# 2.94  2.08  3.80
qprob(gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] /
      gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ]) %>% round(digits = 2)
# 1

# GSJ
qmean(gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ] /
      gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) %>% round(digits = 2)
# At GSJ, GAP on calves was larger than on mothers by a factor of
# mean lower upper
# 2.75  1.72  4.30
qprob(gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ] /
      gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) %>% round(digits = 2)
# 1

# M-C averaging gulfs
tab_temp <- cbind(
  gsj = (gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ] /
         gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]),
  gn = (gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] /
        gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ])
)
qmean(rowMeans(tab_temp)) %>% round(digits = 2)
# Averaging gulfs, GAP on calves was larger than on mothers by a factor of
# mean lower upper
# 2.85  2.17  3.71
qprob(rowMeans(tab_temp)) %>% round(digits = 2)
# 1

# For table (M and C averaging gulfs)
mtemp1 <- aggregate(gap_means ~ mc, gap_pred_means, mean)
fac <- mtemp1[, 1]
mat <- mtemp1[, -1] %>% as.matrix
names(mtemp1)[1:3]
mtemp1[, 1]
avg_mc <- cbind.data.frame(fac, apply(mat, 1, qmean) %>% t)




# GULF comparisons
# Calves
qmean(gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] /
      gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ]) %>% round(digits = 2)
# At GN, GAP on calves was was larger than at GSJ by a factor of
# mean lower upper
# 1.55  0.90  2.46
qprob(gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] /
      gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ]) %>% round(digits = 2)
# 0.94

# Mothers
qmean(gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ] /
      gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) %>% round(digits = 2)
# At GN, GAP on mothers was was larger than at GSJ by a factor of
# mean lower upper
# 1.44  0.79  2.38
qprob(gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ] /
      gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]) %>% round(digits = 2)
# 0.89

# Gulfs averaging M-C
tab_temp2 <- cbind(
  m = (gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "M", ] /
       gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "M", ]),
  c = (gap_means[gap_pred_means$golfo == "N" & gap_pred_means$mc == "C", ] /
       gap_means[gap_pred_means$golfo == "SJ" & gap_pred_means$mc == "C", ])
)
qmean(rowMeans(tab_temp2)) %>% round(digits = 2)
# Averaging mother and calves, the mean GAP at GN was larger than at GSJ by a factor of
# mean lower upper
# 1.49  0.94  2.28
qprob(rowMeans(tab_temp2)) %>% round(digits = 2)
# 0.96


# For table (gulfs averaging mc)
mtemp2 <- aggregate(gap_means ~ golfo, gap_pred_means, mean)
fac <- mtemp2[, 1]
mat <- mtemp2[, -1] %>% as.matrix
names(mtemp2)[1:3]
mtemp2[, 1]
avg_gulf <- cbind.data.frame(fac, apply(mat, 1, qmean) %>% t)

#write.csv(rbind(avg_mc, avg_gulf), "GAP mean table.csv")

# GAF means and sds comparison (with quotients) ----------------------------

# means - relative difference
fd_agg_g; str(gaf_means_samples)
qmean((gaf_means_samples[fd_agg_g$golfo == "N", ] -
       gaf_means_samples[fd_agg_g$golfo == "SJ", ]) /
       gaf_means_samples[fd_agg_g$golfo == "SJ", ]) * 100
# mean GAF at GN was 22.08877 % [1.06845, 49.51850] larger than at GSJ,
#    mean     lower     upper
# 22.08877 -1.06845 49.51850


# means - factor
qmean(gaf_means_samples[fd_agg_g$golfo == "N", ] /
      gaf_means_samples[fd_agg_g$golfo == "SJ", ]) %>% round(digits = 2)
# GAF was higher at GN than at GSJ by a factor of
# mean lower upper
# 1.22  0.99  1.50 (p = 0.97)
qprob(gaf_means_samples[fd_agg_g$golfo == "N", ] /
      gaf_means_samples[fd_agg_g$golfo == "SJ", ]) %>% round(digits = 2)


# Compute sd among years for observed estimated years
gaf_pred_temp <- apply(gaf_pred_y_samples, 1, qmean) %>% t %>% as.data.frame

gaf_sd_samples <- aggregate(gaf_pred_y_samples ~ golfo, fd_agg_gy, sd)
gaf_sds <- apply(gaf_sd_samples[, -1], 1, qmean) %>% t %>% as.data.frame
names(gaf_sds) <- c("sd_mean", "sd_lower", "sd_upper")
temp <- gaf_pred_avg[, c("golfo", "mean", "lower", "upper")]
names(temp)[2:4] <- paste("mu_", c("mean", "lower", "upper"), sep = "")
gaf_pred_mean_sd <- cbind(temp, gaf_sds)
aggregate(prop ~ golfo, fd_agg_gy, sd); gaf_pred_mean_sd
#write.csv(gaf_pred_mean_sd, "GAF predicted means and sds among years.csv")

# golfo   mu_mean  mu_lower  mu_upper    sd_mean   sd_lower   sd_upper
# 1     N 0.2370230 0.2021499 0.2743933 0.05901500 0.04093806 0.07796790
# 2    SJ 0.1951362 0.1677094 0.2238117 0.04029184 0.02243792 0.05656813


# Data ranges
range(fd$prop, na.rm = TRUE) %>% round(digits = 2) # 0.00 0.79
# Estimated means ranges
range(gaf_pred$mean) %>% round(digits = 2)# 0.13 0.36


# Contrasts among years ---------------------------------------------------

# GAP
gap_means_samples <- as.matrix(m1, "gap_pred") %>% t
# aligned with pd_agg_gmy, the df with all predictors at year level.

# average gulfs
gap_means_temp <- aggregate(gap_means_samples ~
                              mc + year, data = pd_agg_gmy, mean)

gap_means <- data.frame(
  mc = gap_means_temp$mc,
  year = gap_means_temp$year,
  variable = NA
)
gap_means$attack <- as.matrix(gap_means_temp[, grep("V", colnames(gap_means_temp))])
gap_means$variable[gap_means$mc == "M"] <- "GAPM"
gap_means$variable[gap_means$mc == "C"] <- "GAPC"

# GAF
# average gulfs
gaf_means_temp <- aggregate(gaf_pred_y_samples ~ year, data = fd_agg_gy, mean)

gaf_means <- data.frame(
  mc = "MC",
  year = gaf_means_temp$year,
  variable = "GAF"
)
gaf_means$attack <- as.matrix(gaf_means_temp[, grep("V", colnames(gaf_means_temp))])

annual_means <- rbind(gap_means, gaf_means)

# number of samples
ns <- ncol(annual_means$attack)

annual_means$year <- as.numeric(as.character(annual_means$year))
annual_means$period <- "1995"
annual_means$period[annual_means$year > 1995] <- "2004-2010"
annual_means$period[annual_means$year > 2010] <- "2011-2013"
annual_means$period[annual_means$year > 2013] <- "2014-2019"

annual_means$period2 <- "1995"
annual_means$period2[annual_means$year > 1995] <- "2004-onward"

# Constrasts

# 2011-2013 > 1995
# GAPC
filt1 <- annual_means$period == "2011-2013" & annual_means$variable == "GAPC"#
filt2 <- annual_means$period == "1995" & annual_means$variable == "GAPC"#
d <- colMeans(annual_means$attack[filt1, ]) - annual_means$attack[filt2, ]
round(sum(d > 0) / ns, 2) # Pr = 1

# GAPM
filt1 <- annual_means$period == "2011-2013" & annual_means$variable == "GAPM"#
filt2 <- annual_means$period == "1995" & annual_means$variable == "GAPM"#
d <- colMeans(annual_means$attack[filt1, ]) - annual_means$attack[filt2, ]
round(sum(d > 0) / ns, 2) # Pr = 0.97

# GAF
filt1 <- annual_means$period == "2011-2013" & annual_means$variable == "GAF"#
filt2 <- annual_means$period == "1995" & annual_means$variable == "GAF"#
d <- colMeans(annual_means$attack[filt1, ]) - annual_means$attack[filt2, ]
round(sum(d > 0) / ns, 2) # Pr = 1

# Contrast only for GN

# Pr(1995 < all)
# GAPC
filt1 <- pd_agg_gmy$mc == "C" & pd_agg_gmy$year != "1995" & pd_agg_gmy$golfo == "N"
filt2 <- pd_agg_gmy$mc == "C" & pd_agg_gmy$year == "1995" & pd_agg_gmy$golfo == "N"
d <- colMeans(gap_means_samples[filt1, ]) -
  gap_means_samples[filt2, ]
round(sum(d > 0) / ns, 2) # Pr = 0.99

# GAF
filt1 <- fd_agg_gy$year != "1995" & fd_agg_gy$golfo == "N"
filt2 <- fd_agg_gy$year == "1995" & fd_agg_gy$golfo == "N"
d <- colMeans(gaf_pred_y_samples[filt1, ]) -
  gaf_pred_y_samples[filt2, ]
round(sum(d > 0) / ns, 2) # Pr = 0.99

# GAPM
filt1 <- pd_agg_gmy$mc == "M" &
  pd_agg_gmy$year %in% as.character(2014:2019) &
  pd_agg_gmy$golfo == "N"
filt2 <- pd_agg_gmy$mc == "M" &
  pd_agg_gmy$year == "1995" &
  pd_agg_gmy$golfo == "N"
d <- colMeans(gap_means_samples[filt1, ]) -
  gap_means_samples[filt2, ]
round(sum(d < 0) / ns, 2) # Pr = 0.79

filt1 <- pd_agg_gmy$mc == "M" &
  pd_agg_gmy$year %in% as.character(2014:2019) &
  pd_agg_gmy$golfo == "N"
filt2 <- pd_agg_gmy$mc == "M" &
  pd_agg_gmy$year %in% as.character(2004:2010) &
  pd_agg_gmy$golfo == "N"
d <- colMeans(gap_means_samples[filt1, ]) -
  colMeans(gap_means_samples[filt2, ])
round(sum(d < 0) / ns, 2) # Pr = 1


# Posterior predictive check ----------------------------------------------

# GAF
gaf_hat_logit <- as.matrix(m1, "gaf_hat_logit") %>% t

y_sim <- sapply(1:ncol(gaf_hat_logit), function(i) {
  rbinom(nrow(fd), size = fd$intt,
         prob = invlogit(gaf_hat_logit[, i] + rnorm(nrow(fd), 0, sd_obs[i])))
}) %>% as.numeric()

par(mfrow = c(1, 2))
hist(y_sim, breaks = 60, main = "simulation", xlim = c(0, 120))
hist(fd$intw, breaks = 15, main = "data", xlim = c(0, 120))
par(mfrow = c(1, 1))


# GAP
gap_hat <- as.matrix(m1, "gap_hat") %>% t
dim(gap_hat)
scale_hat <- as.matrix(m1, "scale")

gap_sim <- sapply(1:ncol(gap_hat), function(i) {
  rnbinom(nrow(pd), mu = gap_hat[, i], size = 1 / scale_hat[i])
}) %>% as.numeric()

par(mfrow = c(1, 2))
hist(y_sim[gap_sim <= 120], breaks = 60, main = "simulation", xlim = c(0, 120))
hist(pd$pa, breaks = 30, main = "data")
par(mfrow = c(1, 1))
# parece subestimar un poco los ceros, pero anda bien.



# Bayesian R2 --------------------------------------------------------------

# Following
# Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019).
# R-squared for Bayesian regression models. The American Statistician.

# gaf_pred_y_samples contains the predictions, without multiplying by N.
str(gaf_pred_y_samples) # 33 year-gulfs, 10000 samples

# define n
hist(fd$intt, breaks = 30)
median(fd$intt) # 48
size_gaf <- median(fd$intt)

# compute var_mu
var_mu_gaf <- apply(gaf_pred_y_samples * size_gaf, 2, var)
# compute residual variance using simulation
nsim <- 1000
var_res_gaf <- matrix(NA, nrow(gaf_pred_mu), ncol(gaf_pred_mu))
for(i in 1:nrow(var_res_gaf)) {
  print(i)
  for(j in 1:ncol(var_res_gaf)) {
    p_samples <- invlogit(rnorm(nsim, gaf_pred_mu[i, j], sd_obs[j]))
    with_sim <- rbinom(nsim, size = size_gaf, p_samples)
    var_res_gaf[i, j] <- var(with_sim)
  }
}
var_res_gaf_mean <- colMeans(var_res_gaf)
plot(density(var_res_gaf_mean, from = 0))

plot(density(R2_gaf <- var_mu_gaf / (var_mu_gaf + var_res_gaf_mean),
             from = 0, to = 1), main = "GAF R2 posterior distribution")
summary(R2_gaf); qmean(R2_gaf)

# This R2 includes variability asociated to not sampling Inf elements from
# each binomial event. If the binomial size were large (say 10000), R2 would be
# a little higher:

# Compute R2 without considering the binomial variance (only variance in p)
var_mu_gaf_p <- apply(gaf_pred_y_samples, 2, var)
var_res_gaf_p <- matrix(NA, nrow(gaf_pred_mu), ncol(gaf_pred_mu))
for(i in 1:nrow(var_res_gaf)) {
  print(i)
  for(j in 1:ncol(var_res_gaf)) {
    var_res_gaf_p[i, j] <- momentsLogitnorm(gaf_pred_mu[i, j], sd_obs[j])["var"]
  }
}
var_res_gaf_mean_p <- colMeans(var_res_gaf_p)
plot(density(var_res_gaf_mean_p, from = 0))
plot(density(R2_gaf_p <- var_mu_gaf_p / (var_mu_gaf_p + var_res_gaf_mean_p),
             from = 0, to = 1), main = "GAF R2 posterior distribution certain")
summary(R2_gaf_p); qmean(R2_gaf_p)
# It approaches 1 as sd_obs approaches 0, which is the binomial R2 when size
# approaches Inf.

# Compare both R2:
plot(density(R2_gaf_p, from = 0, to = 1), main = "GAF R2 posterior distribution",
     ylim = c(0, 13))
lines(density(R2_gaf, from = 0, to = 1), col = 2)
text(0.4, 2, "for p")
text(0.6, 4, "for number of intervals with attacks", col = 2)

# use R2_gaf_p
qmean(R2_gaf_p)
# mean     lower     upper
# 0.1868213 0.1128324 0.2644865

# We choose R2_gaf_p because it represents explained varibility related to the
# process, while R2_gaf also considers variability related to the observation:
# with lower observation intervals (size in the binomial distribution), the
# unexplained variability would be larger, but only because of our poor sampling.


##### R2 for GAP

gap_pred_samples <- as.matrix(m1, pars = "gap_pred") %>% t # year-level predictions
phi_samples <- as.matrix(m1, pars = "phi")

# compute var_mu
var_mu_gap <- apply(gap_pred_samples, 2, var)
# From
# https://mc-stan.org/docs/2_20/functions-reference/nbalt.html
# var(Y) = mu + mu ^ 2 / phi
# E(Y) = mu
var_res_gap <- matrix(NA, nrow(gap_pred_samples), ncol(gap_pred_samples))
for(i in 1:nrow(var_res_gap)) {
  var_res_gap[i, ] <- gap_pred_samples[i, ] + gap_pred_samples[i, ] ^ 2 / phi_samples
}
var_res_gap_mean <- colMeans(var_res_gap)

plot(density(R2_gap <- var_mu_gap / (var_mu_gap + var_res_gap_mean),
             from = 0, to = 1), main = "GAP R2 posterior distribution")
summary(R2_gap)
qmean(R2_gap)
# mean     lower     upper
# 0.2732754 0.2355769 0.3156198

mean(var_mu_gap); mean(var_res_gap_mean)
plot(var_mu_gap ~ var_res_gap_mean)

## R2 for GAP separating mother and calf

nrow(gap_pred_y)      # data with mother-calf variable (mc)
str(gap_pred_samples)
mothers <- gap_pred_y$mc == "M"
calves <- gap_pred_y$mc == "C"


## MOTHERS

gap_pred_samples_m <- gap_pred_samples[mothers, ]

# compute var_mu
var_mu_gap_m <- apply(gap_pred_samples[mothers, ], 2, var)
# residual variance
var_res_gap_m <- matrix(NA,
                        nrow(gap_pred_samples_m),
                        ncol(gap_pred_samples_m))
for(i in 1:nrow(var_res_gap_m)) {
  var_res_gap_m[i, ] <- gap_pred_samples_m[i, ] + gap_pred_samples_m[i, ] ^ 2 / phi_samples
}
var_res_gap_mean_m <- colMeans(var_res_gap_m)

plot(density(R2_gap_m <- var_mu_gap_m / (var_mu_gap_m + var_res_gap_mean_m),
             from = 0, to = 1), main = "GAPM R2 posterior distribution")
summary(R2_gap_m)
qmean(R2_gap_m)
mean(var_mu_gap_m); mean(var_res_gap_mean_m)
# mean     lower     upper
# 0.1537670 0.1038135 0.2161866


## CALVES

gap_pred_samples_c <- gap_pred_samples[calves, ]

# compute var_mu
var_mu_gap_c <- apply(gap_pred_samples[calves, ], 2, var)
# residual variance
var_res_gap_c <- matrix(NA,
                        nrow(gap_pred_samples_c),
                        ncol(gap_pred_samples_c))
for(i in 1:nrow(var_res_gap_c)) {
  var_res_gap_c[i, ] <- gap_pred_samples_c[i, ] + gap_pred_samples_c[i, ] ^ 2 / phi_samples
}
var_res_gap_mean_c <- colMeans(var_res_gap_c)

plot(density(R2_gap_c <- var_mu_gap_c / (var_mu_gap_c + var_res_gap_mean_c),
             from = 0, to = 1), main = "GAPC R2 posterior distribution")
summary(R2_gap_c)
qmean(R2_gap_c)
mean(var_mu_gap_c); mean(var_res_gap_mean_c)

# mean     lower     upper
# 0.1788322 0.1309413 0.2351871


# Plots for posterior predictive distribution -----------------------------

# GAF

fa_y_lims <- matrix(NA, nrow(gaf_pred), 2)
npost <- ncol(gaf_means_logit)
for(i in 1:nrow(gaf_pred_mu)) {
  print(i)
  sim <- rnorm(npost, gaf_pred_mu[i, ], gaf_sds_logit) %>% invlogit
  fa_y_lims[i, ] <- hdint(sim)
}
colnames(fa_y_lims) <- c("sim_lower", "sim_upper")
gaf_pred_sim <- cbind(gaf_pred, fa_y_lims)

# plot
ggplot(gaf_pred_sim, aes(x = year, y = mean, ymin = lower, ymax = upper,
                     colour = predtype)) +
  geom_point() +
  geom_errorbar() +
  geom_linerange(mapping = aes(x = year, y = mean, ymin = sim_lower,
                               ymax = sim_upper), col = "red",
                 size = 2, alpha = 0.3) +
  facet_wrap(vars(golfo)) +
  geom_point(gaf_pred, mapping = aes(x = year, y = prop),
             inherit.aes = FALSE, colour = "black") +
  geom_point(fd, mapping = aes(x = year, y = prop),
             inherit.aes = FALSE, colour = "black", shape = 1,
             alpha = 0.4)
# write.csv(gaf_pred_sim, "GAF prediction table.csv")

# In these plots we add the data (open points) and the posterior predictive
# distribution for every combination of predictor variables (red bands). These are
# 95 % highest density intervals, and should contain approximately 95 % of the data
# points if the model is good.

# GAP

gap_y_lims <- matrix(NA, nrow(gap_pred), 2)
npost <- ncol(gap_pred_samples)
phi_samples <- as.matrix(m1, pars = "phi")

for(i in 1:nrow(gap_pred_samples)) {
  print(i)
  sim <- rnbinom(npost, mu = gap_pred_samples[i, ], size = phi_samples)
  gap_y_lims[i, ] <- hdint(sim)
}
colnames(gap_y_lims) <- c("sim_lower", "sim_upper")
gap_pred_sim <- cbind(gap_pred, gap_y_lims)

ggplot(gap_pred_sim, aes(x = year, y = mean, ymin = lower, ymax = upper,
                         colour = predtype)) +
  geom_point() +
  geom_errorbar() +
  geom_linerange(mapping = aes(x = year, y = mean, ymin = sim_lower,
                               ymax = sim_upper), col = "red",
                 size = 2, alpha = 0.3) +
  facet_grid(rows = vars(golfo), cols = vars(mc)) +
  geom_point(gap_pred_sim, mapping = aes(x = year, y = rate),
             inherit.aes = FALSE, colour = "black") +
  geom_point(pd, mapping = aes(x = year, y = rate),
             inherit.aes = FALSE, colour = "black", shape = 1,
             alpha = 0.4)

# write.csv(gap_pred_sim, "GAP prediction table.csv")