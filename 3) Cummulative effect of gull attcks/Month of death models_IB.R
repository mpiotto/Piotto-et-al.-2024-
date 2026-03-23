                ###### MOnth of death and gull attacks #####

# Article title: Calf mortality in a whale population varies with seabird attacks
# Authors: María Piotto, Iván Barberá, Mariano Sironi, Victoria J. Rowntree,
# Marcela M. Uhart, Macarena Agrelo, Alejandro A. Fernández Ajó, Jon Seger,
#Carina F. Marón.

# Scrip author: Iván Barberá
# Questions can be addressed to ivanbarbera93@gmail.com or to mpiotto@mi.unc.edu.ar


# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayestestR)
library(brms)
library(DHARMa)
library(ggrepel)
library(ggnewscale)

# Gull attack data belong to the Southern Right Whale Research Program (Instituto
# de Conservación de Ballenas, Argentina). Calf mortality data belong to the
# Southern Right Whale Health Monitoring Program of Península Valdés, Argentina
# (University of California, Davis and Instituto de Conservación de Ballenas).
# Requests to use the data should be addressed to Marcela Uhart muhart@ucdavis.edu
# and Mariano Sironi mariano.sironi@icb.org.ar

# intra-annual mortality data
mintra <- read.csv(".....\\Deaths_month_gulf_year.csv", sep = ";")

# In this database you will find the number of dead calves by month (dmonth) and
# by year (dyear) in each gulf. The columns "year", "golfo" (= gulf), "month" and
# "nmonth" indicate the month (nmonth: numeric, month: categorical), the year (year),
# and the gulf (golfo) in which the deaths occured.


# attack variables data
datt <- read.csv("....\\Calf mortality data.csv", sep = ";")


# Functions  --------------------------------------------------------------

hdi_lower <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_low"]])
}

hdi_upper <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_high"]])
}

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

hdint <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, ci$CI_high)
  names(result) <- paste(rep(name, 2), c("lower", "upper"), sep = "_")
  result
}

# equal-tailed intervals
eti_mean <- function(x, name = "mu") {
  ci <- quantile(x, probs = c(0.025, 0.975), method = 8)
  result <- c(ci[1], mean(x), ci[2])
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

eti <- function(x, name = "mu") {
  result <- quantile(x, probs = c(0.025, 0.975), method = 8)
  names(result) <- paste(rep(name, 2), c("lower", "upper"), sep = "_")
  result
}

# inverse cumulative logit
inv_cumlogit <- function(logit_thresholds) {
  cum_probs <- plogis(logit_thresholds)
  p <- c(0, cum_probs, 1) %>% diff
  return(p)
}

mintra$month %>% unique
mintra$month_num <- as.numeric(factor(mintra$month, levels = c(
  "jun", "jul", "ago", "sep", "oct", "nov", "dic"
)))

ord <- with(mintra, order(
  year, golfo, month_num
))
mintra <- mintra[ord, ]


# Longanize intra annual data
d <- do.call("rbind", lapply(1:nrow(mintra), function(i) {
  #i = 2
  N <- mintra[i, "dmonth"]
  reps <- mintra[rep(i, N), c("year", "golfo",
                              "month", "month_num",
                              "dmonth")]
  rownames(reps) <- NULL
  return(reps)
}))


# Prepare data for plots

data_long <- left_join(d, datt[, c("year", "golfo", "gaf", "gapm", "gapc")],
                       by = c("year", "golfo"))

data_long <- data_long[complete.cases(data_long), ]

data_agg <- aggregate(month_num ~ year + golfo, data = d,
                      FUN = mean)
data_agg <- left_join(data_agg, datt[, c("year", "golfo", "gaf", "gapm", "gapc", "dead", "born")],
                      by = c("year", "golfo"))
# get N by year
agg_n <- aggregate(dmonth ~ year + golfo, data = mintra, FUN = sum)
data_agg <- left_join(data_agg, agg_n, by = c("year", "golfo"))

data_agg <- data_agg[complete.cases(data_agg), ]
data_agg$mortality <- data_agg$dead / data_agg$born

agg_long <- pivot_longer(data_agg,
                         cols = which(colnames(data_agg) %in% c("gaf", "gapm", "gapc")),
                         names_to = "attack_variable",
                         values_to = "attack")


# Descriptive plot with lm -------------------------------------------------

ggplot(agg_long, aes(x = attack, y = month_num + 5,
                     color = mortality, size = dmonth)) +

  geom_smooth(method = "lm", se = TRUE, alpha = 0.4, size = 0.4,
              linetype = "dashed", color = "black") +

  geom_point() +

  scale_color_viridis() +
  facet_grid(rows = vars(golfo), cols = vars(attack_variable),
             scales = "free_x") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("Average month of death") +
  xlab("Attack") +
  labs(color = "Mortality", size = "Number of deaths") +
  ylim(7.5, 11)



# Models fit ---------------------------------------------------------------

d2 <- left_join(d, datt[, c("year", "golfo", "gaf", "gapm", "gapc", "dead", "born")],
                by = c("year", "golfo"))

d2 <- d2[complete.cases(d2), ]
# head(d2, 10)

d2$yg <- paste(d2$year, d2$golfo, sep = "_")


# MODEL FIT
# the year-gulf effect variance varies between gulfs.

# For GAPM
model_pam <- brm(
  month_num ~ (1 | gr(yg, by = golfo)) + golfo * gapm,
  data = d2, family = cumulative("logit"),
  cores = 4, chains = 4, warmup = 500, iter = 2500
) # 2669; 1.00 [worst neff and rhat]
# For GAPC
model_pac <- brm(
  month_num ~ (1 | gr(yg, by = golfo)) + golfo * gapc,
  data = d2, family = cumulative("logit"),
  cores = 4, chains = 4, warmup = 500, iter = 2500
) # 2682; 1.00
# For GAF
model_fa <- brm(
  month_num ~ (1 | gr(yg, by = golfo)) + golfo * gaf,
  data = d2, family = cumulative("logit"),
  cores = 4, chains = 4, warmup = 500, iter = 2500
) # 3256; 1.00
models_list <- list("gapm" = model_pam, "gapc" = model_pac, "gaf" = model_fa)

# saveRDS(models_list, "mortalidad_intraanual_models_samples.rds")
# models_list <- readRDS("mortalidad_intraanual_models_samples.rds")

# Predictions for the mean and for hypothetical years ----------------------

at_names <- c("gapm", "gapc", "gaf")

### Data to predict

# edit: include the full range of each predictor marginal to the gulf, so the
# curves don't appear trimmed and the x axis extends in the full range.
# Then, plot ranges without data with a lighter ribbon.

pdata_list <- vector(mode = "list", length = 3)
names(pdata_list) <- at_names
nr <- 150

for(v in at_names) {
  # v = "fa"

  # create initial attack variable sequence from marginal min to marginal max
  at_init <- seq(min(d2[, v]), max(d2[, v]), length.out = nr)

  # add minimum and a little less for GN
  min_gn <- min(d2[d2$golfo == "gn", v])
  at_extra_gn <- c(min_gn, min_gn - 0.0001)

  # add maximum and a little more for GSJ
  max_gsj <- max(d2[d2$golfo == "gsj", v])
  at_extra_gsj <- c(max_gsj, max_gsj + 0.0001)

  # merge extra values with sequence
  at_seq <- c(at_init, at_extra_gn, at_extra_gsj)
  at_seq <- at_seq[order(at_seq)]

  # create data frame to predict
  pdata <- expand.grid(at = at_seq,
                       golfo = c("gn", "gsj"))

  pdata$row <- 1:nrow(pdata)
  names(pdata) <- c(v, "golfo", "row")
  pdata$attack_variable <- v

  # add in_range column to indicate whether there is data or not for the x value
  pdata$in_range <- factor(rep("in", nrow(pdata)), levels = c("in", "out"))

  # get rows out of range:
  out_filter <- which(
    # out golfo nuevo
    (pdata$golfo == "gn" & pdata[, v] < min_gn) |
    # out golfo san josé
    (pdata$golfo == "gsj" & pdata[, v] > max_gsj)
  )

  pdata$in_range[out_filter] <- "out"

  pdata_list[[v]] <- pdata
}

# list to fill with predictions
preds_list <- vector(mode = "list", length = 3)

# extract param names
par_samples <- as.matrix(models_list$pam)
par_names <- colnames(par_samples)


# loop over models
for(m in 1:3) {

# get parameters
# m = 1
# lp is the linear combination of covariates that is added to the
# logit-thresholds
lp_samples <- fitted(models_list[[m]], pdata_list[[m]], scale = "linear", summary = F,
                     re_formula = NA) %>% t

sigmas <- as.matrix(models_list[[m]],
                    variable = c("sd_yg__Intercept:golfogn",
                                 "sd_yg__Intercept:golfogsj")) %>% t
dmgulf <- model.matrix(~ golfo - 1, data = pdata_list[[m]]) # gn primero
sigmas_mat <- dmgulf %*% sigmas
# chequear que dmgulf y sigmas_mat ordenen igual los golfos.


# get intercepts (baseline logit-thresholds)
intercepts <- as.matrix(models_list[[m]],
                        variable = par_names[grep("b_Intercept", par_names)])

# Loop over predictions not to overload RAM

preds_mat <- matrix(NA, nrow(pdata_list[[m]]), 5)
colnames(preds_mat) <- c("mean_mean", "mean_lower", "mean_upper",
                         "year_lower", "year_upper")
# mean_ is the prediction for an average year, year_ is prediction inconditional to years

nsim <- 200 # years to sample
npost <- ncol(lp_samples)
npred <- nrow(lp_samples)

# standardized random years
z_samples <- matrix(rnorm(nsim * npost), npost, nsim)

# to compute inverse cumulative logit
m1 <- matrix(1, npost, nsim)

for(i in 1:npred) {
  print(paste("model ", m, "; row ", i, sep = ""))
  # simulate linear predictors in every posterior sample, corresponding to new
  # years
  #i = 1

  # get simulated years, always using the same standardized normal samples
  matsim <- z_samples * sigmas_mat[i, ] + lp_samples[i, ]

  # create array to compute thresholds' cumulative probabilities
  pcum_arr <- array(NA, dim = c(npost, nsim, 6)) # 6 thresholds
  for(k in 1:6) {
    pcum_arr[, , k] <- plogis(intercepts[, k] - matsim) # importante, se resta.
  }

  # Compute inv_cumlogit

  # (This is the heavy step)
  # faster way:
  enlarged <- abind::abind(pcum_arr, m1, along = 3)
  # using apply is really slow:
  # probs <- apply(enlarged, 1:2, diff)
  # so, compute diff by hand:
  probs <- enlarged[, , 1:7] # make a copy
  for(k in 2:7) probs[, , k] <- enlarged[, , k] - enlarged[, , k-1]

  # Compute means:
  # here I turn classes (1:7) into months (6:12)
  # Apply is slow, use matrix multiplication
  means_raw <- probs[, , 1] # just create the space.
  for(s in 1:nsim) {
    means_raw[, s] <- as.numeric(probs[, s, ] %*% (6:12))
  }
  # npost X nsim matrix.
  # str(means_raw)

  # prediction interval for new years
  preds_mat[i, c("year_lower", "year_upper")] <- eti(as.numeric(means_raw))

  # credible interval for the average year
  preds_mat[i, c("mean_lower", "mean_mean", "mean_upper")] <- eti_mean(rowMeans(means_raw))

}

preds_list[[m]] <- preds_mat
}

# set equal names in pdata
for(i in 1:3) names(pdata_list[[i]])[1] <- "attack"

pdata_table <- do.call("rbind", pdata_list)
preds_table <- do.call("rbind", preds_list)
preds <- cbind(pdata_table, preds_table)

# saveRDS(preds, "mortalidad_intraanual_models_prediction_table.rds")
# preds <- readRDS("mortalidad_intraanual_models_prediction_table.rds")


# Plot --------------------------------------------------------------------

ggplot(agg_long,
       aes(x = attack, y = month_num + 5,
           color = mortality, size = dmonth)) +

  # estimation
  #   ci unconditional to year
  geom_ribbon(data = preds,
              mapping = aes(x = attack, y = mean_mean,
                            ymin = year_lower, ymax = year_upper),
              inherit.aes = F, alpha = 0.0, color = "black",
              linetype = "dashed", size = 0.3) +
  #   ci for average year
  geom_ribbon(data = preds,
              mapping = aes(x = attack, y = mean_mean,
                            ymin = mean_lower, ymax = mean_upper,
                            alpha = in_range),
              inherit.aes = F, color = NA,
              show.legend = FALSE) + # remove legend for in_range

  # set alpha scale for ribbon
  scale_alpha_manual(values = c(0.3, 0.1)) +

  # set new alpha scale for line
  new_scale("alpha") +
  scale_alpha_manual(values = c(1, 0.6)) +

  #   mean
  geom_line(data = preds,
            mapping = aes(x = attack, y = mean_mean, alpha = in_range),
            inherit.aes = F, size = 0.4,
            show.legend = FALSE) + # remove legend for in_range

  # data
  geom_point() +

  # year
  geom_text_repel(data = agg_long,
                  mapping = aes(x = attack, y = month_num + 5, label = year),
                  inherit.aes = FALSE, size = 2,
                  box.padding = 0.1, max.overlaps = Inf) +


  # settings
  scale_color_viridis() +
  facet_grid(rows = vars(golfo), cols = vars(attack_variable),
             scales = "free_x") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("Average month of death") +
  xlab("Attack variable") +
  labs(color = "Mortality", size = "Number of deaths")


# Residuals (posterior predictive check) ----------------------------------

sim_list <- list(mode = "vector", length = 3)
for(i in 1:3) {
  sim_list[[i]] <- predict(models_list[[i]], summary = F) %>% t
}

# get DHARMa residuals
res_list <- list(mode = "vector", length = 3)
for(i in 1:3) {
  res_list[[i]] <- createDHARMa(simulatedResponse = sim_list[[i]],
                                observedResponse = d2$month_num,
                                integerResponse = T)
  plot(res_list[[i]])
}

plot(res_list[[1]])
plot(res_list[[2]])
plot(res_list[[3]])

names(res_list) <- at_names

# add data to compare.
res_list$data <- d2
res_list$notes <- "Data matching the residuals vector is provided to check patterns."

# saveRDS(res_list, "mortalidad_intraanual_models_residuals.rds")



# R-squared -----------------------------------------------------------------
# r2 = var(mu) / (var(mu) + var(y))

# Compute r2 for the annual averages, not the raw observations.

# predict annual averages
v = "gapc"
data_agg$yg <- paste(data_agg$year, data_agg$golfo, sep = "_")
ann_pred_p <- fitted(models_list[[v]], newdata = data_agg, summary = FALSE,
                     scale = "response")
str(ann_pred_p)
ann_pred_m_temp <- ann_pred_p
for(i in 1:dim(ann_pred_p)[3]) ann_pred_m_temp[, , i] <- ann_pred_p[, , i] * ((6:12)[i])
ann_pred_m <- apply(ann_pred_m_temp, 1:2, sum) %>% t
str(ann_pred_m)
# get residual variances
var_res <- apply(ann_pred_m, 2, var)

# predict average of many years
mean_pred_p <- fitted(models_list[[v]], newdata = data_agg, summary = FALSE,
                      scale = "response", re_formula = NA) # do not include raneffs
str(mean_pred_p)
mean_pred_m_temp <- mean_pred_p
for(i in 1:dim(mean_pred_p)[3]) {
  mean_pred_m_temp[, , i] <- mean_pred_p[, , i] * ((6:12)[i])
}
mean_pred_m <- apply(mean_pred_m_temp, 1:2, sum) %>% t
str(mean_pred_m)
# get residual variances
var_mu <- apply(mean_pred_m, 2, var)

r2 <- var_mu / (var_mu + var_res)
plot(density(r2, from = 0, to = 1))
(eti_mean(r2))

# R2 for mortality month models:
#        lower (2.5%)  mean          upper (97.5%)
# gaf     0.03241505    0.24185704    0.50286431
# gapm    0.05812999    0.28145975    0.52149156
# gapc    0.02039316    0.19306263    0.45852898


# R2 for raw observations? ------------------------------------------------

y_pred <- predict(models_list[[v]], newdata = data_agg, summary = FALSE,
                  scale = "response") %>% t
y_pred <- y_pred + 5 # turn into months
var_res_obs <- apply(y_pred, 2, var)
r2_obs <- var_mu / (var_mu + var_res_obs)
plot(density(r2_obs, from = 0, to = 1))
(eti_mean(r2_obs)) # much smaller, makes sense


# Probabilities of positive slopes ----------------------------------------

slope_prob <- lapply(at_names, function(a) {
  # a = "pam"
  d <- as.data.frame(models_list[[a]],
                     variable = c(paste("b_", a, sep = ""),
                                  paste("b_golfogsj:", a, sep = "")))
  colnames(d) <- c("GN", "GSJ_diff")
  d$GSJ <- d$GN + d$GSJ_diff
  d$gulf_avg <- rowMeans(as.matrix(d[, c("GN", "GSJ")]))

  probs <- apply(
    as.matrix(d[, c("GN", "GSJ", "gulf_avg")]),
    2,
    function(x) sum(x > 0) / length(x)
  )

  return(round(probs, 2))
})
names(slope_prob) <- at_names

slope_prob

# $pam
# GN       GSJ      gulf_avg
# 0.98     0.77     0.97

# $pac
# GN       GSJ      gulf_avg
# 0.81     0.30     0.60

# $fa
# GN       GSJ      gulf_avg
# 0.95     0.40     0.83