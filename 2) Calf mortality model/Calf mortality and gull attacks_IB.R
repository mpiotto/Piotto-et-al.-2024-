
                     ###### Calf mortality models #####

# Article title: Calf mortality in a whale population varies with seabird attacks
# Authors: María Piotto, Iván Barberá, Mariano Sironi, Victoria J. Rowntree,
# Marcela M. Uhart, Macarena Agrelo, Alejandro A. Fernández Ajó, Jon Seger,
#Carina F. Marón.

# Scrip author: María Piotto
# Questions can be addressed to mpiotto@mi.unc.edu.ar or to ivanbarbera93@gmail.com

# About this code:

# This script was organized so you can run each model and each prediction
# regardless of the other models or predictions. I.e., you do not need
# to run the whole code to test just a single prediction; you only have to run
# the correspondent model that you would like to test and the prediction in which
# you are interested in.
#
# You will find 5 predictions done with each model. Predictions were presented
# following the order in which they appeared in the manuscript:
#     (1) the predicted absolute increase in calf mortality (the difference between
#         calf mortality when gull attack indexes were at its minimum and maximum
#         observed values in each gulf). Inside this section, you will also find
#         the predicted absolute increase in calf mortality when GAPC is equal
#         to the maximum and minimum value of GAPM (considering and without
#         considering the gulfs).
#     (2) partial predictions for the gull attack index involved in the model;
#     (3) the estimated calves' probability of dying when gull attacks would have
#         been equal to zero;
#     (4) the estimated calves' probability of dying in an average attack year;
#    (4.1) the quotients between (3) and (4); and,
#     (5) partial predictions of the SST anomalies.
#
# Finally, for each model you will find its posterior predictive check and the
# code to calculate its associated R2.
#
# Note that predictions averaging the three models are at the end of this script,
# and that if you like to run them you will need parts of the code that is before
# them. We recommend you to run the whole code before testing these predictions,
# instead of looking for the parts of the code you will need to run just them.


# Library, data and some settings --------------------------------------------

library(jagsUI)
library(arm)
library(tidyverse)
library(cowplot)
library(plyr)
library(extraDistr)
library(colorspace)
theme_set(theme_bw())

# Gull attack data belong to the Southern Right Whale Research Program (Instituto
# de Conservación de Ballenas, Argentina). Calf mortality data belong to the
# Southern Right Whale Health Monitoring Program of Península Valdés, Argentina
# (University of California, Davis and Instituto de Conservación de Ballenas).
# Requests to use the data should be addressed to Marcela Uhart muhart@ucdavis.edu
# and Mariano Sironi mariano.sironi@icb.org.ar


cmd <- read.csv("Calf mortality data.csv",
                header=TRUE, sep = ";")

borns <- cmd$born # number of newborns by year and gulf
deaths <- cmd$dead # number of dead calves by year and gulf
pac <- cmd$egapc # standardized GAPC
pam <- cmd$egapm # standardized GAPM
fa <- cmd$egaf # standardized GAF
n <- 31 # 31 rows in our data
sst <- cmd$esst # standardized SST anomalies
g <- as.factor(cmd$golfo)
gn_dummy <- as.numeric(g == "gn") # Golfo Nuevo
gsj_dummy <- as.numeric(cmd$golfo == "gsj") # Golfo San José
obs_gn <- subset(cmd, golfo == "gn") # observed values in Golfo Nuevo
obs_gsj <- subset(cmd, golfo == "gsj") # observed values in Golfo San José


# For the R2 calculation
# function to summarize posterior distributions:
qmean <- function(x) {
  quants <- quantile(x, probs = c(0.025, 0.975), method = 8)
  res <- c(mean(x), quants)
  names(res) <- c("mean", "lower", "upper")
  return(res)
}

# GAPC model -------------------------------------------------------------------

cat(file = "rp.bug",
    "
model{
# likelihood
for (i in 1:n) {
  deaths[i] ~ dbinom(p[i], borns[i])
  p[i] ~ dbeta(alpha[i], beta[i])
  alpha[i] <- theta[i] * phi[i]
  beta[i] <- (1 - theta[i]) * phi[i]
  logit(theta[i]) <- b0_gn * gn_dummy[i] + b0_gsj * gsj_dummy[i] +
                     b1_gn * gn_dummy[i] * pac[i] + b1_gsj * gsj_dummy[i] * pac[i] +
                     b2 * sst[i]

  phi[i] <-  1 / scale_gn * gn_dummy[i] + 1 / scale_gsj * gsj_dummy[i]
}


# priors

scale_gn ~ dgamma(1, 1 / 0.1) # medida de dispersi?n
scale_gsj ~ dgamma(1, 1 / 0.1)
b0_gn ~ dnorm(0, 1 / 1.5 ^ 2) # nodo estoc?stico (dsps de virulilla)
b0_gsj ~ dnorm(0, 1 / 1.5 ^ 2)
b1_gn ~ dnorm(0, 1 / 1.5 ^ 2)
b1_gsj ~ dnorm(0, 1 / 1.5 ^ 2)
b2 ~ dnorm(0, 1 / 1.5 ^ 2)
}"
)

m.data <- list(deaths = deaths, borns = borns, pac = pac, n = n, sst = sst,
               gn_dummy = gn_dummy, gsj_dummy = gsj_dummy)

inits <- function() list(b0_gn = runif(1, -1, 1),
                         b0_gsj = runif(1, -1, 1),
                         b1_gn = runif(1, -1, 1),
                         b1_gsj = runif(1, -1, 1),
                         b2 = runif(1, -1, 1),
                         scale_gn = runif(1, 0, 12),
                         scale_gsj = runif(1, 0, 12),
                         p = runif(31, 0, 1))

params <- c("b0_gn","b0_gsj", "b1_gn","b1_gsj", "b2", "scale_gn",
            "scale_gsj", "p", "phi", "theta")

ni <- 20000 # number of iterations
nt <- 1 # thining rate
nb <- 5000 # iteration for "burn in"
nc <- 2 # number of chains

jags.sim <- jags(data = m.data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "rp.bug",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb)

print(jags.sim) # Posterior distributions - summary

# Taking the posterior distribution of each parameter
b0_gn_hat <- jags.sim$sims.list$b0_gn
b0_gsj_hat <- jags.sim$sims.list$b0_gsj
b1_gn_hat <- jags.sim$sims.list$b1_gn
b1_gsj_hat <- jags.sim$sims.list$b1_gsj
b2_hat <- jags.sim$sims.list$b2
phi_gn <- 1 / jags.sim$sims.list$scale_gn
phi_gsj <- 1 / jags.sim$sims.list$scale_gsj
phi_hat <- jags.sim$sims.list$phi
theta_hat <- jags.sim$sims.list$theta


# Probability of slopes > 0 -----------------------------------------------

round(sum(b1_gn_hat > 0) / length(b1_gn_hat), 2)   # 0.89
round(sum(b1_gsj_hat > 0) / length(b1_gsj_hat), 2) # 0.78
round(sum(b2_hat < 0) / length(b2_hat), 2)         # 0.75

# Calves' probability of dying - absolute increase-----------------------------

# Golfo Nuevo
in_pac_gn <- expand.grid(pac = c(min(pac[1:15]), max(pac[1:15])),
                         sst = mean(sst))

in_matrix_pac_gn <- matrix(0, 2, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  in_matrix_pac_gn[,i] <- invlogit(b0_gn_hat[i] +
                                   b1_gn_hat[i] * in_pac_gn$pac +
                                   b2_hat[i] * mean(sst))
}

pacmin_gn <- c(in_matrix_pac_gn[1, ]) # post distribution of calf mortality when
# GAPC was equal to its minimum observed value in Golfo Nuevo
pacmax_gn <- c(in_matrix_pac_gn[2, ]) # post distribution of calf mortality when
# GAPC was equal to its maximum observed value in Golfo Nuevo

dif_pac_gn <- pacmax_gn - pacmin_gn # post distribution of the difference (i.e.,
# absolute increase)

# Summarazing the posterior distribution
mean(dif_pac_gn) # 0.16
quantile(dif_pac_gn, probs = c(0.025, 00.975)) # 2.5% = -0.1, 97.5 % = 0.42

# Golfo San José
in_pac_gsj <- expand.grid(pac = c(min(pac[16:31]), max(pac[16:31])),
                            sst = mean(sst))

in_matrix_pac_gsj <- matrix(0, 2, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  in_matrix_pac_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                    b1_gsj_hat[i] * in_pac_gsj$pac +
                                    b2_hat[i] * mean(sst))
}

pacmin_gsj <- c(in_matrix_pac_gsj[1, ]) # post distribution of calf mortality when
# GAPC was equal to its minimum observed value in Golfo San José
pacmax_gsj <- c(in_matrix_pac_gsj[2, ])# post distribution of calf mortality when
# GAPC was equal to its maximum observed value in Golfo San José

dif_pac_gsj <- pacmax_gsj - pacmin_gsj # post distribution of the difference (i.e.,
# absolute increase)

mean(dif_pac_gsj) # 0.05
quantile(dif_pac_gsj, probs = c(0.025, 00.975)) # 2.5% = -0.09, 07.5 % = 0.20
# 2.5% CI is negative because we are working with a difference. Some samples of
# pacmax_gn has lower values than some samples of pacmin_gn. Hence, the difference
# might be negative in some cases.


# Calves' probability of dying - absolute increase considering GAPM values

# This was done to later compare absolute increases of calf mortality when GAPC
# and GAPM increased in the same amount (i.e., to compare GAPC and GAPM effects
# on calf mortality).

# First, standardize GAPM values with the mean and sd of GAPC values, to get
# GAPM extremes in the scale of GAPC.
pam_in_pac_scale <- (cmd$gapm - mean(cmd$gapc)) / sd(cmd$gapc)

pred_pac_gn <- expand.grid(pac = c(min(pam_in_pac_scale[1:15]),
                                    max(pam_in_pac_scale[1:15])),
                           sst = mean(sst))

theta_matrix_pac_gn <- matrix(0, 2, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  theta_matrix_pac_gn[,i] <- invlogit(b0_gn_hat[i] +
                                        b1_gn_hat[i] * pred_pac_gn$pac +
                                        b2_hat[i] * mean(sst))
}

pacmin_gn <- c(theta_matrix_pac_gn[1, ])# post distribution of calf mortality when
# GAPC was equal to GAPM minimum observed value in Golfo Nuevo
pacmax_gn <- c(theta_matrix_pac_gn[2, ]) # post distribution of calf mortality when
# GAPC was equal to GAPM maximum observed value in Golfo Nuevo

dif_pac_gn <- pacmax_gn - pacmin_gn # posterior distribution of the difference

mean(dif_pac_gn) # 0.06
quantile(dif_pac_gn, probs = c(0.025, 00.975)) # 2.5% = -0.06, 07.5 % = 0.12
# 2.5% CI is negative because we are working with a difference. Some samples of
# pacmax_gn has lower values than some samples of pacmin_gn. Hence, the difference
# might be negative in some cases.

# Golfo San José
pred_pac_gsj <- expand.grid(pac = c(min(pam_in_pac_scale[16:31]),
                                    max(pam_in_pac_scale[16:31])),
                            sst = mean(sst))

theta_matrix_pac_gsj <- matrix(0, 2, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  theta_matrix_pac_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                         b1_gsj_hat[i] * pred_pac_gsj$pac +
                                         b2_hat[i] * mean(sst))
}

pacmin_gsj <- c(theta_matrix_pac_gsj[1, ])# post distribution of calf mortality when
# GAPC was equal to GAPM minimum observed value in Golfo San José
pacmax_gsj <- c(theta_matrix_pac_gsj[2, ])# post distribution of calf mortality when
# GAPC was equal to GAPM maximum observed value in Golfo San José

dif_pac_gsj <- pacmax_gsj - pacmin_gsj # posterior distribution of the difference

mean(dif_pac_gsj) # 0.02
quantile(dif_pac_gsj, probs = c(0.025, 00.975)) # 2.5% = -0.06, 97.5 % = 0.09
# 2.5% CI is negative because we are working with a difference. Some samples of
# pacmax_gsj has lower values than some samples of pacmin_gsj. Hence, the difference
# might be negative in some cases

# Without considering the gulf
dif_pac <- rowMeans(as.data.frame(t(rbind(dif_pac_gsj, dif_pac_gn))))

mean(dif_pac) #0.04
quantile(dif_pac, probs = c(0.025, 00.975)) # 2.5% = -0.03 97.5% = 0.09


# GAPC partial predictions ------------------------------------------------

gapc_z <- (-mean(cmd$gapc))/sd(cmd$gapc) # standardizing GAPC = 0

# Golfo Nuevo

pred_pac_gn <- expand.grid(pac = seq(gapc_z, max(pac), length.out = 200),
                           sst = mean(sst))

theta_matrix_pac_gn <- matrix(0, 200, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  theta_matrix_pac_gn[,i] <- invlogit(b0_gn_hat[i] +
                                        b1_gn_hat[i] * pred_pac_gn$pac +
                                        b2_hat[i] * pred_pac_gn$sst)
}
# Now theta_matrix_fa_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of GAPC inside the range [0, GAPC=max].
# Each row represent one posterior distribution (200 rows), and has 30,000
# samples (cols).

# Summarazing post distributions
theta_quan_pac_gn <- apply(theta_matrix_pac_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_pac_gn$mean_t <- rowMeans(theta_matrix_pac_gn)
theta_quan_pac_gn$pac <- pred_pac_gn$pac
names(theta_quan_pac_gn)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_pac_gn <- matrix(0, 200, length(phi_gn))

for (i in 1:length(phi_gn)){
  p_matrix_pac_gn[,i] <- rbeta(200, theta_matrix_pac_gn[,i] * phi_gn[i],
                               (1 - theta_matrix_pac_gn[,i]) * phi_gn[i])
}

p_quan_pac_gn <- apply(p_matrix_pac_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_pac_gn$mean_p <- rowMeans(p_matrix_pac_gn)
names(p_quan_pac_gn)[1:3] <- c("p_lower", "p_median", "p_upper") # naming
t_p_pac_gn <-cbind.data.frame(theta_quan_pac_gn, p_quan_pac_gn)

t_p_pac_gn$epac <- t_p_pac_gn$pac * sd(cmd$gapc) + mean(cmd$gapc) # destandardizing to plot
p_obs_gn <- subset(t_p_pac_gn, epac>2.87) # to work the britness of the ribbons

ggplot(t_p_pac_gn, aes(x = epac, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#C7E020FF", alpha = 0.3)+
  geom_ribbon(p_obs_gn, mapping = (aes(ymin=t_lower, ymax=t_upper)), fill = "#C7E020FF", alpha = 0.3)+
  geom_line(p_obs_gn, mapping = (aes(x = epac, y = mean_t)), colour = "black", size = 0.35)+
  geom_line(size = 0.35, alpha = 0.7, colour = "black")+
  geom_point(obs_gn, mapping = aes(x = gapc, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  scale_x_continuous(breaks = seq(0, 10, by=2.5),
                     labels = scales::number_format(accuracy = 0.01))+
  ylim(0, 1)+
  xlab("GAPC (attacks/h)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Golfo San José

pred_pac_gsj <- expand.grid(pac = seq(gapc_z, max(pac), length.out = 200),
                            sst = mean(sst))

theta_matrix_pac_gsj <- matrix(0, 200, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  theta_matrix_pac_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                         b1_gsj_hat[i] * pred_pac_gsj$pac +
                                         b2_hat[i] * pred_pac_gsj$sst)
}
# Now theta_matrix_fa_gsj has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of GAPC inside the range [0, GAPC=max].
# Each row represent one posterior distribution (200 rows), and has 30,000
# samples (cols).

theta_quan_pac_gsj <- apply(theta_matrix_pac_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_pac_gsj$mean_t <- rowMeans(theta_matrix_pac_gsj)
theta_quan_pac_gsj$pac <- pred_pac_gsj$pac # including to plot
names(theta_quan_pac_gsj)[1:3] <- c("t_lower", "t_median", "t_upper") #naming

p_matrix_pac_gsj <- matrix(0, 200, length(phi_gsj))

for (i in 1:length(phi_gsj)){
  p_matrix_pac_gsj[,i] <- rbeta(200, theta_matrix_pac_gsj[,i] * phi_gsj[i],
                                (1 - theta_matrix_pac_gsj[,i]) * phi_gsj[i])
}

p_quan_pac_gsj <- apply(p_matrix_pac_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_pac_gsj$mean_p <- rowMeans(p_matrix_pac_gsj)
names(p_quan_pac_gsj)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_pac_gsj <-cbind.data.frame(theta_quan_pac_gsj, p_quan_pac_gsj)

t_p_pac_gsj$epac <- t_p_pac_gsj$pac * sd(cmd$gapc) + mean(cmd$gapc) # destandardizing to plot
p_obs_gsj <- subset(t_p_pac_gsj, epac<5.6) # to work the brithness of the ribbons
p_obs_gsj <- subset(p_obs_gsj, epac>0.18) # to work the brithness of the ribbons

ggplot(t_p_pac_gsj, aes(x = epac, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#C7E020FF", alpha = 0.3)+
  geom_ribbon(p_obs_gsj, mapping = (aes(ymin=t_lower, ymax=t_upper)), fill = "#C7E020FF", alpha = 0.3)+
  geom_line(p_obs_gsj, mapping = (aes(x = epac, y = mean_t)), colour = "black", size = 0.35)+
  geom_line(size = 0.35, alpha = 0.7, colour = "black")+
  geom_point(obs_gsj, mapping = aes(x = gapc, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  scale_x_continuous(breaks = seq(0, 10, by=2.5),
                     labels = scales::number_format(accuracy = 0.01))+
  ylim(0, 1)+
  xlab("GAPC (attacks/h)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Probability of dying when GAPC = 0 -----------------------------------

gapc_z <- (-mean(cmd$gapc))/sd(cmd$gapc) # standardizing GAPC = 0

# Golfo Nuevo
pred_pac_gn_z <- expand.grid(pac = gapc_z, sst = mean(sst))
mort_pac_gn_z <- matrix(0, 1, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  mort_pac_gn_z[,i] <- invlogit(b0_gn_hat[i] +
                                b1_gn_hat[i] * pred_pac_gn_z$pac +
                                b2_hat[i] * pred_pac_gn_z$sst)
}
# mort_pac_gn_z is the posterior distribution of the estimated calf mortality
# when GAPC would have been equal to zero in Golfo Nuevo.

# Summarzing the post distribution
apply(mort_pac_gn_z, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
#  2.5%     50%    97.5%
# 0.045    0.13    0.33
rowMeans(mort_pac_gn_z) # 0.15

# Golfo San José
pred_pac_gsj_z <- expand.grid(pac = gapc_z, sst = mean(sst))
mort_pac_gsj_z <- matrix(0, 1, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  mort_pac_gsj_z[,i] <- invlogit(b0_gsj_hat[i] +
                                 b1_gsj_hat[i] * pred_pac_gsj_z$pac +
                                 b2_hat[i] * pred_pac_gsj_z$sst)
}
# mort_pac_gsj_z is the posterior distribution of the estimated calf mortality
# when GAPC would have been equal to zero in Golfo San José.

# Summarazing the post distribution
apply(mort_pac_gsj_z, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
# 2.5%        50%     97.5%
# 1 0.0305299 0.08744052 0.1999709
rowMeans(mort_pac_gsj_z) # 0.09443735


# Probability of dying in an average year ---------------------------------

# b0_gn_hat and b0_gsj_hat are the posterior distributions of the intercepts. As
# the covariables were standardized for model fitting, the intercepts represent
# calf mortality when gull attack (GAPC) and SST anomalies were at their means

# Golfo Nuevo
m_gn_pac <- invlogit(b0_gn_hat)
mean(m_gn_pac) # 0.21
quantile(m_gn_pac, probs = c(0.025, 0.975))
# 2.5%   97.5%
# 0.14   0.299

# Golfo San José
m_gsj_pac <- invlogit(b0_gsj_hat)
mean(m_gsj_pac) # 0.13
quantile(m_gsj_pac, probs = c(0.025, 0.975))
#  2.5%  97.5%
# 0.09   0.19


# Calf mortality in two attack scenarios - quotient (attack = average / attack = 0)

# Golfo Nuevo
in_gn_gapc <- m_gn_pac/mort_pac_gn_z # posterior distribution of the quotient
mean(in_gn_gapc) # 1.716645
quantile(in_gn_gapc, probs = c(0.025, 0.975))
# 2.5%     97.5%
#  0.7730051 3.4135291

# Golfo San José
in_gsj_gapc <- m_gsj_pac/mort_pac_gsj_z # posterior distribution of the quotient
mean(in_gsj_gapc) # 1.728279
quantile(in_gsj_gapc, probs = c(0.025, 0.975))
# 2.5%     97.5%
#  0.5661407 4.4855252


# SST anomalies - partial predictions -----------------------------------------

# Golfo Nuevo
pred_sst_gn <- expand.grid(sst = seq(min(sst), max(sst), length.out = 200),
                            pac = mean(pac))

theta_matrix_sst_gn <- matrix(0, 200, length(b0_gn_hat))
for (i in 1:length(b0_gn_hat)){
  theta_matrix_sst_gn[,i] <- invlogit(b0_gn_hat[i] +
                                      b1_gn_hat[i] * pred_sst_gn$pac +
                                      b2_hat[i] * pred_sst_gn$sst)
}
# Now theta_matrix_sst_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of SST anomalies inside the range
# [SST=min, SST=max]. Each row represent one posterior distribution (200 rows),
# and has 30,000 samples (cols).

# Summarazing posterior distributions:
theta_quan_sst_gn <- apply(theta_matrix_sst_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_sst_gn$mean_t <- rowMeans(theta_quan_sst_gn)
theta_quan_sst_gn$sst <- pred_sst_gn$sst
names(theta_quan_sst_gn)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_sst_gn <- matrix(0, 200, length(phi_gn))

for (i in 1:length(phi_gn)){
  p_matrix_sst_gn[,i] <- rbeta(200, theta_matrix_sst_gn[,i] * phi_gn[i],
                                (1 - theta_matrix_sst_gn[,i]) * phi_gn[i])
}

p_quan_sst_gn <- apply(p_matrix_sst_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_sst_gn$mean_p <- rowMeans(p_quan_sst_gn)
names(p_quan_sst_gn)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_sst_gn <-cbind.data.frame(theta_quan_sst_gn, p_quan_sst_gn)

t_p_sst_gn$esst <- t_p_sst_gn$sst * sd(cmd$sst) + mean(cmd$sst) # destandardizing to plot
obs_gn <- subset(cmd, golfo == "gn") # observed values in Golfo Nuevo

ggplot(t_p_sst_gn, aes(x = esst, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#C7E020FF", alpha = 0.45)+
  geom_line(colour = "black", size = 0.35)+
  ylim(0, 1)+
  geom_point(obs_gn, mapping = aes(x = sst, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  xlab("SST anomalies (°C)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Golfo San José
pred_sst_gsj <- expand.grid(sst = seq(min(sst), max(sst), length.out = 200),
                             pac = mean(pac))

theta_matrix_sst_gsj <- matrix(0, 200, length(b0_gsj_hat))
for (i in 1:length(b0_gsj_hat)){
  theta_matrix_sst_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                          b1_gsj_hat[i] * pred_sst_gsj$pac +
                                          b2_hat[i] * pred_sst_gsj$sst)
}
# Now theta_matrix_sst_gsj has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of SST anomalies inside the range
# [SST=min, SST=max]. Each row represent one posterior distribution (200 rows),
# and has 30,000 samples (cols).

# Summarazing posterior distributions:
theta_quan_sst_gsj <- apply(theta_matrix_sst_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_sst_gsj$mean_t <- rowMeans(theta_quan_sst_gsj)
theta_quan_sst_gsj$sst <- pred_sst_gsj$sst
names(theta_quan_sst_gsj)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_sst_gsj <- matrix(0, 200, length(phi_gsj))
for (i in 1:length(phi_gsj)){
  p_matrix_sst_gsj[,i] <- rbeta(200, theta_matrix_sst_gsj[,i] * phi_gsj[i],
                                 (1 - theta_matrix_sst_gsj[,i]) * phi_gsj[i])
}

p_quan_sst_gsj <- apply(p_matrix_sst_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
names(p_quan_sst_gsj)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_sst_gsj <-cbind.data.frame(theta_quan_sst_gsj, p_quan_sst_gsj)

t_p_sst_gsj$esst <- t_p_sst_gn$sst * sd(cmd$sst) + mean(cmd$sst) # destandardizing to plot
obs_gsj <- subset(cmd, golfo == "gsj") # observed values in Golfo San José

ggplot(t_p_sst_gsj, aes(x = esst, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#C7E020FF", alpha = 0.45)+
  geom_line(colour = "black", size = 0.35)+
  ylim(0, 1)+
  geom_point(obs_gsj, mapping = aes(x = sst, y = mort),
             inherit.aes = FALSE, shape = 1, size = 0.95)+
  xlab("SST anomalies (°C)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Posterior Predictive Check ----------------------------------------------

phi_hat <- jags.sim$sims.list$phi %>% t
theta_hat <- jags.sim$sims.list$theta %>% t
alpha_hat <- theta_hat * phi_hat
beta_hat <- (1 - theta_hat) * phi_hat

m <- matrix(NA, length(deaths), ncol(theta_hat))

npost <- ncol(theta_hat)
for (i in 1:npost) {
  m[,i] <-  rbbinom(n=length(deaths),size=borns, alpha=alpha_hat[,i], beta=beta_hat[,i])
}

plot(density(as.numeric(m))) # marginal posterior distribution
plot(density(deaths)) # to compare marginal with marginal posterior distribution

post_p <- apply(m, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
names(post_p)[1:3] <- c("lower", "median", "upper")
post_p$mean <- (theta_hat %>% rowMeans) * borns

ppc_p <- cbind(cmd,post_p)

ppc_p$golfo <- revalue(ppc_p$golfo,
                       replace = c("gn"="Golfo Nuevo",
                                   "gsj"="Golfo San José"))

ggplot(ppc_p, aes(x=year, y = mean))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(mapping = aes(x = year,
                              ymin = lower, ymax = upper), size = 0.4,  width = 0.5) +
  facet_wrap(~golfo)+
  xlab("Year") +
  ylab("Number of dead calves")+
  geom_point(ppc_p, mapping = aes(x = year, y = dead),
             inherit.aes = FALSE, shape = 1, size = 1.7)+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        legend.text=element_text(size=8),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 8,  colour = "black"),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        strip.text = element_text(size = 10),
        panel.grid.minor = element_blank())


# GAPC model - R2  -----------------------------------------------------------

theta_hat <- jags.sim$sims.list$theta %>% t
phi_hat <- jags.sim$sims.list$phi %>% t

# theta is the expectation for the response (proportion scale), so
# var_theta = var(E(y))

var_theta <- apply(theta_hat, 2, var)

# To compute var(y) we need to get alpha and beta.

alpha_hat <- theta_hat * phi_hat
beta_hat <- (1 - theta_hat) * phi_hat

var_y_samples <- (alpha_hat * beta_hat) /
  ( (alpha_hat + beta_hat) ^ 2 * (alpha_hat + beta_hat + 1) )
# var_y_samples has the predicted variance for every observation and for
# every posterior sample. For posterior sample s, we consider the residual variance
# the average of the observation-wise variances computed with sample s.
# In this case, we take the column means:

var_y <- colMeans(var_y_samples)

pac_model_R2 <- var_theta / (var_theta + var_y)
plot(density(pac_model_R2, from = 0, to = 1), xlim = c(0, 1))

qmean(pac_model_R2) # posterior distribution - summary
# mean  2.5% CI   97.5% CI
# 0.30   0.09      0.54


# GAPM model -----------------------------------------------------------------

cat(file = "rp.bug",
    "
model{
# likelihood
for (i in 1:n) {
  deaths[i] ~ dbinom(p[i], borns[i])
  p[i] ~ dbeta(alpha[i], beta[i])
  alpha[i] <- theta[i] * phi[i]
  beta[i] <- (1 - theta[i]) * phi[i]
  logit(theta[i]) <- b0_gn * gn_dummy[i] + b0_gsj * gsj_dummy[i] +
                     b1_gn * gn_dummy[i] * pam[i] + b1_gsj * gsj_dummy[i] * pam[i] +
                     b2 * sst[i]

  phi[i] <-  1 / scale_gn * gn_dummy[i] + 1 / scale_gsj * gsj_dummy[i]
}

# priors

scale_gn ~ dgamma(1, 1 / 0.1) # medida de dispersi?n
scale_gsj ~ dgamma(1, 1 / 0.1)
b0_gn ~ dnorm(0, 1 / 1.5 ^ 2) # nodo estoc?stico (dsps de virulilla)
b0_gsj ~ dnorm(0, 1 / 1.5 ^ 2)
b1_gn ~ dnorm(0, 1 / 1.5 ^ 2)
b1_gsj ~ dnorm(0, 1 / 1.5 ^ 2)
b2 ~ dnorm(0, 1 / 1.5 ^ 2)
}"
)

m.data <- list(deaths = deaths, borns = borns, pam = pam, n = n, sst = sst,
               gn_dummy = gn_dummy, gsj_dummy = gsj_dummy)

inits <- function() list(b0_gn = runif(1, -1, 1),
                         b0_gsj = runif(1, -1, 1),
                         b1_gn = runif(1, -1, 1),
                         b1_gsj = runif(1, -1, 1),
                         b2 = runif(1, -1, 1),
                         scale_gn = runif(1, 0, 12),
                         scale_gsj = runif(1, 0, 12),
                         p = runif(31, 0, 1))

params <- c("b0_gn","b0_gsj", "b1_gn","b1_gsj", "b2", "scale_gn",
            "scale_gsj", "p", "phi", "theta")

ni <- 20000 # number of iterations
nt <- 1 # thining rate
nb <- 5000 # number of iterations for burn in
nc <- 2 # number of chains

jags.sim <- jags(data = m.data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "rp.bug",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb)

print(jags.sim) # posterior distributions - summary

# Taking the posterior distributions of each parameter
b0_gn_hat <- jags.sim$sims.list$b0_gn
b0_gsj_hat <- jags.sim$sims.list$b0_gsj
b1_gn_hat <- jags.sim$sims.list$b1_gn
b1_gsj_hat <- jags.sim$sims.list$b1_gsj
b2_hat <- jags.sim$sims.list$b2
phi_gn <- 1 / jags.sim$sims.list$scale_gn
phi_gsj <- 1 / jags.sim$sims.list$scale_gsj
phi_hat <- jags.sim$sims.list$phi
theta_hat <- jags.sim$sims.list$theta

# Probability of slopes > 0 -----------------------------------------------

round(sum(b1_gn_hat > 0) / length(b1_gn_hat), 2)   # 0.98
round(sum(b1_gsj_hat > 0) / length(b1_gsj_hat), 2) # 0.97
round(sum(b2_hat < 0) / length(b2_hat), 2)         # 0.73

# Calves' probability of dying - absolute increase  ---------------------------

in_pam_gn <- expand.grid(pam = c(min(pam[1:15]), max(pam[1:15])),
                         sst = mean(sst))
in_matrix_pam_gn <- matrix(0, 2, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  in_matrix_pam_gn[,i] <- invlogit(b0_gn_hat[i] +
                                   b1_gn_hat[i] * in_pam_gn$pam +
                                   b2_hat[i] * mean(sst))
}

pammin_gn <- c(in_matrix_pam_gn[1, ]) # post distribution of calf mortality when
# GAPM was equal to its minimum observed value in Golfo Nuevo
pammax_gn <- c(in_matrix_pam_gn[2, ]) # post distribution of calf mortality when
# GAPM was equal to its maximum observed value in Golfo Nuevo

dif_pam_gn <- pammax_gn - pammin_gn # posterior distribution of the difference

mean(dif_pam_gn) # 0.32
quantile(dif_pam_gn, probs = c(0.025, 00.975)) # 2.5% = 0.03, 97.5 % = 0.59

# Golfo San José
in_pam_gsj <- expand.grid(pam = c(min(pam[16:31]), max(pam[16:31])),
                          sst = mean(sst))

in_matrix_pam_gsj <- matrix(0, 2, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  in_matrix_pam_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                    b1_gsj_hat[i] * in_pam_gsj$pam +
                                    b2_hat[i] * mean(sst))
}

pammin_gsj <- c(in_matrix_pam_gsj[1, ])# post distribution of calf mortality when
# GAPC was equal to its minimum observed value in Golfo San José
pammax_gsj <- c(in_matrix_pam_gsj[2, ])# post distribution of calf mortality when
# GAPC was equal to its maximum observed value in Golfo San José

dif_pam_gsj <- pammax_gsj - pammin_gsj # posterior distribution of the difference

mean(dif_pam_gsj) # 0.17
quantile(dif_pam_gsj, probs = c(0.025, 00.975)) # 2.5% = -0.01, 07.5 % = 0.37
# 2.5% CI is negative because we are working with a difference. Some samples of
# pammax_gn has lower values than some samples of pammin_gn. Hence, the difference
# might be negative in some cases.


# To later compare this absolute increase with GAPC...

dif_pam <- rowMeans(as.data.frame(t(rbind(dif_pam_gsj, dif_pam_gn))))

mean(dif_pam) #0.24
quantile(dif_pam, probs = c(0.025, 00.975)) # 2.5% = 0.08 97.5% = 0.41


# GAPM partial predictions ------------------------------------------------

gapm_z <- (-mean(cmd$gapm))/sd(cmd$gapm) # standardizing GAPM = 0

# Golfo Nuevo
pred_pam_gn <- expand.grid(pam = seq(gapm_z, max(pam), length.out = 200),
                           sst = mean(sst))

theta_matrix_pam_gn <- matrix(0, 200, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  theta_matrix_pam_gn[,i] <- invlogit(b0_gn_hat[i] +
                                      b1_gn_hat[i] * pred_pam_gn$pam +
                                      b2_hat[i] * pred_pam_gn$sst)
}
# Now theta_matrix_pam_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of GAPMC inside the range [0, GAPM=max].
# Each row represent one posterior distribution (200 rows), and has 30,000
# samples (cols).

# Summarazing posterior distributions:
theta_quan_pam_gn <- apply(theta_matrix_pam_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_pam_gn$mean_t <- rowMeans(theta_matrix_pam_gn)
theta_quan_pam_gn$pam <- pred_pam_gn$pam
names(theta_quan_pam_gn)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_pam_gn <- matrix(0, 200, length(phi_gn))

for (i in 1:length(phi_gn)){
  p_matrix_pam_gn[,i] <- rbeta(200, theta_matrix_pam_gn[,i] * phi_gn[i],
                               (1 - theta_matrix_pam_gn[,i]) * phi_gn[i])
}

p_quan_pam_gn <- apply(p_matrix_pam_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_pam_gn$mean_p <- rowMeans(p_matrix_pam_gn)
names(p_quan_pam_gn)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_pam_gn <-cbind.data.frame(theta_quan_pam_gn, p_quan_pam_gn)

t_p_pam_gn$epam <- t_p_pam_gn$pam * sd(cmd$gapm) + mean(cmd$gapm) # destandardizing to plot
p_obs_gn <- subset(t_p_pam_gn, epam>0.39) # to work the britness of the ribbons

ggplot(t_p_pam_gn, aes(x = epam, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#47C16EFF", alpha = 0.2)+
  geom_ribbon(p_obs_gn, mapping = (aes(ymin=t_lower, ymax=t_upper)), fill = "#47C16EFF", alpha = 0.3)+
  geom_line(p_obs_gn, mapping = (aes(x = epam, y = mean_t)), colour = "black", size = 0.35)+
  geom_line(size = 0.35, alpha = 0.7, colour = "black")+
  geom_point(obs_gn, mapping = aes(x = gapm, y = mort),
             inherit.aes = FALSE, colour = "black", shape = 1, size = 0.95)+
  xlab("GAPM (attacks/h)") +
  ylab("Calves' probability of dying")+
  scale_x_continuous(breaks = seq(0, 4, by = 1))+
  ylim(0, 1)+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Golfo San José
pred_pam_gsj <- expand.grid(pam = seq(gapm_z, max(pam), length.out = 200),
                            sst = mean(sst))

theta_matrix_pam_gsj <- matrix(0, 200, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  theta_matrix_pam_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                         b1_gsj_hat[i] * pred_pam_gsj$pam +
                                         b2_hat[i] * pred_pam_gsj$sst)
}
# Now theta_matrix_pam_gsj has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of GAPMC inside the range [0, GAPM=max].
# Each row represent one posterior distribution (200 rows), and has 30,000
# samples (cols).

# Summarazing posterior distributions:
theta_quan_pam_gsj <- apply(theta_matrix_pam_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_pam_gsj$mean_t <- rowMeans(theta_matrix_pam_gsj)
theta_quan_pam_gsj$pam <- pred_pam_gsj$pam
names(theta_quan_pam_gsj)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_pam_gsj <- matrix(0, 200, length(phi_gsj))

for (i in 1:length(phi_gsj)){
  p_matrix_pam_gsj[,i] <- rbeta(200, theta_matrix_pam_gsj[,i] * phi_gsj[i],
                                (1 - theta_matrix_pam_gsj[,i]) * phi_gsj[i])
}

p_quan_pam_gsj <- apply(p_matrix_pam_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_pam_gsj$mean_p <- rowMeans(p_matrix_pam_gsj)
names(p_quan_pam_gsj)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_pam_gsj <-cbind.data.frame(theta_quan_pam_gsj, p_quan_pam_gsj)

t_p_pam_gsj$epam <- t_p_pam_gsj$pam * sd(cmd$gapm) + mean(cmd$gapm) #destandardizing to plot
p_obs_gsj <- subset(t_p_pam_gsj, epam<3.5) # to work the britness of the ribbons
p_obs_gsj <- subset(p_obs_gsj, epam>0.32) # to work the britness of the ribbons

ggplot(t_p_pam_gsj, aes(x = epam, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#47C16EFF", alpha = 0.2)+
  geom_ribbon(p_obs_gsj, mapping = (aes(ymin=t_lower, ymax=t_upper)), fill = "#47C16EFF", alpha = 0.3)+
  geom_line(p_obs_gsj, mapping = (aes(x = epam, y = mean_t)), colour = "black", size = 0.35)+
  geom_line(size = 0.35, alpha = 0.7, colour = "black")+
  ylim(0, 1)+
  geom_point(obs_gsj, mapping = aes(x = gapm, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  xlab("GAPM (attacks/h)") +
  ylab("Calves' probability of dying")+
  scale_x_continuous(breaks = seq(0, 4, by = 1))+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Probability of dying when GAPM = 0 -----------------------------------

gapm_z <- (-mean(cmd$gapm))/sd(cmd$gapm) # standardizing value of GAPM = 0

# Golfo Nuevo
pred_pam_gn_z <- expand.grid(pam = gapm_z, sst = mean(sst))
mort_pam_gn_z <- matrix(0, 1, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  mort_pam_gn_z[,i] <- invlogit(b0_gn_hat[i] +
                                b1_gn_hat[i] * pred_pam_gn_z$pam +
                                b2_hat[i] * pred_pam_gn_z$sst)
}
# mort_pam_gn_z is the posterior distribution of the estimated calf mortality
# when GAPM would have been equal to zero in Golfo Nuevo.

# Sumarazing the post distribution
apply(mort_pam_gn_z, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
# 2.5%       50%     97.5%
# 0.05179664 0.1147066 0.2343343
rowMeans(mort_pam_gn_z) # 0.1218916
round(qmean(mort_pam_gn_z), 2)


# Golfo San José
pred_pam_gsj_z <- expand.grid(pam = gapm_z, sst = mean(sst))
mort_pam_gsj_z <- matrix(0, 1, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  mort_pam_gsj_z[,i] <- invlogit(b0_gsj_hat[i] +
                                   b1_gsj_hat[i] * pred_pam_gsj_z$pam +
                                   b2_hat[i] * pred_pam_gsj_z$sst)
}
# mort_pam_gn_z is the posterior distribution of the estimated calf mortality
# when GAPM would have been equal to zero in Golfo San José.

# Summarizing the post distribution
apply(mort_pam_gsj_z, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
# 2.5%        50%     97.5%
# 0.0386852 0.07272147 0.1367461
rowMeans(mort_pam_gsj_z) # 0.07655123
round(qmean(mort_pam_gsj_z), 2)

# Comparison between gulfs
# Pr(gn > gsj) =
round(sum((mort_pam_gn_z - mort_pam_gsj_z) > 0) / length(mort_pam_gn_z), 2)

# Probability of dying in an average year ---------------------------------

# b0_gn_hat and b0_gsj_hat are the posterior distributions of the intercepts. As
# the covariables were standardized for model fitting, the intercepts represent
# calf mortality when gull attack (GAPC) and SST anomalies were at their means

m_gn_pam <- invlogit(b0_gn_hat)
mean(m_gn_pam) # 0.2075782
quantile(m_gn_pam, probs = c(0.025, 0.975))
#      2.5%     97.5%
# 0.1458938 0.2814253
round(qmean(m_gn_pam), 2)


m_gsj_pam <- invlogit(b0_gsj_hat)
mean(m_gsj_pam) # 0.1266908
quantile(m_gsj_pam, probs = c(0.025, 0.975))
#      2.5%      97.5%
# 0.09016484 0.17519113
round(qmean(m_gsj_pam), 2)


# Calf mortality in two attack scenarios - quotient (attack = average / attack = 0)

# Golfo Nuevo
in_gn_gapm <- m_gn_pam/mort_pam_gn_z # posterior distribution of the quotient
mean(in_gn_gapm) # 1.87
quantile(in_gn_gapm, probs = c(0.025, 0.975))
#    2.5%    97.5%
# 1.053900 3.191544
round(qmean(in_gn_gapm), 2)

# Golfo San José
in_gsj_gapm <- m_gsj_pam/mort_pam_gsj_z # posterior distribution of the quotient
mean(in_gsj_gapm) # 1.79
quantile(in_gsj_gapm, probs = c(0.025, 0.975))
#     2.5%    97.5%
# 0.958834 2.984400
round(qmean(in_gsj_gapm), 2)

## Average the quotients across gulfs
in_pam_avg <- colMeans(rbind(in_gsj_gapm, in_gn_gapm))
round(qmean(in_pam_avg), 2)
# mean lower upper
# 1.84  1.22  2.70
sum(in_pam_avg > 0) / length(in_pam_avg)
# 1

# SST - partial predictions ----------------------------------------------

# Golfo Nuevo
pred_sst_gn <- expand.grid(sst = seq(min(sst), max(sst), length.out = 200),
                            pam = mean(pam))

theta_matrix_sst_gn <- matrix(0, 200, length(b0_gn_hat))
for (i in 1:length(b0_gn_hat)){
  theta_matrix_sst_gn[,i] <- invlogit(b0_gn_hat[i] +
                                         b1_gn_hat[i] * pred_sst_gn$pam +
                                         b2_hat[i] * pred_sst_gn$sst)
}
# Now theta_matrix_sst_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of SST anomalies inside the range
# [SST=min, SST=max]. Each row represent one posterior distribution (200 rows),
# and has 30,000 samples (cols).

# Summarazing posterior distributions:
theta_quan_sst_gn <- apply(theta_matrix_sst_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_sst_gn$mean_t <- rowMeans(theta_quan_sst_gn)
theta_quan_sst_gn$sst <- pred_sst_gn$sst
names(theta_quan_sst_gn)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_sst_gn <- matrix(0, 200, length(phi_gn))
for (i in 1:length(phi_gn)){
  p_matrix_sst_gn[,i] <- rbeta(200, theta_matrix_sst_gn[,i] * phi_gn[i],
                                (1 - theta_matrix_sst_gn[,i]) * phi_gn[i])
}

p_quan_sst_gn <- apply(p_matrix_sst_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_sst_gn$mean_p <- rowMeans(p_quan_sst_gn)
names(p_quan_sst_gn)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_sst_gn <-cbind.data.frame(theta_quan_sst_gn, p_quan_sst_gn)

t_p_sst_gn$esst <- t_p_sst_gn$sst * sd(cmd$sst) + mean(cmd$sst) # destandardizing to plot

ggplot(t_p_sst_gn, aes(x = esst, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#47C16EFF", alpha = 0.45)+
  geom_line(colour = "black", size = 0.35)+
  ylim(0, 1)+
  geom_point(obs_gn, mapping = aes(x = sst, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  xlab("SST anomalies (°C)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Golfo San José
pred_sst_gsj <- expand.grid(sst = seq(min(sst), max(sst), length.out = 200),
                             pam = mean(pam))

theta_matrix_sst_gsj <- matrix(0, 200, length(b0_gsj_hat))
for (i in 1:length(b0_gsj_hat)){
  theta_matrix_sst_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                          b1_gsj_hat[i] * pred_sst_gsj$pam +
                                          b2_hat[i] * pred_sst_gsj$sst)
}

# Now theta_matrix_sst_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of SST anomalies inside the range
# [SST=min, SST=max]. Each row represent one posterior distribution (200 rows),
# and has 30,000 samples (cols).

# Summarazing posterior distributions:
theta_quan_sst_gsj <- apply(theta_matrix_sst_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_sst_gsj$mean_t <- rowMeans(theta_quan_sst_gsj)
theta_quan_sst_gsj$sst <- pred_sst_gsj$sst
names(theta_quan_sst_gsj)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_sst_gsj <- matrix(0, 200, length(phi_gsj))
for (i in 1:length(phi_gsj)){
  p_matrix_sst_gsj[,i] <- rbeta(200, theta_matrix_sst_gsj[,i] * phi_gsj[i],
                                 (1 - theta_matrix_sst_gsj[,i]) * phi_gsj[i])
}

p_quan_sst_gsj <- apply(p_matrix_sst_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
names(p_quan_sst_gsj)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_sst_gsj <-cbind.data.frame(theta_quan_sst_gsj, p_quan_sst_gsj)

t_p_sst_gsj$esst <- t_p_sst_gn$sst * sd(cmd$sst) + mean(cmd$sst) # destandardizing to plot

ggplot(t_p_sst_gsj, aes(x = esst, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#47C16EFF", alpha = 0.45)+
  geom_line(colour = "black", size = 0.35)+
  ylim(0, 1)+
  #scale_x_discrete(breaks = c(0.20, 0.3, 0.4))+
  geom_point(obs_gsj, mapping = aes(x = sst, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  xlab("SST anomalies (°C)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Posterior Predictive Check ----------------------------------------------

phi_hat <- jags.sim$sims.list$phi %>% t
theta_hat <- jags.sim$sims.list$theta %>% t
alpha_hat <- theta_hat * phi_hat
beta_hat <- (1 - theta_hat) * phi_hat

m <- matrix(NA, length(deaths), ncol(theta_hat))

npost <- ncol(theta_hat)
for (i in 1:npost) {
  m[,i] <-  rbbinom(n=length(deaths),size=borns, alpha=alpha_hat[,i], beta=beta_hat[,i])
}

plot(density(as.numeric(m))) # marginal posterior distribution
plot(density(deaths)) # to compare with marginal posterior distribution

post_p <- apply(m, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
names(post_p)[1:3] <- c("lower", "median", "upper")
post_p$mean <- (theta_hat %>% rowMeans) * borns

ppc_p <- cbind(cmd,post_p)

ppc_p$golfo <- revalue(ppc_p$golfo,
                       replace = c("gn"="Golfo Nuevo",
                                   "gsj"="Golfo San José"))

ggplot(ppc_p, aes(x=year, y = mean))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(mapping = aes(x = year,
                              ymin = lower, ymax = upper), size = 0.4,  width = 0.5) +
  facet_wrap(~golfo)+
  xlab("Year") +
  ylab("Number of dead calves")+
  #scale_x_discrete(breaks = c(1995, 2004,2007,2010,2013,2016,2019,2020))+
  geom_point(ppc_p, mapping = aes(x = year, y = dead),
             inherit.aes = FALSE, shape = 1, size = 1.7)+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        legend.text=element_text(size=8),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 8,  colour = "black"),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        strip.text = element_text(size = 10),
        panel.grid.minor = element_blank())


# GAPM model - R2 --------------------------------------------------------------

theta_hat <- jags.sim$sims.list$theta %>% t
phi_hat <- jags.sim$sims.list$phi %>% t

# theta is the expectation for the response (proportion scale), so
# var_theta = var(E(y))

var_theta <- apply(theta_hat, 2, var)

# To compute var(y) we need to get alpha and beta.

alpha_hat <- theta_hat * phi_hat
beta_hat <- (1 - theta_hat) * phi_hat

var_y_samples <- (alpha_hat * beta_hat) /
  ( (alpha_hat + beta_hat) ^ 2 * (alpha_hat + beta_hat + 1) )
# var_y_samples has the predicted variance for every observation and for
# every posterior sample. For posterior sample s, we consider the residual variance
# the average of the observation-wise variances computed with sample s.
# In this case, we take the column means:

var_y <- colMeans(var_y_samples)

pam_model_R2 <- var_theta / (var_theta + var_y)
plot(density(pam_model_R2, from = 0, to = 1), xlim = c(0, 1))

qmean(pam_model_R2)
# mean   2.5% CI   97.5% CI
# 0.44    0.18     0.67


# GAF model---------------------------------------------------------------------

cat(file = "rp.bug",
    "
model{
# likelihood
for (i in 1:n) {
  deaths[i] ~ dbinom(p[i], borns[i])
  p[i] ~ dbeta(alpha[i], beta[i])
  alpha[i] <- theta[i] * phi[i]
  beta[i] <- (1 - theta[i]) * phi[i]
  logit(theta[i]) <- b0_gn * gn_dummy[i] + b0_gsj * gsj_dummy[i] +
                     b1_gn * gn_dummy[i] * fa[i] + b1_gsj * gsj_dummy[i] * fa[i] +
                     b2 * sst[i]

  phi[i] <-  1 / scale_gn * gn_dummy[i] + 1 / scale_gsj * gsj_dummy[i]
}


# priors

scale_gn ~ dgamma(1, 1 / 0.1) # medida de dispersi?n
scale_gsj ~ dgamma(1, 1 / 0.1)
b0_gn ~ dnorm(0, 1 / 1.5 ^ 2) # nodo estoc?stico (dsps de virulilla)
b0_gsj ~ dnorm(0, 1 / 1.5 ^ 2)
b1_gn ~ dnorm(0, 1 / 1.5 ^ 2)
b1_gsj ~ dnorm(0, 1 / 1.5 ^ 2)
b2 ~ dnorm(0, 1 / 1.5 ^ 2)
}"
)

m.data <- list(deaths = deaths, borns = borns, fa = fa, n = n, sst = sst,
               gn_dummy = gn_dummy, gsj_dummy = gsj_dummy)

inits <- function() list(b0_gn = runif(1, -1, 1),
                         b0_gsj = runif(1, -1, 1),
                         b1_gn = runif(1, -1, 1),
                         b1_gsj = runif(1, -1, 1),
                         b2 = runif(1, -1, 1),
                         scale_gn = runif(1, 0, 12),
                         scale_gsj = runif(1, 0, 12),
                         p = runif(31, 0, 1))

params <- c("b0_gn","b0_gsj", "b1_gn","b1_gsj", "b2", "scale_gn", "scale_gsj",
            "p", "theta", "phi")

ni <- 20000 # number of iterations
nt <- 1 # thining rate
nb <- 5000 # number of iterations for "burn in"
nc <- 2 # number of chains

jags.sim <- jags(data = m.data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "rp.bug",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb)

print(jags.sim) # Posterior distributions - summary

# Taking the posterior distributions of each parameter
b0_gn_hat <- jags.sim$sims.list$b0_gn
b0_gsj_hat <- jags.sim$sims.list$b0_gsj
b1_gn_hat <- jags.sim$sims.list$b1_gn
b1_gsj_hat <- jags.sim$sims.list$b1_gsj
b2_hat <- jags.sim$sims.list$b2
phi_gn <- 1 / jags.sim$sims.list$scale_gn
phi_gsj <- 1 / jags.sim$sims.list$scale_gsj
phi_hat <- jags.sim$sims.list$phi
theta_hat <- jags.sim$sims.list$theta

# Probability of slopes > 0 -----------------------------------------------

round(sum(b1_gn_hat > 0) / length(b1_gn_hat), 2)   # 0.99
round(sum(b1_gsj_hat > 0) / length(b1_gsj_hat), 2) # 0.92
round(sum(b2_hat < 0) / length(b2_hat), 2)         # 0.75

# Calves' probability of dying - absolute increase ----------------------------

# Golfo Nuevo
in_fa_gn <- expand.grid(fa = c(min(fa[1:15]), max(fa[1:15])),
                           sst = mean(sst))

in_matrix_fa_gn <- matrix(0, 2, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  in_matrix_fa_gn[,i] <- invlogit(b0_gn_hat[i] +
                                  b1_gn_hat[i] * in_fa_gn$fa +
                                  b2_hat[i] * mean(sst))
}

famin_gn <- c(in_matrix_fa_gn[1, ]) # post distribution of calf mortality when
# GAF is equal to its minimum observed value in Golfo Nuevo
famax_gn <- c(in_matrix_fa_gn[2, ]) # post distribution of calf mortality when
# GAF is equal to its maximum observed value in Golfo Nuevo

dif_fa_gn <- famax_gn - famin_gn # post distribution of the difference

mean(dif_fa_gn) # 0.3
quantile(dif_fa_gn, probs = c(0.025, 00.975)) # 2.5% = 0.03, 97.5 % = 0.57

# Golfo Nuevo
in_fa_gsj <- expand.grid(fa = c(min(fa[16:31]), max(fa[16:31])),
                            sst = mean(sst))

in_matrix_fa_gsj <- matrix(0, 2, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  in_matrix_fa_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                   b1_gsj_hat[i] * in_fa_gsj$fa +
                                   b2_hat[i] * mean(sst))
}

famin_gsj <- c(in_matrix_fa_gsj[1, ])# post distribution of calf mortality when
# GAF is equal to its minimum observed value in Golfo San José
famax_gsj <- c(in_matrix_fa_gsj[2, ])# post distribution of calf mortality when
# GAF is equal to its maximum observed value in Golfo San José

dif_fa_gsj <- famax_gsj - famin_gsj # post distribution of the difference

mean(dif_fa_gsj) # 0.1
quantile(dif_fa_gsj, probs = c(0.025, 00.975)) # 2.5% = -0.04, 07.5 % = 0.26
# 2.5% CI is negative because we are working with a difference. Some samples of
# famax_gn has lower values than some samples of famin_gn. Hence, the difference
# might be negative in some cases.



# GAF - Partial predictions -----------------------------------------------

# Golfo Nuevo

gaf_z <- (-mean(cmd$gaf))/sd(cmd$gaf) # standardizing GAF = 0

pred_fa_gn <- expand.grid(fa = seq(gaf_z, max(fa), length.out = 200),
                          sst = mean(sst))

theta_matrix_fa_gn <- matrix(0, 200, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
  theta_matrix_fa_gn[,i] <- invlogit(b0_gn_hat[i] +
                                     b1_gn_hat[i] * pred_fa_gn$fa +
                                     b2_hat[i] * pred_fa_gn$sst)
}
# Now theta_matrix_pam_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of GAPMC inside the range [0, GAPM=max].
# Each row represent one posterior distribution (200 rows), and has 30,000
# samples (cols).

# Summarazing posterior distributions:
theta_quan_fa_gn <- apply(theta_matrix_fa_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_fa_gn$mean_t <- rowMeans(theta_matrix_fa_gn)
theta_quan_fa_gn$fa <- pred_fa_gn$fa
names(theta_quan_fa_gn)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_fa_gn <- matrix(0, 200, length(phi_gn))

for (i in 1:length(phi_gn)){
  p_matrix_fa_gn[,i] <- rbeta(200, theta_matrix_fa_gn[,i] * phi_gn[i],
                             (1 - theta_matrix_fa_gn[,i]) * phi_gn[i])
}
# Now theta_matrix_pam_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of GAPMC inside the range [0, GAPM=max].
# Each row represent one posterior distribution (200 rows), and has 30,000
# samples (cols).

# Summarazing post distribution...
p_quan_fa_gn <- apply(p_matrix_fa_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_fa_gn$mean_p_gn <- rowMeans(p_matrix_fa_gn)
names(p_quan_fa_gn)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_fa_gn <-cbind.data.frame(theta_quan_fa_gn, p_quan_fa_gn)

t_p_fa_gn$efa <- t_p_fa_gn$fa * sd(cmd$gaf) + mean(cmd$gaf) # de-standardize values to plot
obs_gn <- subset(cmd, golfo == "gn") # observed values for Golfo Nuevo
p_obs_gn <- subset(t_p_fa_gn, efa>0.156) # to change the britness of the ribbons

ggplot(t_p_fa_gn, aes(x = efa, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#287C8EFF", alpha = 0.2)+
  geom_ribbon(p_obs_gn, mapping = (aes(ymin=t_lower, ymax=t_upper)), fill = "#287C8EFF", alpha = 0.3)+
  geom_line(p_obs_gn, mapping = (aes(x = efa, y = mean_t)), colour = "black", size = 0.35)+
  geom_line(size = 0.35, alpha = 0.7, colour = "black")+
  ylim(0, 1)+
  geom_point(obs_gn, mapping = aes(x = gaf, y = mort),
             inherit.aes = FALSE, colour = "black", shape= 1, size = 0.95)+
  xlab("GAF") +
  ylab("Calves' probability of dying")+
  scale_x_continuous(breaks = seq(0, 0.40, by=0.10))+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Golfo San José

pred_fa_gsj <- expand.grid(fa = seq(gaf_z, max(fa), length.out = 200),
                           sst = mean(sst))

theta_matrix_fa_gsj <- matrix(0, 200, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
  theta_matrix_fa_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                      b1_gsj_hat[i] * pred_fa_gsj$fa +
                                      b2_hat[i] * pred_fa_gsj$sst)
}

# Summarazing posterior distributions:
theta_summ_fa_gsj <- apply(theta_matrix_fa_gsj, 1, summary) %>% t
theta_quan_fa_gsj <- apply(theta_matrix_fa_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_fa_gsj$mean_t <- rowMeans(theta_matrix_fa_gsj)
theta_quan_fa_gsj$fa <- pred_fa_gsj$fa
names(theta_quan_fa_gsj)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_fa_gsj <- matrix(0, 200, length(phi_gsj))

for (i in 1:length(phi_gsj)){
  p_matrix_fa_gsj[,i] <- rbeta(200, theta_matrix_fa_gsj[,i] * phi_gsj[i],
                               (1 - theta_matrix_fa_gsj[,i]) * phi_gsj[i])
}

# Summarazing post distribution:
p_quan_fa_gsj <- apply(p_matrix_fa_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_fa_gsj$mean_p <- rowMeans(p_matrix_fa_gsj)
names(p_quan_fa_gsj)[1:3] <- c("p_lower", "p_median", "p_upper")
t_p_fa_gsj <-cbind.data.frame(theta_quan_fa_gsj, p_quan_fa_gsj)

t_p_fa_gsj$efa <- t_p_fa_gsj$fa * sd(cmd$gaf) + mean(cmd$gaf) # de-standardize values to plot
obs_gsj <- subset(cmd, golfo == "gsj") # observed values for Golgo San José
p_obs_gsj <- subset(t_p_fa_gsj, efa<0.275) # to change the britness of the ribbons
p_obs_gsj <- subset(p_obs_gsj, efa>0.103) # to change the britness of the ribbons

ggplot(t_p_fa_gsj, aes(x = efa, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#287C8EFF", alpha = 0.2)+
  geom_ribbon(p_obs_gsj, mapping = (aes(ymin=t_lower, ymax=t_upper)), fill = "#287C8EFF", alpha = 0.3)+
  geom_line(p_obs_gsj, mapping = (aes(x = efa, y = mean_t)), colour = "black", size = 0.35)+
  geom_line(size = 0.35, alpha = 0.7, colour = "black")+
  geom_point(obs_gsj, mapping = aes(x = gaf, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  xlab("GAF") +
  ylab("Calves' probability of dying")+
  scale_x_continuous(breaks = seq(0, 0.40, by=0.10))+
  ylim(0, 1)+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Probability of dying when GAF = 0 -----------------------------------

gaf_z <- (-mean(cmd$gaf))/sd(cmd$gaf) # standardizing GAF = 0

# Golfo Nuevo
pred_fa_gn_z <- expand.grid(fa = gaf_z, sst = mean(sst))
mort_fa_gn_z <- matrix(0, 1, length(b0_gn_hat))

for (i in 1:length(b0_gn_hat)){
     mort_fa_gn_z[,i] <- invlogit(b0_gn_hat[i] +
                                          b1_gn_hat[i] * pred_fa_gn_z$fa +
                                          b2_hat[i] * pred_fa_gn_z$sst)
}
# mort_fa_gn_z has the samples of the posterior distribution of calf mortality
# when GAF would have been equal to zero.

# Summarazing the post distribution
apply(mort_fa_gn_z, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
#       2.5%        50%     97.5%
# 0.0115986 0.05206731 0.2160542

rowMeans(mort_fa_gn_z) # 0.06733213


# Golfo San José
pred_fa_gsj_z <- expand.grid(fa = gaf_z, sst = mean(sst))
mort_fa_gsj_z <- matrix(0, 1, length(b0_gsj_hat))

for (i in 1:length(b0_gsj_hat)){
     mort_fa_gsj_z[,i] <- invlogit(b0_gsj_hat[i] +
                                   b1_gsj_hat[i] * pred_fa_gsj_z$fa +
                                   b2_hat[i] * pred_fa_gsj_z$sst)
}
# mort_fa_gsj_z has the samples of the posterior distribution of calf mortality
# when GAF would have been equal to zero.

# Summarazing the post distribution
apply(mort_fa_gsj_z, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
#  2.5%        50%     97.5%
# 0.00476748 0.03839686 0.1866785
rowMeans(mort_fa_gsj_z) #  0.05294464


# Probability of dying in an average year (of SST and attacks) -----------------

# Golfo Nuevo
m_gn_fa <- invlogit(b0_gn_hat)
# b0_gn_hat and b0_gsj_hat are the posterior distribution of the intercepts. As
# covariables were standardized for model fitting, the intercepts represet
# calf mortality when gull attack (GAF) and SST were at their means

mean(m_gn_fa) # 0.1998302
quantile(m_gn_fa, probs = c(0.025, 0.975))
#      2.5%     97.5%
# 0.1373664 0.2779778

# Golfo San José
m_gsj_fa <- invlogit(b0_gsj_hat)
mean(m_gsj_fa) # 0.1330073
quantile(m_gsj_fa, probs = c(0.025, 0.975))
#       2.5%      97.5%
# 0.09072718 0.18829009


# Calf mortality in two attack scenarios - quotient (attack = average / attack = 0)

# Golfo Nuevo
in_gn_gaf <- m_gn_fa/mort_fa_gn_z # posterior distribution of the quotient
mean(in_gn_gaf) # 4.73
quantile(in_gn_gaf, probs = c(0.025, 0.975))
#      2.5%     97.5%
#  1.099072 13.869955

# Golfo San José
in_gsj_gaf <- m_gsj_fa/mort_fa_gsj_z # posterior distribution of the quotient
mean(in_gsj_gaf) # 6.22
quantile(in_gsj_gaf, probs = c(0.025, 0.975))
#  2.5%   97.5%
# 0.63   28.57


# SST - Partial predictions -------------------------------------------------

# Golfo Nuevo

pred_sst_gn <- expand.grid(sst = seq(min(sst), max(sst), length.out = 200),
                          fa = mean(fa))

theta_matrix_sst_gn <- matrix(0, 200, length(b0_gn_hat))
for (i in 1:length(b0_gn_hat)){
  theta_matrix_sst_gn[,i] <- invlogit(b0_gn_hat[i] +
                                      b1_gn_hat[i] * pred_sst_gn$fa +
                                      b2_hat[i] * pred_sst_gn$sst)
}
# Now theta_matrix_sst_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of SST anomalies inside the range
# [SST=min, SST=max]. Each row represent one posterior distribution (200 rows),
# and has 30,000 samples (cols).

# Summarazing posterior distributions:
theta_quan_sst_gn <- apply(theta_matrix_sst_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_sst_gn$mean_t <- rowMeans(theta_quan_sst_gn)
theta_quan_sst_gn$sst <- pred_sst_gn$sst
names(theta_quan_sst_gn)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_sst_gn <- matrix(0, 200, length(phi_gn))
for (i in 1:length(phi_gn)){
  p_matrix_sst_gn[,i] <- rbeta(200, theta_matrix_sst_gn[,i] * phi_gn[i],
                              (1 - theta_matrix_sst_gn[,i]) * phi_gn[i])
}

p_quan_sst_gn <- apply(p_matrix_sst_gn, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
p_quan_sst_gn$mean_p <- rowMeans(p_quan_sst_gn)
names(p_quan_sst_gn)[1:3] <- c("p_lower", "p_median", "p_upper")

t_p_sst_gn <-cbind.data.frame(theta_quan_sst_gn, p_quan_sst_gn)
t_p_sst_gn$esst <- t_p_sst_gn$sst * sd(cmd$sst) + mean(cmd$sst)
obs_gn <- subset(cmd, golfo == "gn")

ggplot(t_p_sst_gn, aes(x = esst, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#287C8EFF", alpha = 0.45)+
  geom_line(colour = "black", size = 0.35)+
  ylim(0, 1)+
  geom_point(obs_gn, mapping = aes(x = sst, y = mort),
             inherit.aes = FALSE, shape = 1, size = 0.95)+
  xlab("SST anomalies (°C)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Golfo San Jose

pred_sst_gsj <- expand.grid(sst = seq(min(sst), max(sst), length.out = 200),
                            fa = mean(fa))

theta_matrix_sst_gsj <- matrix(0, 200, length(b0_gsj_hat))
for (i in 1:length(b0_gsj_hat)){
  theta_matrix_sst_gsj[,i] <- invlogit(b0_gsj_hat[i] +
                                       b1_gsj_hat[i] * pred_sst_gsj$fa +
                                       b2_hat[i] * pred_sst_gsj$sst)
}
# Now theta_matrix_sst_gn has the posterior distributions of calf mortality in an
# average year (that's theta) for 200 values of SST anomalies inside the range
# [SST=min, SST=max]. Each row represent one posterior distribution (200 rows),
# and has 30,000 samples (cols).

# Summarazing posterior distributions:
theta_quan_sst_gsj <- apply(theta_matrix_sst_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
theta_quan_sst_gsj$mean_t <- rowMeans(theta_quan_sst_gsj)
theta_quan_sst_gsj$sst <- pred_sst_gsj$sst
names(theta_quan_sst_gsj)[1:3] <- c("t_lower", "t_median", "t_upper")

p_matrix_sst_gsj <- matrix(0, 200, length(phi_gsj))
for (i in 1:length(phi_gsj)){
  p_matrix_sst_gsj[,i] <- rbeta(200, theta_matrix_sst_gsj[,i] * phi_gsj[i],
                                 (1 - theta_matrix_sst_gsj[,i]) * phi_gsj[i])
}

p_quan_sst_gsj <- apply(p_matrix_sst_gsj, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
names(p_quan_sst_gsj)[1:3] <- c("p_lower", "p_median", "p_upper")

t_p_sst_gsj <-cbind.data.frame(theta_quan_sst_gsj, p_quan_sst_gsj)
t_p_sst_gsj$esst <- t_p_sst_gsj$sst * sd(cmd$sst) + mean(cmd$sst)
obs_gsj <- subset(cmd, golfo == "gsj")

anom_fa_gsj<-ggplot(t_p_sst_gsj, aes(x = esst, y = mean_t))+
  geom_ribbon(aes(ymin=p_lower, ymax=p_upper),color = "grey15", fill = NA, size = 0.4, linetype = 3)+
  geom_ribbon(aes(ymin=t_lower, ymax=t_upper), fill = "#287C8EFF", alpha = 0.45)+
  geom_line(colour = "black", size = 0.35)+
  ylim(0, 1)+
  #scale_x_discrete(breaks = c(0.20, 0.3, 0.4))+
  geom_point(obs_gsj, mapping = aes(x = sst, y = mort),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 0.95)+
  xlab("SST anomalies (°C)") +
  ylab("Calves' probability of dying")+
  theme(axis.title.x = element_text(size = 12,  colour = "black"),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12,  colour = "black"),
        axis.text.y = element_text(size = 11,  colour = "black"),
        panel.grid.minor = element_blank())


# Posterior Predictive Checks ---------------------------------------------

theta_hat <- jags.sim$sims.list$theta %>% t #transpongo para que me queden las
#muestras de la posterior como columnas y las filas como las filas de mi dataframe
phi_hat <- jags.sim$sims.list$phi %>% t
alpha_hat <- theta_hat * phi_hat
beta_hat <- (1 - theta_hat) * phi_hat

m <- matrix(NA, length(deaths), ncol(theta_hat))

# theta, alpha y beta son cantidades derivadas, no parámetros per se, entonces
# lo que me devuelve JAGS son matrices cuyas columnas tienen muestras de la pos-
# terior y cuyas filas son las filas de la base de datos (que varian por año y
# por golfo).
#  Alpha, beta, phi son cantidades derivadas. Parten de theta que tmb es una can-
# tidad derivada (sin ~, con <-), ésta se obtiene de los datos y de las poste-
# riores de los parámetros. Tengo tantos thetas como muestras de los parámetros
# (fila) y combinacón de datos (col). Alpha y beta se calculan con theta, así que
# surgen de la misma lógica

npost <- ncol(theta_hat)
for (i in 1:npost) {
  m[,i] <-  rbbinom(n=length(deaths),size=borns, alpha=alpha_hat[,i], beta=beta_hat[,i])
}

plot(density(as.numeric(m))) # marginal posterior distribution
plot(density(deaths)) # to compare with marginal posterior distribution

post_p <- apply(m, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
names(post_p)[1:3] <- c("lower", "median", "upper")
post_p$mean <- (theta_hat %>% rowMeans) * borns

ppc_p <- cbind(cmd,post_p)

ppc_p$golfo <- revalue(ppc_p$golfo,
                       replace = c("gn"="Golfo Nuevo",
                                   "gsj"="Golfo San José"))

ggplot(ppc_p, aes(x=year, y = mean))+
  geom_point(col="black", size = 1.5)+
  geom_errorbar(mapping = aes(x = year,
                              ymin = lower, ymax = upper), size = 0.4,  width = 0.5) +
  facet_wrap(~golfo)+
  xlab("Year") +
  ylab("Number dead calves")+
  geom_point(ppc_p, mapping = aes(x = year, y = dead),
             inherit.aes = FALSE, size = 1.7, shape = 1)+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        legend.text=element_text(size=8),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 8,  colour = "black"),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 9, colour = "black"),
        strip.text = element_text(size = 10),
        panel.grid.minor = element_blank())


# GAF model - R2  -----------------------------------------------------------

theta_hat <- jags.sim$sims.list$theta %>% t
phi_hat <- jags.sim$sims.list$phi %>% t

# theta is the expectation for the response (proportion scale), so
# var_theta = var(E(y))

var_theta <- apply(theta_hat, 2, var)

# To compute var(y) we need to get alpha and beta.

alpha_hat <- theta_hat * phi_hat
beta_hat <- (1 - theta_hat) * phi_hat

var_y_samples <- (alpha_hat * beta_hat) /
  ( (alpha_hat + beta_hat) ^ 2 * (alpha_hat + beta_hat + 1) )
# var_y_samples has the predicted variance for every observation and for
# every posterior sample. For posterior sample s, we consider the residual variance
# the average of the observation-wise variances computed with sample s.
# In this case, we take the column means:

var_y <- colMeans(var_y_samples)

fa_model_R2 <- var_theta / (var_theta + var_y)
plot(density(fa_model_R2, from = 0, to = 1), xlim = c(0, 1))

qmean(fa_model_R2)
# mean   2.5% CI   97.5% CI
# 0.40    0.15      0.64


# Predictions averaging the three models ----------------------------------


# Probability of deaths when attacks = 0 (averaging models) ---------------

# Golfo Nuevo
mort_pac_gn_z <- t(mort_pac_gn_z)
mort_pam_gn_z <- t(mort_pam_gn_z)
mort_fa_gn_z <- t(mort_fa_gn_z)
mort_gn_z <- cbind(mort_pac_gn_z,mort_pam_gn_z,mort_fa_gn_z)
avg_gn_z <- rowMeans(mort_gn_z) # posterior distribution averaging models

round(qmean(avg_gn_z), 2)
# mean lower upper
# 0.11  0.06  0.20

# Golfo San José
mort_pac_gsj_z <- t(mort_pac_gsj_z)
mort_pam_gsj_z <- t(mort_pam_gsj_z)
mort_fa_gsj_z <- t(mort_fa_gsj_z)
mort_gsj_z <- cbind(mort_pac_gsj_z, mort_pam_gsj_z, mort_fa_gsj_z)
avg_gsj_z <- rowMeans(mort_gsj_z) # posterior distribution averaging models

round(qmean(avg_gsj_z), 2)
# mean lower upper
# 0.07  0.04  0.13

## First quotient, then average
q_gulfs_dist <- cbind(
  mort_pac_gn_z / mort_pac_gsj_z,
  mort_pam_gn_z / mort_pam_gsj_z,
  mort_fa_gn_z / mort_fa_gsj_z
) %>% rowMeans

round(qmean(q_gulfs_dist), 2)
# mean lower upper
# 2.21  0.76  6.46
round(sum(q_gulfs_dist > 1) / length(q_gulfs_dist), 2)
# 0.9

# Probability of death in an average year (averaging models) ---------------

# Golfo Nuevo
mort_gn <- cbind(m_gn_fa, m_gn_pam, m_gn_pac)
avg_samples_gn <- rowMeans(mort_gn) # posterior distribution averaging models
round(qmean(avg_samples_gn), 2)
# mean lower upper
# 0.21  0.17  0.25

# Golfo San José
mort_gsj <- cbind(m_gsj_fa, m_gsj_pam, m_gsj_pac)
avg_samples_gsj <- rowMeans(mort_gsj) # posterior distribution averaging models
round(qmean(avg_samples_gsj), 2)
# mean lower upper
# 0.13  0.11  0.16

## First quotient, then average
q_gulfs_dist <- rowMeans(mort_gn / mort_gsj)
round(qmean(q_gulfs_dist), 2)
# mean lower upper
# 1.64  1.21  2.18
round(sum(q_gulfs_dist > 1) / length(q_gulfs_dist), 2)
# 1

# Calf mortality in two attack scenarios - quotient (attack = average / attack = 0) (averaging models) ---------------

# First compute quotient (i_), then average

# Without considering the gulf
imort <- rowMeans(t(rbind(in_gn_gapc, in_gn_gapm, in_gn_gaf,
                          in_gsj_gapc, in_gsj_gapc, in_gsj_gapc)))
round(qmean(imort), 2)
# mean lower upper
# 2.26  1.22  4.35
round(sum(imort > 1) / length(imort), 2)
# Pr(q > 1) = 1


## Change the order of average: first average mortality across models, then
## compute quotient. The previous way is preferred.

mort_attack <- rowMeans(cbind(mort_gsj, mort_gn))
mort_z <- rowMeans(cbind(mort_gsj_z, mort_gn_z))
qmort2 <- mort_attack / mort_z
round(qmean(qmort2), 2)
# mean lower upper
# 1.88  1.23  2.71
round(sum(qmort2 > 1) / length(qmort2), 2)
# 1

