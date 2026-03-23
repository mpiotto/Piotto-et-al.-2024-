# Plots for the 3 ordinal regressions (3 attack variables.)

# Library and data --------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayestestR)
library(brms)
library(ggrepel)  # separate text label in plot
library(ggnewscale) # more than one scale for aesthetics
library(DHARMa)

# Gull attack data belong to the Southern Right Whale Research Program (Instituto
# de Conservación de Ballenas, Argentina). Calf mortality data belong to the 
# Southern Right Whale Health Monitoring Program of Península Valdés, Argentina 
# (University of California, Davis and Instituto de Conservación de Ballenas). 
# Requests to use the data should be addressed to Marcela Uhart muhart@ucdavis.edu 
# and Mariano Sironi mariano.sironi@icb.org.ar 

# intra-annual mortality
mintra <- read.csv("...Deaths_month_gulf_year.csv", sep = ";")

# attack variables data
datt <- read.csv("...Caf mortality data.csv", sep = ";")


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

# Editing the data...
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


# Ordinal regression plot --------------------------------------------------

# Predictions for the mean and for hipothetical years 

means_pred <- readRDS("....month_of_death_models_prediction_table.R")


agg_long$golfo <- revalue(agg_long$golfo, 
                          replace = c("gn"="Golfo Nuevo",
                                      "gsj"="Golfo San José"))

agg_long$attack_variable = revalue(agg_long$attack_variable,
                                   replace = c("gaf"="GAF",
                                               "gapc"="GAPC",
                                               "gapm" = "GAPM"))


means_pred$attack_variable = revalue(means_pred$attack_variable,
                                     replace = c("fa"="GAF",
                                                 "pac"="GAPC",
                                                 "pam" = "GAPM"))

means_pred$golfo <- revalue(means_pred$golfo, 
                            replace = c("gn"="Golfo Nuevo",
                                        "gsj"="Golfo San José"))

agg_long$attack_variable = factor(agg_long$attack_variable, 
                                  levels=c("GAPC", "GAPM","GAF"))

means_pred$attack_variable = factor(means_pred$attack_variable, 
                                    levels=c("GAPC", "GAPM","GAF"))

# Plot

ggplot(agg_long,
       aes(x = attack, y = month_num + 5,
           color = mortality, size = dmonth)) + 
  
  # estimation
  #   ci unconditional to year
  geom_ribbon(data = means_pred, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = year_lower, ymax = year_upper),
              inherit.aes = F, alpha = 0.0, color = "black",
              linetype = "dashed", size = 0.3) + 
  #   ci for average year
  geom_ribbon(data = means_pred, 
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
  geom_line(data = means_pred, 
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
  theme(panel.grid.minor = element_blank(),
        legend.position = "botton") + 
  ylab("Average month of death") + 
  xlab("Attack variable") + 
  labs(color = "Mortality", size = "Number of deaths")


# Doing individual plots so I can control the x-axis scales and values

# GAPC --------------------------------------------------------------------

agg_long_gn_pac <- subset(agg_long, attack_variable == "GAPC")

means_pred_gn_pac <- subset(means_pred, attack_variable == "GAPC")


ggplot(agg_long_gn_pac,
       aes(x = attack, y = month_num + 5,
           color = mortality, size = dmonth)) + 
  
  # estimation
  #   ci unconditional to year
  geom_ribbon(data = means_pred_gn_pac, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = year_lower, ymax = year_upper),
              inherit.aes = F, alpha = 0.0, color = "black",
              linetype = 3, size = 0.3) + 
  #   ci for average year
  geom_ribbon(data = means_pred_gn_pac, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = mean_lower, ymax = mean_upper,
                            alpha = in_range),
              inherit.aes = F, color = NA,
              show.legend = FALSE) + # remove legend for in_range
  
  # set alpha scale for ribbon
  scale_alpha_manual(values = c(0.2, 0.1)) + 
  
  # set new alpha scale for line
  new_scale("alpha") +
  scale_alpha_manual(values = c(1, 0.6)) + 
  
  #   mean
  geom_line(data = means_pred_gn_pac, 
            mapping = aes(x = attack, y = mean_mean, alpha = in_range),
            inherit.aes = F, size = 0.4, 
            show.legend = FALSE) + # remove legend for in_range
  
  # data
  geom_point() + 
  
  # year
  geom_text_repel(data = agg_long_gn_pac, 
                  mapping = aes(x = attack, y = month_num + 5, label = year),
                  inherit.aes = FALSE, size = 3,
                  box.padding = 0.1, max.overlaps = Inf) + 
  
  
  # settings
  scale_color_viridis() +
  facet_grid(rows = vars(golfo), cols = vars(attack_variable),
             scales = "free_x") + 
  scale_x_continuous(breaks = seq(0.2 ,9.4, by=3),
                     labels = scales::number_format(accuracy = 0.01))+
  scale_y_continuous(limits = c(7, 12), breaks = (seq(7, 12, by = 2)))+
  ylab("Average month of death") + 
  xlab("GAPC (attacks/h)") + 
  labs(color = "Mortality", size = "Number of deaths")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.text=element_text(size=12),
        legend.position = "bottom",
        legend.title = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_blank(),
        panel.grid.minor.x = element_blank())


# GAPM --------------------------------------------------------------------

agg_long_gn_pam <- subset(agg_long, attack_variable == "GAPM")

means_pred_gn_pam <- subset(means_pred, attack_variable == "GAPM")


ggplot(agg_long_gn_pam,
       aes(x = attack, y = month_num + 5,
           color = mortality, size = dmonth)) + 
  
  # estimation
  #   ci unconditional to year
  geom_ribbon(data = means_pred_gn_pam, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = year_lower, ymax = year_upper),
              inherit.aes = F, alpha = 0.0, color = "black",
              linetype = 3, size = 0.3) + 
  #   ci for average year
  geom_ribbon(data = means_pred_gn_pam, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = mean_lower, ymax = mean_upper,
                            alpha = in_range),
              inherit.aes = F, color = NA,
              show.legend = FALSE) + # remove legend for in_range
  
  # set alpha scale for ribbon
  scale_alpha_manual(values = c(0.16, 0.09)) + 
  
  # set new alpha scale for line
  new_scale("alpha") +
  scale_alpha_manual(values = c(1, 0.6)) + 
  
  #   mean
  geom_line(data = means_pred_gn_pam, 
            mapping = aes(x = attack, y = mean_mean, alpha = in_range),
            inherit.aes = F, size = 0.4, 
            show.legend = FALSE) + # remove legend for in_range
  
  # data
  geom_point() + 
  
  # year
  geom_text_repel(data = agg_long_gn_pam, 
                  mapping = aes(x = attack, y = month_num + 5, label = year),
                  inherit.aes = FALSE, size = 3,
                  box.padding = 0.1, max.overlaps = Inf) + 
  
  # settings
  scale_color_viridis() +
  facet_grid(rows = vars(golfo), cols = vars(attack_variable),
             scales = "free_x") + 
  scale_x_continuous(breaks = seq(0.4, 4, by = 1.1))+
  scale_y_continuous(limits = c(7, 12), breaks = (seq(7, 12, by = 2)))+
  xlab("GAPM (attacks/h)") +
  ylab("Average month of death") +
  labs(color = "Mortality", size = "Number of deaths")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.text=element_text(size=12),
        legend.position = "bottom",
        legend.title = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_blank(),
        panel.grid.minor.x = element_blank())


# GAF --------------------------------------------------------------------

agg_long_gn_fa <- subset(agg_long, attack_variable == "GAF")

means_pred_gn_fa <- subset(means_pred, attack_variable == "GAF")


ggplot(agg_long_gn_fa,
       aes(x = attack, y = month_num + 5,
           color = mortality, size = dmonth)) + 
  
  # estimation
  #   ci unconditional to year
  geom_ribbon(data = means_pred_gn_fa, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = year_lower, ymax = year_upper),
              inherit.aes = F, alpha = 0.0, color = "black",
              linetype = 3, size = 0.3) + 
  #   ci for average year
  geom_ribbon(data = means_pred_gn_fa, 
              mapping = aes(x = attack, y = mean_mean, 
                            ymin = mean_lower, ymax = mean_upper,
                            alpha = in_range),
              inherit.aes = F, color = NA,
              show.legend = FALSE) + # remove legend for in_range
  
  # set alpha scale for ribbon
  scale_alpha_manual(values = c(0.16, 0.09)) + 
  
  # set new alpha scale for line
  new_scale("alpha") +
  scale_alpha_manual(values = c(1, 0.6)) + 
  
  #   mean
  geom_line(data = means_pred_gn_fa, 
            mapping = aes(x = attack, y = mean_mean, alpha = in_range),
            inherit.aes = F, size = 0.4, 
            show.legend = FALSE) + # remove legend for in_range
  
  # data
  geom_point() + 
  
  # year
  geom_text_repel(data = agg_long_gn_fa, 
                  mapping = aes(x = attack, y = month_num + 5, label = year),
                  inherit.aes = FALSE, size = 3,
                  box.padding = 0.1, max.overlaps = Inf) + 
  
  # settings
  scale_color_viridis() +
  facet_grid(rows = vars(golfo), cols = vars(attack_variable),
             scales = "free_x") + 
  scale_x_continuous(breaks = seq(0.10, 0.40, by = 0.10))+
  scale_y_continuous(limits = c(7, 12), breaks = (seq(7, 12, by = 2)))+
  xlab("GAF") +
  ylab("Average month of death") + 
  labs(color = "Mortality", size = "Number of deaths")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.text=element_text(size=12),
        legend.position = "bottom",
        legend.title = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_blank(),
        panel.grid.minor.x = element_blank())



# Residuals analyses ------------------------------------------------------

res_list <- readRDS("....month_of_death_models_residuals.rds")
res_table <- data.frame("res_pam" = res_list$pam$scaledResiduals, 
                        "res_pac" = res_list$pac$scaledResiduals,  
                        "res_fa" = res_list$fa$scaledResiduals)

for(v in names(res_list)[1:3]) plot(res_list[[v]])

plot(res_list[["pac"]])
plot(res_list[["pam"]])
plot(res_list[["fa"]])

