                   ###### Figure 2, Fig. S1 and Fig. S2 #####

# Article title: Calf mortality in a whale population varies with seabird attacks                    
# Authors: María Piotto, Iván Barberá, Mariano Sironi, Victoria J. Rowntree, 
# Marcela M. Uhart, Macarena Agrelo, Alejandro A. Fernández Ajó, Jon Seger, 
#Carina F. Marón.

# Scrip author: María Piotto
# Questions can be addressed to mpiotto@mi.unc.edu.ar or to ivanbarbera93@gmail.com


# Library and data ----------------------------------------------------------

library(viridis)
library(ggplot2)
library(plyr)
library(tidyverse)
theme_set(theme_bw())

# Gull attack data belong to the Southern Right Whale Research Program (Instituto
# de Conservación de Ballenas, Argentina). Requests to use the data should be 
# addressed to Mariano Sironi mariano.sironi@icb.org.ar

pf_pred <- read.csv("...\\GAF and GAP prediction table.csv", 
                    header=TRUE, sep = ";")
fd <- read.csv("...GAF data.csv", sep = ",")
pd <- read.csv("...\\GAP data.csv", sep = ",")
model_summary <- read.csv("...\\GAP-GAF joint model_summary.csv", 
                          header=TRUE, sep = ",")

# About the data contained in pf_pred:
# variable: each gull attack index
# gulf: GN (Golfo Nuevo) vs. GSJ (Golfo San José), where de data was collected.
# year: the year when the data was collected (factor, as includes "avg" for
# the average value of each gull attack index). It has as many levels as years
# in which data was collected plus the "avg" level.
# nyear: the year when the data was collected (numeric).
# observed: observed values of each gull attack index.
# ngaf: observed values for GAF in 1996-2003 (not used for model fitting)
# mean: the mean of the posterior distribution of each gull attack index by year
# and by gulf (or averaging years when it is for the overall mean).
# lower and upper: 95 % credible intervals for the posterior distribution mean.
# sim_lower and sim_ upper: 95% highest density intervals for the posterior 
# predictive distribution in each year.
# mort: observed calf mortality.
# ymort: class of the year (high or low calf mortality years).

# About the data contained in fd:
# year: year when observations were made
# golfo = gulf: where the data was collected (GN: Golfo Nuevo, GJS: Golfo San José)
# intt:total number of intervals of the day
# intw: number of daily 5-min observation intervals with at least one attack on
# either mother or calf
# Obs: day of observations

# Then in pd:

# year: year when observations were made
# golfo = gulf: where the data was collected (GN: Golfo Nuevo, GJS: Golfo San José)
# mc: mother-calf factor. M: mother, C: calf
# pa: daily number of attacks on mothers or calves
# interv: total number of intervals of the day
# h: number of observation hours of the day
# HighLow: class of year (High or low mortality years)

# Finally, model_summary contains the summary of the posterior distributions of
# each parameter of the joint model


# Fig. 2 ------------------------------------------------------------------

# Some edits to plot...

pf_pred$variable = factor(pf_pred$variable, levels=c("GAPC",
                                                     "GAPM",
                                                     "GAF"))

# Edits to plot
pf_pred$gulf <- revalue(pf_pred$gulf, 
                            replace = c("GN" = "Golfo Nuevo",
                                        "GSJ" = "Golfo San José"))
                        
# predtype so we can differ between years and the overall average
pf_pred$predtype <- revalue(pf_pred$year, 
                            replace = c("1995"="year",
                                        "1996"="year",
                                        "1997"="year",
                                        "1998"="year",
                                        "1999"="year",
                                        "2000"="year",
                                        "2001"="year",
                                        "2002"="year",
                                        "2003"="year",
                                        "2004"="year",
                                        "2005"="year",
                                        "2006"="year",
                                        "2007"="year",
                                        "2008"="year",
                                        "2009"="year",
                                        "2010"="year",
                                        "2011"="year",
                                        "2012"="year",
                                        "2013"="year",
                                        "2014"="year",
                                        "2015"="year",
                                        "2016"="year",
                                        "2017"="year",
                                        "2018"="year",
                                        "2019"="year",
                                        "avg"="avg"))

# Adding 2020 so the annual average of each gull attack index will be separated 
# from they overall mean.
df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(df) <-c("variable", "gulf",	"year",	"nyear", "observed","ngaf",	"mean",	"lower",	
                 "upper",	"sim_lower",	"sim_upper", "mort", "ymort")
df[1,] <- c("GAPC", "Golfo Nuevo", 2020, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA) 


# Gull attack pressure on calves

gapc <- subset (pf_pred, variable == "GAPC")

ggplot(gapc, aes(x = year, y = mean, ymin = lower, ymax = upper,
                 shape = predtype, color = predtype)) +
  facet_grid(cols = vars(gulf), scale = "free_y") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#443A83FF", "#C7E020FF"), 
                     name = NULL)+
  geom_errorbar(mapping = aes(x = year, 
                              ymin = lower, ymax = upper), size = 0.6,  width = 0.5) +
  xlab("Year") + 
  ylab("GAPC (attacks/h)") +
  scale_shape_manual(values=c(17,19))+
  scale_x_discrete(breaks = c(1995, 1998, 2001, 2004, 2007, 2010, 2013, 2016, 2019))+
  geom_point(gapc, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 2)+
  geom_point(df, mapping = aes(x = year, y = as.numeric(mean)),
             inherit.aes = FALSE, alpha=0)+ # to add an space between years and avg
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        panel.grid.minor.y = element_blank())


# Gull attack pressure on mothers

gapm <- subset (pf_pred, variable == "GAPM")

ggplot(gapm, aes(x = year, y = mean, ymin = lower, ymax = upper,
                 colour = predtype)) +
  facet_grid(rows = vars(variable), cols = vars(gulf), scale = "free_y") +
  geom_point(size = 2.5) +
  scale_shape_manual(values=c(17,19))+
  scale_color_manual(values = c("#47C16EFF", "#47C16EFF"), 
                     name = NULL)+
  geom_errorbar(mapping = aes(x = year, 
                              ymin = lower, ymax = upper), size = 0.6,  width = 0.5) +
  ylab("Mean") +
  scale_x_discrete(breaks = c(1995, 1998, 2001, 2004, 2007, 2010, 2013, 2016, 2019))+
  geom_point(gapm, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 2)+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(),
        axis.text.y = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 12))

# Now I invent a new point so that both GAPC and GAPM would be presented in 
# the same scale and you would be able to compare them. This point will be included 
# with alpha = 0 as it is not part of a real data (you will not see the point 
# in the graph).

df[2,] <- c("GAPM", "Golfo Nuevo", 2020, NA, NA, NA, 13.7, #13.7 is the upper CI 95 of GAPC in 2011
            NA, NA, NA, NA, NA, NA) 

ggplot(gapm, aes(x = year, y = mean, ymin = lower, ymax = upper,
                 colour = predtype, shape = predtype)) +
  facet_grid(cols = vars(gulf), scale = "free_y") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#443A83FF", "#47C16EFF"), 
                     name = NULL)+
  ylab("GAPM (attacks/h)")+
  xlab("Year")+
  geom_errorbar(mapping = aes(x = year, 
                              ymin = lower, ymax = upper), size = 0.6,  width = 0.5) +
  scale_x_discrete(breaks = c(1995, 1999, 2003,  2007, 2011,  2015, 2019))+
  geom_point(gapm, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 2)+
  geom_point(df, mapping = aes(x = year, y = as.numeric(mean)),
             inherit.aes = FALSE, alpha=0)+
  scale_shape_manual(values=c(17,19))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        panel.grid.minor.y = element_blank())

# Again, adding 2020 so the annual average will be separated from they overall mean
df[3,] <- c("GAF", "Golfo Nuevo", 2020, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA) 

gaf <- subset (pf_pred, variable == "GAF")

ggplot(gaf, aes(x = year, y = mean, ymin = lower, ymax = upper,
                 colour = predtype, shape = predtype)) +
  facet_grid(cols = vars(gulf), scale = "free_y") +
  geom_point(size = 2.5) +
  scale_shape_manual(values=c(17,19))+
  scale_color_manual(values = c("#443A83FF","#287C8EFF"), 
                     name = NULL)+
  geom_errorbar(mapping = aes(x = year, 
                              ymin = lower, ymax = upper), size = 0.6,  width = 0.5) +
  geom_point(gaf, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 2)+
  geom_point(gaf, mapping = aes(x = year, y = ngaf),
             inherit.aes = FALSE, shape = 0, colour = "black", size = 2)+
  scale_x_discrete(breaks = c(1995, 2001,  2007, 2013,  2019))+
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2))+
  xlab("Year") + 
  ylab("GAF") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11, colour = "black"),
        panel.grid.minor.y = element_blank())


# Observed calf mortality

mgn<-subset(pf_pred, gulf == "Golfo Nuevo")
mean(na.omit(mgn$mort)) # this is the average mortality in Golfo Nuevo, you could
# check it looking for the mort value when year is equal to avg and gulf is equal
# to GN
mgsj<-subset(pf_pred, gulf == "Golfo San José")
mean(na.omit(mgsj$mort)) # this is the average mortality in Golfo San José
# As we need a continuous x-axis so as the dots are connected by a line, we could not
# add the average mortality by golf in the same plot. To do so, we will include
# those mortality values as if they would have corresponded to 2021. In this way, 
# means could be treated differently from the other values and would be shown
# separated of the other points as well.
# Obviously, the "2021" label will not be included, as is not correct. This is
# just a little trick to include means in this plots while been able to connect
# the other point with a line.

df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <-c("gulf", "year", "nyear","mort")

df[1,] <- c("Golfo Nuevo", NA, 2021, 0.219)
df[2,] <- c("Golfo San José", NA, 2021, 0.118)

df0 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df0) <-c("gulf",	"year", "nyear","mort")
df0[1,] <- c("Golfo San José", NA, 2001, 0.55) # to modify the scale so
# it is the same as the one used in GAF-plot

hy <- subset(pf_pred, ymort == "high") # to differentiate years of high and
# low calf mortality

ggplot(pf_pred, aes(x = nyear, y = mort)) +
  facet_grid(cols = vars(gulf)) +
  geom_smooth(method = "loess", se = FALSE, 
              span = 0.18, size = 0.7, linetype = 3, color = "black")+
  geom_point(size = 2.2, color = "black", shape = 21, stroke = 0.7, fill = "white") +
  xlab("Year") + 
  ylab("Observed calf mortality") +
  labs(color="Mortality")+
  scale_x_continuous(breaks = seq(1995,2021, by=6))+
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2))+
  geom_point(df, mapping = aes(x = as.numeric(nyear), y = as.numeric(mort)),
             inherit.aes = FALSE, size = 2, color = "black", shape = 2,
             stroke = 0.9)+ # I add the mean of mortality per golf in 2021, so 
  #it is separeted from the others values
  geom_point(hy, mapping = aes(x = as.numeric(nyear), y = as.numeric(mort)),
             inherit.aes = FALSE, shape = 1, size = 0.9, color = "black",
             stroke = 1.5)+ # high mortality data
  geom_point(df0, mapping = aes(x = as.numeric(nyear), y = as.numeric(mort)),
             inherit.aes = FALSE, alpha = 0)+ #just to move the scale
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11,  colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11, colour = "black"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())


# Posterior predictive checks ---------------------------------------------

pd <- pd[pd$interv >0, ] # remove rows with 0 hours of observation

# set equal levels for factors and complete names

fd$golfo <- revalue(fd$golfo, replace = c("GSJ" = "Golfo San José", "GN" = "Golfo Nuevo"))
pd$golfo <- revalue(pd$golfo, replace = c("SJ" = "Golfo San José", "N" = "Golfo Nuevo"))
pd$mc <- revalue(pd$mc, 
                 replace = c("M"="Mother",
                             "C"="Calf"))

years <- unique(fd$year)
#years == unique(pd$year) # OK

fd$year <- factor(fd$year, levels = as.character(years))
pd$year <- factor(pd$year, levels = as.character(years))

fd$golfo <- factor(fd$golfo, levels = c("Golfo Nuevo", "Golfo San José"))
pd$golfo <- factor(pd$golfo, levels = c("Golfo Nuevo", "Golfo San José"))
pd$mc <- factor(pd$mc, levels = c("Mother", "Calf"))

# order

fd <- fd[with(fd, order(year, golfo)), ]
pd <- pd[with(pd, order(year, golfo, mc)), ]


# observed proportion and rate
pd$rate <- pd$pa / pd$h
fd$prop <- fd$intw / fd$intt


gap_pred <- subset(pf_pred, variable != "GAF") 
gaf_pred <- subset(pf_pred, variable == "GAF")


# Plots

ggplot(gaf_pred, aes(x = year, y = mean, ymin = lower, ymax = upper,
                     colour = predtype)) +
  geom_linerange(mapping = aes(x = year, y = mean, ymin = sim_lower,
                               ymax = sim_upper), col = "grey60", 
                 size = 2, alpha = 0.3) +
  geom_point(fd, mapping = aes(x = year, y = prop),
            inherit.aes = FALSE, colour = "black", shape = 1, alpha = 0.6)+
  geom_point(size = 2.5) + 
  geom_errorbar(size = 0.6,  width = 0.5) +
  geom_point(gaf_pred, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 0, colour = "black", size = 2) +
  scale_color_manual(values = c("#443A83FF", "#287C8EFF"), 
                     name = NULL)+
  facet_wrap(vars(gulf)) +
  xlab("Year") + 
  ylab("GAF") +
  scale_x_discrete(breaks = c(1995, 2001, 2007, 2013, 2019))+
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.3))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text = element_text(size = 14),
        panel.grid.minor = element_blank())

gapm_pred <- subset(gap_pred, variable == "GAPM")
pdm <- subset(pd, mc == "Mother")

ggplot(gapm_pred, aes(x = year, y = mean, ymin = lower, ymax = upper,
                      colour = predtype)) +
  geom_linerange(mapping = aes(x = year, y = mean, ymin = sim_lower,
                               ymax = sim_upper), col = "grey60", 
                 size = 2, alpha = 0.3) +
  geom_point(pdm, mapping = aes(x = year, y = rate),
             inherit.aes = FALSE, colour = "black", shape = 1,
             alpha = 0.6)+
  geom_point(size = 2.5) + 
  geom_errorbar(size = 0.6,  width = 0.5) +
  geom_point(gapm_pred, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 0, size = 2) +
  scale_color_manual(values = c("#443A83FF", "#47C16EFF"), 
                     name = NULL)+ 
  xlab("Year") + 
  ylab("GAPM (attacks/h)") +
  scale_x_discrete(breaks = c(1995, 2001, 2007, 2013, 2019))+
  scale_y_continuous(breaks = seq(0, 10, by = 5))+
  facet_grid(cols = vars(gulf)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text = element_text(size = 14),
        panel.grid.minor = element_blank())

gapc_pred <- subset(gap_pred, variable == "GAPC")
pdc <- subset(pd, mc == "Calf")

ggplot(gapc_pred, aes(x = year, y = mean, ymin = lower, ymax = upper,
                      colour = predtype)) +
  geom_linerange(mapping = aes(x = year, y = mean, ymin = sim_lower,
                               ymax = sim_upper), col = "grey60", 
                 size = 2, alpha = 0.3) +
  geom_point(pdc, mapping = aes(x = year, y = rate),
             inherit.aes = FALSE, colour = "black", shape = 1,
             alpha = 0.6)+
  geom_point(size = 2.5) + 
  geom_errorbar(size = 0.6,  width = 0.5) +
  geom_point(gapc_pred, mapping = aes(x = year, y = observed),
             inherit.aes = FALSE, shape = 0, size = 2) +
  scale_color_manual(values = c("#443A83FF", "#C7E020FF"), 
                     name = NULL)+ 
  xlab("Year") + 
  ylab("GAPC (attacks/h)") +
  scale_x_discrete(breaks = c(1995, 2001, 2007, 2013, 2019))+
  facet_grid(cols = vars(gulf)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12, colour = "black"),
        strip.text = element_text(size = 14),
        panel.grid.minor = element_blank())


# Filled points are observed means,
# open ones are observed daily data (proportion or rate),
# blue points are posterior means for the means, 
# blue bars are the 95 % equal tails credible intervals for those means,
# and red bars are 95 % highest density intervals for the posterior predictive
# distribution in each year. 
# When the model has good coverage, approximately 95 % of the data should fall 
# within the red bars.


# N eff and R hat ---------------------------------------------------------

min(model_summary$n_eff, na.rm = TRUE) # 1260.376
max(model_summary$Rhat, na.rm = TRUE)  # 1.008666

hist(model_summary$n_eff, xlab = "Effective sample size", main = NULL,
     breaks = 50)
text(3500, 225, "min: 1260.376")

hist(model_summary$Rhat, xlab = "R_hat", main = NULL,
     breaks = 35)
text(1.006, 500, "max: 1.008666")
