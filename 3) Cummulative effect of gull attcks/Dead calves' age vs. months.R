
                      ##### Dead calf age over months #####

# Article title: Calf mortality in a whale population varies with seabird attacks                    
# Authors: María Piotto, Iván Barberá, Mariano Sironi, Victoria J. Rowntree, 
# Marcela M. Uhart, Macarena Agrelo, Alejandro A. Fernández Ajó, Jon Seger, 
#Carina F. Marón.
                     
# Scrip author: María Piotto
# Questions can be addressed to mpiotto@mi.unc.edu.ar or to ivanbarbera93@gmail.com
               
                             
# Library and dataset ---------------------------------------------------------
                      
library(tidyverse)
library(plyr)
theme_set(theme_bw())


# Calf mortality data belong to the Southern Right Whale Health Monitoring Program 
# of Península Valdés, Argentina (University of California, Davis and Instituto 
# de Conservación de Ballenas). Requests to use the data should be addressed to 
# Marcela Uhart muhart@ucdavis.edu and Mariano Sironi mariano.sironi@icb.org.ar 


d <- read.csv("...Dead calves age.csv", header=TRUE, sep = ";")

# About the database:
#
# You will find the monthly average length of dead calves and the monthly
# absolute number of dead calves with healed and open umbilical cord. The column
# named "value" has these measures. "sc" indicates if the measures are on calves'
# "size" (length) or "cords", and "condition" indicates the umbilical cord 
# condition associated to "cords" measures. Finally, the columns names "month" 
# and "year" indicate the month and the class of year in which the data were 
# collected. "high": high mortality years; "low": low mortality years;
# and "all": the sum between low and high mortality years.


# Dead calves' length vs. month -----------------------------------------------

m <- subset(d, sc == "size")
mt <- subset(m, year == "all")
cor(mt$month, mt$value, method = "pearson")

m$year <- revalue(m$year, 
                  replace = c("high"="high mortality",
                             "low"="low mortality"))

m$Mortality <- m$year

sh <- subset(m, m$year == "high mortality")
cor(sh$month, sh$value, method = "pearson")

sl <- subset(m, year = "low mortality")
cor(sl$month, sl$value, method = "pearson")

ggplot(subset(m, year != "all"), aes(x = as.numeric(month), y = as.numeric(value)))+
  geom_point(shape = 1, colour = "black", size = 2.5)+
  geom_point(sh, mapping = aes(x = as.numeric(month), y = value),
             inherit.aes = FALSE, shape = 1, colour = "black", size = 2.5, stroke = 2)+
  geom_hline(yintercept = 4.9, linetype="dashed", size = 0.85)+ # estimated average length of SRW neonates (Christiansen et al., 2022)
  #scale_y_continuous(breaks = seq(2, 8, by = 3))+
  scale_x_continuous(breaks = seq(6, 12, by = 2))+
  ylab("Length of dead calves (m)")+
  xlab("Month")+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.text=element_text(size=13),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 13,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 13, colour = "black"),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=13))



# Umbilical cord condition vs. month ------------------------------------------

c <- subset(d, sc == "cord")

c$year <- revalue(c$year, 
                replace = c("high"="High mortality years",
                            "low"="Low mortality years"))

ca <- subset(c, c$year == "all") # independently on high-low calf mortality years
ct <- as.data.frame(ca %>% group_by(month) %>% dplyr::summarize(value = sum(value)))
# monthly absolute number of calves whose umbilical cord condition was recorded

co <- subset(ca, ca$condition == "open") # independently on high-low mortality years,
# how many calves have their cord open each month
co[7,1] <- 12 #adding December
co[7,3] <- 0 # adding 0 to December as no dead calf with open cord was recorded in that month
co$rel <- co$value / ct$value # relative frequency of calves with open cords
cor(co$month, co$rel, method = "pearson")

ch <- subset(ca, ca$condition == "healed") # independently on high-low mortality years,
# how many calves have their cord healed each month
ch_r <- c(6, NA, 0, "cord", "healed") # adding June with 0 dead calves with healed umbilical cord
ch <- rbind(ch_r, ch)
ch$rel <- as.numeric(ch$value) / ct$value # relative frequency of calves with healed cords
cor(as.numeric(ch$month), as.numeric(ch$rel), method = "pearson")

p <- subset(c, year != "all")

ggplot(p, aes(x = month, y = value, color = condition, fill = condition)) +
  geom_bar(stat="identity", alpha = 0.9)+
  facet_grid(cols=vars(year))+
  ylab("Absolute frequency")+
  xlab("Month")+
  scale_fill_manual(values = c("#C7E020FF","#47C16EFF"), 
                    name = NULL)+
  scale_color_manual(values = c("#C7E020FF","#47C16EFF"), 
                     name = NULL)+
  scale_y_continuous(breaks = seq(0, 70, by = 30))+
  scale_x_continuous(breaks = c(6, 8, 10, 12))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 13,  colour = "black"),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 13, colour = "black"),
        strip.text = element_text(size = 14),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

