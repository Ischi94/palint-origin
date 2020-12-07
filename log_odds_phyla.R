# loading packages
library(tidyverse) # for visualization
library(here) # for project tidiness

# Functions ---------------------------------------------------------------


# this functions calculates the log Odds ratio with Wald Confidence Intervals.
# With a default value for alpha at 0.05, which yields the 95% confidence intervals 
# for the computed odds ratio, based on the Wald approximation.  

calculate_log_odds <- function(n00, n01, n10, n11, alpha = 0.05) {
  ###
  # Compute the odds ratio between two binary outcomes, x and y
  # n00 = number of cases where x = 0 and y = 0
  # n01 = number of cases where x = 0 and y = 1
  # n10 = number of cases where x = 1 and y = 0
  # n11 = number of cases where x = 1 and y = 1
  ###
  
  # the odds ratio
  OR <- (n00 * n11) / (n01 * n10)
  
  # Compute the Wald confidence intervals:
  siglog <- sqrt((1 / n00) + (1 / n01) + (1 / n10) + (1 / n11))
  zalph <- qnorm(1 - alpha / 2)
  logOR <- log(OR)
  loglo <- logOR - zalph * siglog
  loghi <- logOR + zalph * siglog
  # produce output
  log_odds_output <-
    tibble(
      LowerCI = round(loglo, 3),
      OR = round(logOR, 3),
      UpperCI = round(loghi, 3)
    )
  log_odds_output
}



# function to calculate log Odds that taxa show higher probability of origination
# after cooling-ccoling palaeclimate interaction, per phyla

interaction_log_odds <- function(my_phylum, my_bins = unique(dat_final$bins)) {
  ###  warming-warming
  ww_raw <-
    subset(dat_final, trend.st7 >= 0 &
             warming >= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for warming-warming palaeclimate interaction with log odds as output
  ww_pred <-
    predict(warming_pal_int, newdata = ww_raw, type = "link")
  # how many warming-warming cases show higher probability
  h_ww <-
    length(subset(ww_pred, ww_pred > warming_sum$coefficients[1]))
  # how many warming-warming cases show lower probability
  l_ww <-
    length(subset(ww_pred, ww_pred < warming_sum$coefficients[1]))
  
  ###  cooling-warming
  cw_raw <-
    subset(dat_final, trend.st7 <= 0 &
             warming >= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for cooling-warming palaeclimate interaction with log odds as output
  cw_pred <-
    predict(warming_pal_int, newdata = cw_raw, type = "link")
  # how many cooling-warming cases show higher probability
  h_cw <-
    length(subset(cw_pred, cw_pred > warming_sum$coefficients[1]))
  # how many cooling-warming cases show lower probability
  l_cw <-
    length(subset(cw_pred, cw_pred < warming_sum$coefficients[1]))
  
  ###  warming cooling
  wc_raw <-
    subset(dat_final, trend.st6 >= 0 &
             cooling <= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for warming-cooling palaeclimate interaction with log odds as output
  wc_pred <- predict(cooling_pal_int, newdata = wc_raw,
                     type = "link")
  # how many warming-cooling cases show higher probability
  h_wc <- length(subset(wc_pred,
                        wc_pred > cooling_sum$coefficients[1]))
  # how many warming-cooling cases show lower probability
  l_wc <- length(subset(wc_pred,
                        wc_pred < cooling_sum$coefficients[1]))
  
  ###  cooling cooling
  cc_raw <-
    subset(dat_final, trend.st6 <= 0 &
             cooling <= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for cooling-cooling palaeclimate interaction with log odds as output
  cc_pred <- predict(cooling_pal_int, newdata = cc_raw,
                     type = "link")
  # how many cooling-cooling cases show higher probability
  h_cc <- length(subset(cc_pred,
                        cc_pred > cooling_sum$coefficients[1]))
  # how many cooling-cooling cases show lower probability
  l_cc <- length(subset(cc_pred,
                        cc_pred < cooling_sum$coefficients[1]))
  
  ##############
  # log odds
  
  # calculate log odds for cooling-cooling
  my_output <- calculate_log_odds(h_cc, l_cc, h_wc, l_wc)
  my_output
}

# Data --------------------------------------------------------------------

# load data
load(here("data/final_data.RData"))

# Split short term temperature change into warming and cooling, to calculate the results for each:
dat_final$cooling<-ifelse(dat_final$change.prev<0, dat_final$change.prev, NA)
dat_final$warming<-ifelse(dat_final$change.prev>0, dat_final$change.prev, NA)

# GLMM's
load(here("data/pal_int_models.RData"))

# final glmm for paleoclimate interaction and warming
warming_pal_int <- pal_int_models[[1]]
warming_sum <- summary(warming_pal_int)

# final glmm for paleoclimate interaction and cooling
cooling_pal_int <- pal_int_models[[2]]
cooling_sum <- summary(cooling_pal_int)



# Per phyla ---------------------------------------------------------------

# get phyla names
phyl_names <- dat_final %>% 
  distinct(phylum) %>%
  drop_na() %>% 
  # remove hemichordata and nematoda due to few data
  filter(phylum != "Hemichordata" & phylum != "Nematoda") %>% 
  pull() %>% 
  sort()

# build dataframe
log_odds <- tibble(name = c("Total", phyl_names, "Stage 14:29", "Stage 30:45", 
                            "Stage 46:61", "Stage 62:77", "Stage 78:94"),
                     lower_CI = numeric(length(name)), estimate = numeric(length(name)), 
                   upper_CI = numeric(length(name)))


# calculate log odds that taxa have higher origination rates after cooling-cooling
# palaeoclimate interaction

# Total
all_phyla <- unique(dat_final$phylum)
Total <- interaction_log_odds(all_phyla)
log_odds[1, 2:4] <- Total

# Annelida
Annelida <- interaction_log_odds("Annelida")
log_odds[2, 2:4] <- Annelida

# Arthropoda
Arthropoda <- interaction_log_odds("Arthropoda")
log_odds[3, 2:4] <- Arthropoda

# Brachiopoda
Brachiopoda <- interaction_log_odds("Brachiopoda")
log_odds[4, 2:4] <- Brachiopoda

# Bryozoa
Bryozoa <- interaction_log_odds("Bryozoa")
log_odds[5, 2:4] <- Bryozoa

# Chordata
Chordata <- interaction_log_odds("Chordata")
log_odds[6, 2:4] <- Chordata

# Cnidaria
Cnidaria <- interaction_log_odds("Cnidaria")
log_odds[7, 2:4] <- Cnidaria

# Echinodermata
Echinodermata <- interaction_log_odds("Echinodermata")
log_odds[8, 2:4] <- Echinodermata

# Foraminifera
Foraminifera <- interaction_log_odds("Foraminifera")
log_odds[9, 2:4] <- Foraminifera

# Mollusca
Mollusca <- interaction_log_odds("Mollusca")
log_odds[10, 2:4] <- Mollusca

# Porifera
Porifera <- interaction_log_odds("Porifera")
log_odds[11, 2:4] <- Porifera




# Through time ------------------------------------------------------------

# take equaly spaced sequences through time
age_seq <- dat_final %>% 
  as_tibble() %>% 
  distinct(bins) %>% 
  arrange(bins) %>% 
  pull(bins)

# Tremadocian to Lochkovian
trem_loch <- interaction_log_odds(all_phyla, age_seq[1:16])
log_odds[12, 2:4] <- trem_loch

# Pragian to Artinskian
prag_art <- interaction_log_odds(all_phyla, age_seq[17:32]) 
log_odds[13, 2:4] <- prag_art

# Kungurian to Pliensbachian
kun_pli <- interaction_log_odds(all_phyla, age_seq[33:48]) 
log_odds[14, 2:4] <- kun_pli

# Toarcian to Turonian
toa_tur <- interaction_log_odds(all_phyla, age_seq[49:64]) 
log_odds[15, 2:4] <- toa_tur

# Coniacian to Pleistocene
coni_plei <- interaction_log_odds(all_phyla, age_seq[65:81]) 
log_odds[16, 2:4] <- coni_plei

# add type column
log_odds <- log_odds %>% add_column(type = c("total", rep("phyla", 10), rep("stages", 5)))

# save data
# save(log_odds, file= here("data/log_odds_phyla.RData"))


# Plotting ----------------------------------------------------------------


#plot it
###forest plot
forest_plot <- log_odds %>% 
  # reorder labels
  mutate(name = as_factor(name)) %>% 
  mutate(name = fct_reorder(name, desc(name))) %>% 
  ggplot(aes(x = estimate, y = name)) +
  geom_linerange(aes(xmin = lower_CI, xmax = upper_CI)) +
  geom_hline(yintercept = 15.5, colour = "grey50") +
  geom_hline(yintercept = 5.5, colour = "grey50") +
  annotate(geom = "segment", x = 1.64, xend = 1.64, y = 0, yend = 16, 
           size = 2.5, colour = "darkred", alpha = 0.2) +
  geom_linerange(aes(xmin = lower_CI, xmax = upper_CI),
                 size = 1, colour = "grey45") +
  geom_point(aes(fill = type), size = 2, shape = 21, stroke = 0.5, colour = "grey25" ) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major.x=element_line(colour = "grey", linetype = "dotted"),
        text = element_text(family = "sans"), 
        legend.position = "none", 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_blank()) + 
  labs(x= "Log Odds ratio \n (Origination | Cooling-Cooling)", y = NULL)+
  scale_x_continuous(breaks = seq(0, 4, by = 1)) +
  # add annotations
  # box
  annotate(geom = "rect", xmin = -0.1, xmax = 1.15, 
           ymin = - 0.6, ymax = 0.4, fill = "white") +
  # text
  annotate(geom = "text", x = 0.4, y = -0.1,
           colour = "grey30", label = "increasing likelihood", size = 2.5) +
  # arrow
  annotate(geom = "segment", x = 0.95, y = -0.1,  
           xend = 1.2, yend = -0.1, arrow = arrow(length = unit(2, "mm")), 
           colour = "grey40", size = 0.35) +
  coord_cartesian(xlim = c(0, 4), ylim = c(-0.05, 16)) +
  scale_fill_manual(values = c("grey40", "#d5a069", "indianred"))

forest_plot

# save plot
ggsave(plot = forest_plot, filename = here("figures/Log_Odds.png"), 
       width = 12.7, height = 9, units = "cm")

ggsave(plot = forest_plot, filename = here("figures/Log_Odds.pdf"), 
       width = 12.7, height = 9, units = "cm", dpi = 500)
