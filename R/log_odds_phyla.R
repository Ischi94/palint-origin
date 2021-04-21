# loading packages
library(tidyverse) # for visualization
library(here) # for project tidiness
library(divDyn) # for stage data

# load self-defined functions 
source(here("R/functions.R"))

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
  # remove hemichordata, nematoda, and hyolitha due to few data
  filter(!phylum %in% c("Hemichordata", "Hyolitha", "Nematoda")) %>% 
  pull() %>% 
  sort()

# build dataframe
log_odds <- tibble(name = c("Total", phyl_names, 
                            "Tremadocian-Lochkovian", 
                            "Pragian-Artinskian", 
                            "Kungurian-Pliensbachian", 
                            "Toarcian-Turonian", 
                            "Coniacian-Pleistocene"),
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

# take equally spaced sequences through time
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
  geom_vline(xintercept = 0, colour = "darkred") +
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
  annotate(geom = "rect", xmin = 0.4, xmax = 1.65, 
           ymin = - 0.6, ymax = 0.4, fill = "white") +
  # text
  annotate(geom = "text", x = 0.9, y = -0.1,
           colour = "grey30", label = "increasing likelihood", size = 2.5) +
  # arrow
  annotate(geom = "segment", x = 1.45, y = -0.1,  
           xend = 1.7, yend = -0.1, arrow = arrow(length = unit(2, "mm")), 
           colour = "grey40", size = 0.35) +
  coord_cartesian(xlim = c(0, 4), ylim = c(-0.05, 16)) +
  scale_fill_manual(values = c("grey40", "#d5a069", "indianred"))

forest_plot

# save plot
ggsave(plot = forest_plot, filename = here("figures/Log_Odds.png"), 
       width = 12.7, height = 9, units = "cm")

ggsave(plot = forest_plot, filename = here("figures/Log_Odds.pdf"), 
       width = 12.7, height = 9, units = "cm", dpi = 500)



# improve grouping --------------------------------------------------------

log_odds_impr <- log_odds[1:11, ] 

# Paleozoic
pal <- interaction_log_odds(all_phyla, 14:51) %>% 
  add_column(name = "Paleozoic", .before = "LowerCI", 
             type = "stages") %>% 
  select(name, lower_CI = LowerCI, estimate = OR, upper_CI = UpperCI, type)

log_odds_impr <- log_odds_impr %>% 
  add_row(pal)

# Mesozoic
mes <- interaction_log_odds(all_phyla, 52:81) %>% 
  add_column(name = "Mesozoic", .before = "LowerCI", 
             type = "stages") %>% 
  select(name, lower_CI = LowerCI, estimate = OR, upper_CI = UpperCI, type)

log_odds_impr <- log_odds_impr %>% 
  add_row(mes)

# Cenozoic
cen <- interaction_log_odds(all_phyla, 79:92) %>% 
  add_column(name = "Cenozoic", .before = "LowerCI", 
             type = "stages") %>% 
  select(name, lower_CI = LowerCI, estimate = OR, upper_CI = UpperCI, type)

log_odds_impr <- log_odds_impr %>% 
  add_row(cen)

# save data
# save(log_odds_impr, file= here("data/log_odds_impr.RData"))

#plot it
###forest plot
forest_plot_impr <- log_odds_impr %>% 
  # reorder labels
  mutate(name = as_factor(name)) %>% 
  mutate(name = fct_reorder(name, desc(name))) %>% 
  ggplot(aes(x = estimate, y = name)) +
  geom_linerange(aes(xmin = lower_CI, xmax = upper_CI)) +
  geom_vline(xintercept = 0, colour = "darkred") +
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
  annotate(geom = "rect", xmin = 0.4, xmax = 1.65, 
           ymin = - 0.6, ymax = 0.4, fill = "white") +
  # text
  annotate(geom = "text", x = 0.9, y = -0.1,
           colour = "grey30", label = "increasing likelihood", size = 2.5) +
  # arrow
  annotate(geom = "segment", x = 1.45, y = -0.1,  
           xend = 1.7, yend = -0.1, arrow = arrow(length = unit(2, "mm")), 
           colour = "grey40", size = 0.35) +
  coord_cartesian(xlim = c(0, 4), ylim = c(-0.05, 14)) +
  scale_fill_manual(values = c("grey40", "#d5a069", "indianred"))

forest_plot_impr

# save plot
ggsave(plot = forest_plot_impr, filename = here("figures/Log_Odds_impr.png"), 
       width = 12.7, height = 9, units = "cm")

ggsave(plot = forest_plot_impr, filename = here("figures/Log_Odds_impr.pdf"), 
       width = 12.7, height = 9, units = "cm", dpi = 500)


