# loading packages
library(tidyverse) # for visualization
library(here) # for clean file storage
library(deeptime) # for a nice plotting of geologic stages, see 
                  # https://github.com/willgearty/deeptime
library(divDyn) # for origination rate calculation
library(patchwork) # for putting plots together


# load data 
data(stages)

load(file = here("data/occurrence_sqs_data.RData"))

# convert it to tibble
datsqs <- as_tibble(datsqs)

# load theme
source(here("R/functions.R"))

# calculate rates
# interesting rates: 
# tOri: Number of originating taxa, taxa that have first occurrences in the focal bin, 
#           and last occurrences after it.
# oriProp: Proportional originations including single-interval taxa: 
#          (tOri + tSing) / (tThrough + tOri + tExt + tSing).
# oriPC: Per capita origination rates of Foote (1999). -log(tOri/(tOri + tThrough)). 
#           Values are not normalized with bin lengths
# ori3t: Three-timer origination rates of Alroy (2008). log(t2u/t3).
#
# oriC3t: Corrected three-timer origination rates of Alroy (2008). 
#           ori3t[i] + log(samp3t[i-1]).
# oriGF: Gap-filler origination rates of Alroy(2014). 
#           log((t2u + tPart)/(t3+tPart+tGFd))
# O2f3: Second-for-third origination proportions of Alroy (2015).
#
# ori2f3: Second-for-third origination rates (based on Alroy, 2015). 
#           Transformed to the usual rate form with log(1/(1-O2f3)).



# Per phyla ---------------------------------------------------------------

phyla <- c("Annelida", "Arthropoda", 
  "Brachiopoda", "Bryozoa", "Chordata", 
  "Cnidaria", "Echinodermata", "Hemichordata", 
  "Mollusca", "Porifera")

# preallocate empty list
metrics <- vector("list", length = length(phyla))

# give names to the list for subsetting
names(metrics) <- phyla

# calculate metric per phyla and save it in the list as a tibble
# note the additional phylum column, which is necessary for collapsing of the list later
for (i in phyla){
  metrics[[i]] <- datsqs %>% 
  # go through phyla
  filter(phylum == i) %>% 
  # calculate metrics
  divDyn(., tax = "genus", bin = "stg") %>% 
  as_tibble() %>% 
  # add phylum names and geologic ages calculated with stages file
  mutate(phylum = i, age = stages$mid[1:nrow(.)]) %>% 
  # order it properly
  select(phylum, age, stg, everything())  
}

# collapse or flatten the list to a dataframe for facetting  
metrics <- bind_rows(metrics)


# if you want to showcase a different metric in the plot, just change
# oriPC to whatever metric you fancy

origination_phyla <- ggplot(metrics) +
  geom_line(aes(x = age, y = oriPC), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(480, 0)) +
  coord_geo(size = 2.5, height = unit(0.75, "lines"), alpha = 2/3) +
  labs(x = "age [myr]", y = "Per-capita origination rate") +
  facet_wrap( ~ phylum, scales = "free", ncol = 2) +
  my_theme +
  theme(strip.text.x = element_text(size = 14, family = "sans"), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_blank())

# save it
ggsave(plot = origination_phyla, file=here("figures/per_capita_ori_rate.png"), 
       width = 210, height = 297, units = "mm")


# same for diversity

diversity_phyla <- ggplot(metrics) +
  geom_line(aes(x = age, y = divCSIB), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(480, 0)) +
  coord_geo(size = 2.5, height = unit(0.75, "lines"), alpha = 2/3) +
  labs(x = "age [myr]", y = "Corrected sampled-in-bin diversity") +
  facet_wrap( ~ phylum, scales = "free", ncol = 2) +
  my_theme +
  theme(strip.text.x = element_text(size = 14, family = "sans"), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.x = element_blank())

# save it
ggsave(plot = diversity_phyla, file=here("figures/corrected_sib_diversity.png"), 
       width = 210, height = 297, units = "mm")


# Total -------------------------------------------------------------------

total_metrics <- datsqs %>% 
  divDyn(., tax = "genus", bin = "stg") %>% 
  as_tibble() %>% 
  mutate(age = stages$mid[1:nrow(.)]) %>% 
  select(age, stg, everything())  



# Diversity (richness) ----------------------------------------------------


# Corrected sampled-in-bin diversity
divCSIB <- ggplot(total_metrics) +
  geom_step(aes(x = age, y = divCSIB), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(475, 0)) +
  coord_geo(size = list(2.5,3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2/3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"), 
            dat = list("epochs", "periods"), abbrv = TRUE) +
  labs(x = "age [myr]", y = "Number of genera", 
       title = "Corrected sampled-in-bin diversity") +
  my_theme

# Boundary-crosser diversity
divBC <- ggplot(total_metrics) +
  geom_step(aes(x = age, y = divBC), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(475, 0)) +
  coord_geo(size = list(2.5,3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2/3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"), 
            dat = list("epochs", "periods"), abbrv = TRUE) +
  labs(x = "age [myr]", y = "Number of genera", title = "Boundary-crosser diversity") +
  my_theme 

# Range-through diversity
divRT <- ggplot(total_metrics) +
  geom_step(aes(x = age, y = divRT), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(475, 0)) +
  coord_geo(size = list(2.5,3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2/3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"), 
            dat = list("epochs", "periods"), abbrv = TRUE) +
  labs(x = "age [myr]", y = "Number of genera", title = "Range-through diversity") +
  my_theme 

# arrange them on one plot
diversity <- divCSIB / divBC / divRT + plot_annotation(tag_levels = 'A')

# save it
ggsave(plot = diversity, file=here("figures/diversity_rates.png"), 
       width = 210, height = 297, units = "mm")



# Origination -------------------------------------------------------------

# Per capita origination rate
oriPC <- ggplot(total_metrics) +
  geom_step(aes(x = age, y = oriPC), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(475, 0)) +
  coord_geo(size = list(2.5,3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2/3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"), 
            dat = list("epochs", "periods"), abbrv = TRUE) +
  labs(x = "age [myr]", y = "Origination proportion", 
       title = "Per capita origination rate") +
  my_theme

# Three-timer origination rate
ori3t <- ggplot(total_metrics) +
  geom_step(aes(x = age, y = ori3t), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(475, 0)) +
  coord_geo(size = list(2.5,3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2/3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"), 
            dat = list("epochs", "periods"), abbrv = TRUE) +
  labs(x = "age [myr]", y = "Origination rate", title = "Three-timer origination rate") +
  my_theme 

# Second-for-third origination rate
ori2f3 <- ggplot(total_metrics) +
  geom_step(aes(x = age, y = ori2f3), colour = "grey15", size = 0.9) +
  scale_x_reverse(limits = c(475, 0)) +
  coord_geo(size = list(2.5,3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2/3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"), 
            dat = list("epochs", "periods"), abbrv = TRUE) +
  labs(x = "age [myr]", y = "Origination proportion", title = "Second-for-third origination rate") +
  my_theme 

# arrange them on one plot
origination <- oriPC / ori3t / ori2f3 + plot_annotation(tag_levels = 'A')

# save it
ggsave(plot = origination, file=here("figures/origination_rates.png"), 
       width = 210, height = 297, units = "mm")



# Number of observations --------------------------------------------------

observations <- datsqs %>%
  # count the observations per phyla
  count(phylum, name = "observations") %>%
  # plot it, but reorder levels
  ggplot(aes(x = reorder(phylum, observations), y = observations)) +
  geom_segment(aes(xend = reorder(phylum, observations), yend = 0), 
               size = 2, colour = "#354E71", alpha = 0.95) +
  geom_point(size = 4, shape = 21, colour = "grey20", fill = "#841F27") +  
  coord_flip() +
  my_theme +
  labs(x = "", y = "Number of observations")


# per stage, but only those with a sufficient number of observations
# 

# define colours
my_colours <- c(wesanderson::wes_palettes$Darjeeling1, 
                wesanderson::wes_palettes$Moonrise2,
                wesanderson::wes_palettes$GrandBudapest1[1:2])

stage_observations <- datsqs %>% 
  count(phylum, stg) %>% 
  left_join(stages) %>% 
  select(phylum, stg, n, mid) %>% 
  filter(phylum %in% phyla)  %>% 
  ggplot() +
  geom_col(aes(stg, n, fill = phylum), position = "stack") +
  my_theme +
  labs(x = "Geologic stages", y = "Number of observations") +
  scale_fill_manual(values = my_colours, guide = guide_legend(nrow = 2), 
                    name = NULL) +
  scale_x_continuous(breaks = seq(15, 95, by = 5))


# now proportional
prop_observations <- datsqs %>% 
  count(phylum, stg) %>% 
  left_join(stages) %>% 
  select(phylum, stg, n, mid) %>% 
  filter(phylum %in% phyla)  %>% 
  ggplot() +
  geom_col(aes(stg, n, fill = phylum), position = "fill") +
  my_theme +
  labs(x = "Geologic stages", y = "Percentage of observations") +
  scale_fill_manual(values = my_colours, guide = guide_legend(nrow = 2), 
                    name = NULL) +
  scale_x_continuous(breaks = seq(15, 95, by = 5))
  
    
# stack them and add a common legend
combined_observations <- stage_observations / guide_area() / prop_observations + 
         plot_annotation(tag_levels = 'A') +
         plot_layout(guides = "collect", heights = c(6,1,6)) 
          
  


ggsave(plot = combined_observations, file=here("figures/combined_observations.png"), 
       width = 210, height = 297, units = "mm") 



# isotope data ------------------------------------------------------------

# load data
veizer <- read_csv(file = here("data/raw_isotope_data.csv"))

# temperature data
isotemp <- read_csv(file=here("data/TimeSeriesUsed.csv")) 

# plot veizer isotope data
isotope_plot <- veizer %>%
  drop_na(d18O) %>% 
  mutate(fossil = as.factor(fossil), 
         fossil = fct_collapse(fossil, 
                               Brachiopod = c("brachiopod", "Brachiopod"), 
                               Belemnite = c("belemnite", "Belemnite"), 
                               Bivalve = c("bivalve", "Bivalve"), 
                               Foraminifera = c("Plankticf", "PlankticF", "BenthicF"))) %>% 
  filter(fossil %in% c("Brachiopod", "Belemnite", "Bivalve", "Foraminifera"), 
         gts2012 <= 500) %>% 
  ggplot(aes(gts2012, d18O, fill = fossil)) +
  geom_hline(yintercept = 0, linetype = "dotted", 
             colour = "grey20") +
  geom_point(shape = 21, colour = "grey30", 
             stroke = 1.2, alpha = 0.1, 
             show.legend = FALSE) +
  scale_fill_manual(values = my_colours[c(1, 2, 5, 7)]) +
  scale_x_reverse(limits = c(500, -5)) +
  coord_cartesian(expand = FALSE) +
  facet_wrap(~fossil, ncol = 1) +
  labs(x = "age [myr]", 
       y = expression(delta^18~O)) +
  my_theme +
  theme(strip.text.x = element_text(size = 12, 
                                    colour = "grey20"))

# plot average temperature
temp_plot <- isotemp %>% 
  rename(stg = Stage) %>% 
  left_join(stages) %>% 
  ggplot(aes(mid, Temp)) +
  geom_line(colour = "grey40", 
            size = 1.2) +
  coord_geo(size = list(2.5, 3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2 / 3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"),
            dat = list("epochs", "periods"), abbrv = TRUE, 
            xlim = c(500, -5)) +
  scale_x_reverse() +
  ylim(c(14, 40)) +
  labs(y = "Temperature [°C]", 
       x = "age [myr]") +
  my_theme

# combine plots 
combined_isotope <- isotope_plot / temp_plot +
  plot_layout(heights = c(2,1))

ggsave(plot = combined_isotope, file = here("figures/combined_isotopes.png"), 
       width = 210, height = 297, units = "mm") 



# time series of all parameter ----------------------------------------------

# origination signal
signal_plot <- 
  dat_final %>% 
  mutate(stg = as.numeric(as.character(bins))) %>% 
  left_join(stages) %>% 
  group_by(mid, systemCol) %>% 
  summarise(ori_signal = sum(origination)) %>% 
  ungroup() %>% 
  ggplot(aes(mid, ori_signal, fill = systemCol)) +
  geom_line(aes(group = 1)) +
  geom_point(shape = 21, size = 3, 
             alpha = 0.8, show.legend = FALSE) +
  scale_fill_identity() +
  scale_x_reverse(limits = c(500, -5)) +
  coord_cartesian(expand = FALSE, ylim = c(-10, 225)) +
  labs(x = "age [myr]", 
       y = "Origination events") +
  my_theme

# continental fragmentation
load(file = here("data/continental_fragmentation_data.RData"))
cont_plot <- 
  zaffos %>% 
  ggplot(aes(age_seq, fragmentation_index)) +
  geom_line() +
  scale_x_reverse(limits = c(500, -5)) +
  labs(x = "age [myr]", 
       y = "Fragmentation index") +
  coord_cartesian(expand = FALSE, ylim = c(0.26, 0.6)) +
  my_theme

# averaged temperature
temp_plot <- isotemp %>% 
  rename(stg = Stage) %>% 
  left_join(stages) %>% 
  ggplot(aes(mid, Temp)) +
  geom_line() +
  coord_geo(size = list(2.5, 3), height = list(unit(1.5, "lines"), unit(0.75, "lines")),
            alpha = 2 / 3, rot = list(90, 0),
            pos = list("bottom", "bottom"),  skip = c("Holocene", "Pleistocene"),
            dat = list("epochs", "periods"), abbrv = TRUE, 
            xlim = c(500, -5)) +
  scale_x_reverse() +
  ylim(c(14, 40)) +
  labs(y = "Temperature [°C]", 
       x = "age [myr]") +
  my_theme  

# combine plots
combined_parameter <- signal_plot / cont_plot / temp_plot + 
  plot_annotation(tag_levels = 'A')

ggsave(plot = combined_parameter, file = here("figures/combined_parameter.png"), 
       width = 210, height = 297, units = "mm")

