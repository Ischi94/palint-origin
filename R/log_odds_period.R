# loading packages
library(tidyverse) # for visualization
library(here) # for clean file storage
library(divDyn) # for origination rate calculation

# load custom functions
source(here("R/functions.R"))

# load data ---------------------------------------------------------------


# load stage data 
data(stages)

# load log_odds
load(here("data/log_odds_phyla.RData"))

# load occurrence data
load(file = here("data/occurrence_sqs_data.RData"))

# convert it to tibble
datsqs <- as_tibble(datsqs) %>% 
  count(phylum, stg) %>% 
  rename(name = phylum)



# calculate weighted log odds ---------------------------------------------


# get the log-odds per period by multiplying the 
# number of taxa within a period with their corresponding log-odds
taxa <- list(c(14:29), # Tremadocian - Lochkovian
     c(30:45), # Pragian - Artinskian 
     c(46:61), # Kungurian - Pliensbachian
     c(62:77), # Toarcian - Turonian
     c(78:94)) %>% # Coniacian - Pleistocene
  map_dfr(separate_groups) %>% 
  # weighing changed the actual scale
  # rescale between 0 and 1 to facilitate comparison
  # to observed data
  mutate(estimate = range_scale(estimate), 
         lower_CI = estimate - sd_est*1.96, 
         upper_CI = estimate + sd_est*1.96) 



# get the actual data and normalize between 0 and 1 
# to facilitate comparison
period <- log_odds[12:16,] %>%  # select stages
  select(-c(name, type)) %>% 
  # normalize between 0 and 1
  pivot_longer(cols = c(estimate, lower_CI, upper_CI)) %>% 
  mutate(value = range_scale(value)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  unnest(cols = c(estimate, lower_CI, upper_CI))
  


# visualize ---------------------------------------------------------------

# comparison
tibble(y = taxa$estimate, 
       ymin = taxa$lower_CI, 
       ymax = taxa$upper_CI, 
       x = period$estimate, 
       xmin = period$lower_CI, 
       xmax = period$upper_CI, 
       period = c(
         "Tremadocian - Lochkovian",
         "Pragian - Artinskian", 
         "Kungurian - Pliensbachian",
         "Toarcian - Turonian",
         "Coniacian - Pleistocene"
       )) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed") +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), 
                  size = 0.8, colour = "grey20") +
  geom_pointrange(aes(xmin = xmin, xmax = xmax), 
                  size = 0.8, colour = "grey20") +
  geom_point(size = 2, colour = "white") +
  geom_label(aes(label = period),
             nudge_y = 0.1) +
  labs(y = "Log-odds of taxa per period [std]", 
       x = "Log-odds per period [std]") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), 
                     limits = c(0, 1.05)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  my_theme +
  theme(panel.grid.major.x = element_blank())

# absolute differece
tibble(y = taxa$estimate, 
       ymin = taxa$lower_CI, 
       ymax = taxa$upper_CI, 
       x = period$estimate, 
       xmin = period$lower_CI, 
       xmax = period$upper_CI, 
       period = c(
         "Tremadocian - Lochkovian",
         "Pragian - Artinskian", 
         "Kungurian - Pliensbachian",
         "Toarcian - Turonian",
         "Coniacian - Pleistocene"
       )) %>% 
  mutate(diff_est = abs(y - x)) %>% 
  select(period, diff_est) %>% 
  bind_cols(log_odds[12:16,]) %>% 
  ggplot(aes(x = estimate, y = diff_est)) +
  geom_smooth(method = "lm", 
              colour = "#354E71", 
              fill = "grey80", 
              size = 2) +
  geom_point(size = 4.5, colour = "grey20") +
  geom_point(size = 2.8, colour = "white") +
  geom_label(aes(label = period),
             nudge_y = 0.03, 
             nudge_x = c(-0.33, 0.25, 0, 0, 0)) +
  labs(y = "Absolute difference", 
       x = "Log-odds per period") +
  my_theme +
  theme(panel.grid.major.x = element_blank())
