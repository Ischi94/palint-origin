# loading packages
library(tidyverse) # for visualization
library(here) # for clean file storage
library(divDyn) # for origination rate calculation
library(patchwork) # for combination of plots

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


# changes through time ----------------------------------------------------

datsqs %>% 
  full_join(log_odds) %>% 
  drop_na() %>% 
  mutate(weighted_estimate = n*estimate) %>% 
  group_by(name, stg) %>% 
  summarise(est = mean(weighted_estimate, na.rm = TRUE)) %>% 
  mutate(log_odds = (est - min(est)) / (max(est) - min(est))) %>% 
  full_join(log_odds) %>% 
  drop_na(stg) %>% 
  ggplot(aes(stg, fct_reorder(name, estimate), fill = log_odds)) +
  geom_tile() +
  labs(fill = "Log-odds [std]", x = "age [myr]", 
       y = NULL) +
  scale_x_continuous(breaks = c(25, 42, 59, 75, 95),
                     labels = c(400, 300, 200, 100, 0)) +
  scale_fill_viridis_c(option = "B") +
  my_theme +
  theme(panel.grid.major.x = element_blank())


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
dat_comp <- tibble(x = taxa$estimate, 
       xmin = taxa$lower_CI, 
       xmax = taxa$upper_CI, 
       y = period$estimate, 
       ymin = period$lower_CI, 
       ymax = period$upper_CI, 
       period = c(
         "Tremadocian - Lochkovian",
         "Pragian - Artinskian", 
         "Kungurian - Pliensbachian",
         "Toarcian - Turonian",
         "Coniacian - Pleistocene"
       )) 


# visualize
dat_comp_plot <- dat_comp %>% 
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed") +
  geom_pointrange(aes(ymin = ymin, ymax = ymax), 
                  size = 0.8, colour = "grey20") +
  geom_pointrange(aes(xmin = xmin, xmax = xmax), 
                  size = 0.8, colour = "grey20") +
  geom_point(size = 2, colour = "white") +
  geom_label(aes(label = period),
             nudge_y = 0.13, 
             nudge_x = c(-0.04, 0, 0, 0, 0.03)) +
  labs(x = "Log-odds of taxa per period [std]", 
       y = "Log-odds per period [std]") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  my_theme +
  theme(panel.grid.major.x = element_blank())




# calculate absolute difference
dat_diff <- dat_comp %>% 
  mutate(diff_est = abs(y - x)) %>% 
  select(period, diff_est) %>% 
  bind_cols(log_odds[12:16,])

# calculate correlation coefficient
diff_r <- cor.test(dat_diff$diff_est, dat_diff$estimate) %>% 
  broom::tidy() %>% 
  mutate(across(c(estimate, conf.low, conf.high), ~round(.x, 2)))

diff_r <- paste0(
  "r = ",
  diff_r$estimate, " [",
  diff_r$conf.low, ", ", 
  diff_r$conf.high, "]")

# visualize
dat_diff_plot <- dat_diff %>% 
  ggplot(aes(x = estimate, y = diff_est)) +
  geom_smooth(method = "lm", 
              colour = "#354E71", 
              fill = "grey80", 
              size = 2) +
  geom_point(size = 4.5, colour = "grey20") +
  geom_point(size = 2.8, colour = "white") +
  geom_label(aes(label = period),
             nudge_y = 0.08, 
             nudge_x = c(-0.25, 0.2, 0, 0, -0.1)) +
  annotate("text", x = 3.2, y = 0.58,
           label = diff_r, size = 4) +
  labs(y = "Absolute difference", 
       x = "Log Odds ratio \n (Origination | Cooling-Cooling)") +
  my_theme +
  theme(panel.grid.major.x = element_blank())

log_odds_explained <- dat_comp_plot / dat_diff_plot + 
  plot_annotation(tag_levels = 'A')


# save
ggsave(log_odds_explained, filename = here("figures/log_odds_by_taxa.png"))  

