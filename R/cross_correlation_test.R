library(tidyverse)
library(here)
library(broom)

# load in custom functions
source(here("R/functions.R"))

# load data ---------------------------------------------------------------

# temperature data
temperature <- read_csv(here("data/TimeSeriesUsed.csv"))

# continental fragmentation index binned to stages
cont_frag <- read_csv(here("data/continental_fragmentation_data_binned.csv"))



# combine data ------------------------------------------------------------

# all data
all_data <- temperature %>% 
  select(stg = Stage, temp = Temp) %>% 
  left_join(cont_frag)


# select only those data cooler than the average
cool_data <- all_data %>% 
  filter(temp <= mean(temp))



# calculate cross correlation ---------------------------------------------

# for all data
ccf_all <- all_data %>% 
  summarise(cor_test = list(tidy(cross_corr(.))), 
          n.used = length(stg)) %>% 
  mutate(lag = map(cor_test, pluck, 1), 
         acf = map(cor_test, pluck, 2),
         ci = cross_corr_ci(n.used)) %>% 
  select(lag, acf, ci) %>% 
  unnest(c(lag, acf)) %>% 
  add_column(type = "All data")

# for cool data
ccf_cool <- cool_data %>% 
  summarise(cor_test = list(tidy(cross_corr(.))), 
            n.used = length(stg)) %>% 
  mutate(lag = map(cor_test, pluck, 1), 
         acf = map(cor_test, pluck, 2),
         ci = cross_corr_ci(n.used)) %>% 
  select(lag, acf, ci) %>% 
  unnest(c(lag, acf)) %>% 
  add_column(type = "Cool subset")





# pearson r at lag one ----------------------------------------------------


# all data
corr_all <- cor.test(diff(all_data$temp), diff(all_data$fragmentation_index), 
           alternative = "two.sided", 
           method = "pearson") %>% 
  tidy() %>% 
  select(estimate, conf.low, conf.high) %>% 
  mutate(across(everything(), round, 2)) %>% 
  transmute(
    estimate = str_c(estimate, "[", 
                     sep = " "),
    estimate = str_c(estimate, conf.low, 
                          sep = ""), 
    estimate = str_c(estimate, conf.high, 
                     sep = ", "),
    estimate = str_c(estimate, "]", 
                     sep = "")) %>% 
  pull()

# cool data
corr_cool <- cor.test(diff(cool_data$temp), diff(cool_data$fragmentation_index), 
                     alternative = "two.sided", 
                     method = "pearson") %>% 
  tidy() %>% 
  select(estimate, conf.low, conf.high) %>% 
  mutate(across(everything(), round, 2)) %>% 
  transmute(
    estimate = str_c(estimate, "[", 
                     sep = " "),
    estimate = str_c(estimate, conf.low, 
                     sep = ""), 
    estimate = str_c(estimate, conf.high, 
                     sep = ", "),
    estimate = str_c(estimate, "]", 
                     sep = "")) %>% 
  pull()




# visualize correlation ---------------------------------------------------

  
# merge data and plot
correlation_plot <- ccf_cool %>% 
  add_column(cor_estimate = corr_cool) %>% 
  full_join(ccf_all %>% 
              add_column(cor_estimate = corr_all)) %>% 
  mutate(cor_estimate = paste0("r = ", cor_estimate)) %>% 
  ggplot(aes(lag, acf, colour = type)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_segment(aes(xend = lag, 
                   y = 0, yend = acf)) +
  geom_line(aes(lag, ci), 
            linetype = "dashed") +
  geom_line(aes(lag, -ci), 
            linetype = "dashed") + 
  annotate(geom = "rect", xmin = -9.5, xmax = -2, 
           ymin = 0.295, ymax = 0.33, 
           colour = "white", fill = "white") +
  geom_text(aes(y = 0.31, x = -6, 
                label = cor_estimate),
            size = 4, colour = "grey20", 
            fontface = "plain", family = "serif") +
  scale_color_manual(values = c("#841F27", "#354E71")) +
  facet_wrap(~type) +
  labs(y = "Cross-correlation", x = "Lag") +
  my_theme +
  theme(legend.position = "none", 
        strip.text.x = element_text(size = 11, 
                                    colour = "grey20")) 

ggsave(correlation_plot, filename = here("figures/correlation_plot.png"))  




# overall cross correlation -----------------------------------------------


# get subsets of data
dat_subset <- list(
  all_data, # all data
  cool_data, # cool data only
  all_data %>%
    filter(stg <= 51), # Paleozoic
  all_data %>%
    filter(stg > 51 & stg < 82), # Mesozoic
  all_data %>%
    filter(stg >= 82) # Cenozoic
) 

# cross-correlation between detrended temperature and detrended continental
# fragmentation index
pearson_r <- dat_subset %>% 
  map_dfr(cor_test_r) %>% 
  add_column(.before = "estimate", 
             type = c("All Data", 
                      "Cool subset", 
                      "Paleozoic", 
                      "Mesozoic", 
                      "Cenozoic"))

# check for auto-correllation in residuals
pearson_r <- dat_subset %>% 
  map_dfr(~ lmtest::dwtest(diff(.x$temp) ~ diff(.x$fragmentation_index)) %>% 
            tidy() %>% 
            select(statistic, p.value)) %>% 
  add_column(.before = "statistic", 
             type = c("All Data", 
                      "Cool subset", 
                      "Paleozoic", 
                      "Mesozoic", 
                      "Cenozoic")) %>% 
  rename(DW_p_value = p.value) %>% 
  full_join(pearson_r) %>% 
  select(type, estimate, conf.low,
         conf.high, statistic, DW_p_value)

# save data
# write_csv(pearson_r, here("data/pearson_r.csv"))

# visualize autocorrelation
pearson_r_plot <- pearson_r %>% 
  ggplot(aes(estimate, fct_reorder(type, c(5, 4, 3, 2, 1)), 
             xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, alpha = 0.8,
             colour = "#841F27", size = 1.2) +
  geom_pointrange(colour = "grey30", size = 1.2) +
  geom_point(size = 3, colour = "white") +
  labs(y = NULL, 
       x = "Pearson correlation coefficient") +
  my_theme

ggsave(pearson_r_plot, filename = here("figures/pearson_r_plot.png")) 



