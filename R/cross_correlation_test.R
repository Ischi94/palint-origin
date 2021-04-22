library(tidyverse)
library(here)
library(broom)

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




# calculate r squared -----------------------------------------------------


# all data
corr_all <- cor.test(all_data$temp, all_data$fragmentation_index, 
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
corr_cool <- cor.test(cool_data$temp, cool_data$fragmentation_index, 
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


