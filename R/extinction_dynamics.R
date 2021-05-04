# loading packages
library(tidyverse) # for visualization
library(here) # for clean file storage
library(divDyn) # for origination rate calculation
library(broom) # for dataframe tidying
library(patchwork) # to combine images

# source function for correlations
source(here("R/functions.R"))

# load data ---------------------------------------------------------------


data("stages")

# load occurrence data
load(file = here("data/occurrence_sqs_data.RData"))

# load data
load(here("data/final_data.RData"))


# calculate extinction metrics from sqs data
ext_rate <- divDyn(datsqs, tax = "genus", bin = "stg") %>% 
  as_tibble() %>% 
  filter(stg %in% c(15:94))

# calculate origination signal
ori_signal <- dat_final %>% 
  mutate(stg = as.numeric(as.character(bins))) %>% 
  group_by(stg) %>% 
  summarise(ori_signal = sum(origination)) %>% 
  ungroup() %>% 
  filter(stg %in% c(15:94))


cor_ori_ext <- function(ext_rate_col, 
                        lagged = FALSE,
                        detrended = FALSE, 
                        ext_input = NULL) {
  
  ext_input <- ext_rate[, ext_rate_col] %>% 
    pull(ext_rate_col)
  
  ori_input <- ori_signal$ori_signal
  
  if(lagged == TRUE) {
    ext_input <- lag(ext_input)
  }

  if(detrended == TRUE) {
    ext_input <- diff(ext_input)
    ori_input <- diff(ori_input)
  }

  cor.test(ext_input, ori_input)  %>%
    tidy() %>%
    select(estimate, conf.low, conf.high) %>%
    add_column(metric = ext_rate_col)
}




# direct correlation ------------------------------------------------------

# get metric names
metrics <- c("tExt", "extProp", "extPC", 
             "ext3t", "extC3t", "extGF",
             "E2f3", "ext2f3")

# without detrending
plot1 <- metrics %>% 
  map_dfr(cor_ori_ext) %>% 
  ggplot(aes(x = estimate, y = metric, 
             xmin = conf.low, 
             xmax = conf.high)) +
  geom_vline(xintercept = 0, 
             colour = "darkred", 
             linetype = "dashed", 
             size = 1) +
  geom_pointrange(size = 1.2, colour = "grey20") +
  geom_point(colour = "white", size = 4) +
  labs(y = NULL, x = "Pearson correlation coefficient") +
  my_theme


# with detrending
plot2 <- metrics %>% 
  map_dfr(cor_ori_ext, detrended = TRUE) %>% 
  ggplot(aes(x = estimate, y = metric, 
             xmin = conf.low, 
             xmax = conf.high)) +
  geom_vline(xintercept = 0, 
             colour = "darkred", 
             linetype = "dashed", 
             size = 1) +
  geom_pointrange(size = 1.2, colour = "grey20") +
  geom_point(colour = "white", size = 4) +
  labs(y = NULL, x = "Pearson correlation coefficient") +
  my_theme



# laged correlation -------------------------------------------------------

# without detrending
plot3 <- metrics %>% 
  map_dfr(cor_ori_ext, lagged = TRUE) %>% 
  ggplot(aes(x = estimate, y = metric, 
             xmin = conf.low, 
             xmax = conf.high)) +
  geom_vline(xintercept = 0, 
             colour = "darkred", 
             linetype = "dashed", 
             size = 1) +
  geom_pointrange(size = 1.2, colour = "grey20") +
  geom_point(colour = "white", size = 4) +
  labs(y = NULL, x = "Pearson correlation coefficient") +
  my_theme


# with detrending
plot4 <- metrics %>% 
  map_dfr(cor_ori_ext, lagged = TRUE, detrended = TRUE) %>% 
  ggplot(aes(x = estimate, y = metric, 
             xmin = conf.low, 
             xmax = conf.high)) +
  geom_vline(xintercept = 0, 
             colour = "darkred", 
             linetype = "dashed", 
             size = 1) +
  geom_pointrange(size = 1.2, colour = "grey20") +
  geom_point(colour = "white", size = 4) +
  labs(y = NULL, x = "Pearson correlation coefficient") +
  my_theme



# combine images ----------------------------------------------------------

# for direct correlation
direct_correlation <- plot1 + plot2 +
  plot_annotation(tag_levels = "A")

ggsave(plot = direct_correlation, 
       filename = here("figures/ori_ext_direct.png"))

# for lagged correlation
lagged_correlation <- plot3 + plot4 +
  plot_annotation(tag_levels = "A")

ggsave(plot = lagged_correlation, 
       filename = here("figures/ori_ext_lagged.png"))
