library(brms)
library(bayesplot)
library(here)
library(tidyverse)

# load self-defined functions 
source(here("R/functions.R"))

# set up colour scheme
color_scheme_set("red")

# read in data
brms_best <- readRDS(here("data/brms_best.rds"))


# model convergence -------------------------------------------------------


# effective number of samples
neff_plot <- neff_ratio(brms_best) %>% 
  enframe(name = "estimate", value = "neff_ratio") %>% 
  mutate(estimate = c("Cooling-Cooling", 
                      "All Other", 
                      "Sigma Cooling-Cooling", 
                      "Sigma All Other", 
                      "nu", 
                      "LP"), 
         estimate = as_factor(estimate), 
         estimate = fct_reorder(estimate, neff_ratio)) %>% 
  filter(estimate != "LP") %>% 
  ggplot(aes(neff_ratio, estimate)) +
  geom_vline(xintercept = 0.1, 
             colour = "#841F27", alpha = 0.8, 
             linetype = "dashed", size = 1) +
  geom_segment(aes(x = 0, xend = neff_ratio, 
                   y = estimate, yend = estimate), 
               colour = "#841F27", size = 0.8) +
  geom_point(shape = 21, size = 6, 
             fill = "#841F27", colour = "grey30", 
             stroke = 1) +
  labs(y = NULL, x = "Effective sample size ratio") +
  my_theme


# rhat values
rhat_plot <- rhat(brms_best) %>% 
  enframe(name = "estimate", value = "rhat") %>% 
  filter(!str_detect(estimate, 'prior')) %>% 
  mutate(estimate = c("Cooling-Cooling", 
                      "All Other", 
                      "Sigma Cooling-Cooling", 
                      "Sigma All Other", 
                      "nu", 
                      "LP")) %>% 
  ggplot(aes(rhat, estimate)) +
  geom_point(shape = 21, size = 3, 
             fill = "#841F27", colour = "grey30", 
             stroke = 1) +
  geom_vline(xintercept = 1.05,
             colour = "#841F27", alpha = 0.4, 
             linetype = "dashed") +
  labs(y = NULL, x = "Rhat value") +
  my_theme



# posterior predictive checks ---------------------------------------------


# posterior predictive check
pp_plot <- pp_check(brms_best) +
  my_theme


# set up labeller
my_labeller <- as_labeller(
  c("b_pal.intcooling_cooling" = "Cooling-Cooling", 
    "b_pal.intother" = "All Other", 
    "b_sigma_pal.intcooling_cooling" = "Sigma Cooling-Cooling", 
    "b_sigma_pal.intother" =  "Sigma All Other",
    "nu" = "nu",
    "lp__" = "lp"
  )
)

# trace plot for convergence
trace_plot <-  mcmc_trace(brms_best,
                          facet_args = list(labeller = my_labeller)) +
  my_theme 



# trace plot for convergence
trank_plot <- mcmc_rank_overlay(brms_best, 
                                facet_args = list(labeller = my_labeller)) +
  my_theme


# save images -------------------------------------------------------------

# select all plots from environment
gg_charts <- ls()[str_detect(ls(),"plot")][2:5]

all_plots <- mget(gg_charts)

# set image names
gg_names <- paste0(gg_charts, ".png")

# save all images at once
walk2(gg_names, all_plots, ~ggsave(filename = .x, plot = .y, 
                                   device = "png", 
                                   path = here("figures")))
