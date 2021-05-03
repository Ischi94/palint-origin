# load libraries
library(tidyverse)
library(here)
library(divDyn)
library(brms)
library(patchwork)

# source ggplot theme
source(here("R/functions.R"))



# load data ---------------------------------------------------------------


# final data 
load(file = here("data/final_data.RData"))

# Split short term temperature change into warming and cooling, 
# to calculate the results for each:
dat_final <- dat_final %>% 
  mutate(cooling = ifelse(change.prev < 0, change.prev, NA), 
         warming = ifelse(change.prev > 0, change.prev, NA))



# calculate onshore-offshore affinity -------------------------------------

# load occurrence data
load(file = here("data/occurrence_sqs_data.RData"))

# load keys for environmental binding
data("keys")

# calculate affinity on genus level
genus_affinity <- datsqs %>% 
  as_tibble() %>% 
  mutate(depenv = if_else(environment %in% keys$depenv$onshore, "onshore", 
                          if_else(environment %in% keys$depenv$offshore, "offshore", 
                                  "unknown"))) %>% 
  select(phylum, class, order, family, genus, depenv) %>% 
  filter(depenv != "unknown") %>% 
  mutate(depenv = as_factor(depenv)) %>% 
  count(phylum, class, order, family, genus, depenv, .drop = FALSE) %>% 
  pivot_wider(names_from = depenv, values_from = n) %>% 
  mutate(affin = if_else(offshore > onshore, "offshore", "onshore")) %>% 
  select(-c(offshore, onshore))


# load models -------------------------------------------------------------

# load fitted GLMMS for predictions
load(here("data/pal_int_models.RData"))

warm_interaction_final <- pal_int_models[[1]]
cool_interaction_final <- pal_int_models[[2]]


# split data --------------------------------------------------------------

offshore <- dat_final %>% 
  left_join(genus_affinity) %>% 
  drop_na(affin) %>% 
  filter(affin == "offshore")

onshore <- dat_final %>% 
  left_join(genus_affinity) %>% 
  drop_na(affin) %>% 
  filter(affin == "onshore")


# Make predictions based on GLMM's ----------------------------------------


### offshore

#  warming warming
ww_raw <- subset(offshore, trend.st7 >=0 & warming >= 0)
ww_pred <- predict(warm_interaction_final, newdata = ww_raw,
                   type = "response")

#  cooling warming
cw_raw <- subset(offshore, trend.st7 <=0 & warming >= 0)
cw_pred <- predict(warm_interaction_final, newdata = cw_raw,
                   type = "response")

#  warming cooling 
wc_raw <- subset(offshore, trend.st4 >=0 & cooling <= 0)
wc_pred <- predict(cool_interaction_final, newdata = wc_raw,
                   type = "response")

#  cooling cooling 
cc_raw <- subset(offshore, trend.st4 <=0 & cooling <= 0)
cc_pred <- predict(cool_interaction_final, newdata = cc_raw,
                   type = "response")

# make a dataframe with the output
offshore_prob <- tibble(ori.prob = c(cc_pred, wc_pred, cw_pred,ww_pred), 
               pal.int = c(rep("CC", length(cc_pred)),
                           rep("WC", length(wc_pred)),
                           rep("CW", length(cw_pred)),
                           rep("WW", length(ww_pred)))) %>% 
  mutate(
    # absolute values
    ori.prob = ori.prob * 100,
    # combine the rest
    pal.int = fct_collapse(pal.int, 
                           cooling_cooling = "CC", 
                           other = c("CW", "WC", "WW")))


### same for onshore

#  warming warming
ww_raw <- subset(onshore, trend.st7 >=0 & warming >= 0)
ww_pred <- predict(warm_interaction_final, newdata = ww_raw,
                   type = "response")

#  cooling warming
cw_raw <- subset(onshore, trend.st7 <=0 & warming >= 0)
cw_pred <- predict(warm_interaction_final, newdata = cw_raw,
                   type = "response")

#  warming cooling 
wc_raw <- subset(onshore, trend.st4 >=0 & cooling <= 0)
wc_pred <- predict(cool_interaction_final, newdata = wc_raw,
                   type = "response")

#  cooling cooling 
cc_raw <- subset(onshore, trend.st4 <=0 & cooling <= 0)
cc_pred <- predict(cool_interaction_final, newdata = cc_raw,
                   type = "response")

# make a dataframe with the output
onshore_prob <- tibble(ori.prob = c(cc_pred, wc_pred, cw_pred,ww_pred), 
                        pal.int = c(rep("CC", length(cc_pred)),
                                    rep("WC", length(wc_pred)),
                                    rep("CW", length(cw_pred)),
                                    rep("WW", length(ww_pred)))) %>% 
  mutate(
    # absolute values
    ori.prob = ori.prob * 100,
    # combine the rest
    pal.int = fct_collapse(pal.int, 
                           cooling_cooling = "CC", 
                           other = c("CW", "WC", "WW")))



# Bayesian Regression -----------------------------------------------------


# set monte carlo parameters
CHAINS <- 4
ITER <- 1000
WARMUP <- 500
BAYES_SEED <- 1234
options(mc.cores = parallel::detectCores())  # Use all cores


### offshore

# run the model usin brms and rcpp
brms_best_offshore <- brm(
  # we suppress the intercept by setting ori.prob ~ 0 + pal.int, 
  # brms returns coefficients for each of the groups, 
  # and these coefficients represent group means.
  bf(ori.prob ~ 0 + pal.int, sigma ~ 0 + pal.int), 
  family = student,
  data = offshore_prob,
  prior = c(
    # Set group mean prior as visualized above
    # these are fairly informative prior capturing realistic expectations
    set_prior("normal(15, 30)", class = "b", lb = 0, ub = 100),
    # Set group variance prior using cauchy(0, 1)
    set_prior("cauchy(0, 1)", class = "b", dpar = "sigma"),
    set_prior("exponential(1.0/29)", class = "nu")),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED)



# extract the posterior samples for each of the groups, subtract them from each other,
# and then calculate the credible interval
offshore_post <- posterior_samples(brms_best_offshore) %>% 
  # calculate differences
  transmute(diff_medians = b_pal.intcooling_cooling - b_pal.intother) %>% 
  as_tibble()


### same for onshore

# run the model usin brms and rcpp
brms_best_onshore <- brm(
  # we suppress the intercept by setting ori.prob ~ 0 + pal.int, 
  # brms returns coefficients for each of the groups, 
  # and these coefficients represent group means.
  bf(ori.prob ~ 0 + pal.int, sigma ~ 0 + pal.int), 
  family = student,
  data = onshore_prob,
  prior = c(
    # Set group mean prior as visualized above
    # these are fairly informative prior capturing realistic expectations
    set_prior("normal(15, 30)", class = "b", lb = 0, ub = 100),
    # Set group variance prior using cauchy(0, 1)
    set_prior("cauchy(0, 1)", class = "b", dpar = "sigma"),
    set_prior("exponential(1.0/29)", class = "nu")),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED)



# extract the posterior samples for each of the groups, subtract them from each other,
# and then calculate the credible interval
onshore_post <- posterior_samples(brms_best_onshore) %>% 
  # calculate differences
  transmute(diff_medians = b_pal.intcooling_cooling - b_pal.intother) %>% 
  as_tibble()



# combine posteriors ------------------------------------------------------

combined_post <- offshore_post %>% 
  add_column(affin = "offshore") %>% 
  full_join(onshore_post %>% 
              add_column(affin = "onshore")) 

# visualize raw posteriors
raw_post <- combined_post %>% 
  ggplot(aes(diff_medians, fill = affin)) +
  geom_density(alpha = 0.9, 
               colour = "white") +
  scale_fill_manual(values = c("#354E71", "#841F27"), 
                    labels = c("Onshore", "Offshore")) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "Difference in medians", 
       y = NULL,
       fill = "Environment") +
  my_theme +
  theme(legend.position = c(0.85, 0.9)) +
  ggtitle('A')


# calculate the contrast
contrast_post <- combined_post %>% 
  pivot_wider(names_from = affin, values_from = diff_medians, 
              values_fn = list) %>% 
  unnest(cols = c(offshore, onshore)) %>% 
  transmute(diff_cont = offshore - onshore) 

# calculate the mean and ci
ci_contrast_post <- contrast_post %>% 
  summarise(mean_val = mean(diff_cont),
            hpdi_89 = list(bayestestR::ci(diff_cont, ci = 0.89)), 
            hpdi_95 = list(bayestestR::ci(diff_cont, ci = 0.95))) %>% 
  unnest()

# visualize the contrast
raw_contrast <- contrast_post %>% 
  ggplot(aes(diff_cont)) +
  geom_density(fill = "grey60", colour = "grey50", 
               alpha = 0.8) +
  geom_linerange(aes(y = 0, x = mean_val,
                     xmin = CI_low, xmax = CI_high), 
                 size = 2, colour = "grey20",
                 data = ci_contrast_post) +
  geom_linerange(aes(y = 0, x = mean_val,
                     xmin = CI_low1, xmax = CI_high1), 
                 size = 1, colour = "grey20",
                 data = ci_contrast_post) +
  geom_point(aes(y = 0, x = mean_val), 
             size = 5, colour = "grey10",
             data = ci_contrast_post) +
  geom_point(aes(y = 0, x = mean_val), 
             size = 4, colour = "white",
             data = ci_contrast_post) +
  scale_y_continuous(breaks = NULL) +
  labs(x = "Contrast", 
       y = NULL) +
  my_theme +
  ggtitle('B')

## convert to percentage change

# get average odds
odds <- exp(summary(warm_interaction_final)$coefficients[1])
# transform to probability
prob_warm <- odds / (1 + odds)

# same for cool
odds <- exp(summary(cool_interaction_final)$coefficients[1])
prob_cool <- odds / (1 + odds)

# calculate average origination likelihood
av <- (prob_cool + prob_warm)/2

# convert credible intervals to percentage change
ci_contrast_perc <- ci_contrast_post %>% 
  mutate(across(everything(), ~(.x/(av*100))*100))


perc_contrast <- contrast_post %>% 
  mutate(diff_cont = (diff_cont/(av*100))*100) %>% 
  ggplot(aes(diff_cont)) +
  geom_density(fill = "grey60", colour = "grey50", 
               alpha = 0.8) +
  geom_linerange(aes(y = 0, x = mean_val,
                     xmin = CI_low, xmax = CI_high), 
                 size = 2, colour = "grey20",
                 data = ci_contrast_perc) +
  geom_linerange(aes(y = 0, x = mean_val,
                     xmin = CI_low1, xmax = CI_high1), 
                 size = 1, colour = "grey20",
                 data = ci_contrast_perc) +
  geom_point(aes(y = 0, x = mean_val), 
             size = 5, colour = "grey10",
             data = ci_contrast_perc) +
  geom_point(aes(y = 0, x = mean_val), 
             size = 4, colour = "white",
             data = ci_contrast_perc) +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, '%')) +
  labs(x = "Contrast as percentage change", 
       y = NULL) +
  my_theme + 
  ggtitle('C')


## combine plots
combined_environment <- raw_post & (raw_contrast / perc_contrast) +
  plot_annotation(tag_levels = 'A')

# save plot
ggsave(plot = combined_environment, 
       filename = here("figures/environmental_differences.png"))
