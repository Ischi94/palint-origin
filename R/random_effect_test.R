# load libraries
library(tidyverse)
library(lme4)
library(geiger)
library(here)
library(brms)
library(patchwork)

# load self-defined functions 
source(here("R/functions.R"))


# load and prepare data ---------------------------------------------------

# load data
load(here("data/final_data.RData"))

# add constant random effect, see N. W. Galwey in 'Introduction to Mixed
# Modelling: Beyond Regression and Analysis of Variance' on pages 213-214
dat_final <- dat_final %>%
  add_column(const_eff = 1)

# Split short term temperature change into warming and cooling, 
# to calculate the results for each:
dat_final$cooling<-ifelse(dat_final$change.prev<0, dat_final$change.prev, NA)
dat_final$warming<-ifelse(dat_final$change.prev>0, dat_final$change.prev, NA)


# fit models with constant second random effect ---------------------------


# Iterate through each warming
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
warm_interaction1 = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warming:", var, "+(1|genus)", "+(1|const_eff)")
  glmer(form, data=dat_final, family="binomial", 
        control=glmerControl(check.nlev.gtr.1="ignore"))
})


# Make data frame for model output
warm_interaction_df1 <- model_df(warm_interaction1)

# choose final model
warm_interaction_final1 <- warm_interaction1[[which(warm_interaction_df1$dAIC==0)]]



# Iterate through each cooling
cool_interaction1 = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooling:", var, "+(1|genus)", "+(1|const_eff)")
  glmer(form, data=dat_final, family="binomial", 
        control=glmerControl(check.nlev.gtr.1="ignore"))
})


# Make data frame for model output
cool_interaction_df1 <- model_df(cool_interaction1)

# choose final model
cool_interaction_final1 <- cool_interaction1[[which(cool_interaction_df1$dAIC==0)]]



# fit models with actual random effects --------------------------------------



# model taking  both short-term and long-term temperature at each stage into account/ 
# for warming
# Iterate through each warming
warm_interaction2 = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warming:", var, "+(1|genus)", "+(1|bins)")
  glmer(form, data=dat_final, family="binomial")
})


# Make data frame for model output
warm_interaction_df2 <- model_df(warm_interaction2)

# choose final model
warm_interaction_final2 <- warm_interaction2[[which(warm_interaction_df2$dAIC==0)]]



# Iterate through each cooling
cool_interaction2 = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooling:", var, "+(1|genus)", "+(1|bins)")
  glmer(form, data=dat_final, family="binomial")
})


# Make data frame for model output
cool_interaction_df2 <- model_df(cool_interaction2)

# choose final model
cool_interaction_final2 <- cool_interaction2[[which(cool_interaction_df2$dAIC==0)]]



# compare models via information criterion --------------------------------

# set up names for table
model_name <- c("Added effect", "Null model")

# for warming
comparison_warm <- anova(warm_interaction_final1, warm_interaction_final2) %>% 
  as_tibble() %>% 
  add_column(Type = "Warming", 
             Model = model_name, 
             .before = "npar") %>% 
  select(Type, Model, AIC:logLik)
  
  
# for cooling
comparison_cool <- anova(cool_interaction_final1, cool_interaction_final2) %>% 
  as_tibble() %>% 
  add_column(Type = "Cooling",
             Model = model_name, 
             .before = "npar") %>% 
  select(Type, Model, AIC:logLik)


# combine to one tibble and save as csv
comparison_warm %>% 
  full_join(comparison_cool) %>% 
  write_csv(file = here("data/random_effect_test.csv"))



# Make predictions based on GLMM's ----------------------------------------

# make informed predictions based on a subset (palaeoclimate interaction)
# type = response gives us probability instead of log Odds

#  warming warming
ww_raw <- subset(dat_final, trend.st8 >=0 & warming >= 0)
ww_pred <- predict(warm_interaction_final2, newdata = ww_raw,
                   type = "response")

#  cooling warming
cw_raw <- subset(dat_final, trend.st8 <=0 & warming >= 0)
cw_pred <- predict(warm_interaction_final2, newdata = cw_raw,
                   type = "response")

#  warming cooling 
wc_raw <- subset(dat_final, trend.st6 >=0 & cooling <= 0)
wc_pred <- predict(cool_interaction_final2, newdata = wc_raw,
                   type = "response")

#  cooling cooling 
cc_raw <- subset(dat_final, trend.st6 <=0 & cooling <= 0)
cc_pred <- predict(cool_interaction_final2, newdata = cc_raw,
                   type = "response")

# make a dataframe with the output
prob_comparison <- tibble(
  ori.prob = c(cc_pred, wc_pred, cw_pred, ww_pred),
  pal.int = c(
    rep("CC", length(cc_pred)),
    rep("WC", length(wc_pred)),
    rep("CW", length(cw_pred)),
    rep("WW", length(ww_pred)))) %>%
  mutate(
    # absolute values
    ori.prob = ori.prob * 100,
    # combine the rest
    pal.int = fct_collapse(
      pal.int,
      cooling_cooling = "CC",
      other = c("CW", "WC", "WW")))


# calculate average origination response ----------------------------------

# for warm
odds <- exp(summary(warm_interaction_final2)$coefficients[1])
prob_warm <- odds / (1 + odds)

# for cool
odds <- exp(summary(cool_interaction_final2)$coefficients[1])
prob_cool <- odds / (1 + odds)

# mean
av <- (prob_cool + prob_warm)/2


# bayesian estimation -----------------------------------------------------

# set monte carlo parameters
CHAINS <- 4
ITER <- 1000
WARMUP <- 500
BAYES_SEED <- 1234
options(mc.cores = parallel::detectCores())  # Use all cores

# run the model usin brms and rcpp
brms_best <- brm(
  # we suppress the intercept by setting ori.prob ~ 0 + pal.int, 
  # brms returns coefficients for each of the groups, 
  # and these coefficients represent group means.
  bf(ori.prob ~ 0 + pal.int, sigma ~ 0 + pal.int), 
  family = student,
  data = prob_comparison,
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
brms_best_post <- posterior_samples(brms_best) %>% 
  # Rescale sigmas (sigma terms are on a log scale, 
  # so we need to exponentiate them back to the scale of the data)
  mutate_at(vars(contains("sigma")), exp) %>% 
  # we need to log nu
  mutate(nu = log10(nu)) %>% 
  # calculate differences
  mutate(diff_medians = b_pal.intcooling_cooling - b_pal.intother,
         diff_sigma = b_sigma_pal.intcooling_cooling - b_sigma_pal.intother)


# we can use tidyMCMC from broom to calculate the difference and 
# create confidence intervals
brms_best_tidy <- 
  broom::tidy(brms_best_post, conf.int = TRUE, conf.level = 0.89, 
              estimate.method = "median", conf.method = "HPDinterval")

# calculate effect
best_results <- brms_best_tidy %>% 
  filter(column == "diff_medians") %>% 
  transmute(estimate = median, conf.low = min, conf.high = max) %>% 
  # percentage change
  mutate(across(everything(), 
                ~./(av*100)*100, 
                .names = "{col}.perc")) %>% 
  pivot_longer(cols = everything()) %>% 
  add_column(type = rep(c("raw", "percentage"), each = 3)) %>% 
  mutate(name = rep(c("estimate", "conf.low", "conf.high"), 2)) %>% 
  pivot_wider(names_from = "name") %>% 
  add_column(random_effect = "two")



# comparison with actual data ---------------------------------------------

load(here("data/effect_plot_data.RData"))

# bring results in similar format
actual_results <- effect_plot_data %>% 
  map_dfr(filter, method == "best") %>%  
  mutate(type = c("raw", "percentage")) %>% 
  add_column(random_effect = "one")

# merge data
raw_diff <- best_results %>% 
  full_join(actual_results) %>% 
  # visualize difference in medians
  filter(type == "raw") %>% 
  ggplot(aes(estimate, 
             xmin = conf.low,
             xmax = conf.high, 
             random_effect)) +
  geom_pointrange(size = 1.3, colour = "grey30") +
  geom_point(size = 3.5, colour = "white") +
  labs(x = "Difference in medians", 
       y = NULL) +
  scale_y_discrete(labels = c("Reported results", "Additional random effect")) +
  my_theme

  

# merge data
perc_diff <- best_results %>% 
  full_join(actual_results) %>% 
  # visualize difference in medians
  filter(type == "percentage") %>% 
  ggplot(aes(estimate, 
             xmin = conf.low,
             xmax = conf.high, 
             random_effect)) +
  geom_pointrange(size = 1.3, colour = "grey30") +
  geom_point(size = 3.5, colour = "white") +
  labs(x = "Percentage change", 
       y = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, '%')) +
  scale_y_discrete(labels = c("Reported results", "Additional random effect")) +
  my_theme

# combine plots
combined_plot <- raw_diff / perc_diff +
  plot_annotation(tag_levels = 'A')

# save it
ggsave(plot = combined_plot, filename = here("figures/second_effect_test.png"), 
       units = "cm")
