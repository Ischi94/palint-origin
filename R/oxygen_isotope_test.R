# load libraries
library(divDyn)
library(tidyverse)
library(here)
library(lme4)
library(geiger)
library(brms)

# load self-defined functions 
source(here("R/functions.R"))



# load data ---------------------------------------------------------------

# fossil data already cleaned with Data_preparation.R
load(here("data/cleaned_data.RData"))


# alternative paleotemperature data 
isotemp <- read_csv(here("data/song_et_al.csv"))
         
# stage data
data("stages")


# bin temperature data ----------------------------------------------------

# build stage data which contains the same million year id as song et al. for binning
stages_seq <- stages %>% 
  as_tibble() %>% 
  select(bottom, top, stg) %>% 
  rowwise() %>% 
  mutate(age_seq = list(seq(bottom, top, -0.001))) %>% 
  ungroup() %>% 
  unnest(age_seq) %>% 
  select(age_seq, stg) %>% 
  ungroup()

# bin temperature data
isotemp <- isotemp %>% 
  select(age_seq = age_ma, 
         Temp = mean_temperature_C) %>% 
  right_join(stages_seq) %>% 
  drop_na(Temp) %>% 
  group_by(age_seq) %>% 
  top_n(n = -1) %>% 
  ungroup() %>% 
  select(stg, Temp)




# calculate trends --------------------------------------------------------


# calculate stage difference in temperature (short-term change)

# prepare zaffos data
isotemp2 <- isotemp %>% 
  # combine fragmentation_index in list column
  group_by(stg) %>% 
  nest() %>%
  # apply function to each column
  mutate(Temp = map(data, "Temp")) %>% 
  ungroup() %>% 
  select(-data) %>% 
  # calculate short-term change in fragmentation index
  mutate(change.prev = map_dbl(unique(isotemp$stg), short_term_temp, j = 1))


# now move gradually back in time for each trend
isotemp2 <- isotemp2 %>% 
  mutate(trend.st1 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 2),
         trend.st2 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 3), 
         trend.st3 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 4),
         trend.st4 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 5), 
         trend.st5 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 6), 
         trend.st6 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 7), 
         trend.st7 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 8), 
         trend.st8 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 9), 
         trend.st9 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 10), 
         trend.st10 = map_dbl(unique(isotemp2$stg), long_term_temp, j = 11)) 


# remove everything older than the ordovician
isotemp2 <- isotemp2 %>% 
  filter(stg >= 14) %>% 
  select(bins = stg, change.prev:trend.st10)



# 6 Calculate multiple long-term temperature trends -------------------------


# build dummies for calculation
dumbo <- dat_range
dat_safe <- dat_range[,c("genus","FAD","LAD")]

# add bins with 0 and fill in with 0 for survival and 1 origination
namevector <- as.character(c(1:94))
dat_safe[ , namevector] <- NA
dat_safe <- dat_safe[, c(4:length(dat_safe[0,]))]

for (i in 1:length(dat_safe[,1])) {
  dat_safe[i,dumbo[i,"FAD"]] <- 1 
  dat_safe[i,dumbo[i,"LAD"]] <- 0
  ifelse(dumbo[i,"FAD"]!= dumbo[i,"LAD"]-1 & dumbo[i,"FAD"]!= dumbo[i,"LAD"],
         dat_safe[i,(dumbo[i,"FAD"]+1):(dumbo[i,"LAD"])] <- 0, NA)
}

# bind it again 
rownames(dat_safe) <-  make.names(dumbo[,"genus"], unique=TRUE)


# Reorganise ranges through time data so we can bind it to temperature data
dat_safe <- dat_safe %>%
  as_tibble() %>%
  mutate(genus = rownames(dat_safe)) %>%
  group_by(genus) %>%
  gather(-genus, key="bins", value="origination", na.rm = T, factor_key = T) %>%
  arrange(desc(bins), .by_group = T) %>% 
  mutate(bins = as.double(as.numeric(bins)))


# Now bind the two
dat_temp <- full_join(dat_safe, isotemp2)

#add taxonomical levels using the dummy and environmental preference
dumbo <- dumbo[, 1:5]

# Now bind the two
dat_temp <- full_join(dat_temp, dumbo)

# order it properly
dat_temp <- dat_temp %>%
  dplyr::select(phylum:family, genus:trend.st10) 


# Calculate GLMM's ----------------------------------------------------------------


# Split short term temperature change into warming and cooling, 
# to calculate the results for each:
dat_final <- dat_temp %>% 
  mutate(cooling = ifelse(change.prev < 0, change.prev, NA), 
         warming = ifelse(change.prev > 0, change.prev, NA)) %>% 
  ungroup()


# for warming
# model taking  both short-term and long-term temperature at each stage into account/ 
# for warming
# Iterate through each warming
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
warm_interaction = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warming:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})


# Make data frame for model output
warm_interaction_df <- model_df(warm_interaction)

# choose final model
warm_interaction_final <- warm_interaction[[which(warm_interaction_df$dAIC==0)]]

# summary
sum_warm_interaction <- summary(warm_interaction_final) 

# for cooling
# model taking both short-term and long-term temperature at each stage into account/
# for warming
# Iterate through each cooling
cool_interaction = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooling:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})


# Make data frame for model output
cool_interaction_df <- model_df(cool_interaction)

# choose final model
cool_interaction_final <- cool_interaction[[which(cool_interaction_df$dAIC==0)]]

# summary
sum_cool_interaction <- summary(cool_interaction_final) 

# Make predictions based on GLMM's ----------------------------------------

# make informed predictions based on a subset (palaeoclimate interaction)
# type = response gives us probability instead of log Odds

#  warming warming
ww_raw <- subset(dat_final, trend.st1 >=0 & warming >= 0)
ww_pred <- predict(warm_interaction_final, newdata = ww_raw,
                   type = "response")

#  cooling warming
cw_raw <- subset(dat_final, trend.st1 <=0 & warming >= 0)
cw_pred <- predict(warm_interaction_final, newdata = cw_raw,
                   type = "response")

#  warming cooling 
wc_raw <- subset(dat_final, trend.st6 >=0 & cooling <= 0)
wc_pred <- predict(cool_interaction_final, newdata = wc_raw,
                   type = "response")

#  cooling cooling 
cc_raw <- subset(dat_final, trend.st6 <=0 & cooling <= 0)
cc_pred <- predict(cool_interaction_final, newdata = cc_raw,
                   type = "response")

# make a dataframe with the output
prob_comparison <- tibble(ori.prob = c(cc_pred, wc_pred, cw_pred,ww_pred), 
               pal.int = c(rep("WW", length(ww_pred)),
                           rep("WC", length(wc_pred)),
                           rep("CW", length(cw_pred)),
                           rep("CC", length(cc_pred)))) %>%
  mutate(
    # absolute values
    ori.prob = ori.prob * 100,
    # combine the rest
    pal.int = fct_collapse(
      pal.int,
      cooling_cooling = "CC",
      other = c("CW", "WC", "WW")))



# average response --------------------------------------------------------


# the predictions are now in percentage (probability), but the intercept, to which we 
# want to compare our predictions with, is still in log Odds.
# transform logit intercept into probability
# for warm
odds <- exp(sum_warm_interaction$coefficients[1])
prob_warm <- odds / (1 + odds)

# for cool
odds <- exp(sum_cool_interaction$coefficients[1])
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
  add_column(type = rep(c("Difference in medians", 
                          "Percentage Change"), each = 3)) %>% 
  mutate(name = rep(c("estimate", "conf.low", "conf.high"), 2)) %>% 
  pivot_wider(names_from = "name") %>% 
  add_column(isotope_data = "Song et al. 2019")



# compare results ---------------------------------------------------------

# load data fitted with Veizer and Prokoph isotope data
load(here("data/effect_plot_data.RData"))

# bring results in similar format
actual_results <- effect_plot_data %>% 
  map_dfr(filter, method == "best") %>%  
  mutate(type = c("Difference in medians", "Percentage Change")) %>% 
  add_column(isotope_data = "Veizer & Prokoph 2015")

# merge data
best_results %>% 
  full_join(actual_results) %>% 
  select(-method) %>% 
  write_csv(here("data/oxygen_isotope_test.csv"))
