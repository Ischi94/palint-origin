library(tidyverse)
library(chronosphere)
library(divDyn)
library(here)
library(lme4)
library(geiger) 

# select is conflicted, define to use dplyr's select
select <- dplyr::select

# load self-defined functions 
source("functions.R")


# load data ---------------------------------------------------------------

# load zaffos fragmentation index
zaffos <- fetch(dat = "som", var = "zaffos-fragment") %>% 
  as_tibble() %>% 
  select(age_seq = X, fragmentation_index = fragmentation.index)

# load stage data for binning
data(stages)


# load origination data
load(here("data/final_data.RData"))



# bin fragmentation index -------------------------------------------------


# build stage data which contains the same million year id as staffos for binning
stages_seq <- stages %>% 
  as_tibble() %>% 
  select(bottom, top, stg) %>% 
  rowwise() %>% 
  mutate(age_seq = list(seq(bottom, top, -0.001))) %>% 
  ungroup() %>% 
  unnest(age_seq) %>% 
  select(age_seq, stg) %>% 
  ungroup()
  
# bin fragmentation index
zaffos_binned <- zaffos %>% 
  right_join(stages_seq) %>% 
  drop_na(fragmentation_index) %>% 
  group_by(age_seq) %>% 
  top_n(n = -1) %>% 
  ungroup() %>% 
  select(stg, fragmentation_index)


# calculate stage difference in fragmentation (short-term change)

# prepare zaffos data
zaffos_binned2 <- zaffos_binned %>% 
  # combine fragmentation_index in list column
  group_by(stg) %>% 
  nest() %>%
  # apply function to each column
  mutate(fragmentation_index = map(data, "fragmentation_index")) %>% 
  ungroup() %>% 
  # add dummy column 
  add_column(change.prev = double(length = length(unique(zaffos_binned$stg))))


# calculate short-term change in fragmentation index
zaffos_binned2 <- zaffos_binned2 %>% 
  select(-data) %>% 
  mutate(change.prev = map_dbl(unique(zaffos_binned2$stg), short_term, j = 1))


# now move gradually back in time for each trend
zaffos_trends <- zaffos_binned2 %>% 
  mutate(trend.st1 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 2),
    trend.st2 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 3), 
    trend.st3 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 4),
    trend.st4 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 5), 
    trend.st5 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 6), 
    trend.st6 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 7), 
    trend.st7 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 8), 
    trend.st8 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 9), 
    trend.st9 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 10), 
    trend.st10 = map_dbl(unique(zaffos_binned2$stg), long_term, j = 11)) 



# combine data ------------------------------------------------------------



# remove temperature trends and add fragmentation index trends to the origination data
dat_final_frag <- dat_final %>% 
  as_tibble() %>% 
  select(phylum:origination, stg = bins) %>% 
  mutate(stg = as.numeric(stg)) %>% 
  left_join(zaffos_trends)

# save data
# save(dat_final_frag, file = here("data/final_fragmentation_data.RData"))


# Calculate GLMM's ----------------------------------------------------------------

# Split short term fragmentation change into decrease and increase, 
# to calculate the results for each:
dat_final_frag <- dat_final_frag %>% 
  mutate(decrease = if_else(change.prev < 0, change.prev, NA_real_), 
         increase = if_else(change.prev > 0, change.prev, NA_real_))

# for increase
# model taking  both short-term and long-term fragmentation at each stage into account/ 
# Iterate through each increase
vars = names(dplyr::select(dat_final_frag, trend.st1:trend.st10)) 
incr_interaction = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~increase:", var, "+(1|genus)")
  glmer(form, 
        data = dat_final_frag, 
        family = "binomial",
        nAGQ = 25)
})


# Make data frame for model output
incr_interaction_df <- model_df(incr_interaction)

# choose final model
incr_interaction_final <- incr_interaction[[which(incr_interaction_df$dAIC==0)]]

# summary
sum_incr_interaction <- summary(incr_interaction_final) 

# for decrease
# model taking  both short-term and long-term fragmentation at each stage into account/ 
# Iterate through each decrease
decr_interaction = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~decrease:", var, "+(1|genus)")
  glmer(form, 
        data = dat_final_frag, 
        family = "binomial", 
        nAGQ= 25)
})


# Make data frame for model output
decr_interaction_df <- model_df(decr_interaction)

# choose final model
decr_interaction_final <- decr_interaction[[which(decr_interaction_df$dAIC==0)]]

# summary
sum_decr_interaction <- summary(decr_interaction_final) 


# save it as a list
# fragm_int_df <- list(incr_interaction_df, decr_interaction_df)
# save(fragm_int_df, file = here("data/fragm_int_df.RData"))
# fragm_int_models <- list(incr_interaction_final, decr_interaction_final)
# save(fragm_int_models, file = here("data/fragm_int_models.RData"))


# Make predictions based on GLMM's ----------------------------------------

# make informed predictions based on a subset (palaeoclimate interaction)
# type = response gives us probability instead of log Odds

#  warming warming
ww_raw <- subset(dat_final, trend.st7 >=0 & warming >= 0)
ww_pred <- predict(warm_interaction_final, newdata = ww_raw,
                   type = "response")

#  cooling warming
cw_raw <- subset(dat_final, trend.st7 <=0 & warming >= 0)
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
prob <- tibble(ori.prob = c(cc_pred, wc_pred, cw_pred,ww_pred), 
               pal.int = c(rep("CC", length(cc_pred)),
                           rep("WC", length(wc_pred)),
                           rep("CW", length(cw_pred)),
                           rep("WW", length(ww_pred))))

# save it
# save(prob, file = here("data/violin_plot_data.RData"))

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


# define theme
my_theme <- theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                  panel.grid.major.y=element_line(colour = "grey", linetype = "dotted"),
                  text = element_text(family = "sans"), 
                  panel.grid.major.x = element_blank())