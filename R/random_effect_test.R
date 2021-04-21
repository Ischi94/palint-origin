# load libraries
library(tidyverse)
library(lme4)
library(geiger)
library(here)

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




