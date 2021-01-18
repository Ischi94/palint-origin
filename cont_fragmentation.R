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

#  increase increase
ww_raw <- subset(dat_final_frag, trend.st7 >=0 & increase >= 0)
ww_pred <- predict(incr_interaction_final, newdata = ww_raw,
                   type = "response")

#  decrease increase
cw_raw <- subset(dat_final_frag, trend.st7 <=0 & increase >= 0)
cw_pred <- predict(incr_interaction_final, newdata = cw_raw,
                   type = "response")

#  increase decrease 
wc_raw <- subset(dat_final_frag, trend.st6 >=0 & decrease <= 0)
wc_pred <- predict(decr_interaction_final, newdata = wc_raw,
                   type = "response")

#  decrease decrease 
cc_raw <- subset(dat_final_frag, trend.st6 <=0 & decrease <= 0)
cc_pred <- predict(decr_interaction_final, newdata = cc_raw,
                   type = "response")

# make a dataframe with the output
prob_fragm <- tibble(ori.prob = c(cc_pred, wc_pred, cw_pred,ww_pred), 
               fragm.int = c(rep("Decrease-Decrease", length(cc_pred)),
                           rep("Increase-Decrease", length(wc_pred)),
                           rep("Decrease-Increase", length(cw_pred)),
                           rep("Increase-Increase", length(ww_pred))))

# save it
# save(prob_fragm, file = here("data/fragmentation_plot_data.RData"))

# the predictions are now in percentage (probability), but the intercept, to which we 
# want to compare our predictions with, is still in log Odds.
# transform logit intercept into probability
# for increase
odds <- exp(sum_incr_interaction$coefficients[1])
prob_incr <- odds / (1 + odds)

# for decrease
odds <- exp(sum_decr_interaction$coefficients[1])
prob_decr <- odds / (1 + odds)


# mean
av <- mean(c(prob_incr, prob_decr))
# 0.134989

# calculate summaries
prob_fragm_sum <- prob_fragm %>% 
  group_by(fragm.int) %>% 
  summarise(ci = list(mean_cl_boot(ori.prob) %>% 
                        rename(mean=y, lwr=ymin, upr=ymax))) %>% 
  unnest(cols = c(ci)) 




# plot data ---------------------------------------------------------------

cont_fragm_plot <- ggplot(aes(x = mean, y = fragm.int), data = prob_fragm_sum) +
  geom_vline(xintercept = av, colour = "grey5") +
  geom_linerange(aes(xmin = lwr, xmax = upr), size = 1.5, colour = "grey50") +
  geom_point(shape = 21, size = 3, fill = "grey50", colour = "grey10") +
  coord_cartesian(xlim = c(0.115, 0.141)) +
  labs(x = "Origination probability", y = NULL) +
  scale_x_continuous(labels = function(x) paste0(x*100, '%')) +
  annotate(geom = "curve", x = 0.1326, y = 1.25,
           xend = 0.1347, yend = 1.54, 
           curvature = -.35, arrow = arrow(length = unit(1.5, "mm")), 
           colour = "grey30", size = 0.4) +
  annotate(geom = "segment", x = 0.1326, y = 1.25,
           xend = 0.1326, yend = 1.58, colour = "grey30", size = 0.4) +
  annotate(geom = "rect", xmin = 0.126, xmax = 0.1323,
           ymin = 1.225, ymax = 1.66, fill = "white") +
  annotate(geom = "text", x = 0.1307, y = c(1.51, 1.32),
           colour = "grey30", label = c("Overall", "mean"), size = 3) +
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank())

ggsave(plot = cont_fragm_plot, filename = here("figures/cont_fragm.png"), 
       width = 12.7, height = 9, units = "cm")

ggsave(plot = cont_fragm_plot, filename = here("figures/cont_fragm.pdf"), 
       width = 12.7, height = 9, units = "cm", dpi = 500)
