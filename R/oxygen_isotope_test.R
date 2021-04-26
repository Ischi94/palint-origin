# load libraries
library(divDyn)
library(tidyverse)
library(here)
library(lme4)
library(geiger)

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
prob <- tibble(ori.prob = c(cc_pred, wc_pred, cw_pred,ww_pred), 
               pal.int = c(rep("WW", length(cc_pred)),
                           rep("WC", length(wc_pred)),
                           rep("CW", length(cw_pred)),
                           rep("CC", length(ww_pred))))


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


# visualize ---------------------------------------------------------------

ggplot(prob, aes(x=pal.int, y=ori.prob, fill=pal.int)) + 
  geom_hline(yintercept = av)+
  geom_violin() +
  scale_fill_manual(values = c("#354E71", "#354E71","#841F27", "#841F27"))+
  scale_y_continuous(name = "Origination Probability", limits=c(0, 0.3), 
                     breaks = seq(0, 0.3, by= 0.05), 
                     labels = scales::percent_format(accuracy = 1)) +
  xlab(NULL) +
  my_theme +
  theme(legend.position = "none", 
        axis.ticks.length.x.bottom = unit(0.25, "cm")) +
  # add half a violin to visualise palaeoclimate interactions
  ggnewscale::new_scale_fill() +
  see::geom_violinhalf(data = filter(prob, pal.int == "CW" | pal.int == "WC"), 
                       aes(colour = pal.int, fill = pal.int)) +
  scale_fill_manual(values = c("#841F27", "#354E71")) +
  scale_colour_manual(values=c("#841F27", "#354E71")) +
  # add grey lines to visualise medians per group 
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = c(0.915, 0.35, 0.53, 0.42), 
               colour = "grey60") +
  # add outer layer
  geom_violin(fill = NA) +
  # add arrows
  annotate(geom = "curve", x = 4, y = 0.22,  # overall
           xend = 4.45, yend = 0.138, 
           curvature = -.325, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  annotate(geom = "curve", x = 1.25, y = 0.055, # per group
           xend = 0.68, yend = 0.13, 
           curvature = -.5, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  # add lines
  annotate(geom = "segment", x = c(4, 1.25), y = c(0.22, 0.055), # overall
           xend = c(4.55, 0.575), yend = c(0.22, 0.055), colour = "grey40") +
  # add background box
  annotate(geom = "rect", xmin = 4, xmax = 4.5, # overall
           ymin = 0.235, ymax = 0.26, fill = "white") +
  # add text
  annotate(geom = "text", x = 4.275, y = c(0.25, 0.235), # overall
           colour = "grey30", label = c("Overall", "median"), size = 3.5) +
  annotate(geom = "text", x = 0.9, y = c(0.045, 0.03, 0.015), # per group
           colour = "grey30", label = c("Median", "per", "group"), size = 3.5) +
  scale_x_discrete(labels = c("Cooling-Cooling", "Cooling-Warming", 
                              "Warming-Cooling", "Warming-Warming"), 
                   guide = guide_axis(n.dodge=2))

