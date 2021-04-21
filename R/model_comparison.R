# loading packages
library(tidyverse) # for cleaning and visualization
library(lme4) # for GLMM's
library(geiger) # for AICw
library(here) # for project tidying

# load self-defined functions 
source(here("R/functions.R"))

# Analysis ----------------------------------------------------------------


# load data
load(here("data/final_data.RData"))

# Split short term temperature change into warming and cooling, to calculate the results for each:
dat_final$cooling<-ifelse(dat_final$change.prev<0, dat_final$change.prev, NA)
dat_final$warming<-ifelse(dat_final$change.prev>0, dat_final$change.prev, NA)


dat_final$cooler <- ifelse(dat_final$Temp.mean < dat_final$lag1, dat_final$Temp.mean, NA )
dat_final$warmer <- ifelse(dat_final$Temp.mean > dat_final$lag1, dat_final$Temp.mean, NA )

###  


# Make data frame for model output
model_comparison <- tibble(models = c("Temp.warm","Temp.warm.PI",  
                                          "Temp.cool" , "Temp.cool.PI",
                                          "Temp.lag.warm", "Temp.lag.warm.PI",
                                          "Temp.lag.cool", "Temp.lag.cool.PI",
                                          "Temp.lag.int.warm","Temp.lag.int.warm.PI",
                                          "Temp.lag.int.cool", "Temp.lag.int.cool.PI") ,
                               traditional= c(rep(c("yes", "no"), 6)),
                               type= rep(c("warming", "warming", "cooling", "cooling"), 3),
                               overdispersed=NA,
                               AIC=NA, BIC=NA)


# Temp_warm ----------------------------------------------------------------

# model taking only the short-term temperature at each stage into account/ for warming
Temp.warm <- glmer(formula = "origination ~ warmer + (1|genus)", 
                   family = "binomial", data = dat_final, nAGQ= 25)


summarise_model(Temp.warm)

# Temp_warm.PI  ------------------------------------------------------------


# model taking  both short-term and long-term temperature at each stage into account/ 
# for warming
# Iterate through each warming
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
Temp.warm.PI = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warmer+warming:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})

# select final model from the list of paleoclimate interactions
Temp.warm.PI <- select_model(Temp.warm.PI)

# summarise it and write it into model_comparison
summarise_model(Temp.warm.PI)


# Temp_cool ----------------------------------------------------------------

# model taking only the short-term temperature at each stage into account/ for cooling
Temp.cool <- glmer(formula = "origination ~ cooler + (1|genus)", 
                   family = "binomial", data = dat_final, nAGQ= 25)

# summarise it and write it into model_comparison
summarise_model(Temp.cool)


# Temp_cool.PI  ------------------------------------------------------------


# model taking  both short-term and long-term temperature at each stage into account/ for cooling
# Iterate through each cooling
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
Temp.cool.PI = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooler+cooling:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})

# select final model from the list of paleoclimate interactions
Temp.cool.PI <- select_model(Temp.cool.PI)

# summarise it and write it into model_comparison
summarise_model(Temp.cool.PI)



# Temp.lag.warm  ------------------------------------------------------------


# model taking  both short-term and long-term temperature at each stage into account/ for warming
# Iterate through each warming
vars = names(dplyr::select(dat_final, lag1:lag10)) 
Temp.lag.warm = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warmer+", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})


# select final model from the list of paleoclimate interactions
Temp.lag.warm <- select_model(Temp.lag.warm)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.warm)


# Temp.lag.warm.PI  ------------------------------------------------------------


# model taking  both short-term and long-term temperature at each stage into account/ for warming
# Iterate through each warming
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
Temp.lag.warm.PI = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warmer+lag10+warming:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})

# select final model from the list of paleoclimate interactions
Temp.lag.warm.PI <- select_model(Temp.lag.warm.PI)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.warm.PI)


# Temp.lag.cool ------------------------------------------------------------

# model taking  both short-term and long-term temperature at each stage into account/ for cooling
# Iterate through each cooling
vars = names(dplyr::select(dat_final, lag1:lag10)) 
Temp.lag.cool = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooler+", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})

# select final model from the list of paleoclimate interactions
Temp.lag.cool <- select_model(Temp.lag.cool)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.cool)


# Temp.lag.cool.PI  ------------------------------------------------------------


# model taking  both short-term and long-term temperature at each stage into account/ for cooling
# Iterate through each cooling
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
Temp.lag.cool.PI = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooler+lag10+cooling:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})


# select final model from the list of paleoclimate interactions
Temp.lag.cool.PI <- select_model(Temp.lag.cool.PI)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.cool.PI)


# Temp.lag.int.warm -----------------------------------------------------------------

# model taking paleoclimate interaction only into account/ for warming
# Iterate through each warming interaction
vars = names(dplyr::select(dat_final, lag1:lag10)) 
Temp.lag.int.warm = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warmer*", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})

# select final model from the list of paleoclimate interactions
Temp.lag.int.warm <- select_model(Temp.lag.int.warm)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.int.warm)


# Temp.lag.int.warm.PI -----------------------------------------------------------------

# model taking paleoclimate interaction only into account/ for warming
# Iterate through each warming interaction
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
Temp.lag.int.warm.PI = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~warmer*lag10+warming:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})


# select final model from the list of paleoclimate interactions
Temp.lag.int.warm.PI <- select_model(Temp.lag.int.warm.PI)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.int.warm.PI)



# Temp.lag.int.cool -----------------------------------------------------------------

# model taking paleoclimate interaction only into account/ for cooling
# Iterate through each cooling interaction
vars = names(dplyr::select(dat_final, lag1:lag10)) 
Temp.lag.int.cool = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooler*", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})


# select final model from the list of paleoclimate interactions
Temp.lag.int.cool <- select_model(Temp.lag.int.cool)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.int.cool)


# Temp.lag.int.cool.PI -----------------------------------------------------------------

# model taking paleoclimate interaction only into account/ for cooling
# Iterate through each cooling interaction
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
Temp.lag.int.cool.PI = lapply(setNames(vars, vars), function(var) {
  form = paste("origination~cooler*lag10+cooling:", var, "+(1|genus)")
  glmer(form, data=dat_final, family="binomial", nAGQ= 25)
})

# select final model from the list of paleoclimate interactions
Temp.lag.int.cool.PI <- select_model(Temp.lag.int.cool.PI)

# summarise it and write it into model_comparison
summarise_model(Temp.lag.int.cool.PI)


# save it
# save(model_comparison,file = here("data/model_comparison.RData"))

# combine all (final) models in a list and save it as RData. Make sure that you
# use the glmms and not the function output from final model, they have the same names.
# final_models <- list(Temp.warm, Temp.warm.PI,
#      Temp.cool, Temp.cool.PI,
#      Temp.lag.warm, Temp.lag.warm.PI,
#      Temp.lag.cool, Temp.lag.cool.PI,
#      Temp.lag.int.warm, Temp.lag.int.warm.PI,
#      Temp.lag.int.cool, Temp.lag.int.cool.PI)

# save it
# save(final_models,file = here("data/final_models.RData"))


# plotting ----------------------------------------------------------------


# tidy it for plotting

# tidy data
tidy_comparison <- model_comparison %>%
  # sum AIC of each type (warming, cooling) up, distinguished by whether they include
  # paleolicmate interaction or not
  group_by(type, traditional) %>% 
  summarise(meanAIC = sum(AIC), meanBIC = sum(BIC)) %>% 
  # take the average AIC for each and calculate deltaAIC
  group_by(type) %>%
  mutate(meandAIC_TR = abs(meanAIC - min(meanAIC)), 
         meandBIC_TR = abs(meanBIC - min(meanBIC))) %>% 
  # cut it up into two colums for plotting
  filter(meandAIC_TR != 0) %>% 
  mutate(meandAIC_PI = 0) %>% 
  # choose only relevant columns, 
  # meandAIC_TR = traditional models
  # meandAIC_PI = paleoclimate models
  select(type, meandAIC_TR, meandAIC_PI)
  
# plot it

# first make it long for legend plotting
tidy_comparison_long <- model_comparison %>%
  group_by(type, traditional) %>% 
  summarise(meanAIC = sum(AIC)) %>% 
  group_by(type) %>%
  mutate(meandAIC_TR = abs(meanAIC - min(meanAIC)))


model_comparison_plot <- ggplot(tidy_comparison, aes(x=meandAIC_PI, xend=meandAIC_TR, y=type)) + 
  geom_segment(aes(yend = type), colour = c("#354E71", "#841F27"),  size = 1, show.legend = FALSE) +
  geom_point(data = tidy_comparison_long, aes(meandAIC_TR, fill = traditional), 
             shape = 21, colour = "grey40", stroke = 0.3, size = 3) +
  scale_fill_manual(name = NULL, 
                     values = c("grey20", "grey60"),
                     labels = c("Palaeoclimate\nInteraction added", "Traditional Model")) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE, override.aes = list(size = c(2, 2)))) +
  xlim(700,0) +
  labs(x = expression(paste("Mean ",Delta, "AIC")), y = NULL) +
  scale_y_discrete(labels = c("Cooling", "Warming")) +
  # add annotations
  # box
  annotate(geom = "rect", xmin = 480, xmax = 190, 
           ymin = 0.41, ymax = 0.6, fill = "white") +
  # text
  annotate(geom = "text", x = 350, y = 0.5,
           colour = "grey30", label = "increasing model performance", size = 2.5) +
  # arrow
  annotate(geom = "segment", x = 150, y = 0.49,  
           xend = 90, yend = 0.49, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"), 
        legend.position = c(0.3, 0.9), 
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "grey50"), 
        legend.key = element_blank(), 
        legend.box.margin = margin(-8,-4,-5,-6), 
        legend.spacing.x = unit(0, units = "mm"), 
        legend.title = element_text(size = 5), 
        legend.text=element_text(size = 5)) 


# save it
ggsave(plot = model_comparison_plot, here("figures/model_comparison.png"),
       width = 9, height = 6, units = "cm")

ggsave(plot = model_comparison_plot, here("figures/model_comparison.pdf"),
       width = 9, height = 6, units = "cm", dpi = 500)


# per model comparison
model_comparison_grouped <- model_comparison %>% 
  # add group colum
  select(models, type, aic = AIC) %>% 
  # add group colum
  add_column(group_dummy = rep(c("Temp", "Temp+Lag", "Temp+Lag+Int"), each = 4)) %>%
  unite("models", c(group_dummy, type), remove = FALSE) %>% 
  mutate(with_pi = lead(aic)) %>%
  filter(row_number() %% 2 == 1) %>% 
  select(models, type, traditional = aic, with_pi) %>% 
  pivot_longer(cols = c(traditional, with_pi), values_to = "aic", names_to = "aic_type") %>% 
  mutate(dummy_col = rep(c(1, 1.5, 2, 2.5, 3, 3.5), each = 2)) %>% 
  ggplot(aes(aic, dummy_col)) +
  geom_point(aes(group = type, fill = aic_type), 
             shape = 21, colour = "grey40", size = 3, stroke = 0.4) +
  geom_line(aes(group = models, colour = type), size = 0.8) +
  geom_hline(yintercept = c(1.75, 2.75), colour = "grey") +
  scale_fill_manual(values = c("grey60", "grey20"), 
                    name = "Model", 
                    labels = c("Traditional", "Palaeoclimate\nInteraction")) +
  scale_color_manual(values = c("#354E71", "#841F27"), name = "Short-term\nChange") +
  scale_y_continuous(breaks = c(1.25, 2.25, 3.25), 
                     labels = c("Temperature", "Temperature\n+\nLag", "Temperature\n+\nLag\n+\nInteraction")) +
  labs(y = NULL, x = "AIC") +
  # add annotations
  # box
  annotate(geom = "rect", xmin = 18000, xmax = 19500, 
           ymin = 3.4, ymax = 3.5, fill = "white") +
  # text
  annotate(geom = "text", x = 19150, y = 3.457,
           colour = "grey30", label = "increasing model performance", size = 3.5) +
  # arrow
  annotate(geom = "segment", x = 18550, y = 3.457,  
           xend = 18200, yend = 3.457, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  # box
  annotate(geom = "rect", xmin = 19700, xmax = 19550, 
           ymin = 1.9, ymax = 3.4, fill = "white") +
  # text
  annotate(geom = "text", x = 19650, y = 2.9, angle = 270, 
           colour = "grey30", label = "decreasing model complexity", size = 3.5) +
  # arrow
  annotate(geom = "segment", x = 19650, y = 2.38,  
           xend = 19650, yend = 2.15, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(family = "sans"), 
        legend.position = c(0.8, 0.7), 
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_blank(),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8)
  )
  

# save it
ggsave(plot = model_comparison_grouped, here("figures/model_comparison_grouped.png"))
