# loading packages
library(tidyverse) # for visualization
library(lme4) # for GLMM's
library(geiger) # for AICw
library(here) # for project tidiness
library(patchwork) # for combining of plots
library(infer) # for bootstrapping
library(brms) # for Stan-based models with standard R syntax (Bayesian Regression)


# Functions ---------------------------------------------------------------

# take a model and produce a Data frame (tibble for nicer printing) with model output
# and AIC, BIC, deltaAIC, AICweights and deltaBIC
model_df <- function(my_model) {
  # make data frame for model output
  df <-
    tibble(
      model = names(dplyr::select(dat_final, trend.st1:trend.st10)),
      intercept = character(10),
      interaction = character(10),
      AIC = numeric(10),
      BIC = numeric(10)
    )
  # Run loop to fill interaction_warm (coefficients and p-values)
  for (i in df$model) {
    sum <- summary(my_model[[i]])
    df[df$model == i, "intercept"] <-
      paste(
        round(sum$coefficients[1, 1], 2),
        sep = " ",
        "+-",
        round(sum$coefficients[1, 2], 2),
        ifelse(
          sum$coefficients[1, 4] < 0.001,
          "***",
          ifelse(
            sum$coefficients[1, 4] < 0.01,
            "**",
            ifelse(sum$coefficients[1, 4] < 0.05, "*", "")
          )
        )
      )
    df[df$model == i, "interaction"] <-
      paste(
        round(sum$coefficients[2, 1], 2),
        sep = " ",
        "+-",
        round(sum$coefficients[2, 2], 2),
        ifelse(
          sum$coefficients[2, 4] < 0.001,
          "***",
          ifelse(
            sum$coefficients[2, 4] < 0.01,
            "**",
            ifelse(sum$coefficients[2, 4] < 0.05, "*", "")
          )
        )
      )
    df[df$model == i, "AIC"] <- as.numeric(round(sum$AICtab[[1]], 1))
    df[df$model == i, "BIC"] <- as.numeric(round(sum$AICtab[[2]], 1))
  }
  
  # Add column weith AIC weights and deltaBIC
  df$dAIC <- as.numeric(round(aicw(df$AIC)$delta, 1))
  df$AICweights <- as.numeric(signif(aicw(df$AIC)$w, 3))
  df$dBIC <- as.numeric(round(aicw(df$BIC)$delta, 1))
  df
}

# take the predictions for a palaeoclimate interaction and produce a ggplot to
# (visually) check for normality
my_qqplot <- function(PI, title) {
  ggpubr::ggqqplot(prob$ori.prob[prob$pal.int == PI], title = title, ggtheme = my_theme)
}

# take the predictions and test whether they significantly are above or below the
# baseline. This is a one sided Wilcoxon rank sum test 
# (equivalent to the Mann-Whitney test)
my_wilcoxtest <- function(PI, alternative) {
  wilcox.test(prob$ori.prob[prob$pal.int == PI],
              overall,
              paired = FALSE,
              alternative = alternative)
}

# Calculate GLMM's ----------------------------------------------------------------


# load data
load(here("data/final_data.RData"))

# Split short term temperature change into warming and cooling, 
# to calculate the results for each:
dat_final$cooling<-ifelse(dat_final$change.prev<0, dat_final$change.prev, NA)
dat_final$warming<-ifelse(dat_final$change.prev>0, dat_final$change.prev, NA)


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
# Iterate through each warming
vars = names(dplyr::select(dat_final, trend.st1:trend.st10)) 
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


# save it as a list
# pal_int_df <- list(warm_interaction_df, cool_interaction_df)
# save(pal_int_df, file = here("data/pal_int_df.RData"))
# pal_int_models <- list(warm_interaction_final, cool_interaction_final)
# save(pal_int_models, file = here("data/pal_int_models.RData"))


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

# produce a violin plot
# black line shows overall mean
violin <- ggplot(prob, aes(x=pal.int, y=ori.prob, fill=pal.int)) + 
  geom_hline(yintercept = av)+
  geom_violin() +
  scale_fill_manual(values = c("#354E71", "#354E71","#841F27", "#841F27"))+
  scale_y_continuous(name = "Origination Probability", limits=c(0, 0.3), 
                     breaks = seq(0, 0.3, by= 0.05), 
                     labels = scales::percent_format(accuracy = 1)) +
  xlab(NULL) +
  my_theme+
  theme(legend.position = "none", 
        axis.ticks.length.x.bottom = unit(0.25, "cm")) +
  # add half a violin to visualise palaeoclimate interactions
  ggnewscale::new_scale_fill() +
  see::geom_violinhalf(data = filter(prob, pal.int == "CW" | pal.int == "WC"), 
                       aes(colour = pal.int, fill = pal.int)) +
  scale_fill_manual(values = c("#841F27", "#354E71")) +
  scale_colour_manual(values=c("#841F27", "#354E71")) +
  # add grey lines to visualise means per group 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = c(0.7, 0.37, 0.58, 0.52), 
               colour = "grey60") +
  # add outer layer
  geom_violin(fill = NA) +
  # add arrows
  annotate(geom = "curve", x = 4, y = 0.22,  # overall
           xend = 4.45, yend = 0.138, 
           curvature = -.325, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  annotate(geom = "curve", x = 1.05, y = 0.055, # per group
           xend = 0.58, yend = 0.146, 
           curvature = -.7, arrow = arrow(length = unit(2.5, "mm")), 
           colour = "grey40") +
  # add lines
  annotate(geom = "segment", x = c(4, 1.05), y = c(0.22, 0.055), # overall
           xend = c(4.55, 0.475), yend = c(0.22, 0.055), colour = "grey40") +
  # add background box
  annotate(geom = "rect", xmin = 4, xmax = 4.5, # overall
           ymin = 0.235, ymax = 0.26, fill = "white") +
  # add text
  annotate(geom = "text", x = 4.275, y = c(0.25, 0.235), # overall
           colour = "grey30", label = c("Overall", "mean"), size = 3.5) +
  annotate(geom = "text", x = 0.75, y = c(0.045, 0.03, 0.015), # per group
           colour = "grey30", label = c("Mean", "per", "group"), size = 3.5) +
  # add stars for significance
  annotate("text", x = c("CC", "CW", "WC", "WW"), 
           y = c(0.27, 0.283, 0.172, 0.185), colour = "grey20",
           label= c("***", "***", "***", "***")) +
  scale_x_discrete(labels = c("Cooling-Cooling", "Cooling-Warming", 
                              "Warming-Cooling", "Warming-Warming"), 
                   guide = guide_axis(n.dodge=2))


ggsave(plot = violin, filename = here("figures/violin_plot.png"), 
       width = 12.7, height = 9, units = "cm")

ggsave(plot = violin, filename = here("figures/violin_plot.pdf"), 
       width = 12.7, height = 9, units = "cm", dpi = 500)



# Significance testing ----------------------------------------------------


# testing whether the mean of a palaeoclimate Interaction
# significantly lies above or below the overall mean. 

# test for normality, visual
CCqq <- my_qqplot("CC", "Cooling-Cooling") 
CWqq <- my_qqplot("CW", "Cooling-Warming")  
WCqq <- my_qqplot("WC", "Warming-Cooling")
WWqq <- my_qqplot("WW", "Warming-Warming")

# combine it
qqplots <- CCqq + CWqq + WCqq + WWqq + plot_annotation(tag_levels = 'A')

# save it
ggsave(plot = qqplots, filename = here("figures/qqplots.png"), 
       width = 12.7, height = 9, units = "cm")


# performing Anova even though no perfect fit to normal data based on qq-plots
# anova is not very sensitive to moderate deviations from normality; 
# simulation studies, using a variety of non-normal distributions, have 
# shown that the false positive rate is not affected very much by this 
# violation of the assumption (Glass et al. 1972, Harwell et al. 1992, 
# Lix et al. 1996). 

# Compute the analysis of variance
res.aov <- aov(ori.prob ~ pal.int, data = prob)

# Summary of the analysis
summary(res.aov)

# As the ANOVA test is significant, we can compute Tukey HSD 
# (Tukey Honest Significant Differences) for performing multiple
# pairwise-comparison between the means of groups
TukeyHSD(res.aov) # all significant

# test wether the group mean significantly lies above or below the overall mean
# baseline is the overall response
overall <- prob$ori.prob[prob$pal.int=="CC" | prob$pal.int=="CW" | 
                           prob$pal.int=="WC" | prob$pal.int=="WW"]

# Cooling-Cooling is significantly greater
my_wilcoxtest("CC", "greater") 

# Cooling-Warming is significantly lower
my_wilcoxtest("CW", "less") 

# Warming-Cooling is significantly lower
my_wilcoxtest("WC", "less") 

# Warming-Warming is significantly lower
my_wilcoxtest("WW", "less") 



# Simulation based tests --------------------------------------------------

# Significance testing with classical methods (Anova, Tukey HSD, Wilcox) relies on 
# various assumptions. Instead of trying to fit them all, we can use bootstrapping, 
# which is more flexible and equally robust. 

# we use the following workflow:
# Step 1: Calculate the difference of means
# Step 2: Use simulation to model a world where the difference is null
# Step 3: Look at the difference in the null world
# Step 4: Calculate the probability that the difference could exist in null world
# Step 5: Decide if the difference is statistically significant based on probability

# first we need to define the groups we want to compare. We are interested in the 
# increase of origination probability after cooling-cooling, compared to all other 
# palaeoclimate interactions
prob_comparison <- prob %>%
  mutate(pal.int = fct_collapse(pal.int, 
                                cooling_cooling = "CC", 
                                other = c("CW", "WC", "WW")))   

# Step 1: Calculate difference of means
diff_prob <- prob_comparison %>% 
  specify(ori.prob ~ pal.int) %>% 
  calculate("diff in means", order = c("cooling_cooling", "other")) 


# Generate a bootstrapped distribution of the difference in means 
boot_means <- prob_comparison %>% 
  specify(ori.prob ~ pal.int) %>% 
  generate(reps = 2000, type = "bootstrap") %>% 
  calculate("diff in means", order = c("cooling_cooling", "other"))

# calculate confidence interval
boot_confint <- boot_means %>% get_confidence_interval()

# plot bootstrapped distribution of differences in means
# Red line shows observed difference; shaded area shows 95% confidence interval
boot_distr <- boot_means %>% 
  visualize() + 
  shade_confidence_interval(boot_confint,
                            color = "#354E71", fill = "lightblue") +
  geom_vline(xintercept = diff_prob$stat, size = 1, color = "#841F27") +
  labs(title = NULL, 
       x = "Cooling-Cooling vs. all other", y = "Count") +
  my_theme


# Step 2: Use simulation to model a world where the difference is null
# computationally demanding (2000 repetitions)
diff_prob_null <- prob_comparison %>% 
  specify(ori.prob ~ pal.int) %>% 
  hypothesize(null = "independence") %>% 
  generate(reps = 2000, type = "permute") %>% 
  calculate("diff in means", order = c("cooling_cooling", "other"))

# Step 3: Put actual observed difference in the null world and see if it fits
# Simulation-based null distribution of difference in means
boot_distr_null <- diff_prob_null %>% 
  visualize() + 
  geom_vline(xintercept = diff_prob$stat, size = 1, color = "#841F27") +
  labs(x = "Cooling-Cooling vs. all other",
       y = "Count", title = NULL) +
  my_theme

# Step 4: Calculate probability that observed difference could exist in null world
diff_prob_null %>% 
  get_p_value(obs_stat = diff_prob, direction = "both") %>% 
  mutate(p_value_clean = scales::pvalue(p_value))


# combine to one plot 
bootstrapping_plot <- boot_distr / boot_distr_null + plot_annotation(tag_levels = 'A')

# save it
ggsave(plot = bootstrapping_plot, filename = here("figures/bootstrapping.png"), 
       width = 12.7, height = 9, units = "cm")


# Bayesian Regression -----------------------------------------------------

# Instead of frequentist methods and null hypothesis testing, we can use bayesian
# inference.
# For comparison of group differences, the BEST regression (Bayesian estimation 
# supersedes the t test) is suitable


# set monte carlo parameters
CHAINS <- 4
ITER <- 1000
WARMUP <- 500
BAYES_SEED <- 1234
options(mc.cores = parallel::detectCores())  # Use all cores

# # specify the distribution of origination probability
mean(prob$ori.prob) # 0.125
sd(prob$ori.prob) # 0.034
min(prob$ori.prob) # 0.040
max(prob$ori.prob) # 0.278


# run the model usin brm and rcpp
brms_best <- brm(
    # we suppress the intercept by setting ori.prob ~ 0 + genre, 
    # brms returns coefficients for each of the groups, 
    # and these coefficients represent group means.
    bf(ori.prob ~ 0 + pal.int, sigma ~ 0 + pal.int), 
    family = student,
    data = prob_comparison,
    prior = c(
      # Set group mean prior, see distribution characteristics
      set_prior("normal(0.124, 0353)", class = "b", lb = 0.029, ub = 0.288),
      # Ser group variance priors
      # we use the default exponential prior with a rate of 1/29
      set_prior("exponential(1.0/29)", class = "nu")),
    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED, 
    # save it as brms_best.rds
    file = here("data/brms_best"))



# extract the posterior samples for each of the groups, subtract them from each other,
# and then calculate the credible interval
brms_best_post <- posterior_samples(brms_best) %>% 
  # Rescale sigmas (sigma terms are on a log scale, 
  # so we need to exponentiate them back to the scale of the data)
  mutate_at(vars(contains("sigma")), exp) %>% 
  # we need to log nu
  mutate(nu = log10(nu)) %>% 
  # calculate differences
  mutate(diff_means = b_pal.intcooling_cooling - b_pal.intother,
         diff_sigma = b_sigma_pal.intcooling_cooling - b_sigma_pal.intother)


# we can use tidyMCMC from broom to calculate the difference and 
# create confidence intervals
brms_best_tidy <- 
  broom::tidy(brms_best_post, conf.int = TRUE, conf.level = 0.95, 
           estimate.method = "median", conf.method = "HPDinterval") 


# difference in means ------------------------------------------------

# Extract simulation results
boot_results <- tibble(estimate = diff_prob$stat,
                           conf.low = boot_confint$lower_ci,
                           conf.high = boot_confint$upper_ci)

# Extract bayesian best results
best_results <- brms_best_tidy %>% 
  filter(column == "diff_means") %>% 
  select(estimate = mean, conf.low = min, conf.high = max)

# combine
diff_in_means <- 
full_join(boot_results, best_results) %>% 
  add_column(method = c("bootstrapping", "best")) %>% 
  # transform to percentage
  mutate(estimate = estimate*100, conf.low = conf.low*100, conf.high = conf.high*100)



# Percentage change -------------------------------------------------------


percent_change <- diff_in_means %>% 
  mutate(estimate = (estimate/(av*100))*100,
         conf.low = (conf.low/(av*100))*100,
         conf.high = (conf.high/(av*100))*100)

# Effect size -------------------------------------------------------------

# safe both cases (origination vs. all other) as vectors
treatment <- filter(prob_comparison, pal.int == "cooling_cooling") %>% pull(ori.prob)
control <- filter(prob_comparison, pal.int == "other") %>% pull(ori.prob)

# cohens d
effsize_raw <- effsize::cohen.d(treatment, control)

  
# the bayesian way
# calculate effect size following BEST Kruschke 
# (u_treat - u_contr)/ sqr((sd_treat^2+sd_treat^2)/2)
effsize_bayes <- brms_best_post %>% 
  # transform to tibble
  as_tibble() %>% 
  # calculate effect size per iteration
  transmute(numerator = b_pal.intcooling_cooling - b_pal.intother, 
            sigma_sum = b_sigma_pal.intcooling_cooling^2 + b_sigma_pal.intother^2, 
            sigma_div = sigma_sum/2, 
            denumerator = sqrt(sigma_div), 
            effect_size = numerator/denumerator) %>% 
  # calculate median effect size and 95 % CI using Hmisc
  summarise(ci = list(enframe(Hmisc::smedian.hilow(effect_size)))) %>% 
  unnest(cols = c(ci)) %>% 
  spread(name, value)

# combine to dataframe
effect_size <- tribble(
  ~conf.low,                 ~estimate,             ~conf.high,             ~method,
# ------------------------/---------------------/--------------------------/---------- 
effsize_raw$conf.int[[1]], effsize_raw$estimate, effsize_raw$conf.int[[2]],  "raw",
effsize_bayes$Lower,       effsize_bayes$Median, effsize_bayes$Upper,        "bayes"
)


# Combined plot -----------------------------------------------------------

# first save all the necessary data sets as list
# effect_plot_data <- list(diff_in_means, percent_change, effect_size)
# save(effect_plot_data, file = here("data/effect_plot_data.RData"))


# plot it

# define colours 
my_colours <- c("grey40", "#d5a069", "indianred")

# difference in means
mean_diff <- ggplot(diff_in_means, aes(x = estimate, y = method)) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high),
                 size = 1, colour = "grey60") +
  geom_point(aes(fill = method), size = 3, colour = "grey25", 
             shape = 21, stroke = 1) +
  labs(y = NULL, x = "Difference in means") +
  my_theme +
  theme(panel.grid.major.x=element_line(colour = "grey", linetype = "dotted"), 
        panel.grid.major.y = element_blank(), 
        legend.position = "none", 
        panel.grid.minor.x = element_blank()) +
  scale_y_discrete(labels = c("Bayesian estimation", "Bootstrapping")) +
  scale_fill_manual(values = my_colours[2:3])


# percentage change caused by cooling-cooling
perc_change <- ggplot(percent_change, aes(x = estimate, y = method)) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high),
                 size = 1, colour = "grey60") +
  geom_point(aes(fill = method), size = 3, colour = "grey25",
             shape = 21, stroke = 1) +
  labs(y = NULL, x = "Percentage change") +
  my_theme +
  theme(panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"), 
        panel.grid.major.y = element_blank(), 
        legend.position = "none", 
        panel.grid.minor.x = element_blank()) +
  scale_y_discrete(labels = c("Bayesian estimation", "Bootstrapping")) + 
  scale_x_continuous(labels = function(x) paste0(x, '%')) +
  scale_fill_manual(values = my_colours[2:3])


# effect size cohens d
cohens_d <- ggplot(effect_size, aes(x = estimate, y = method)) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high),
                 size = 1, colour = "grey60") +
  geom_point(aes(fill = method), size = 3, colour = "grey25", 
             shape = 21, stroke = 1) +
  labs(y = NULL, x = "Effect size (Cohen's d)") +
  my_theme +
  theme(panel.grid.major.x=element_line(colour = "grey", linetype = "dotted"), 
        panel.grid.major.y = element_blank(), 
        legend.position = "none",
        panel.grid.minor.x = element_blank()) +
  scale_y_discrete(labels = c("Bayesian estimation", "Raw data")) +
  scale_fill_manual(values = my_colours[2:1])


# arrange plot and annotate
combined_effect_size <- mean_diff / perc_change / cohens_d + 
  plot_annotation(tag_levels = 'A')

# save it
ggsave(plot = combined_effect_size, filename = here("figures/combined_effect_sizes.png"), 
       width = 11, height = 11, units = "cm")

ggsave(plot = combined_effect_size, filename = here("figures/combined_effect_sizes.pdf"), 
       width = 11, height = 11, units = "cm", dpi = 500)
