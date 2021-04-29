# this script contains all self-defined functions used throughout the data processing

### 
# this functions calculates the log Odds ratio with Wald Confidence Intervals.
# With a default value for alpha at 0.05, which yields the 95% confidence intervals 
# for the computed odds ratio, based on the Wald approximation.  
calculate_log_odds <- function(n00, n01, n10, n11, alpha = 0.05) {
  ###
  # Compute the odds ratio between two binary outcomes, x and y
  # n00 = number of cases where x = 0 and y = 0
  # n01 = number of cases where x = 0 and y = 1
  # n10 = number of cases where x = 1 and y = 0
  # n11 = number of cases where x = 1 and y = 1
  ###
  
  # the odds ratio
  OR <- (n00 * n11) / (n01 * n10)
  
  # Compute the Wald confidence intervals:
  siglog <- sqrt((1 / n00) + (1 / n01) + (1 / n10) + (1 / n11))
  zalph <- qnorm(1 - alpha / 2)
  logOR <- log(OR)
  loglo <- logOR - zalph * siglog
  loghi <- logOR + zalph * siglog
  # produce output
  log_odds_output <-
    tibble(
      LowerCI = round(loglo, 3),
      OR = round(logOR, 3),
      UpperCI = round(loghi, 3)
    )
  log_odds_output
}
###

###
# testing whether a model is overdispersed
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
###

###
# function for adding AIC, BIC and overdispersion to a predefined data frame. 
summarise_model <- function(model){
  # convert to character for subsetting
  model_chr <- deparse(substitute(model))
  # summary
  sum <- summary(model) 
  # AIC
  model_comparison$AIC[model_comparison$models == model_chr] <<- 
    as.numeric(round(sum$AICtab[[1]], 1)) 
  # BIC
  model_comparison$BIC[model_comparison$models == model_chr] <<- 
    as.numeric(round(sum$AICtab[[2]], 1)) 
  # test for overdispersion
  model_comparison$overdispersed[model_comparison$models == model_chr] <<- 
    ifelse(overdisp_fun(model)[4]< 0.05, "yes", "no") 
  filter(model_comparison, models == model_chr)
}
###

###
# function for selecting the best model from a list with paleoclimate interactions 
select_model <- function(model){
  # model summary for each model in the list
  model_sum <- map(model, summary)
  # get the AIC values for each model, substitute with BIC in the second call to get 
  # BIC values
  model_sum_aic <- map(model_sum, "AICtab")
  model_sum_aic <- map(model_sum_aic, "AIC")
  # change the list to a vector and select the model with the lowest AIC value
  best_model <- model_sum_aic %>% unlist() %>% 
    enframe() %>% 
    filter(value == min(value))
  # choose final model
  model[[best_model[[1]]]]
}
###


###
# function to calculate log Odds that taxa show higher probability of origination
# after cooling-ccoling palaeclimate interaction, per phyla
interaction_log_odds <- function(my_phylum, my_bins = unique(dat_final$bins)) {
  ###  warming-warming
  ww_raw <-
    subset(dat_final, trend.st7 >= 0 &
             warming >= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for warming-warming palaeclimate interaction with log odds as output
  ww_pred <-
    predict(warming_pal_int, newdata = ww_raw, type = "link")
  # how many warming-warming cases show higher probability
  h_ww <-
    length(subset(ww_pred, ww_pred > warming_sum$coefficients[1]))
  # how many warming-warming cases show lower probability
  l_ww <-
    length(subset(ww_pred, ww_pred < warming_sum$coefficients[1]))
  
  ###  cooling-warming
  cw_raw <-
    subset(dat_final, trend.st7 <= 0 &
             warming >= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for cooling-warming palaeclimate interaction with log odds as output
  cw_pred <-
    predict(warming_pal_int, newdata = cw_raw, type = "link")
  # how many cooling-warming cases show higher probability
  h_cw <-
    length(subset(cw_pred, cw_pred > warming_sum$coefficients[1]))
  # how many cooling-warming cases show lower probability
  l_cw <-
    length(subset(cw_pred, cw_pred < warming_sum$coefficients[1]))
  
  ###  warming cooling
  wc_raw <-
    subset(dat_final, trend.st6 >= 0 &
             cooling <= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for warming-cooling palaeclimate interaction with log odds as output
  wc_pred <- predict(cooling_pal_int, newdata = wc_raw,
                     type = "link")
  # how many warming-cooling cases show higher probability
  h_wc <- length(subset(wc_pred,
                        wc_pred > cooling_sum$coefficients[1]))
  # how many warming-cooling cases show lower probability
  l_wc <- length(subset(wc_pred,
                        wc_pred < cooling_sum$coefficients[1]))
  
  ###  cooling cooling
  cc_raw <-
    subset(dat_final, trend.st6 <= 0 &
             cooling <= 0 & phylum %in% my_phylum & bins %in% my_bins)
  # predictions for cooling-cooling palaeclimate interaction with log odds as output
  cc_pred <- predict(cooling_pal_int, newdata = cc_raw,
                     type = "link")
  # how many cooling-cooling cases show higher probability
  h_cc <- length(subset(cc_pred,
                        cc_pred > cooling_sum$coefficients[1]))
  # how many cooling-cooling cases show lower probability
  l_cc <- length(subset(cc_pred,
                        cc_pred < cooling_sum$coefficients[1]))
  
  ##############
  # log odds
  
  # calculate log odds for cooling-cooling
  my_output <- calculate_log_odds(h_cc, l_cc, h_wc, l_wc)
  my_output
}
###


###
# set up function to calculate a regression between two rows/ stages for fragmentation index
# this function calculates the short-term change in the fragmentation index
short_term <- function(i, j) {
  dum1 <- filter(zaffos_binned,
                 stg %in% stg[between(stg, i-j, i)])
  
  dum2 <-  dum1 %>% 
    unnest(fragmentation_index) %>% 
    lm(formula = fragmentation_index ~ stg)
  
  -dum2$coefficients[2]
}
### 

###
# same for temperature
short_term_temp <- function(i, j) {
  dum1 <- filter(isotemp,
                 stg %in% stg[between(stg, i-j, i)])
  
  dum2 <-  dum1 %>% 
    unnest(Temp) %>% 
    lm(formula = Temp ~ stg)
  
  -dum2$coefficients[2]
}
### 

###
# function that calculates the long-term trends for the fragmenation index
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term, and thus shift the long 
# term lm result up by i-1.
long_term <- function(i, j) {
  dum1 <- filter(zaffos_binned2,
                 stg %in% stg[between(stg, i-j, i-1)])
  
  if(nrow(dum1) != 0) {
    dum2 <-  dum1 %>% 
      unnest(fragmentation_index) %>% 
      lm(formula = fragmentation_index ~ stg)
    
    -dum2$coefficients[2]
  } else NA
}
###

###
# same for temperatures
long_term_temp <- function(i, j) {
  dum1 <- filter(isotemp2,
                 stg %in% stg[between(stg, i-j, i-1)])
  
  if(nrow(dum1) != 0) {
    dum2 <-  dum1 %>% 
      unnest(Temp) %>% 
      lm(formula = Temp ~ stg)
    
    -dum2$coefficients[2]
  } else NA
}
###

### 
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
###

###
# take the predictions for a palaeoclimate interaction and produce a ggplot to
# (visually) check for normality
my_qqplot <- function(PI, title) {
  ggpubr::ggqqplot(prob$ori.prob[prob$pal.int == PI], title = title, ggtheme = my_theme)
}
###

###
# take the predictions and test whether they significantly are above or below the
# baseline. This is a one sided Wilcoxon rank sum test 
# (equivalent to the Mann-Whitney test)
my_wilcoxtest <- function(PI, alternative) {
  wilcox.test(prob$ori.prob[prob$pal.int == PI],
              overall,
              paired = FALSE,
              alternative = alternative)
}
###

###
# calculate the cross correlation without plottting and in a tidy way
cross_corr <- function(.data) {
  ccf(.data$temp, .data$fragmentation_index, 
      plot = FALSE, lag.max = 10)
}
###

###
# calculate the ci for the cross correlation
cross_corr_ci <- function(n.used) {
  qnorm((1 + 0.95)/2)/sqrt(n.used)
}

###

### calculate pearsons r correlation and return a tibble
cor_test_r <- function(.data) {
  cor.test(.data$temp, .data$fragmentation_index) %>% 
    tidy() %>% 
    select(estimate, conf.low, conf.high, p.value)
}
###

###
# separate sqs data into temporal groups
separate_groups <- function(time_interval) {
  datsqs %>%
    filter(stg %in% time_interval) %>% 
    full_join(log_odds) %>% 
    drop_na() %>% 
    mutate(weighted_estimate = n*estimate) %>% 
    summarise(mean_w_est = mean(weighted_estimate, na.rm = TRUE), 
              sd_est = sd(weighted_estimate, na.rm = TRUE), 
              sum_n = sum(n, na.rm = TRUE)) %>% 
    transmute(estimate = mean_w_est/sum_n, 
              sd_est = sd_est/sum_n)
}

###

### 
# scale/ normalize column between 0 and 1
range_scale <- function(col_name) {
  (col_name - min(col_name)) /
    (max(col_name) - min(col_name))
}

###

###
# custom ggplot theme
# define theme
my_theme <- theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                  panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor = element_blank(),
                  text = element_text(family = "sans"),
                  strip.background = element_blank())

###