# set working directory to where files are stored using the "rstudioapi" package

# Getting the path of this script
current_path = rstudioapi::getActiveDocumentContext()$path 

# setting it as working directory 
setwd(dirname(current_path ))

# load libraries
library(tidyverse) # for data processing
library(here) # for project tidyness
library(flextable) # for creating nice tables
library(officer) # to render tables to word file



# # number of taxa within phyla -------------------------------------------

# load cleaned pbdb data
load(here("data/final_data.RData"))

# set up rows for phyla
my_rows <- dat_final %>% distinct(phylum) %>% na.omit(.) %>% pull()

# set upt empty data frame
pbdb_smry <- tibble(Phylum = my_rows, Class = numeric(12), Order = numeric(12), 
       Family = numeric(12), Genus = numeric(12))

# set up taxonomic resolution
taxon <- c("class", "order", "family", "genus")

# go through data and count number of taxa, and fill in dataframa
for (i in seq_along(my_rows)) {
  for (j in seq_along(taxon)) {
    pbdb_smry[i, j+1] <- dat_final %>%
      filter(phylum == my_rows[i]) %>%
      select(taxon[j]) %>% 
      distinct() %>%
      tally() %>%
      pull() 
  }
}

# order alphabetically
pbdb_smry <- pbdb_smry %>% arrange(Phylum)

# make a flextable
pbdb_fxt <- flextable(pbdb_smry) %>% theme_zebra() %>% autofit()

# open docx-file and add flextable
my_doc <- read_docx() %>% body_add_flextable(pbdb_fxt)


# Model comparison --------------------------------------------------------

# load table
load(here("data/model_comparison.RData"))

# make model structure more clear
model_comparison$models <- 
  c("~ Warming", "~ Warming + Pal. Int.",
    "~ Cooling + Lag", "~ Cooling + Pal. Int.",
    "~ Warming + Lag", "~ Warming + Lag + Pal. Int.",
    "~ Cooling + Lag", "~ Cooling + Lag + Pal. Int.",
    "~ Warming + Lag + Warming:Lag", "~ Warming + Lag + Warming:Lag + Pal. Int.",
    "~ Cooling + Lag + Cooling:Lag", "~ Cooling + Lag + Cooling:Lag + Pal. Int.")

# indicate significance
model_comparison$overdispersed <- rep("no ***", 12)

# make flextable
model_comparison_fxt <- model_comparison %>% 
  flextable(col_keys = c("models", "overdispersed", "AIC", "BIC")) %>% 
  set_header_labels(models = "Model Structure", overdispersed = "Overdispersed") %>% 
  theme_zebra() %>% 
  autofit() %>% 
  hline(i = c(2, 4, 6, 8, 10), border = officer::fp_border(width = 1)) 
  

# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(model_comparison_fxt, pos = "after")



# All palaeoclimate interactions ------------------------------------------

# load list of pal_int_df
load(here("data/pal_int_df.RData"))

# first table for intercept and interaction
pal_int_df_fxt1 <- pal_int_df %>% 
  bind_rows() %>% 
  add_column(Type = c(rep("Warming", 10), rep("Cooling", 10)), 
             .before = "model") %>% 
  select(Type, model, intercept, interaction) %>% 
  flextable() %>% 
  set_header_labels(model = "Model", intercept = "Intercept", 
                    interaction = "Interaction") %>% 
  theme_zebra() %>% 
  autofit() %>% 
  hline(i = 10, border = officer::fp_border(width = 1)) %>% 
  merge_v(j = ~ Type)

# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(pal_int_df_fxt1, pos = "after")


# second table for AIC and BIC
pal_int_df_fxt2 <- pal_int_df %>% 
  bind_rows() %>% 
  add_column(Type = c(rep("Warming", 10), rep("Cooling", 10)), 
             .before = "model") %>% 
  select(Type, model, AIC, BIC, dAIC, dBIC) %>% 
  flextable() %>% 
  set_header_labels(model = "Model", dAIC = "\u0394AIC",
                    dBIC = "\u0394BIC") %>% 
  theme_zebra() %>% 
  autofit() %>% 
  hline(i = 10, border = officer::fp_border(width = 1)) %>% 
  merge_v(j = ~ Type)

# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(pal_int_df_fxt2, pos = "after")



# palaeoclimate interactions final -----------------------------------------------

# load list of palaeoclimate interactions
load(here("data/pal_int_models.RData"))

# get the models from the list
warming <- pal_int_models[[1]]
cooling <- pal_int_models[[2]]

# produce model summaries
warming_tidy <- warming %>% broom::tidy()
cooling_tidy <- cooling %>% broom::tidy()

model_tidy <- warming_tidy %>% full_join(cooling_tidy)

model_tidy_fxt <- model_tidy %>% 
  add_column(Model = c(rep("Warming", 3), rep("Cooling", 3)), 
                          .before = "term") %>% 
  mutate(p.value = gtools::stars.pval(p.value)) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(Term = c("Intercept", "Warming:Trend.st7", "Random effect", 
                  "Intercept", "Cooling:Trend.st6", "Random effect")) %>% 
  select(Model, Term, everything(), -term) %>% 
  flextable() %>% 
  set_header_labels(term = "Term", estimate = "Estimate", std.error = "Std.error", 
                    statistic = "Z value", p.value = "P value", group = "Group") %>% 
  theme_zebra() %>% 
  autofit() %>% 
  hline(i = 3, border = officer::fp_border(width = 1)) %>% 
  merge_v(j = ~ Model) 
  
# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(model_tidy_fxt, pos = "after")



# model performance
warming_glance <- warming %>% broom::glance()
cooling_glance <- cooling %>% broom::glance()

model_perf <- warming_glance %>% full_join(cooling_glance)

model_perf_fxt <- model_perf %>% 
  add_column(Model = c("Warming", "Cooling"), .before = "sigma") %>% 
  add_column(Overdispersion = rep("no ***", 2)) %>% 
  mutate_if(is.numeric, round, 2) %>%
  flextable() %>% 
  set_header_labels(sigma = "Sigma", logLik = "LogLik", deviance = "Deviance", 
                    df.residual = "DF residual") %>% 
  theme_zebra() %>% 
  autofit()

# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(model_perf_fxt, pos = "after")




# Violin plot -------------------------------------------------------------

# load prob
load(here("data/violin_plot_data.RData"))

# summarise each palaeclimate response using quantiles
prob_fxt <- prob %>% group_by(pal.int) %>% 
  summarise(my_quant = list(quantile(ori.prob, c(0.25, 0.5, 0.75))), 
            q = list(c(0.25, 0.5, 0.75))) %>% 
  unnest(c(my_quant, q)) %>% 
  pivot_wider(names_from = q, values_from = my_quant) %>%
  mutate(pal.int = replace(pal.int, 
                           pal.int == c("CC", "CW", "WC", "WW"),
                           c("Cooling-Cooling", "Cooling-Warming", 
                             "Warming-Cooling", "Warming-Warming"))) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  flextable() %>% 
  set_header_labels(pal.int = "Palaeoclimate Interaction", '0.25' = "Lower Quartile",
                    '0.5' = "Median", '0.75' = "Upper Quartile") %>% 
  theme_zebra() %>% 
  autofit()
  
# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(prob_fxt, pos = "after")





# effect size -------------------------------------------------------------

# load effect_plot_data
load(here("data/effect_plot_data.RData"))


effect_data_fxt <- effect_plot_data %>% 
  bind_rows() %>%
  add_column(Parameter = c(rep("Difference in means", 2),
                           rep("Percentage change", 2), 
                           rep("Cohen's d", 2)), .before = "estimate") %>% 
  mutate_if(is.numeric, round, 2) %>%
  mutate(Method = c("Bootstrapping", "Bayesian Estimate", 
                    "Bootstrapping", "Bayesian Estimate", 
                    "Raw Data", "Bayesian Estimate")) %>% 
  select(Parameter, conf.low, estimate, conf.high, Method, -method) %>% 
  flextable() %>% 
  set_header_labels(estimate = "Estimate", conf.low = "Lower CI",
                    conf.high = "Upper CI", method = "Method") %>% 
  theme_zebra() %>% 
  autofit() %>% 
  merge_v(j = ~ Parameter) %>% 
  hline(i = c(2,4), border = officer::fp_border(width = 1)) 

# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(effect_data_fxt, pos = "after")




# log Odds ratio ----------------------------------------------------------

# load log_odds
load(here("data/log_odds_phyla.RData"))

log_odds_fxt <- log_odds %>% 
  select(-type) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  flextable() %>% 
  set_header_labels(name = "Group", lower_CI = "Lower CI", 
                    estimate = "Log Odds ratio",
                    upper_CI = "Upper CI", method = "Method") %>% 
  theme_zebra() %>% 
  autofit() %>% 
  hline(i = 11, border = officer::fp_border(width = 1)) 

  
# add to word file
my_doc <- my_doc %>% 
  body_add_break() %>% 
  body_add_flextable(log_odds_fxt, pos = "after")


# Make word file ----------------------------------------------------------

# convert to word file/ add input to empty docx
print(my_doc, target = here("figures/supplemental_tables.docx"))

