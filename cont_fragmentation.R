library(tidyverse)
library(chronosphere)
library(divDyn)

select <- dplyr::select

# load zaffos fragmentation index
zaffos <- fetch(dat = "som", var = "zaffos-fragment") %>% 
  as_tibble() %>% 
  select(age_seq = X, fragmentation_index = fragmentation.index)

# load stage data for binning
data(stages)

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

# set up function to calculate a regression between two rows/ stages
short_term <- function(i, j) {
  dum1 <- filter(zaffos_binned2,
                 zaffos_binned2$stg %in% zaffos_binned2$stg[between(zaffos_binned2$stg, i-j, i)])
  
  dum2 <-  dum1 %>% 
    unnest(fragmentation_index) %>% 
    lm(formula = fragmentation_index ~ stg)
  
  -dum2$coefficients[2]
}

# calculate short-term change in fragmentation index
zaffos_binned2 <- zaffos_binned2 %>% 
  select(-data) %>% 
  mutate(change.prev = map_dbl(unique(zaffos_binned2$stg), short_term, j = 1))


# Use lm() to determine slope of the long-term fragmentation trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term, and thus shift the long 
# term lm result up by i-1.
long_term <- function(i, j) {
  dum1 <- filter(zaffos_binned2,
                 zaffos_binned2$stg %in% zaffos_binned2$stg[between(zaffos_binned2$stg, i-j, i-1)])
  
  if(nrow(dum1) != 0) {
    dum2 <-  dum1 %>% 
    unnest(fragmentation_index) %>% 
    lm(formula = fragmentation_index ~ stg)
  
  -dum2$coefficients[2]
  } else NA
}

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



