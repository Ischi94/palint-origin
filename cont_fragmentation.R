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

for (i in unique(zaffos_binned2$stg)) {
  dum1 <- filter(zaffos_binned2,
                 zaffos_binned2$stg %in% zaffos_binned2$stg[between(zaffos_binned2$stg, i-1, i)])
  
  dum2 <-  dum1 %>% 
    unnest(fragmentation_index) %>% 
    lm(formula = fragmentation_index ~ stg)
  
  zaffos_binned2[zaffos_binned2$stg==i, "change.prev"] <- -dum2$coefficients[2]
}

# Calculate lags 1 to 10
zaffos_binned <- zaffos_binned2 %>%  
  select(-data) %>% 
  mutate(lag1 = lag(fragmentation_index, order_by = stg), 
         lag2 = lag(lag1, order_by = stg),
         lag3 = lag(lag2, order_by = stg),
         lag4 = lag(lag3, order_by = stg),
         lag5 = lag(lag4, order_by = stg),
         lag6 = lag(lag5, order_by = stg),
         lag7 = lag(lag6, order_by = stg),
         lag8 = lag(lag7, order_by = stg),
         lag9 = lag(lag8, order_by = stg),
         lag10 = lag(lag9, order_by = stg)) 

# Use lm() to determine slope of the long-term fragmentation trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term, and thus shift the long 
# term lm result up by +1.
for (i in unique(zaffos_binned$stg)) {
  sub1 <- filter(zaffos_binned2,
                 zaffos_binned2$stg %in% zaffos_binned2$stg[between(zaffos_binned2$stg, i-1, i)])
  
  lin1 <- lm(Temp ~ age, data = sub1)
  isotemp2[isotemp2$Stage == i + 1, "trend.st1"] <-
    -lin1$coefficients[2]
  
  sub2 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               2, i)])
  lin2 <- lm(Temp ~ age, data = sub2)
  isotemp2[isotemp2$Stage == i + 1, "trend.st2"] <-
    -lin2$coefficients[2]
  
  sub3 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               3, i)])
  lin3 <- lm(Temp ~ age, data = sub3)
  isotemp2[isotemp2$Stage == i + 1, "trend.st3"] <-
    -lin3$coefficients[2]
  
  sub4 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               4, i)])
  lin4 <- lm(Temp ~ age, data = sub4)
  isotemp2[isotemp2$Stage == i + 1, "trend.st4"] <-
    -lin4$coefficients[2]
  
  sub5 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               5, i)])
  lin5 <- lm(Temp ~ age, data = sub5)
  isotemp2[isotemp2$Stage == i + 1, "trend.st5"] <-
    -lin5$coefficients[2]
  
  sub6 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               6, i)])
  lin6 <- lm(Temp ~ age, data = sub6)
  isotemp2[isotemp2$Stage == i + 1, "trend.st6"] <-
    -lin6$coefficients[2]
  
  sub7 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               7, i)])
  lin7 <- lm(Temp ~ age, data = sub7)
  isotemp2[isotemp2$Stage == i + 1, "trend.st7"] <-
    -lin7$coefficients[2]
  
  sub8 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               8, i)])
  lin8 <- lm(Temp ~ age, data = sub8)
  isotemp2[isotemp2$Stage == i + 1, "trend.st8"] <-
    -lin8$coefficients[2]
  
  sub9 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               9, i)])
  lin9 <- lm(Temp ~ age, data = sub9)
  isotemp2[isotemp2$Stage == i + 1, "trend.st9"] <-
    -lin9$coefficients[2]
  
  sub10 <-
    filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i -
                                                               10, i)])
  lin10 <- lm(Temp ~ age, data = sub10)
  isotemp2[isotemp2$Stage == i + 1, "trend.st10"] <-
    -lin10$coefficients[2]
}

