# this script contains all self-defined functions used throughout the data processing

# set up function to calculate a regression between two rows/ stages for fragmentation index
# this function calculates the short-term change in the fragmentation index
short_term <- function(i, j) {
  dum1 <- filter(zaffos_binned2,
                 zaffos_binned2$stg %in% zaffos_binned2$stg[between(zaffos_binned2$stg, i-j, i)])
  
  dum2 <-  dum1 %>% 
    unnest(fragmentation_index) %>% 
    lm(formula = fragmentation_index ~ stg)
  
  -dum2$coefficients[2]
}

# function that calculates the long-term trends for the fragmenation index
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