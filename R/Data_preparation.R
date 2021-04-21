# download libraries
library(divDyn)
library(tidyverse)
library(here)

# load self-defined functions 
source(here("R/functions.R"))



# 1 load in downloaded PBDB data ------------------------------------------


# load in data. The data contains all occurrences and parameters of the pbdb
# data base
load(file=here("data/allData_2021-04-21.RData"))



# 2 Filtering by taxonomic data ----------------------------------------------


# Filtering started with omission of occurrences that were not identified to the level 
# of genus, which is indicated in the accepted_rank field or with 
# non-informative genus names (empty quotes)
pbdb <- pbdb[pbdb$accepted_rank %in% c("genus", "species"), ]
pbdb <- pbdb[pbdb$genus!="", ]

# The filtering process continues with the selection of phyla that have or had recorded marine
# species. Phyla that have considerable terrestrial and marine records 
# (chordates and arthropods) are included based on lower taxonomic ranks.

marineNoPlant <- c("", "Agmata", "Annelida", "Bilateralomorpha", "Brachiopoda", 
                   "Bryozoa", "Calcispongea", "Chaetognatha", "Cnidaria", "Ctenophora",
                   "Echinodermata", "Entoprocta", "Foraminifera", "Hemichordata", 
                   "Hyolitha", "Mollusca", "Nematoda", "Nematomorpha", "Nemertina",
                   "Onychophora", "Petalonamae", "Phoronida",  "Platyhelminthes",
                   "Porifera", "Problematica", "Rhizopodea","Rotifera", 
                   "Sarcomastigophora", "Sipuncula", "Uncertain", "Vetulicolia","" )

# Then a logical vector was defined for each row, suggesting whether the phylum of that
# occurrence is present or not in the above defined set marineNoPlant.
bByPhyla <- pbdb$phylum %in% marineNoPlant

# The rest of the data were saved for further filtering. Classes that are likely to be 
# marine animals were extracted from this subset.
needClass <- c("Acanthodii", "Actinopteri", "Actinopterygii", "Agnatha", 
               "Cephalaspidomorphi", "Chondrichthyes", "Cladistia", 
               "Coelacanthimorpha", "Conodonta", "Galeaspida", "Myxini", 
               "Osteichthyes", "Petromyzontida", "Plagiostomi", "Pteraspidomorphi",
               "Artiopoda",  "Branchiopoda", "Cephalocarida", "Copepoda", 
               "Malacostraca", "Maxillopoda", "Megacheira", "Merostomoidea", 
               "Ostracoda", "Paratrilobita", "Pycnogonida", "Remipedia", 
               "Thylacocephala", "Trilobita", "Xiphosura" )

# logical vector of rows indicating occurrences
bNeedClass <- pbdb$class %in% needClass

# The mammalian orders Sirenia and Cetacea were also included, as well as pinniped 
# families from the order Carnivora.
needMammalOrd <- c("Cetacea", "Sirenia")
bMammalOrder <- pbdb$order %in% needMammalOrd

# the carnivores
needFam <- c("Otariidae", "Phocidae", "Desmatophocidae")
bNeedMamFam <- pbdb$family %in% needFam

# Some reptile orders were also included:
needReptOrd<-c("Eosauropterygia", "Hupehsuchia", "Ichthyosauria", "Placodontia",
               "Sauropterygia", "Thalattosauria")

# the logical vector for the total data
bRept <- pbdb$order %in% needReptOrd

# Families of sea turtles are also added to the analyzed set.
needTurtleFam <- c("Cheloniidae", "Protostegidae", "Dermochelyidae",
                   "Dermochelyoidae", "Toxochelyidae", "Pancheloniidae" )

# the logical vector for the total data
bTurtle <- pbdb$family%in%needTurtleFam

# And then, we subsetted the data with these multiple filters. 
# In this sense, the logical OR | works as a union operator between the taxonomic 
# groups.
pbdb <- pbdb[bByPhyla | bNeedClass | bMammalOrder | bNeedMamFam | bRept | bTurtle, ]

# After the filtering by taxonomy was finished, we can make sure that potential 
# homonymies will not affect the results by combining the class names and genus 
# names to create individual  entries:
pbdb$clgen <- paste(pbdb$class, pbdb$genus)


# 3 Filtering by sedimentary environment ----------------------------------


# Some of the taxa above contain freshwater and/or terrestrial groups, which need to 
# be omitted
omitEnv <- c("\"floodplain\"", "alluvial fan", "cave", "\"channel\"", "channel lag" ,
             "coarse channel fill", "crater lake", "crevasse splay", "dry floodplain",
             "delta plain", "dune", "eolian indet.", "fine channel fill", 
             "fissure fill", "fluvial indet.", "fluvial-lacustrine indet.", 
             "fluvial-deltaic indet.", "glacial", "interdune", "karst indet.", 
             "lacustrine - large", "lacustrine - small", "lacustrine delta front", 
             "lacustrine delta plain", "lacustrine deltaic indet.", "lacustrine indet.",
             "lacustrine interdistributary bay", "lacustrine prodelta", "levee",
             "loess", "mire/swamp", "pond", "sinkhole", "spring", "tar", 
             "terrestrial indet.", "wet floodplain")

# We can omit occurrences with these entries with a similar command as above.
pbdb <- pbdb[!pbdb$environment%in%omitEnv, ]

# Collections that came from unlithified sediments can yield fossils with
# unusually good preservation. As these are more frequent in younger sites and occur 
# heterogeneously, sampling bias can be reduced by omitting such collections from the
# data.
pbdb <- pbdb[pbdb$lithification1!="unlithified", ]



# 4 Stratigraphic binning to stages ---------------------------------------


# load stages and key file
data(stages)
data(keys)

# rename file
dat <- pbdb

# the 'stg' entries (lookup)
stgMin <- categorize(dat[ ,"early_interval"], keys$stgInt)
stgMax <- categorize(dat[ ,"late_interval"], keys$stgInt)

# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)

# empty container
dat$stg <- rep(NA, nrow(dat))

# select entries, where
stgCondition <- c(
  # the early and late interval fields indicate the same stg
  which(stgMax==stgMin),
  # or the late_intervar field is empty
  which(stgMax==-1))

# in these entries, use the stg indicated by the early_interval
dat$stg[stgCondition] <- stgMin[stgCondition]

#Stage-level assignments have been a problem in the earliest Paleozoic 
# (Cambrian and Ordovician periods). Therefore, the assignments above are not 
# perfect in these periods: numerous entries are not properly processed based on 
# the interval lookup table. To make use of these collections in the analyses that 
# follow, we processed these data separately.

# For the Ordovician assignments, Wolfgang Kiessling compiled tables on formations 
# and biozones.
load(url(
  "https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/ordStrat.RData"))

#You can use these with the following script.
source(
  "https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2019-05-31/ordProcess.R")

# omit all entries without a stage assignment and older than Ordovician
dat <- subset(dat, !is.na(stg))


# Sampling-standardized values can be calculated with a single function 
# called subsample().
set.seed(43)
datsqs <- subsample(dat, bin="stg", tax="clgen", coll="collection_no", q=0.8,
                    iter=1, ref="reference_no",singleton="ref", type="sqs", 
                    duplicates=FALSE, excludeDominant=TRUE, largestColl =TRUE,
                    output="list", na.rm=TRUE, FUN = NULL)
# save it
datsqs <- datsqs$results[[1]]


# remove problematica and uncertain phyla
datsqs <- datsqs[datsqs$phylum != "Problematica" & datsqs$phylum != "Uncertain",]

# save(datsqs, file=here("data/occurrence_sqs_data.RData"))

# clean up
rm(list=ls()[! ls() %in% c("datsqs", "stages")])

# order them properly
names(datsqs)
dat <- datsqs %>%
  dplyr::select(phylum, class, order, family, genus, stg)  

# sometimes genus entries have brackets within the characters
# subset the genera with brackets and assign them to the overarching genus
dat$genus <- gsub("\\s*\\([^\\)]+\\)","",as.character(dat$genus))


# remove the NA's
dat <- subset(dat, dat$class != "NO_CLASS_SPECIFIED" &
                dat$order != "NO_ORDER_SPECIFIED" &
                dat$family !=  "NO_FAMILY_SPECIFIED")


# check for NA's
dat <- na.omit(dat)

# save higher taxonomy
tax <- dat[,c(1:5)]
  
# transform to range data
dat_range <- fadlad(dat, tax="genus", bin="stg")

#remove singletons as they add noise 
dat_range <- subset(dat_range, dat_range$FAD != dat_range$LAD)

# when genera range to the recent (bin 95), remove them
dat_range <- dat_range[dat_range$LAD < 95,]

# assign genus names to a column
dat_range$genus <- rownames(dat_range)

# merge them together to get higher taxonomy
dat_range <- merge(dat_range,tax, by  =  "genus")

# remove the duplicates that were generated by merging
dat_range <- dat_range[!duplicated(dat_range),]

# order them properly
dat_range <- dat_range %>%
  dplyr::select(phylum, class, order, family, genus, FAD, LAD)   


# save it as a binary R-file
# no singletons, no NA's, binned to stages, no problematica and uncertain phyla, 
# no unlitified, sampling-standardized, no genera ranging to the modern, 
# no brackets (duplicates) in the genus field
# save(dat_range, file=here("data/cleaned_data.RData"))



# 5 Calculate short-term temperature change ---------------------------------


# prepare weizer & prokoph temperature data. We are loading the file which was 
# already processed as described
# in the methods paragraph and as described by Reddin et al. 2017
isotemp <- read.csv(file=here("data/TimeSeriesUsed.csv"), header=TRUE,row.names = 1) 

# assign age to isotope data
isotemp$age <- stages$mid[14:94]


# Calculate average Temperature per bin starting in the Ordovician (481.55 myr)

isotemp2<- isotemp %>%
  subset(age<=481.55000) %>%
  group_by(Stage) %>%
  summarise(Temp.mean = mean(Temp), n = n()) %>%
  arrange(desc(Stage))
isotemp2$age <- stages$mid[94:14]

# Add a column with short-term change (change.prev) in Temperature 
# (temperature change from the previous to the focal bin)
isotemp2$change.prev <- double(length = 81)

# fill in the values for short-term change using the lm as we do for the long 
# term trend as well
for (i in unique(isotemp$Stage)) {
  dum1 <- filter(isotemp,
                 isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-1, i)])
  dum2 <-  lm(Temp ~ age, data = dum1)
  isotemp2[isotemp2$Stage==i, "change.prev"] <- -dum2$coefficients[2]
}

# Equal the number of temperature time bins 
isotemp2 %<>% subset(Stage > 13 & Stage < 95) %>%
  transform(bins = as.factor(Stage))  
# We need bins as factor later on when binding with fossil data


# Calculate lags 1 to 10
isotemp2 <- isotemp2 %>%  
    mutate(lag1 = lag(Temp.mean, order_by = bins), 
           lag2 = lag(lag1, order_by = bins),
           lag3 = lag(lag2, order_by = bins),
           lag4 = lag(lag3, order_by = bins),
           lag5 = lag(lag4, order_by = bins),
           lag6 = lag(lag5, order_by = bins),
           lag7 = lag(lag6, order_by = bins),
           lag8 = lag(lag7, order_by = bins),
           lag9 = lag(lag8, order_by = bins),
           lag10 = lag(lag9, order_by = bins)) 

# clean up
rm(list=ls()[! ls() %in% c("dat_range", "stages", "isotemp", "isotemp2")])



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
  arrange(desc(bins), .by_group = T) 


# Now bind the two
dat_temp <- full_join(dat_safe, isotemp2)

#add taxonomical levels using the dummy and environmental preference
dumbo <- dumbo[, 1:5]

# Now bind the two
dat_temp <- full_join(dat_temp, dumbo)

# order it properly
dat_temp <- dat_temp %>%
  dplyr::select(phylum, class, order, family, genus, bins, age, origination, Temp.mean, 
                change.prev, lag1:lag10) %>%
  transform(bins = as.numeric(as.character(bins)))

# 7 Calculate multiple long-term temperature trends with varying l --------


# Use lm() to determine slope of the long-termtemperature trends
# If we want to investigate the interaction between long term change and short term,
# we need to exclude the short term bin from the long term, and thus shift the long 
# term lm result up by +1.

for (i in unique(isotemp$Stage)) {
  sub1<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-1, i)])
  lin1<-lm(Temp~age, data=sub1)
  isotemp2[isotemp2$Stage==i+1, "trend.st1"]<- -lin1$coefficients[2]
  sub2<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-2, i)])
  lin2<-lm(Temp~age, data=sub2)
  isotemp2[isotemp2$Stage==i+1, "trend.st2"]<- -lin2$coefficients[2]
  sub3<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-3, i)])
  lin3<-lm(Temp~age, data=sub3)
  isotemp2[isotemp2$Stage==i+1, "trend.st3"]<- -lin3$coefficients[2]
  sub4<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-4, i)])
  lin4<-lm(Temp~age, data=sub4)
  isotemp2[isotemp2$Stage==i+1, "trend.st4"]<- -lin4$coefficients[2]
  sub5<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-5, i)])
  lin5<-lm(Temp~age, data=sub5)
  isotemp2[isotemp2$Stage==i+1, "trend.st5"]<- -lin5$coefficients[2]
  sub6<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-6, i)])
  lin6<-lm(Temp~age, data=sub6)
  isotemp2[isotemp2$Stage==i+1, "trend.st6"]<- -lin6$coefficients[2]
  sub7<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-7, i)])
  lin7<-lm(Temp~age, data=sub7)
  isotemp2[isotemp2$Stage==i+1, "trend.st7"]<- -lin7$coefficients[2]
  sub8<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-8, i)])
  lin8<-lm(Temp~age, data=sub8)
  isotemp2[isotemp2$Stage==i+1, "trend.st8"]<- -lin8$coefficients[2]
  sub9<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-9, i)])
  lin9<-lm(Temp~age, data=sub9)
  isotemp2[isotemp2$Stage==i+1, "trend.st9"]<- -lin9$coefficients[2]
  sub10<-filter(isotemp, isotemp$Stage %in% isotemp$Stage[between(isotemp$Stage, i-10, i)])
  lin10<-lm(Temp~age, data=sub10)
  isotemp2[isotemp2$Stage==i+1, "trend.st10"]<- -lin10$coefficients[2]
}

# clean up
rm(list=ls()[! ls() %in% c("dat_temp", "isotemp2")])

# Set bins back to factor for joining
dat_temp <- dat_temp %>% 
  transform(bins = as.factor(as.character(bins)))

# Now bind the two
dat_final <- full_join(dat_temp, isotemp2)

# remove reduntant columns
dat_final <- dat_final %>% 
  select(-Stage, -n)

# save it to your working directory
save(dat_final, file = here("data/final_data.RData"))



