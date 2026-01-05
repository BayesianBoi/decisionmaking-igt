rm(list=ls())

library(dplyr)
library(reshape2)
library(foreach)

setwd("Data/Raw/")

load("Steingroever_2014/IGTdata.rdata")

# Create long data
all_dat <- melt(choice_100, c("subjID", "trial"), value.name="deck")
all_dat$loss <- melt(lo_100, c("subjID", "trial"), value.name="loss")$loss
all_dat$gain <- melt(wi_100, c("subjID", "trial"), value.name="gain")$gain
all_dat <- all_dat[order(all_dat$subjID),]

# Modify IDs
all_dat$subjID <- as.numeric(gsub(x = all_dat$subjID, pattern = "^Subj_", replacement = ""))

# Join IDs with studies
all_dat <- full_join(all_dat, as.data.frame(index_100) %>%
                       mutate(Subj = as.numeric(as.character(Subj))), 
                     by = c("subjID" = "Subj")) %>%
  mutate(trial = as.numeric(gsub(x = trial, 
                                 pattern = "^Choice_", 
                                 replacement = "")))

# Set study ID to integer
all_dat$group_id <- as.integer(as.factor(all_dat$Study))

# Read in Ahn2014 and Fridberg data
ahn2014_HC <- read.table("Ahn_2014/IGTdata_HC.txt", header=T) %>% 
  mutate(Study = "ahn2014_HC",
         group_id = 8)
ahn2014_Amph <- read.table("Ahn_2014/IGTdata_Amph.txt", header=T) %>% 
  mutate(Study = "ahn2014_Amph",
         group_id = 9)
ahn2014_Hero <- read.table("Ahn_2014/IGTdata_Hero.txt", header=T) %>% 
  mutate(Study = "ahn2014_Hero",
         group_id = 10)
frid_HC <- read.table("Fridberg_2010/IGTdata_HC.txt", header=T) %>% 
  mutate(Study = "frid_HC",
         group_id = 11) %>%
  select(-deckCopy)
frid_Cbis <- read.table("Fridberg_2010/IGTdata_Cbis.txt", header=T) %>% 
  mutate(Study = "frid_Cbis",
         group_id = 12) %>%
  select(-deckCopy)

all_dat <- full_join(all_dat, rbind(ahn2014_HC,
                                    ahn2014_Amph,
                                    ahn2014_Hero,
                                    frid_HC,
                                    frid_Cbis))
saveRDS(all_dat, file = "../Preprocessed/all_dat_12_studies_IGT_Data_stan_ready.rds")

dataList <- foreach(i=seq_along(unique(all_dat$group_id))) %do% {
  # Subset study data
  rawdata <- all_dat %>% filter(group_id==i) 
  
  # Individual Subjects
  subjList <- unique(rawdata[,"subjID"])  # list of subjects x blocks
  numSubjs <- length(subjList)  # number of subjects
  
  Tsubj <- as.vector( rep( 0, numSubjs ) ) # number of trials for each subject
  group_id <- as.vector( rep( 0, numSubjs ) ) # group identifier for each subject
  
  for ( sIdx in 1:numSubjs )  {
    curSubj     <- subjList[ sIdx ]
    Tsubj[sIdx] <- sum( rawdata$subjID == curSubj )  # Tsubj[N]
    group_id[sIdx] <- rawdata$group_id[rawdata$subjID == curSubj][1]
  }
  
  maxTrials <- max(Tsubj)
  
  RLmatrix <- SRLmatrix <- array( 0, c(numSubjs, maxTrials ) )
  Ydata    <- array(1, c(numSubjs, maxTrials) )
  
  for ( subjIdx in 1:numSubjs )   {
    #number of trials for each subj.
    useTrials                      <- Tsubj[subjIdx]
    currID                         <- subjList[ subjIdx ]
    rawdata_curSubj                <- subset( rawdata, rawdata$subjID == currID )
    RLmatrix[subjIdx, 1:useTrials] <- rawdata_curSubj[, "gain"] -1 * abs( rawdata_curSubj[ , "loss" ])
    
    for ( tIdx in 1:useTrials ) {
      Y_t                     <- rawdata_curSubj[ tIdx, "deck" ] # chosen Y on trial "t"
      Ydata[ subjIdx , tIdx ] <- Y_t
      # For binarizing 
      if ( RLmatrix[subjIdx, tIdx] > 0 ) {
        SRLmatrix[subjIdx, tIdx] <- 1
      } else if ( RLmatrix[subjIdx, tIdx] == 0 ) {
        SRLmatrix[subjIdx, tIdx] <- 0
      } else {
        SRLmatrix[subjIdx, tIdx] <- -1
      }
    }
  }
  
  list(
    N       = numSubjs,
    T       = maxTrials,
    Tsubj   = Tsubj ,
    Srewlos = SRLmatrix,     
    rewlos  = RLmatrix / 100,
    ydata   = Ydata
  )
}
names(dataList) <- unique(all_dat$Study)

# Save all data
saveRDS(dataList, file = "../Preprocessed/dataList_12datasets_stan_ready.rds")
