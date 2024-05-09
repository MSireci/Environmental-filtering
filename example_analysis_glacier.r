rm(list=ls())
require(tidyverse)
require(dplyr)
library(tidyr)
source("correlation_analysis.r")

# EXAMPLE FOR CROSSSECTIONAL DATA OF THE GLACIER BIOME
load("./crossscdata.RData")

# From the file select the species present in the desired biome

abddata<- proj %>% filter(classification=="Glacier") %>% mutate( nruns = n_distinct(run_id) ) %>% 
  group_by(otu_id) %>% 
  mutate( mf = mean(count/nreads), occ = n_distinct(run_id)/nruns ) %>% 
  ungroup() %>% filter(nreads > 10^4 )%>% select(otu_id, count, project_id, sample_id, run_id, classification, nreads, nruns, mf, occ)
  rm(proj)

# The information of each dataset is encoded in the project id classification

project_id<- (abddata %>% select(project_id) %>%unique())$project_id

# Construct the phylogenetic distance file of the communities of this datset (project_id) and order them in n=17 bin

dist_bin<-distance_bin(project_id,abddata,17)

# For simplicity and  for example use the distance dataset ready for the Glacier biome charged in the dropbox folder reported in the readme file


