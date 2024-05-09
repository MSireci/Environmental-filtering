rm(list=ls())
require(tidyverse)
require(dplyr)
library(tidyr)
source("correlation_analysis.r")

# EXAMPLE FOR CROSSSECTIONAL DATA OF THE GLACIER BIOME
load("./crossscdata.RData")

# 1-Obtain the phylogenetic distances
# From the file select the OTU present in the desired biome

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

load("./dist_bin_ERP017997_Glacier.RData")

# 2- Obtain the Abundances estimators

d_Glacier<-calculate_abd(abddata_Glacier)
rm(abddata_Glacier)
  
# 3- Construct the pearson correlation for each species couples and link it to their phylogenetic distance

q_Glacier<-calculate_p2(d_Glacier, 17,dist_bin_ERP017997_Glacier, 10, 20, 15, 10)
rm(dist_bin_ERP017997_Glacier)
# see the description in correlation_analysis for an explanation of the last four input parameters. They are just parameters that split the file in subsets easier to apply left.join

save(q_Glacier, file="./observable_Glacier.RData")

# 3- Obtain the average decay

p_Glacier<-av_q(q_Glacier)
p_Glacier<- p_Glacier %>% mutate(env="Glacier")
save(p_Glacier, file="./correlation_Glacier.RData")

# 4- Null/Randomized model= Check if the correlation decay disappears when the distances are randomized

p_rand_Glacier<-randomization(q_Glacier, rand, 1.96)
rm(q_Glacier)

p_rand_Glacier <-p_rand_Glacier %>%mutate(env="Glacier")
save(p_rand_Glacier, file="./correlation_Glacier_rand.RData")

p_Gut<-rbind(p_Gut, p_rand_Gut)
save(p_Gut, file="./correlation_Gut2_new.RData")
