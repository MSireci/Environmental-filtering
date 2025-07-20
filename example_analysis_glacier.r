rm(list=ls())
require(tidyverse)
require(dplyr)
library(tidyr)
source("correlation_analysis.r")

# EXAMPLE FOR CROSSSECTIONAL DATA OF THE GLACIER BIOME

load("./crossscdata.RData")


#1- Obtain the OTU abundances for the desired biome

# the table  will report the count, abundances and occurences of each OTU. Furthermore, we keep the code identifying  the dataset (project_id), the sample and the run id.
abddata_Glacier<- datatax %>% filter(classification=="Glacier") %>% mutate( nruns = n_distinct(run_id) ) %>% 
  group_by(otu_id) %>% 
  mutate( mf = mean(count/nreads), occ = n_distinct(run_id)/nruns ) %>% 
  ungroup() %>% filter(nreads > 10^4 )%>% select(otu_id, count, project_id, sample_id, run_id, classification, nreads, nruns, mf, occ)
  rm(datatax)



# 2-Obtain the phylogenetic distances of the communities in the biomes

# The information of each dataset is encoded in the project id classification

project_id<- (abddata_Glacier %>% select(project_id) %>%unique())$project_id

# Construct the phylogenetic distance file of the communities of this datset (project_id) and order them in n=17 bin

dist_bin<-distance_bin(project_id,abddata,17)

# For simplicity and  for example use the distance dataset ready for the Glacier biome charged in the dropbox folder reported in the readme file

load("./dist_bin_ERP017997_Glacier.RData")

# 3- Obtain the fluctuation estimators

d_Glacier<-calculate_abd(abddata_Glacier)
rm(abddata_Glacier)
  
# 4- Construct the pearson correlation for each species couples and link it to their phylogenetic distance

q_Glacier<-calculate_p2(d_Glacier, 17,dist_bin_ERP017997_Glacier, 10, 20, 15, 10)
rm(dist_bin_ERP017997_Glacier)
# see the description in correlation_analysis.r for an explanation of the last four input parameters. They are just parameters that split the file in subsets easier to apply left.join

save(q_Glacier, file="./observable_Glacier.RData")

# 5- Obtain the average decay

p_Glacier<-av_q(q_Glacier)
p_Glacier<- p_Glacier %>% mutate(env="Glacier")
save(p_Glacier, file="./correlation_Glacier.RData")

# 6- Null/Randomized model= Check if the correlation decay disappears when the distances are randomized

p_rand_Glacier<-randomization(q_Glacier, rand, 1.96)
rm(q_Glacier)

p_rand_Glacier <-p_rand_Glacier %>%mutate(env="Glacier")
save(p_rand_Glacier, file="./correlation_Glacier_rand.RData")

# 7- plots
require(ggplot2)
require(viridis)

cbp2 <- c( "#7D3102", "#23C710", "#E87203",
           "#f9280B","#F3F011", "#8C0fAD","#0b4cf9", "#0ecab0")

#choose environment
p<-p_Glacier
# Choose data type
p1<- p %>% filter(type=="Data")  %>% filter(nbin>10^3) %>% mutate(corr_th=exp(-3.5*exp(ld*1/3)))

#p1$env<-as.factor(p$env)

# Choose observable
p2 <- ggplot(p1 %>% filter(observable=="eta3"))+ theme()+
  
  geom_hline( yintercept = 0, color = "gray", size = 2 ) +
  aes(
    x=exp(ld),
    y=(Corr),
    shape=env,
    color=env
    
    
  )+
  geom_point(,size=3.5, stroke=2)+
  scale_colour_manual(values=cbp2)+
  scale_shape_manual(values=1:8)+
  #geom_errorbar(data=p %>%filter(type=="Null Model") %>% filter(nbin> 10^3), aes(ymin=eta1-error1, ymax=eta1+error1, x=exp(ld)), width=0.1)+
  geom_smooth(data=p1, colour="black", aes(x=exp(ld), y=(corr_th),group = 1, weight = nbin), se = FALSE, size=1.5)+
  scale_x_log10( " " ) +
  scale_y_continuous(" " )+
  mytheme_main+labs(colour=" ", shape=" ")+theme(legend.position = c("none"))+
  labs(colour=" ", shape=" ")
  #facet_wrap(~ observable, nrow = 2,  scale = "free_y")


p2

