rm(list=ls())
require(tidyverse)
require(dplyr)
library(tidyr)
source("correlation_analysis.r")

# EXAMPLE FOR CROSSSECTIONAL DATA OF THE GLACIER BIOME
# UPLOAD THE DATA
load("./crosssecdata.RData")


#1- Obtain the OTU abundances for the desired biome

## Calculate abundances and keep the zeros

abddata0<- datatax %>% filter(classification=="Glacier")  %>% filter(nreads > 10^4 )
sp<- unique(abddata0$otu_id)

sample_tot<- unique(abddata0$sample_id)

abddata0_new<-abddata0[0,]

for(sample_ch in sample_tot ){
  
  abddata0_1<-abddata0 %>% filter(sample_id==sample_ch)
  proj_id<-unique(abddata0_1$project_id)
  run<-unique(abddata0_1$run_id)
  reads<-unique(abddata0_1$nreads)
  env<-unique(abddata0_1$classification)
  sp1<-unique(abddata0_1$otu_id)
  
  absent<- setdiff(sp,sp1)
  absent0<- as.data.frame(absent) %>% mutate(count=0,project_id=proj_id, sample_id=sample_ch, run_id=run,nreads=reads, classification=env)
  colnames(absent0)[1]<-"otu_id"
  abddata0_2<- rbind(abddata0_1,absent0) %>% unique()
  abddata0_new<- rbind(abddata0_new,abddata0_2)
  rm(abddata0_1,abddata0_2,absent,absent0)
}

abddata<- abddata0_new %>% filter(classification=="Glacier") %>% mutate( nruns = n_distinct(run_id) ) %>% 
  group_by(otu_id) %>% mutate( mf = mean(count/nreads), occ = n_distinct(run_id)/nruns ) %>% 
  ungroup() %>% filter(nreads > 10^4 )%>% select(otu_id, count, project_id, sample_id, run_id, classification, nreads, nruns, mf, occ)

save(abddata,file="./abd_Glacier0.RData")
rm(datatax,abddata0_new,sample_tot,sample_ch,env,proj_id,reads,run,sp,sp1)



# 2-Obtain the phylogenetic distances of the communities in the biomes

# The information of each dataset is encoded in the project id classification

project_id<- (abddata %>% select(project_id) %>%unique())$project_id

# Construct the phylogenetic distance file of the communities of this datset (project_id) and order them in n=17 bin

dist_bin<-distance_bin(project_id,abddata,17)

# For simplicity and  for example use the distance dataset ready for the Glacier biome charged in the dropbox folder reported in the readme file

load("./dist_bin_ERP017997_Glacier.RData")

# 3- Obtain the fluctuation estimators

d_Glacier<-calculate_abd(abddata)
save(d_Glacier, file="./estimator_flu_Glacier0.RData")
rm(abddata_Glacier)

load("./estimator_flu_Glacier0.RData")

# 4- Construct the pearson correlation for each species couples and link it to their phylogenetic distance

q_Glacier<-calculate_p2(d_Glacier,17,dist_bin_ERP017997_Glacier)

rm(dist_bin_ERP017997_Glacier)

save(q_Glacier, file="./observable_Glacier.RData")

# 5- Obtain the average decay

p_Glacier<-av_q(q_Glacier)

p_Glacier<- p_Glacier %>% mutate(env="Glacier")
save(p_Glacier, file="./correlation_Glacier.RData")

# 6- Null/Randomized model= Check if the correlation decay disappears when the distances are randomized
# rand is the number of randomizations runs

p_rand_Glacier<-randomization(q_Glacier, 10, 1.96)
rm(q_Glacier)

p_rand_Glacier <-p_rand_Glacier %>%mutate(env="Glacier")
save(p_rand_Glacier, file="./correlation_Glacier_rand.RData")


############################# Decay plots ############################################
mytheme_main <- theme_bw() + theme(
  legend.title  = element_text(family="Helvetica", size=25, color = "black"
  ),
  legend.key = element_blank(),
  legend.text  = element_text(family="Helvetica", size=30, color = "black"),
  panel.background = element_rect(fill="transparent"),
  #plot.background = element_rect(fill="transparent", colour = NA),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=25, color = "black"),
  panel.border = element_blank(),
  axis.title = element_text( family="Helvetica", size=35, color = "black"),
  axis.text = element_text( family="Helvetica", size=12, color = "black"),
  axis.line = element_line(size = 1., color = "black"),
  axis.ticks = element_line(size = 1.,color = "black"),
  legend.background = element_rect(fill="transparent")
)



color<-c("red2","blue2","darkorange2","springgreen4","deeppink2","gold","mediumpurple2","tan4")
cbp2 <- c( "#7D3102", "#23C710", "#E87203",
           "#f9280B","#F3F011", "#8C0fAD","#0b4cf9", "#0ecab0")

p<-p_glacier

p1<- p %>% filter(type=="Data")  %>% filter(nbin>10^3) %>% mutate(corr_th=exp(-3.5*exp(ld*1/3)))

#p1$env<-as.factor(p$env)
p2 <- ggplot(p1 %>% filter(observable=="eta3"))+ theme()+
  
  geom_hline( yintercept = 0, color = "gray", size = 2 ) +
  aes(
    x=exp(ld),
    y=(Corr),
    #shape=env,
   # color=env
    
    
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

