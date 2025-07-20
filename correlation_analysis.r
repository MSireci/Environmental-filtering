
####################### Environmental filtering- Sireci, Muñoz, Grilli ###########

require(tidyverse)
require(dplyr)
require(phyloseq)
require(adephylo)
require(ape)
require(tidytree)

#######################- Functions- #################################################


###########################################-1-#calculate philogenetic distance 


# a) Construct phylogenetic tree

tree<-function( class, proj){
  
  f <- read_tree_greengenes("./97_otus.tree")
  
  
  a<-as_tibble(f$tip.label)
  colnames(a)<-c("otu_id")
  tips_tree<-a %>% mutate(otu_id=as.integer(otu_id))
  rm(a)
  
  com<-unique(proj %>% filter(nreads>10^4, project_id==class)%>% select(sample_id))$sample_id
  run<-unique(proj %>% filter(nreads>10^4, project_id==class)%>% select(run_id,sample_id))
  
  
  for(sid in com){
    
    run1<-c((unique(run %>% filter(sample_id==sid)))$run_id)[1]  
    
    
    #species that I need
    sp_class<-  c((proj%>% filter( project_id ==class, sample_id==sid, run_id==run1) %>% select(otu_id) %>% unique())$otu_id %>% unique()) 
    
    
    
    #Otu_id to drop
    
    to_drop<-c((tips_tree%>% filter(otu_id, !otu_id %in% sp_class) %>% mutate(otu_id=as.character(otu_id)))$otu_id )
    rm( sp_class)
    
    
    f2<-drop.tip(f, to_drop)
    rm( to_drop)
    
    filename <- paste("./trees/tree_", class,"_", sid, "_", run1, ".tsv", sep = "")
    write.tree(f2, filename)
    
    
  }
  return(0)
}

# b) calculate distance all bin together

distance <-function(project_id_name, abddata){
  
  
  
  com<-unique(abddata%>% select(sample_id))$sample_id
  run<-c(rep(0,length(com)))
  
  df_dist<-data.frame()
  
  for(sid in com){
    sid
    print(sid)
    #for each community select one run, the first
    com2<- abddata %>% filter(sample_id==sid) %>% select(otu_id,run_id)
    rid <-(unique(com2$run_id))
    
    
    filename <- paste("./trees/tree_", project_id_name,"_", sid, "_", rid, ".tsv", sep = "")
    f <- read_tree_greengenes(filename)
    system.time(d2<-cophenetic(f))
    
    #calculate distances with cutoff at ten decimals, and not reporting the doubles
    df_dist <-data.frame(otu2=rownames(d2)[row(d2)], otu1=colnames(d2)[col(d2)], dist=c(d2), run_id= rid) %>% filter(as.numeric(otu1) > as.numeric(otu2))   %>% 
      mutate( dist = round(dist,10) )
    rm(d2)
  }
  return(df_dist)
}
# c) Calculate distance one bin at the time

distance_bin <-function(project_id_name, abddata,nbins){
  
  
  df_dist <-data.frame()
  com<-unique(abddata%>% select(sample_id))$sample_id
  run<-c(rep(0,length(com)))
  
  
  df_dist<-data.frame()
  ####Here I divide in three,n because it is too large
  for(sid in com){
    
    print(sid)
    
    #for each community select one run, the first
    com2<- abddata %>% filter(sample_id==sid)%>% select(otu_id,run_id,project_id)
    rid <-(unique(com2$run_id))[1]#<-qui prende solo un run_id per comunità
    
    
  
    filename <- paste("./trees/tree_", project_id_name,"_", sid, "_", rid, ".tsv", sep = "")
    f <- read_tree_greengenes(filename)
    system.time(d2<-cophenetic(f))
    
    #calculate distances with cutoff at ten decimals, and not reporting the doubles
    df_dist <-rbind( data.frame(otu2=rownames(d2)[row(d2)], otu1=colnames(d2)[col(d2)], dist=c(d2))%>% filter(as.numeric(otu1) > as.numeric(otu2)) %>% 
                       mutate( dist = round(dist,10)),df_dist )
  
    rm(d2)
    df_dist<-df_dist  %>% unique()%>%  filter(!is.na(dist))
    #bin the distance
    
  }
  
  print("end")
  

  df_dist_bin<-df_dist  %>% unique() %>%  filter(!is.na(dist))%>% mutate( otu1 = as.integer(as.character(otu1)) ,
                                                                          otu2 = as.integer(as.character(otu2)) ,
                                                                          dist = as.double(dist))
  print("a")  
  
  df_dist_bin<-df_dist_bin %>% filter(otu1>otu2)  %>% unique()
  
  
  print("b")
  
  df_dist_bin<-df_dist_bin %>%   mutate( b = as.integer( (log(dist)-min(log(dist)))/( max(log(dist))-min(log(dist)))*nbins ))        
  
  
  
  print("end2")
  return(df_dist_bin)
}

## d) Tree randomization
tree_randomization<-function(d2)
{
  
  
  ## randomizzazione albero
  otudist <- d2 %>% select( otu1, otu2, dist, b) %>% distinct()
  
  oturatio <- d2 %>% select( otu1, otu2, ratio1, ratio2, sample_couple)
  
  oturndmap <- otudist %>% select(otu1) %>% distinct() %>% mutate( otu1rnd = sample(otu1) ) ## mappa da otu a otu rnd
  
  otudistrnd <- otudist %>% left_join( oturndmap) %>% left_join( oturndmap %>% rename( otu2 = otu1, otu2rnd = otu1rnd ) ) %>% 
    mutate( otu1 = otu1rnd, otu2 = otu2rnd) %>% select( -otu1rnd, -otu2rnd )
  
  dtot1rndtree <- oturatio %>% left_join( otudistrnd ) 
  
  return(dtot1rndtree)
  
}
# e) Distances randomization
distances_rand<-function(dtot1)
{
  otudist <- dtot1 %>% select( otu1, otu2, dist, b) %>% distinct()
  oturatio <- dtot1 %>%   select(otu1, otu2, sample_id, q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b) 
  nd = length(otudist$dist)
  ord = sample(1:nd)
  dv = otudist$dist[ord]
  bv = otudist$b[ord]
  otudistrnd <- otudist %>% mutate( dist = dv, b = bv)  
  dtot1rnddist <- oturatio %>% left_join( otudistrnd )
  
  
  return( dtot1rnddist)
}



############################################################## Correlation Analysis
# The abjective is to have a file relating distances and correlation in each bin.


#a) Calculate statistical estimators
# As a first step, for each species of the considered biome we calculate the five estimators q1,q2, q3, q4, q5.
calculate_abd<-function(abddata)
{
  
  d<-abddata %>% select(otu_id, count, nreads, project_id, sample_id, run_id) %>% mutate(abd=count/nreads, log_abd=log(count/nreads))
  
  rm(abddata)
  
  d<- d %>% group_by(otu_id) %>% mutate( av_abd=mean(abd), std_abd=sqrt(mean(abd*abd)-mean(abd)*mean(abd)), 
                                         av_log_abd=mean(log_abd), std_log_abd=sqrt(mean(log_abd*log_abd)-mean(log_abd)*mean(log_abd)), rank=rank(abd)-1) %>%
    filter(std_abd>0) %>% mutate(rank=rank/max(rank)) %>%ungroup()
  
  
  d_tot<-d %>% mutate(q1=(abd-av_abd)/av_abd, q2=(count-av_abd*nreads)/(av_abd*nreads), q3=(abd-av_abd)/std_abd, q4=(log_abd-av_log_abd)/std_log_abd, q5=2*rank-1) %>%
    select(otu_id, project_id, sample_id, q1, q2, q3,q4, q5)
  rm(d)
  return(d_tot)
}


# b) From species estimators for each coubple, we  construct a data frame relating for each species couuple its correlation quantifiers q1 etc, with the phylogenetic distance. This is done for each bin.
calculate_p<-function(d_tot, dist, bin,h, c,h2, k)
{
  # EXPLICATION OF INPUTS h, c, h2, k. 
  # the dataframes are really heavy, and hence to perform the "left-join " with the distance file  I need to cut them in subsets and perform a loop.
  # The amount of species changes in each bin (17 in total), largest the bin more species it does contation. 
  # Hence, if the considered bin is smaller than h2 (normally h2=15) and larger than h (h=10) , I can subset the file in pieces of c containing lines (normally c=20).
  # For bin larger than h2 I considered a smaller division k, normally k=15.
  
  df<-as.data.frame(dist) %>% filter(b==bin) %>%  mutate( otu1 = as.integer(as.character(otu1)) ,
                                                          otu2 = as.integer(as.character(otu2)) )
  
  rm(dist)
  
  print(bin)

  
  
  
  
  
  df2<-tibble(
    
    
    otu1=as.integer(),otu2=as.integer(), dist=as.double(), b=as.integer(), sample_id=as.integer(), 
    q1_a=as.double(), q2_a=as.double(), q3_a=as.double(), q4_a=as.double(), q5_a=as.double(), 
    q1_b=as.double(), q2_b=as.double(), q3_b=as.double(), q4_b=as.double(), q5_b=as.double())
  
  
  
  if(bin>=h & bin<h2)
  {
    # now we make a for cycle with the samples_coumple
    sample_tot<-d_tot %>% select(sample_id) %>% unique()
    sample_tot<-c(sample_tot$sample_id)
    
    i<-1
    
    
    # c parameters for division of cycles
    c
    
    a<-as.integer(as.integer(length(sample_tot))/c)
    # b<-as.integer(length(couple))-c*a
    
    
    
    
    
    
    while(i<c*a )
    {
      print(i)
      
      d_tot1<-d_tot %>% group_by(sample_id) %>% filter(sample_id %in% sample_tot[i:i+c-1]) %>% ungroup()
      
      
      
      
      
      df1<- df %>%  left_join( d_tot1%>% rename(otu1 = otu_id, q1_a = q1, q2_a = q2, q3_a = q3, q4_a = q4, q5_a = q5) %>% mutate(otu1 = as.integer(otu1) )  ) 
      
      
      df1<- df1 %>%  left_join( d_tot1%>% rename(otu2 = otu_id, q1_b = q1, q2_b = q2, q3_b = q3, q4_b = q4, q5_b = q5) %>% mutate(otu2 = as.integer(otu2) ), 
                                by=c("otu2"="otu2", "sample_id"="sample_id")  ) %>%
        select(otu1, otu2, dist, b, sample_id, q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b) 
      rm(d_tot1)
      
      df1<-df1 %>% filter(!is.na(q1_a), !is.na(q1_b) )  %>% distinct(otu1, otu2, sample_id, .keep_all=TRUE)
      
      df2<-rbind(df2,df1)
      rm(df1)
      i<-i+c
      
    }
    d_tot1<-d_tot %>% group_by(sample_id) %>% filter(sample_id %in% sample_tot[i:length(sample_tot)]) %>% ungroup()
    
    
    
    
    rm(d_tot)
    df1<- df %>%  left_join( d_tot1%>% rename(otu1 = otu_id, q1_a = q1, q2_a = q2, q3_a = q3, q4_a = q4, q5_a = q5) %>% mutate(otu1 = as.integer(otu1) )  ) 
    
    
    df1<- df1 %>%  left_join( d_tot1%>% rename(otu2 = otu_id, q1_b = q1, q2_b = q2, q3_b = q3, q4_b = q4, q5_b = q5) %>% mutate(otu2 = as.integer(otu2) ), 
                              by=c("otu2"="otu2", "sample_id"="sample_id")  ) %>%
      select(otu1, otu2, dist, b, sample_id, q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b) 
    
    rm(d_tot1)
    df1<-df1 %>% filter(!is.na(q1_a), !is.na(q1_b) )  %>% distinct(otu1, otu2, sample_id, .keep_all=TRUE)
    
    df<-rbind(df2,df1) 
    rm(df1,df2)
    
    
    
  }
  else if (bin>=h2){
    # now we make a for cycle with the samples_coumple
    sample_tot<-d_tot %>% select(sample_id) %>% unique()
    sample_tot<-c(sample_tot$sample_id)
    
    i<-1
    
    
    #k parameter for cycle division
    k
    
    a<-as.integer(as.integer(length(sample_tot))/k)
    # b<-as.integer(length(couple))-c*a
    
    
    
    
    
    
    while(i<k*a )
    {
      print(i)
      
      d_tot1<-d_tot %>% group_by(sample_id) %>% filter(sample_id %in% sample_tot[i:i+k-1]) %>% ungroup()
      
      
      
      
      
      df1<- df %>%  left_join( d_tot1%>% rename(otu1 = otu_id, q1_a = q1, q2_a = q2, q3_a = q3, q4_a = q4, q5_a = q5) %>% mutate(otu1 = as.integer(otu1) )  ) 
      
      
      df1<- df1 %>%  left_join( d_tot1%>% rename(otu2 = otu_id, q1_b = q1, q2_b = q2, q3_b = q3, q4_b = q4, q5_b = q5) %>% mutate(otu2 = as.integer(otu2) ), 
                                by=c("otu2"="otu2", "sample_id"="sample_id")  ) %>%
        select(otu1, otu2, dist, b, sample_id, q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b) 
      rm(d_tot1)
      
      df1<-df1 %>% filter(!is.na(q1_a), !is.na(q1_b) )  %>% distinct(otu1, otu2, sample_id, .keep_all=TRUE)
      
      df2<-rbind(df2,df1)
      rm(df1)
      i<-i+k
      
    }
    d_tot1<-d_tot %>% group_by(sample_id) %>% filter(sample_id %in% sample_tot[i:length(sample_tot)]) %>% ungroup()
    
    
    
    
    rm(d_tot)
    df1<- df %>%  left_join( d_tot1%>% rename(otu1 = otu_id, q1_a = q1, q2_a = q2, q3_a = q3, q4_a = q4, q5_a = q5) %>% mutate(otu1 = as.integer(otu1) )  ) 
    
    
    df1<- df1 %>%  left_join( d_tot1%>% rename(otu2 = otu_id, q1_b = q1, q2_b = q2, q3_b = q3, q4_b = q4, q5_b = q5) %>% mutate(otu2 = as.integer(otu2) ), 
                              by=c("otu2"="otu2", "sample_id"="sample_id")  ) %>%
      select(otu1, otu2, dist, b, sample_id, q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b) 
    
    rm(d_tot1)
    df1<-df1 %>% filter(!is.na(q1_a), !is.na(q1_b) )  %>% distinct(otu1, otu2, sample_id, .keep_all=TRUE)
    
    df<-rbind(df2,df1) 
    rm(df1,df2)
    
    
    
    
  }
  else{
    df<- df %>%  left_join( d_tot%>% rename(otu1 = otu_id, q1_a = q1, q2_a = q2, q3_a = q3, q4_a = q4, q5_a = q5) %>% mutate(otu1 = as.integer(otu1) )  ) 
    
    
    df<- df %>%  left_join( d_tot%>% rename(otu2 = otu_id, q1_b = q1, q2_b = q2, q3_b = q3, q4_b = q4, q5_b = q5) %>% mutate(otu2 = as.integer(otu2) ), 
                            by=c("otu2"="otu2", "sample_id"="sample_id")  ) %>%
      select(otu1, otu2, dist, b, sample_id, q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b) 
    rm(d_tot)
    
    df<-df %>% filter(!is.na(q1_a), !is.na(q1_b) )  %>% distinct(otu1, otu2, sample_id, .keep_all=TRUE)
  }
  return(df)
}





###################### c) The operation above is repaeted for each bin.
calculate_p2<-function(dtot,bmax, dist, h, c, h2, k)
{
  
  
  g_q<-calculate_p(dtot, dist,0, h, c, h2,k)
  
  for(b in 1:bmax){
    g1<-calculate_p(dtot,dist, b ,h,c, h2, k)
    g_q<-rbind(g_q,g1)
    rm(g1)
  }
  return(g_q)
}


############################### d) Finally,  we calculated correlations for each couple and  perform averages over species couples in each bin. 


av_q<-function(g_q)
{
  c<-1.96
  
  g_q<-g_q %>%  mutate(p1=q1_a*q1_b,p2=q2_a*q2_b,p3=q3_a*q3_b, p4=q4_a*q4_b, p5=q5_a*q5_b, b=b ) %>% 
    select(otu1, otu2, dist, b, sample_id, p1, p2, p3, p4, p5)
  
  
  g<-g_q %>% group_by(b) %>% summarise(ld=mean(log(dist)), eta1=mean(p1), var_eta1=var(p1), eta2=mean(p2), var_eta2=var(p2),  eta3=mean(p3), var_eta3=var(p3), 
                                       eta4=mean(p4), var_eta4=var(p4),  eta5=mean(p5), var_eta5=var(p5), nbin=n()) %>% ungroup() 
  
  rm(g_q)
  
  
  g<-g %>% mutate(    error1=sqrt(var_eta1)*c/sqrt(nbin),
                      error2=sqrt(var_eta2)*c/sqrt(nbin),
                      error3=sqrt(var_eta3)*c/sqrt(nbin),
                      error4=sqrt(var_eta4)*c/sqrt(nbin),
                      error5=sqrt(var_eta5)*c/sqrt(nbin),
                      type="Data") %>%  select(b,ld,nbin, eta1, error1, eta2, error2, eta3, error3, eta4, error4,eta5, error5, type )
  
  g<-g %>%pivot_longer(c(eta1, eta2,eta3, eta4, eta5), names_to = "observable", values_to = "Corr")
  g1 <- g %>% filter(observable == "eta1" ) %>% rename(error=error1) %>% select(-c(error2, error3, error4, error5))
  g2 <- g %>% filter(observable == "eta2" ) %>% rename(error=error2) %>% select(-c(error1, error3, error4, error5))
  g3 <- g %>% filter(observable == "eta3" ) %>% rename(error=error3) %>% select(-c(error2, error1, error4, error5))
  g4 <- g %>% filter(observable == "eta4" ) %>% rename(error=error4) %>% select(-c(error2, error3, error1, error5))
  g5 <- g %>% filter(observable == "eta5" ) %>% rename(error=error5) %>% select(-c(error2, error3, error4, error1))
  
  g<-rbind(g1, g2) 
  rm(g1, g2)
  g<-rbind(g, g3)
  g<-rbind(g,g4)
  g<-rbind(g, g5)
  rm(g3,g4,g5)
  
  
  
  
  return(g)
}



################### e) Doing the same, but doing correlation inter and intra taxonomic phylia



tax_av<-function(observable, taxonomy, group){
  
  # Attaccare tassonomia al file
  
  q_tax<-observable  %>% left_join( taxonomy %>% rename(otu1=otu_id)) %>% rename(group1=as.character(group))
  rm(observable)
  q_tax<- q_tax %>% left_join( taxonomy %>% rename(otu2=otu_id)) %>% rename(group2=as.character(group))
  
  remove(taxonomy)
  
  #Aggiunger inter/intra
  q_tax<-q_tax %>% mutate(tax_rel=ifelse(group1==group2, "Intra", "Inter"))
  
  #unire le tassonomie
  q_tax<-q_tax %>% mutate(tax_couple=ifelse(group1>group2, paste(group1, group2, sep="-"), paste(group2, group1, sep="-"))) %>%select(-c(group1, group2))
  
  
  
  c<-1.96
  
  q_tax<-q_tax %>%  mutate(p1=q1_a*q1_b,p2=q2_a*q2_b,p3=q3_a*q3_b, p4=q4_a*q4_b, p5=q5_a*q5_b, b=b ) %>% 
    select(otu1, otu2, dist, b, sample_id, p1, p2, p3, p4, p5, tax_couple,tax_rel)
  
  
  g<-q_tax %>% group_by(b, tax_couple, tax_rel)  %>% summarise(ld=mean(log(dist)), eta1=mean(p1), var_eta1=var(p1), eta2=mean(p2), var_eta2=var(p2),  eta3=mean(p3), var_eta3=var(p3), 
                                                               eta4=mean(p4), var_eta4=var(p4),  eta5=mean(p5), var_eta5=var(p5), nbin=n()) %>% ungroup() 
  
  rm(q_tax)
  
  g<-g %>% mutate(    error1=sqrt(var_eta1)*c/sqrt(nbin),
                      error2=sqrt(var_eta2)*c/sqrt(nbin),
                      error3=sqrt(var_eta3)*c/sqrt(nbin),
                      error4=sqrt(var_eta4)*c/sqrt(nbin),
                      error5=sqrt(var_eta5)*c/sqrt(nbin),
                      type="Data") %>%  select(b,ld,nbin, eta1, error1, eta2, error2, eta3, error3, eta4, error4,eta5, error5, type, tax_couple,tax_rel )
  
  
  
  g<-g %>%pivot_longer(c(eta1, eta2,eta3, eta4, eta5), names_to = "observable", values_to = "Corr")
  g1 <- g %>% filter(observable == "eta1" ) %>% rename(error=error1) %>% select(-c(error2, error3, error4, error5))
  g2 <- g %>% filter(observable == "eta2" ) %>% rename(error=error2) %>% select(-c(error1, error3, error4, error5))
  g3 <- g %>% filter(observable == "eta3" ) %>% rename(error=error3) %>% select(-c(error2, error1, error4, error5))
  g4 <- g %>% filter(observable == "eta4" ) %>% rename(error=error4) %>% select(-c(error2, error3, error1, error5))
  g5 <- g %>% filter(observable == "eta5" ) %>% rename(error=error5) %>% select(-c(error2, error3, error4, error1))
  
  g<-rbind(g1, g2) 
  rm(g1, g2)
  g<-rbind(g, g3)
  g<-rbind(g,g4)
  g<-rbind(g, g5)
  rm(g3,g4,g5)
  
  
  
  
  
  
  
  
  return(g)
}


############################# g) Same as before, but for times series instead of cross sectional data


abd_dt<-function(abddata){
  
  
  ##  Calculate abd
  d<-abddata  %>% mutate(abd=count/nreads, log_abd=log(count/nreads))
  
  rm(abddata)
  
  d<- d %>% group_by(otu_id) %>% mutate( av_abd=mean(abd), std_abd=sqrt(mean(abd*abd)-mean(abd)*mean(abd)), 
                                         av_log_abd=mean(log_abd), std_log_abd=sqrt(mean(log_abd*log_abd)-mean(log_abd)*mean(log_abd)), rank=rank(abd)-1) %>%
    filter(std_abd>0) %>% mutate(rank=rank/max(rank)) %>%ungroup()
  
  
  d_tot<-d %>% mutate(q1=(abd-av_abd)/av_abd, q2=(count-av_abd*nreads)/(av_abd*nreads), q3=(abd-av_abd)/std_abd, q4=(log_abd-av_log_abd)/std_log_abd, q5=2*rank-1) %>%
    select(otu_id, project_id, sample_id, q1, q2, q3,q4, q5, classification, host_id, experiment_day)
  rm(d)
  return(d_tot)
}


calc_p_dt<-function(d_tot, dist, bin,dt, c)
{
  df<-as.data.frame(dist) %>% filter(b==bin) %>%  mutate( otu1 = as.integer(as.character(otu1)) ,
                                                          otu2 = as.integer(as.character(otu2)) )
  
  rm(dist)
  
  print(bin)
  
  
  d<-tibble( otu1=as.integer(),otu2=as.integer(), dist=as.double(), b=as.integer(), sample1=as.integer(),  sample2=as.integer(),
             q1_a=as.double(), q2_a=as.double(), q3_a=as.double(), q4_a=as.double(), q5_a=as.double(), 
             q1_b=as.double(), q2_b=as.double(), q3_b=as.double(), q4_b=as.double(), q5_b=as.double(), host_id=as.character(), classification=as.character(),
             dt=as.numeric())
  
  
  d_tot<-d_tot %>% mutate(time= experiment_day)
  
  d_tot2<-  d_tot  %>% transform(experiment_day= experiment_day-dt) %>% filter(!is.na(experiment_day)) 
  
  time<-c((d_tot %>% select(time) %>% unique())$time)
  
  
  time_max<-max(time)
  
  
  
  
  
  df1<- df %>%  left_join( d_tot%>% rename(otu1 = otu_id, q1_a = q1, q2_a = q2, q3_a = q3, q4_a = q4, q5_a = q5, time1=time, sample1=sample_id) %>% 
                             mutate(otu1 = as.integer(otu1) )  ) 
  
  df1<- df1 %>%  left_join( d_tot2%>% rename(otu2 = otu_id, q1_b = q1, q2_b = q2, q3_b = q3, q4_b = q4, q5_b = q5,time2=time, sample2=sample_id) %>% mutate(otu2 = as.integer(otu2) ), 
                            by=c("experiment_day"="experiment_day","otu2"="otu2" , "host_id"= "host_id", "classification"="classification")  )
  
  rm(d_tot2)
  
  
  df1<-df1 %>% filter(abs(time1-time2)==dt) %>% mutate(dt=dt) %>%
    select(otu1, otu2, dist, b,sample1, sample2,  q1_a, q2_a, q3_a, q4_a, q5_a, q1_b, q2_b, q3_b, q4_b, q5_b, host_id,dt, classification) 
  
  
  df1<-df1 %>% filter(!is.na(q1_a), !is.na(q1_b) )  %>% distinct(otu1, otu2, sample1, sample2, .keep_all=TRUE)
  
  d<-rbind(d,df1)
  rm(df1)
  
  
  return(d)
} 

calculate_p2_dt<-function(d_tot, bmax, dist, dt,c)
{
  
  
  g_q<-calc_p_dt(d_tot, dist, 0,dt,c)
  
  for(b in 1:bmax){
    g1<-calc_p_dt(d_tot, dist, b,dt,c)
    g_q<-rbind(g_q,g1)
    rm(g1)
  }
  return(g_q)
}

av_q_dt<-function(g_q)
{
  c<-1.96
  ##Come voglio fare le medie ?
  
  
  g_q<-g_q %>%  mutate(p1=q1_a*q1_b,p2=q2_a*q2_b,p3=q3_a*q3_b, p4=q4_a*q4_b, p5=q5_a*q5_b, b=b ) %>% 
    select(otu1, otu2, dist, b, sample1, sample2, p1, p2, p3, p4, p5, host_id, dt, classification)
  
  
  #media dentro ogni individuo
  g<-g_q %>% group_by(b, host_id,dt) %>% summarise(ld=mean(log(dist)), eta1=mean(p1), var_eta1=var(p1), eta2=mean(p2), var_eta2=var(p2),  eta3=mean(p3), var_eta3=var(p3), 
                                                   eta4=mean(p4), var_eta4=var(p4),  eta5=mean(p5), var_eta5=var(p5), nbin=n() ) %>% ungroup() 
  
  
  # Ora faccio anche media sugli individui:
  # q<-g_q %>% 
  #   group_by(b) %>% summarise(ld=mean(log(dist)), eta1=mean(p1), var_eta1=var(p1), eta2=mean(p2), var_eta2=var(p2),  eta3=mean(p3), var_eta3=var(p3), 
  #                            eta4=mean(p4), var_eta4=var(p4),  eta5=mean(p5), var_eta5=var(p5), nbin=n(), host_id="Average", dt=dt) %>% ungroup() 
  
  
  rm(g_q)
  g<-g %>% mutate(    error1=sqrt(var_eta1)*c/sqrt(nbin),
                      error2=sqrt(var_eta2)*c/sqrt(nbin),
                      error3=sqrt(var_eta3)*c/sqrt(nbin),
                      error4=sqrt(var_eta4)*c/sqrt(nbin),
                      error5=sqrt(var_eta5)*c/sqrt(nbin),
                      type="Data") %>%  select(b, host_id, dt, ld,nbin, eta1, error1, eta2, error2, eta3, error3, eta4, error4,eta5, error5, type )
  
  g<-g %>%pivot_longer(c(eta1, eta2,eta3, eta4, eta5), names_to = "observable", values_to = "Corr")
  g1 <- g %>% filter(observable == "eta1" ) %>% rename(error=error1) %>% select(-c(error2, error3, error4, error5))
  g2 <- g %>% filter(observable == "eta2" ) %>% rename(error=error2) %>% select(-c(error1, error3, error4, error5))
  g3 <- g %>% filter(observable == "eta3" ) %>% rename(error=error3) %>% select(-c(error2, error1, error4, error5))
  g4 <- g %>% filter(observable == "eta4" ) %>% rename(error=error4) %>% select(-c(error2, error3, error1, error5))
  g5 <- g %>% filter(observable == "eta5" ) %>% rename(error=error5) %>% select(-c(error2, error3, error4, error1))
  
  g<-rbind(g1, g2) 
  rm(g1, g2)
  g<-rbind(g, g3)
  g<-rbind(g,g4)
  g<-rbind(g, g5)
  rm(g3,g4,g5)
  
  
  
  
  return(g)
}

 
# f) Randomizations to create null model. 
#It just consist in running the analysis after having randomized the structure.
#There are two ways of randomize: by tree (see function tree_randomization) or distance randomization  (see function distances_rand).  
# By default I am using the distances_rand, but you can change the function in the first function line.
# rand is the number of runs of the randomization underwhich we average.


randomization<-function(gratio, rand, c)
{
  grandtree0<- distances_rand(gratio)
  
  dtotbinned<-av_q(grandtree0)
  rm(grandtree0)
  
  
  
  for(i in 1:rand)
  {
    print(i)
    grandtree<- distances_rand(gratio)
    
    
    
    dtotbinned2<-av_q(grandtree)
    
    dtotbinned<-rbind(dtotbinned, dtotbinned2)
    rm(dtotbinned2)
    
  }
  c<-1.96

  
  dtotbinned_rand_dist<- dtotbinned %>% filter(nbin>1)%>% group_by(b, observable)  %>% summarise(ld=sum(nbin/sum(nbin)*ld), 
                                                                                                 Corr_av=sum(nbin*Corr/sum(nbin)), error=c*sqrt(sum((Corr-Corr_av)*(Corr-Corr_av)*nbin)/(sum(nbin)-1))/(sqrt(sum(nbin))),
                                                                                                 nbin=sum(nbin) ) %>% ungroup()                                          
  
  dtotbinned_rand_dist<-  dtotbinned_rand_dist %>% rename(Corr=Corr_av) %>%
    mutate(type="Null Model") %>%
    mutate(nbin=as.integer(nbin/rand)) %>%
    select( b ,   ld  , nbin,  error ,type ,  observable,    Corr )
  
  
  
  
  
  
  
  return(dtotbinned_rand_dist)
}


