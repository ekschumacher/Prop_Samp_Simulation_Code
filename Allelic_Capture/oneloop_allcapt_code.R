####################################
########## Load Libraries ##########
####################################

library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)

###################################
###### Load Directory #############
###################################

#containing sub-folders
my_dir = "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\Attempt5_allcombinations_100reps\\Simulations"
setwd(my_dir)
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

###function
source("G:/Shared drives/ten_oaks_gen/Fa_sample_funcs.R")
##functions
colMax <- function(data) sapply(data, max, na.rm = TRUE)
sample.pop<-function(genind_obj,vect_pop_ID,vect_samp_sizes){
  p<-length(vect_pop_ID)
  if (p>1) {
    for (p in 1:length(vect_pop_ID))
      alleles[p,]<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)
    alleles<-colSums(alleles)
  } else {alleles<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)}
  
  alleles
}    

#list of combinations
#combination sub-folder directories
combinations = c("\\highMig", "\\lowMig")

#list of scenarios
#simulation sub-folder directories
scenarios = c("\\scen1",
              "\\scen2",
              "\\scen3",
              "\\scen4",
              "\\scen5",
              "\\scen6",
              "\\scen7",
              "\\scen8",
              "\\scen9")

#####################################################
############## Run Sampling Code ####################
#####################################################
#creating list of vectors representing which individuals to sample from genind object
#sampling equally from each population - equal strategy
#5 pops in each scenario
rows_to_samp_equal = list(length = length(scenarios))
rows_to_samp_equal[[1]] = c(sample(1:30,15), sample(31:130,15), sample(131:230,15), sample(231:330,15), sample(331:1500,15))
rows_to_samp_equal[[2]] = c(sample(1:40,15), sample(41:190,15), sample(191:340,15), sample(341:490,15), sample(491:1500,15))
rows_to_samp_equal[[3]] = c(sample(1:50,15), sample(51:250,15), sample(251:450,15), sample(451:650,15), sample(651:1500,15))
rows_to_samp_equal[[4]] = c(sample(1:100,15), sample(101:300,15), sample(301:500,15), sample(501:700,15), sample(701:1500,15))
rows_to_samp_equal[[5]] = c(sample(1:150,15), sample(151:350,15), sample(351:550,15), sample(551:750,15), sample(750:1500,15))
rows_to_samp_equal[[6]] = c(sample(1:200,15), sample(201:450,15), sample(451:700,15), sample(701:950,15), sample(951:1500,15))
rows_to_samp_equal[[7]] = c(sample(1:200,15), sample(201:500,15), sample(501:800,15), sample(801:1100,15), sample(1101:1500,15))
rows_to_samp_equal[[8]] = c(sample(1:290,15), sample(291:590,15), sample(591:890,15), sample(891:1190,15), sample(1191:1500,15))
rows_to_samp_equal[[9]] = c(sample(1:300,15), sample(301:600,15), sample(601:900,15), sample(901:1200,15), sample(1201:1500,15))

#sampling 5% from each population
#5 pops in each scenario
rows_to_samp_prop = list(length = length(scenarios))
rows_to_samp_prop[[1]] = c(sample(1:30,2), sample(31:130,5), sample(131:230,5), sample(231:330,5), sample(331:1500,59))
rows_to_samp_prop[[2]] = c(sample(1:40,2), sample(41:190,8), sample(191:340,8), sample(341:490,8), sample(491:1500,51))
rows_to_samp_prop[[3]] = c(sample(1:50,3), sample(51:250,10), sample(251:450,10), sample(451:650,10), sample(651:1500,43))
rows_to_samp_prop[[4]] = c(sample(1:100,5), sample(101:300,10), sample(301:500,10), sample(501:700,10), sample(701:1500,40))
rows_to_samp_prop[[5]] = c(sample(1:150,8), sample(151:350,10), sample(351:550,10), sample(551:750,10), sample(750:1500,38))
rows_to_samp_prop[[6]] = c(sample(1:200,10), sample(201:450,13), sample(451:700,13), sample(701:950,13), sample(951:1500,28))
rows_to_samp_prop[[7]] = c(sample(1:200,10), sample(201:500,15), sample(501:800,15), sample(801:1100,15), sample(1101:1500,20))
rows_to_samp_prop[[8]] = c(sample(1:290,15), sample(291:590,15), sample(591:890,15), sample(891:1190,15), sample(1191:1500,15))
rows_to_samp_prop[[9]] = c(sample(1:300,15), sample(301:600,15), sample(601:900,15), sample(901:1200,15), sample(1201:1500,15))

##matrices
highmig_all_existing_by_sp_reps <- array(dim = c(100,9,9))
lowmig_all_existing_by_sp_reps <- array(dim = c(100,9,9))

##intermediate
highmig_equal_alleles <- array(dim = c(100,9,9))
highmig_prop_alleles <- array(dim = c(100,9,9))

highmig_equal_per <- array(dim = c(100,9,9))
highmig_prop_per <- array(dim = c(100,9,9))

##outputs
highmig_alleles_existing_by_cat <- matrix(nrow = 9, ncol = 9)
highmig_equal_mean_all_cap <- matrix(nrow = 9, ncol = 9)
highmig_prop_mean_all_cap <- matrix(nrow = 9, ncol = 9)
highmig_equal_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)
highmig_prop_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)

##intermediate
lowmig_equal_alleles <- array(dim = c(100,9,9))
lowmig_prop_alleles <- array(dim = c(100,9,9))

lowmig_equal_per <- array(dim = c(100,9,9))
lowmig_prop_per <- array(dim = c(100,9,9))

##outputs
lowmig_alleles_existing_by_cat <- matrix(nrow = 9, ncol = 9)
lowmig_equal_mean_all_cap <- matrix(nrow = 9, ncol = 9)
lowmig_prop_mean_all_cap <- matrix(nrow = 9, ncol = 9)
lowmig_equal_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)
lowmig_prop_mean_all_cap_per <- matrix(nrow = 9, ncol = 9)

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  #loop going through each scenario for each combination
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    #loop going through every replicate for each scenario
    for(k in 1:length(list_files)) {
      #creating a temporary genind object
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      #calculating the number of individuals per population 
      n_ind <- table(temp_genind@pop)
      
      ##create genepop file
      Spp_tot_genpop <- genind2genpop(temp_genind)
      
      ##separate by population 
      Spp_tot_genind_sep <- seppop(temp_genind)
      
      ##determine number of alleles captured by the equal sampling 
      alleles_cap_equal <- colSums(temp_genind@tab[rows_to_samp_equal[[j]],], na.rm = T)
      ##determine number of alleles captured by the equal sampling 
      alleles_cap_prop <- colSums(temp_genind@tab[rows_to_samp_prop[[j]],], na.rm = T)
      ##determine the total number of alleles in every category 
      allele_cat_tot <- get.allele.cat(Spp_tot_genpop, c(1:5), 2, n_ind)
      
      #if high migration 
      if(i == 1) {
        #saving number of alleles total in high migration scenario
        for (a in 1:length(allele_cat_tot)) highmig_all_existing_by_sp_reps[k,a,j] <- sum((allele_cat_tot[[a]])>0,na.rm=T)
        
        ##all the high migration allele # per category
        for (b in 1:length(allele_cat_tot)) highmig_equal_alleles[k,b,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[b]]]>0))
        for (c in 1:length(allele_cat_tot)) highmig_prop_alleles[k,c,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[c]]]>0))
       
         ##high migration percent alleles captured per category
        for (d in 1:length(allele_cat_tot)) highmig_equal_per[k,d,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[d]]]>0)/length(allele_cat_tot[[d]]),4)
        for (e in 1:length(allele_cat_tot)) highmig_prop_per[k,e,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[e]]]>0)/length(allele_cat_tot[[e]]),4)
      
  } #else { #if low migration
    #saving number of alleles total in low migration scenario
    for (a in 1:length(allele_cat_tot)) lowmig_all_existing_by_sp_reps[k,a,j] <- sum((allele_cat_tot[[a]])>0,na.rm=T)
    
    ##saving # of alleles captured per category   
    for (b in 1:length(allele_cat_tot)) lowmig_equal_alleles[k,b,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[b]]]>0))
    for (c in 1:length(allele_cat_tot)) lowmig_prop_alleles[k,c,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[c]]]>0))
    ##saving percent  of alleles per category
    for (d in 1:length(allele_cat_tot)) lowmig_equal_per[k,d,j] <- round(sum(alleles_cap_equal[allele_cat_tot[[d]]]>0)/length(allele_cat_tot[[d]]),4)
    for (e in 1:length(allele_cat_tot)) lowmig_prop_per[k,e,j] <- round(sum(alleles_cap_prop[allele_cat_tot[[e]]]>0)/length(allele_cat_tot[[e]]),4)
    
      }
    }
  }
#}

##existing alleles by species 
for(j in 1:length(highmig_alleles_existing_by_cat[1,])) {
  for (f in 1:length(allele_cat_tot)) {
    
    highmig_alleles_existing_by_cat[f,j] <-  signif(mean(highmig_all_existing_by_sp_reps[,j,f]), 3)
    highmig_equal_mean_all_cap[f,j] <- signif(mean(highmig_equal_alleles[,j,f]), 3)
    highmig_prop_mean_all_cap[f,j] <- signif(mean(highmig_prop_alleles[,j,f]),3)
    
    ##calculate percent of allele capture
    highmig_equal_mean_all_cap_per[f,j] <- signif(mean(highmig_equal_per[,j,f]),3)*100
    highmig_prop_mean_all_cap_per[f,j] <- signif(mean(highmig_prop_per[,j,f]),3)*100
    
  }
}

##name rows and columns 
rownames(highmig_alleles_existing_by_cat) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                               "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(highmig_alleles_existing_by_cat) <- list_allele_cat

write.csv(highmig_alleles_existing_by_cat, "highmig_alleles_existing_by_cat.csv")

##create a new df = equal sampling
highmig_all_cap_equal_df <- matrix(nrow = length(highmig_equal_mean_all_cap_per[,1]), 
                                   ncol = length(highmig_equal_mean_all_cap_per[,1]))

for(m in 1:length(highmig_equal_mean_all_cap_per[,1])){
  for(n in 1:length(highmig_equal_mean_all_cap_per[,1])){
    highmig_all_cap_equal_df[m,n] <- paste0(highmig_equal_mean_all_cap_per[m,n], "%", " ", "(", highmig_equal_mean_all_cap[m,n], ")")
  }
}  

rownames(highmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                        "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_all_cap_equal_df) <- list_allele_cat

##make new table using proportional sampling 
highmig_all_cap_prop_df <- matrix(nrow = length(highmig_prop_mean_all_cap_per[,1]), 
                                  ncol = length(highmig_prop_mean_all_cap_per[,1]))


for(o in 1:length(highmig_prop_mean_all_cap_per[,1])){
  for(p in 1:length(highmig_prop_mean_all_cap_per[,1])){
    highmig_all_cap_prop_df[o,p] <- paste0(highmig_prop_mean_all_cap_per[o,p], "%", " ", "(", highmig_prop_mean_all_cap[o,p], ")")
  }
}  

rownames(highmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                        "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_all_cap_prop_df) <- list_allele_cat

##write out data frames
write.csv(highmig_all_cap_equal_df, "highmig_all_cap_equal_df.csv")
write.csv(highmig_all_cap_prop_df, "highmig_all_cap_prop_df.csv")

################################
###### Low migration tables ####
################################

##existing alleles by species 
for(j in 1:length(lowmig_alleles_existing_by_cat[1,])) {
  for (f in 1:length(allele_cat_tot)) {
    
    lowmig_alleles_existing_by_cat[f,j] <-  signif(mean(lowmig_all_existing_by_sp_reps[,j,f]), 3)
    lowmig_equal_mean_all_cap[f,j] <- signif(mean(lowmig_equal_alleles[,j,f]), 3)
    lowmig_prop_mean_all_cap[f,j] <- signif(mean(lowmig_prop_alleles[,j,f]), 3)
    
    ##calculate percent of allele capture
    lowmig_equal_mean_all_cap_per[f,j] <- signif(mean(lowmig_equal_per[,j,f]), 3)*100
    lowmig_prop_mean_all_cap_per[f,j] <- signif(mean(lowmig_prop_per[,j,f]), 3)*100 
    
  }
}

##name rows and columns 
rownames(lowmig_alleles_existing_by_cat) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                              "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(lowmig_alleles_existing_by_cat) <- list_allele_cat

##write out files
write.csv(lowmig_alleles_existing_by_cat, "lowmig_alleles_existing_by_cat.csv")
write.csv(lowmig_equal_mean_all_cap, "lowmig_equal_mean_all_cap.csv")
write.csv(lowmig_prop_mean_all_cap, "lowmig_prop_mean_all_cap.csv")
write.csv(lowmig_equal_mean_all_cap_per, "lowmig_equal_mean_all_cap_per.csv")
write.csv(lowmig_prop_mean_all_cap_per, "lowmig_prop_mean_all_cap_per.csv")

##create a new df = equal sampling
lowmig_all_cap_equal_df <- matrix(nrow = length(lowmig_equal_mean_all_cap_per[,1]), 
                                  ncol = length(lowmig_equal_mean_all_cap_per[,1]))

for(m in 1:length(lowmig_equal_mean_all_cap_per[,1])){
  for(n in 1:length(lowmig_equal_mean_all_cap_per[,1])){
    lowmig_all_cap_equal_df[m,n] <- paste0(lowmig_equal_mean_all_cap_per[m,n], "%", " ", "(", lowmig_equal_mean_all_cap[m,n], ")")
  }
}  

rownames(lowmig_all_cap_equal_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                       "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_all_cap_equal_df) <- list_allele_cat

##make new table using proportional sampling 
lowmig_all_cap_prop_df <- matrix(nrow = length(lowmig_prop_mean_all_cap_per[,1]), 
                                 ncol = length(lowmig_prop_mean_all_cap_per[,1]))


for(o in 1:length(lowmig_prop_mean_all_cap_per[,1])){
  for(p in 1:length(lowmig_prop_mean_all_cap_per[,1])){
    lowmig_all_cap_prop_df[o,p] <- paste0(lowmig_prop_mean_all_cap_per[o,p], "%", " ", "(", lowmig_prop_mean_all_cap[o,p], ")")
  }
}  

rownames(lowmig_all_cap_prop_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                      "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_all_cap_prop_df) <- list_allele_cat

##write out data frames
write.csv(lowmig_all_cap_equal_df, "lowmig_all_cap_equal_df.csv")
write.csv(lowmig_all_cap_prop_df, "lowmig_all_cap_prop_df.csv")