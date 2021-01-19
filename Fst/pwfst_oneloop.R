########################################
########### Libraries ##################
########################################

library(adegenet)
library(diveRsity)
library(hierfstat)

########################################
############# Load Documents ###########
########################################

#root directory
#containing sub-folders
my_dir = "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\Attempt5_allcombinations_100reps\\Simulations"
setwd(my_dir)

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

#loop converting .arp to .gen for all combinations and scenarios
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    import_arp2gen_files(paste(my_dir,combinations[i],scenarios[j],sep=""), ".arp$")
  }
}

#****************************************************************************************************************************************************
highmig_hierfstat <- list(list(), list(), list(), list(),  list(),
                          list(), list(), list(), list())

lowmig_hierfstat <- list(list(), list(), list(), list(),  list(),
                         list(), list(), list(), list())

##arrays low and high mig 
highmig_pwfst_array <- array(dim = c(5,5,100,9))
lowmig_pwfst_array <- array(dim = c(5,5,100,9))

##results dfs
highmig_fst_df <- matrix(nrow = 9, ncol = 3)
lowmig_fst_df <- matrix(nrow = 9, ncol = 3)

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    for(k in 1:length(list_files)) {
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      if(i == 1) {
        
        ##convert to hierfstat data frame and store in list
        highmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array - 
        highmig_pwfst_array[,,k,j] <-  pairwise.neifst(highmig_hierfstat[[j]][[k]])
        
        ##calculate max, min, and mean pwfst for every replicate
        highmig_fst_min_mean_max[1,k,j] <- mean(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        highmig_fst_min_mean_max[2,k,j] <- min(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        highmig_fst_min_mean_max[3,k,j] <- max(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        
        ##write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
        for(a in 1:length(scenarios)){
          for(b in 1:length(highmig_fst_min_mean_max[,1,1])){
            highmig_pwfst_output[a,b] <- round(mean(highmig_fst_min_mean_max[b,,a]),3)
          }
        }
        
      } else {
        
        ##convert to hierfstat data frame and store in a list
        lowmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array 
        lowmig_pwfst_array[,,k,j] <-  pairwise.neifst(lowmig_hierfstat[[j]][[k]])
        
        ##calculate max, min, and mean pwfst for every replicate
        lowmig_fst_min_mean_max[1,k,j] <- mean(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        lowmig_fst_min_mean_max[2,k,j] <- min(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        lowmig_fst_min_mean_max[3,k,j] <- max(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        
        ##write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
        for(a in 1:length(scenarios)){
          for(b in 1:length(lowmig_fst_min_mean_max[,1,1])){
            lowmig_pwfst_output[a,b] <- round(mean(lowmig_fst_min_mean_max[b,,a]),3)
          }
        }
      }
    }
  }
  
}

##name data frames - high mig
rownames(highmig_pwfst_output) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                              "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_pwfst_output) <- c("Mean Fst", "Min Fst", "Max Fst")

##name data frames - low mig
rownames(lowmig_pwfst_output) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                              "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_pwfst_output) <- c("Mean Fst", "Min Fst", "Max Fst")



##write out files
write.csv(highmig_pwfst_output,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\highmig_fst_df.csv")
write.csv(lowmig_pwfst_output,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\lowmig_fst_df.csv")

