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
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array 
        highmig_pwfst_array[,,k,j] <-  pairwise.neifst(highmig_hierfstat[[j]][[k]])
        
        ##take the mean pwfst and save it    
        highmig_fst_df[j,1] <- signif(summary(highmig_pwfst_array[,,,j])[4], 2)
        ##save min pwfst
        highmig_fst_df[j,2] <- signif(summary(highmig_pwfst_array[,,,j])[1])
        ##save max pwfst
        highmig_fst_df[j,3] <- signif(summary(highmig_pwfst_array[,,,j])[6])
        
        ##reduce to 3 sigfigs
        highmig_fst_df <- signif(highmig_fst_df, digits = 3)
        
      } else {
        
        ##convert to hierfstat data frame and store in a list
        lowmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array 
        lowmig_pwfst_array[,,k,j] <-  pairwise.neifst(lowmig_hierfstat[[j]][[k]])
        
        ##take mean pwfst
        lowmig_fst_df[j,1] <- summary(lowmig_pwfst_array[,,,j])[4]
        ##save min pwfst
        lowmig_fst_df[j,2] <- summary(lowmig_pwfst_array[,,,j])[1]
        ##save max pwfst
        lowmig_fst_df[j,3] <- summary(lowmig_pwfst_array[,,,j])[6]
        
        ##reduce to 3 sigfigs
        lowmig_fst_df <- signif(lowmig_fst_df, digits = 3)
        
      }
    }
  }
  
}

##name data frames - high mig
rownames(highmig_fst_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                              "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(highmig_fst_df) <- c("Mean Fst", "Min Fst", "Max Fst")

##name data frames - low mig
rownames(lowmig_fst_df) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                              "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")

colnames(lowmig_fst_df) <- c("Mean Fst", "Min Fst", "Max Fst")



##write out files
write.csv(highmig_fst_df,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\highmig_fst_df.csv")
write.csv(lowmig_fst_df,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\lowmig_fst_df.csv")

