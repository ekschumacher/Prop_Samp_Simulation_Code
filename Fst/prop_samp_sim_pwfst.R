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
my_dir = "C:\\Users\\eschumacher\\Documents\\kaylee_code_1_19_21\\samp_pop_sims\\Simulations"
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

#****************************************************************************************************************************************************
##lists to store hierfstat data frames - high mig
highmig_hierfstat <- list(list(), list(), list(), list(),  list(),
                          list(), list(), list(), list())

##lists to store hierfstat data frames - low mig
lowmig_hierfstat <- list(list(), list(), list(), list(),  list(),
                         list(), list(), list(), list())

##4D arrays to save pwfst for each scenario and replicate - high and low migration combinations 
highmig_pwfst_array <- array(dim = c(5,5,100,9))
lowmig_pwfst_array <- array(dim = c(5,5,100,9))

##3D arrays to store mean, max, and minimum pwfst for each replicate and scenario in high and low migration combinations
highmig_fst_min_mean_max <- array(dim = c(3,100,9))
lowmig_fst_min_mean_max <- array(dim = c(3,100,9))

##2D arrays to store averaged results for mean, min, and max Fst in high and low migration combinations
highmig_fst_df <- matrix(nrow = 9, ncol = 3)
lowmig_fst_df <- matrix(nrow = 9, ncol = 3)

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    setwd(paste0(my_dir,combinations[i],"_highSamp",scenarios[j]))
    list_files = list.files(path = paste0(my_dir,combinations[i],"_highSamp",scenarios[j]), pattern = ".gen$")
    for(k in 1:length(list_files)) {
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      if(i == 1) {
        
        ##convert to hierfstat data frame and store in list
        highmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array - high mig
        #dims 1 and 2: Nei's pwfst values between all 5 populations 
        #dim 3: Replicates (100)
        #dim 4: Scenarios 1 - 9 (very different population sizes --> equal population sizes)
        highmig_pwfst_array[,,k,j] <-  pairwise.neifst(highmig_hierfstat[[j]][[k]])
        
        ##calculate max, min, and mean pwfst for every replicate - high mig
        #dim 1: Mean, min, max pwfst
        #dim 2: Replicates
        #dim 3: Scenarios 1 - 9 (very different population sizes --> equal population sizes)
        highmig_fst_min_mean_max[1,k,j] <- mean(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        highmig_fst_min_mean_max[2,k,j] <- min(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        highmig_fst_min_mean_max[3,k,j] <- max(highmig_pwfst_array[,,k,j], na.rm = TRUE)
        
        #Calculate the mean of min/max/mean pwfst and then do it across scenario - high mig
        ##Result df output: 
        #dim 1: Scenarios 1 - 9
        #dim 2: Pwfst 
        for(a in 1:length(scenarios)){
          for(b in 1:length(highmig_fst_min_mean_max[,1,1])){
            highmig_fst_df[a,b] <- round(mean(highmig_fst_min_mean_max[b,,a]),3)
          }
        }
        
      } else {
        
        ##convert to hierfstat data frame and store in a list
        lowmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        ##pairwise fst results for all scenarios and replicates stored in a 4D array - low mig 
        #dims 1 and 2: Nei's pwfst values between all 5 populations 
        #dim 3: Replicates (100)
        #dim 4: Scenarios 1 - 9 (very different population sizes --> equal population sizes)
        lowmig_pwfst_array[,,k,j] <-  pairwise.neifst(lowmig_hierfstat[[j]][[k]])
        
        ##calculate max, min, and mean pwfst for every replicate - low mig
        #dim 1: Mean, min, max pwfst
        #dim 2: Replicates
        #dim 3: Scenarios 1 - 9 (very different population sizes --> equal population sizes)
        lowmig_fst_min_mean_max[1,k,j] <- mean(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        lowmig_fst_min_mean_max[2,k,j] <- min(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        lowmig_fst_min_mean_max[3,k,j] <- max(lowmig_pwfst_array[,,k,j], na.rm = TRUE)
        
        ##write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
        ##Result df output: 
        #dim 1: Scenarios 1 - 9
        #dim 2: Pwfst 
        for(a in 1:length(scenarios)){
          for(b in 1:length(lowmig_fst_min_mean_max[,1,1])){
            lowmig_fst_df[a,b] <- round(mean(lowmig_fst_min_mean_max[b,,a]),3)
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
write.csv(highmig_fst_df,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\highmig_fst_df.csv")
write.csv(lowmig_fst_df,"G:\\Shared drives\\Emily_Schumacher\\simulation_code\\lowmig_fst_df.csv")

