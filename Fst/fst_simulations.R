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

#looping over combinations, scenarios, and replicates
#saving results in 2D arrays
for(i in 1:length(combinations)) {
  for(j in 1:length(scenarios)) {
    setwd(paste(my_dir,combinations[i],scenarios[j],sep=""))
    list_files = list.files(path = paste(my_dir,combinations[i],scenarios[j],sep=""), pattern = ".gen$")
    for(k in 1:length(list_files)) {
      temp_genind = read.genepop(list_files[[k]], ncode=3)
      
      if(i == 1) {
        
        highmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
      
      } else {
      
       lowmig_hierfstat[[j]][[k]] <- genind2hierfstat(temp_genind)
        
        
      }
    }
  }
}
##########high migration code
##scenarios for calculating pwfst arrays - high mig
highmig_pwfst_scen1_array <- array(dim = c(5,5,100))
highmig_pwfst_scen2_array <- array(dim = c(5,5,100))
highmig_pwfst_scen3_array <- array(dim = c(5,5,100))
highmig_pwfst_scen4_array <- array(dim = c(5,5,100))
highmig_pwfst_scen5_array <- array(dim = c(5,5,100))
highmig_pwfst_scen6_array <- array(dim = c(5,5,100))
highmig_pwfst_scen7_array <- array(dim = c(5,5,100))
highmig_pwfst_scen8_array <- array(dim = c(5,5,100))
highmig_pwfst_scen9_array <- array(dim = c(5,5,100))

##final mean pwfst matrices
highmig_meanfst_scen1 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen2 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen3 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen4 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen5 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen6 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen7 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen8 <- matrix(nrow = 5, ncol = 5)
highmig_meanfst_scen9 <- matrix(nrow = 5, ncol = 5)

for(n in 1:length(highmig_hierfstat[[1]])){
    
  highmig_pwfst_scen1_array[,,n] <- pairwise.neifst(highmig_hierfstat[[1]][[n]])
  highmig_pwfst_scen2_array[,,n] <- pairwise.neifst(highmig_hierfstat[[2]][[n]])
  highmig_pwfst_scen3_array[,,n] <- pairwise.neifst(highmig_hierfstat[[3]][[n]])
  highmig_pwfst_scen4_array[,,n] <- pairwise.neifst(highmig_hierfstat[[4]][[n]])
  highmig_pwfst_scen5_array[,,n] <- pairwise.neifst(highmig_hierfstat[[5]][[n]])
  highmig_pwfst_scen6_array[,,n] <- pairwise.neifst(highmig_hierfstat[[6]][[n]])
  highmig_pwfst_scen7_array[,,n] <- pairwise.neifst(highmig_hierfstat[[7]][[n]])
  highmig_pwfst_scen8_array[,,n] <- pairwise.neifst(highmig_hierfstat[[8]][[n]])
  highmig_pwfst_scen9_array[,,n] <- pairwise.neifst(highmig_hierfstat[[9]][[n]])
  
  
  
}

##calculate mean pwfst tables 

for(o in 1:5){
  
  for(p in 1:5){
    
    highmig_meanfst_scen1[o,p] <- mean(highmig_pwfst_scen1_array[o,p,])
    highmig_meanfst_scen2[o,p] <- mean(highmig_pwfst_scen2_array[o,p,])
    highmig_meanfst_scen3[o,p] <- mean(highmig_pwfst_scen3_array[o,p,])
    highmig_meanfst_scen4[o,p] <- mean(highmig_pwfst_scen4_array[o,p,])
    highmig_meanfst_scen5[o,p] <- mean(highmig_pwfst_scen5_array[o,p,])
    highmig_meanfst_scen6[o,p] <- mean(highmig_pwfst_scen6_array[o,p,])
    highmig_meanfst_scen7[o,p] <- mean(highmig_pwfst_scen7_array[o,p,])
    highmig_meanfst_scen8[o,p] <- mean(highmig_pwfst_scen8_array[o,p,])
    highmig_meanfst_scen9[o,p] <- mean(highmig_pwfst_scen9_array[o,p,])
    
    
  }
  
  
}

##########low migration code
##scenarios for calculating pwfst arrays - high mig
lowmig_pwfst_scen1_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen2_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen3_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen4_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen5_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen6_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen7_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen8_array <- array(dim = c(5,5,100))
lowmig_pwfst_scen9_array <- array(dim = c(5,5,100))

##final mean pwfst matrices
lowmig_meanfst_scen1 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen2 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen3 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen4 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen5 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen6 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen7 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen8 <- matrix(nrow = 5, ncol = 5)
lowmig_meanfst_scen9 <- matrix(nrow = 5, ncol = 5)

for(n in 1:length(lowmig_hierfstat[[1]])){
  
  lowmig_pwfst_scen1_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[1]][[n]])
  lowmig_pwfst_scen2_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[2]][[n]])
  lowmig_pwfst_scen3_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[3]][[n]])
  lowmig_pwfst_scen4_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[4]][[n]])
  lowmig_pwfst_scen5_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[5]][[n]])
  lowmig_pwfst_scen6_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[6]][[n]])
  lowmig_pwfst_scen7_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[7]][[n]])
  lowmig_pwfst_scen8_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[8]][[n]])
  lowmig_pwfst_scen9_array[,,n] <- pairwise.neifst(lowmig_hierfstat[[9]][[n]])
  
  
  
}

##calculate mean pwfst tables 

for(o in 1:5){
  
  for(p in 1:5){
    
    lowmig_meanfst_scen1[o,p] <- mean(lowmig_pwfst_scen1_array[o,p,])
    lowmig_meanfst_scen2[o,p] <- mean(lowmig_pwfst_scen2_array[o,p,])
    lowmig_meanfst_scen3[o,p] <- mean(lowmig_pwfst_scen3_array[o,p,])
    lowmig_meanfst_scen4[o,p] <- mean(lowmig_pwfst_scen4_array[o,p,])
    lowmig_meanfst_scen5[o,p] <- mean(lowmig_pwfst_scen5_array[o,p,])
    lowmig_meanfst_scen6[o,p] <- mean(lowmig_pwfst_scen6_array[o,p,])
    lowmig_meanfst_scen7[o,p] <- mean(lowmig_pwfst_scen7_array[o,p,])
    lowmig_meanfst_scen8[o,p] <- mean(lowmig_pwfst_scen8_array[o,p,])
    lowmig_meanfst_scen9[o,p] <- mean(lowmig_pwfst_scen9_array[o,p,])
    
  }
  
  
}

##write out Fsts - high migration 
highmig_mean_max_min <- matrix(nrow = 9, ncol = 3)

highmig_mean_max_min[1,1] <- as.numeric(mean(highmig_meanfst_scen1, na.rm = TRUE))
highmig_mean_max_min[1,2] <- as.numeric(min(highmig_meanfst_scen1, na.rm = TRUE))
highmig_mean_max_min[1,3] <- as.numeric(max(highmig_meanfst_scen1, na.rm = TRUE))

highmig_mean_max_min[2,1] <- as.numeric(mean(highmig_meanfst_scen2, na.rm = TRUE))
highmig_mean_max_min[2,2] <- as.numeric(min(highmig_meanfst_scen2, na.rm = TRUE))
highmig_mean_max_min[2,3] <- as.numeric(max(highmig_meanfst_scen2, na.rm = TRUE))

highmig_mean_max_min[3,1] <- as.numeric(mean(highmig_meanfst_scen3, na.rm = TRUE))
highmig_mean_max_min[3,2] <- as.numeric(min(highmig_meanfst_scen3, na.rm = TRUE))
highmig_mean_max_min[3,3] <- as.numeric(max(highmig_meanfst_scen3, na.rm = TRUE))

highmig_mean_max_min[4,1] <- as.numeric(mean(highmig_meanfst_scen4, na.rm = TRUE))
highmig_mean_max_min[4,2] <- as.numeric(min(highmig_meanfst_scen4, na.rm = TRUE))
highmig_mean_max_min[4,3] <- as.numeric(max(highmig_meanfst_scen4, na.rm = TRUE))

highmig_mean_max_min[5,1] <- as.numeric(mean(highmig_meanfst_scen5, na.rm = TRUE))
highmig_mean_max_min[5,2] <- as.numeric(min(highmig_meanfst_scen5, na.rm = TRUE))
highmig_mean_max_min[5,3] <- as.numeric(max(highmig_meanfst_scen5, na.rm = TRUE))

highmig_mean_max_min[6,1] <- as.numeric(mean(highmig_meanfst_scen6, na.rm = TRUE))
highmig_mean_max_min[6,2] <- as.numeric(min(highmig_meanfst_scen6, na.rm = TRUE))
highmig_mean_max_min[6,3] <- as.numeric(max(highmig_meanfst_scen6, na.rm = TRUE))

highmig_mean_max_min[7,1] <- as.numeric(mean(highmig_meanfst_scen7, na.rm = TRUE))
highmig_mean_max_min[7,2] <- as.numeric(min(highmig_meanfst_scen7, na.rm = TRUE))
highmig_mean_max_min[7,3] <- as.numeric(max(highmig_meanfst_scen7, na.rm = TRUE))

highmig_mean_max_min[8,1] <- as.numeric(mean(highmig_meanfst_scen8, na.rm = TRUE))
highmig_mean_max_min[8,2] <- as.numeric(min(highmig_meanfst_scen8, na.rm = TRUE))
highmig_mean_max_min[8,3] <- as.numeric(max(highmig_meanfst_scen8, na.rm = TRUE))

highmig_mean_max_min[9,1] <- as.numeric(mean(highmig_meanfst_scen9, na.rm = TRUE))
highmig_mean_max_min[9,2] <- as.numeric(min(highmig_meanfst_scen9, na.rm = TRUE))
highmig_mean_max_min[9,3] <- as.numeric(max(highmig_meanfst_scen9, na.rm = TRUE))

##name rows and columns 
rownames(highmig_mean_max_min) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                    "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(highmig_mean_max_min) <- c("Mean_Fst", "Min_Fst", "Max_Fst")

##write out csv
write.csv(highmig_mean_max_min, "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\highmig_mean_max_min.csv")

##write out Fsts - high migration 
lowmig_mean_max_min <- matrix(nrow = 9, ncol = 3)

lowmig_mean_max_min[1,1] <- as.numeric(mean(lowmig_meanfst_scen1, na.rm = TRUE))
lowmig_mean_max_min[1,2] <- as.numeric(min(lowmig_meanfst_scen1, na.rm = TRUE))
lowmig_mean_max_min[1,3] <- as.numeric(max(lowmig_meanfst_scen1, na.rm = TRUE))

lowmig_mean_max_min[2,1] <- as.numeric(mean(lowmig_meanfst_scen2, na.rm = TRUE))
lowmig_mean_max_min[2,2] <- as.numeric(min(lowmig_meanfst_scen2, na.rm = TRUE))
lowmig_mean_max_min[2,3] <- as.numeric(max(lowmig_meanfst_scen2, na.rm = TRUE))

lowmig_mean_max_min[3,1] <- as.numeric(mean(lowmig_meanfst_scen3, na.rm = TRUE))
lowmig_mean_max_min[3,2] <- as.numeric(min(lowmig_meanfst_scen3, na.rm = TRUE))
lowmig_mean_max_min[3,3] <- as.numeric(max(lowmig_meanfst_scen3, na.rm = TRUE))

lowmig_mean_max_min[4,1] <- as.numeric(mean(lowmig_meanfst_scen4, na.rm = TRUE))
lowmig_mean_max_min[4,2] <- as.numeric(min(lowmig_meanfst_scen4, na.rm = TRUE))
lowmig_mean_max_min[4,3] <- as.numeric(max(lowmig_meanfst_scen4, na.rm = TRUE))

lowmig_mean_max_min[5,1] <- as.numeric(mean(lowmig_meanfst_scen5, na.rm = TRUE))
lowmig_mean_max_min[5,2] <- as.numeric(min(lowmig_meanfst_scen5, na.rm = TRUE))
lowmig_mean_max_min[5,3] <- as.numeric(max(lowmig_meanfst_scen5, na.rm = TRUE))

lowmig_mean_max_min[6,1] <- as.numeric(mean(lowmig_meanfst_scen6, na.rm = TRUE))
lowmig_mean_max_min[6,2] <- as.numeric(min(lowmig_meanfst_scen6, na.rm = TRUE))
lowmig_mean_max_min[6,3] <- as.numeric(max(lowmig_meanfst_scen6, na.rm = TRUE))

lowmig_mean_max_min[7,1] <- as.numeric(mean(lowmig_meanfst_scen7, na.rm = TRUE))
lowmig_mean_max_min[7,2] <- as.numeric(min(lowmig_meanfst_scen7, na.rm = TRUE))
lowmig_mean_max_min[7,3] <- as.numeric(max(lowmig_meanfst_scen7, na.rm = TRUE))

lowmig_mean_max_min[8,1] <- as.numeric(mean(lowmig_meanfst_scen8, na.rm = TRUE))
lowmig_mean_max_min[8,2] <- as.numeric(min(lowmig_meanfst_scen8, na.rm = TRUE))
lowmig_mean_max_min[8,3] <- as.numeric(max(lowmig_meanfst_scen8, na.rm = TRUE))

lowmig_mean_max_min[9,1] <- as.numeric(mean(lowmig_meanfst_scen9, na.rm = TRUE))
lowmig_mean_max_min[9,2] <- as.numeric(min(lowmig_meanfst_scen9, na.rm = TRUE))
lowmig_mean_max_min[9,3] <- as.numeric(max(lowmig_meanfst_scen9, na.rm = TRUE))

##name rows and columns 
rownames(lowmig_mean_max_min) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5",
                                    "Scenario 6", "scenario 7", "Scenario 8", "Scenario 9")
colnames(lowmig_mean_max_min) <- c("Mean_Fst", "Min_Fst", "Max_Fst")

##write out csv
write.csv(lowmig_mean_max_min, "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\lowmig_mean_max_min.csv")
