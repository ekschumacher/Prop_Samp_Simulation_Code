########################################
########### Libraries ##################
########################################

library(adegenet)
library(diveRsity)
library(hierfstat)

######################################
########## Calculate PW Fst ##########
######################################
quac_dir <- "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\q_acerifolia"
quen_dir <- "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\q_engelmannii"
quog_dir <- "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\q_oglethorpensis"

species_names <- c("q_acerifolia", "q_engelmannii", "q_oglethorpensis")

##final df with combined species
species_pwfst_df <- matrix(nrow = 3, ncol = 3)

##########QUAC
setwd(quac_dir)
quac_list <- list.files(quac_dir,".gen$")

quac_genind_list <- list() 

quac_hierfstat <- list()

##quac pwfst array
quac_pwfst_array <- array(dim = c(4,4,100))

##min, max, mean of replicates 
quac_mean_max_min_fst <- matrix(nrow = 3, ncol = 100)

#loop going through every replicate for each scenario
for(k in 1:length(quac_list)) {
  #creating a temporary genind object
  quac_genind_list[[k]] <- read.genepop(quac_list[[k]], ncode=3)
  
  quac_hierfstat[[k]] <- genind2hierfstat(quac_genind_list[[k]])
  
  quac_pwfst_array[,,k] <- pairwise.neifst(quac_hierfstat[[k]])
  
  ##calculate statistics for QUAC
  quac_mean_max_min_fst[1,k] <- mean(quac_pwfst_array[,,k], na.rm = TRUE)
  quac_mean_max_min_fst[2,k] <- min(quac_pwfst_array[,,k], na.rm = TRUE)
  quac_mean_max_min_fst[3,k] <- max(quac_pwfst_array[,,k], na.rm = TRUE)
  
}
    
    ##write loops to calculate the mean of min/max/mean pwfst and then do it across scenario
  #  for(a in 1:length(scenarios)){
  #    for(b in 1:length(scenarios)){
  #      species_pwfst_df[a,b] <- round(mean(species_array_pwfst[b,,a]),3)
  #    }
  #  }
    
  
    }
    
}

##name rows and columns 
rownames(species_pwfst_df) <- species_names
colnames(species_pwfst_df) <- c("MeanFst", "MinFst", "MaxFst")

##write table out
write.csv(species_pwfst_df, "case_study_pwfst_df.csv")


