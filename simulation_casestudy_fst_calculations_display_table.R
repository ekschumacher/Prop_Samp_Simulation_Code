########################################
########### Libraries ##################
########################################

library(adegenet)
library(diveRsity)
library(hierfstat)

########################################
############# Load Documents ###########
########################################
quac_path <- "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\q_acerifolia"
setwd("G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\q_acerifolia")

##import files 
quac_list <- list.files(path = quac_path, pattern = ".gen$")

##create list
quac_genind_list <- list()

##run loop
for(i in 1:length(quac_list)) {
  quac_genind_list[[i]] <- read.genepop(quac_list[[i]], ncode=3)
}

######################################
########## Calculate PW Fst ##########
######################################
##create array 
quac_array <- array(dim = c(4,4,100))

##run loop to calculate pwfst
for(i in 1:length(quac_genind_list)) {

    quac_array[,,i] <- as.matrix(pairwise.fst(quac_genind_list[[i]]))
    
}

######################################################
############# Average Across Replicates ##############
######################################################
##mean pwfst data frame
quac_mean_pwfst_df <- matrix(nrow = 4, ncol = 4)

##calulate mean Fsts for all 100 replicates
for(i in 1:4) {
  for(j in 1:4) {
    
    quac_mean_pwfst_df[i,j] <- mean(quac_array[i,j,])
    
  }
}

######################################################
################## Write Out Files ###################
######################################################
##write out matrix 
write.csv(quac_mean_pwfst_df, "quac_mean_pwfst_df.csv")

##create boxplot
pdf("G:\\Shared drives\\Emily_Schumacher\\Morton-REU-master\\case_studies\\quac_mean_pwfst.pdf", width = 8, height = 6)
boxplot(quac_mean_pwfst_df, ylim = c(0,0.35),
        main = "Mean PWFst for Q. acerifolia simulations")
dev.off()

###################################################################
################## Load in oglethorpensis files ###################
###################################################################
quog_case_studies <- "G:\\Shared drives\\Emily_Schumacher\\Morton-REU-master\\case_studies\\Simulations\\q_oglethorpensis"

setwd(quog_case_studies)

import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}

#converting all simulation files from arlequin format to genepop format
import_arp2gen_files(quog_case_studies, ".arp$")

##import files 
quog_list <- list.files(quog_case_studies, ".gen$")

##create list
quog_genind_list <- list()

##run loop
for(i in 1:length(quog_list)) {
  quog_genind_list[[i]] <- read.genepop(quog_list[[i]], ncode=3)
}

######################################
########## Calculate PW Fst ##########
######################################
##create array 
quog_array <- array(dim = c(length(levels(quog_genind_list[[1]]@pop)),
                            length(levels(quog_genind_list[[1]]@pop)),
                            100))

##run loop to calculate pwfst
for(i in 1:length(quog_genind_list)) {
  
  quog_array[,,i] <- as.matrix(pairwise.fst(quog_genind_list[[i]]))
  
}



########################################################
############# Load englemanii Documents ################
########################################################

quen_case_studies <- "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\Simulations\\q_engelmannii\\"

setwd(quen_case_studies)

##import files 
quen_list <- list.files(pattern = ".gen$")

##create list
quen_genind_list <- list()

##run loop
for(i in 1:length(quen_list)) {
  quen_genind_list[[i]] <- read.genepop(quen_list[[i]], ncode=3)
}

######################################
########## Calculate PW Fst ##########
######################################
##create array 
quen_array <- array(dim = c(length(levels(quen_genind_list[[1]]@pop)),length(levels(quen_genind_list[[1]]@pop)),
                            length(quen_genind_list)))

##run loop to calculate pwfst
for(i in 1:length(quen_genind_list)) {
  
  quen_array[,,i] <- as.matrix(pairwise.fst(quen_genind_list[[i]]))
  
}

######################################################
############# Average Across Replicates ##############
######################################################
##mean pwfst data frame
quen_mean_pwfst_df <- matrix(nrow = length(levels(quen_genind_list[[1]]@pop)), 
                             ncol = length(levels(quen_genind_list[[1]]@pop)))

##calulate mean Fsts for all 100 replicates
for(i in 1:length(levels(quen_genind_list[[1]]@pop))) {
  for(j in 1:length(levels(quen_genind_list[[1]]@pop))) {
    
    quen_mean_pwfst_df[i,j] <- mean(quen_array[i,j,])
    
  }
}

rownames(quen_mean_pwfst_df) <- c("Pop1", "Pop2", "Pop3","Pop4")
colnames(quen_mean_pwfst_df) <- c("Pop1", "Pop2", "Pop3","Pop4")

######################################################
################## Write Out Files ###################
######################################################
##write out matrix 
setwd("G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies")
write.csv(quen_mean_pwfst_df, "quen_mean_pwfst_df.csv")

##create boxplot
pdf("G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\quen_mean_pwfst.pdf", width = 8, height = 6)
boxplot(quen_mean_pwfst_df, ylim = c(0,0.15),
        main = "Mean PWFst for Q. englemannii simulations")
dev.off()

#####################################################################
################# Create Data Table for Pub Display #################
#####################################################################
case_studies <- "G:\\Shared drives\\Emily_Schumacher\\simulation_code\\case_studies\\"

##load in df for each species 
setwd(case_studies)

pwfst_df_list <- list.files(pattern = ".csv$")

pwfst_list <- list()

for(d in 1:length(pwfst_df_list)){
  
  pwfst_list[[d]] <- read.csv(pwfst_df_list[[d]])
  
  pwfst_list[[d]] <- pwfst_list[[d]][,-1]
  
  pwfst_list[[d]][pwfst_list[[d]] == 0] <- NA
  
}

##begin input into a dataframe 
#create df
case_study_df <- matrix(nrow = 3, ncol = 3)

##set up dataframe, each row is a case study, each column is a metric 
rownames(case_study_df) <- c("QUAC","QUEN","QUOG")
colnames(case_study_df) <- c("Mean_PWFst", "Min_PWFst", "Max_PWFst")


for(c in 1:length(pwfst_list)){
  
  case_study_df[c,2:3] <- c(min(pwfst_list[[c]], na.rm = TRUE), max(pwfst_list[[c]], na.rm = TRUE))
  
  case_study_df[c, 1] <- mean(colMeans(pwfst_list[[c]], na.rm = TRUE))
  
}

##write table out
write.csv(case_study_df, "case_study_pwfst_df.csv")


