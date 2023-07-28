#####1.0 Installing Packages####
source("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/1_1_simulations_source.R")
library(tidyverse)
library(ggstatsplot)
library(ggsci)
library(parallel)
library(doParallel)
library(Rmpi)
library(doMPI)
library(doSNOW)
library(foreach)


length_df<- unique(simu_df$Expt_No)

#Define variables
var<- c('Expt_No', 'Sample_Type')
var2<- c('Expt_No', 'Sample_Type', 'Timepoint')
input_df = simu_df
output_df = output_df
output_csv = T
output_dir = './simulations'
output_csv_name = 'simu_output_df.csv'

{
  output_df<- data.frame(matrix(ncol = 6, nrow =0))
  colnames(output_df)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
  output_df <- output_df %>%
    mutate_at(c('Sample_Type', 
                'VP_Type'), as.character) %>%
    mutate_at(c('Expt_No', 'VP', 'VP_SE', 
                'VP_R_Squared'), as.numeric)
  
}

if (file.exists('./simulations/logs')){
  print('Logs folder already exist') 
} else{
dir.create('./simulations/logs')
}


#parallel run. 
cl<- parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

foreach::foreach (mtd = 1:12,
                  .combine = rbind) %dopar% {
  
  print(paste("Started analysis using method ", mtd))
  #sink(paste0('./simulations/logs/simu_output_log_method', mtd,'.txt'))
  try(simu_calc_funct_list[[mtd]](simu_df))
  #sink()
}



parallel::stopCluster(cl)



write.csv(output_df, paste0(output_dir, '/', output_csv_name),
          row.names = F)
