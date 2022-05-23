##Loading libraries
setwd("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R")

{
  library(tidyverse)
  
}
#local.min.max from SpatialEco package might be handy

NJ1<- read.csv("NJ1.csv")

plots_vp(NJ1)
