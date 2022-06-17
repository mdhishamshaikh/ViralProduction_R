#A script to install and load required packages, and to source global functions

####1. Installing and loading libraries####

#Installing BiocManager
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

#List of packages to install
packages_to_load<- c("tidyverse", 
                     "flowWorkspace")

#Checking if packages are already present. If absent, then installing packages from BiocManager 
for (pack in packages_to_load){
  if(!requireNamespace(pack))
    BiocManager::install(pack, force = T)
  }

#Loading libraries  
for (pack in packages_to_load){
   library(pack, character.only = T)
}

