### Viral Production Assay: Step 2 ###

# Importing used functions, these can be found in viral_production_step2_source.R
source("viral_production_step2_source.R")

# Importing data from Step 1
data_test <- read.csv("NJ1.csv") # Change name of column => this is so that I can compare if same output is produced
data <- read.csv("NJ1.csv")
names(data)[names(data) == 'Expt_No'] <- 'Station_Number' # Changing column name to something more appropriate => you can also use index of column (5)

