### Viral Production Assay: Step 2 ###

# Importing used functions, these can be found in viral_production_step2_source.R & vp_calc_functions.R
source("viral_production_step2_source.R")
source("vp_calc_functions.R")

# Importing data from Step 1
data <- read.csv("NJ1.csv")
names(data)[names(data) == 'Expt_No'] <- 'Station_Number' # Changing column name to something more appropriate => you can also use index of column (5)

