source("0_vp_source.R")
source("vp_functions.R")
library("tidyverse")


####2.1 VIPCAL - with separate replicate.
#We will average all the VP and VPC replicates, followed by suvtracting them to get Diff


vp_sr<- df_sr_tp(data)

VPCL_SR<- vipcal_sr(vp_sr)



####2.2 AVERAGE replicates VP and VPC. Followed by subtraction for Diff (SE)

vp_avg<- df_avg_tp(data)
vp_avg<- vp_avg[vp_avg$Sample_Type != 'Diff',]

VPCL_AVG<- vipcal_avg_sd(vp_avg)
#none of the VPC had any production


####2.3 Calculate the differenc curve first and the calc avg vipcal (SE)
vp_avg_diff<- df_avg_tp(data)
VPCL_AVG_Diff<- vipcal_avg_sd(vp_avg_diff)
#diff has positive values


####2.4 aVERAGED REPLICATE . NO STANDARD ERROR. 

###2.4.1 FIRST VP, VPC AND THEN SUBTREACT FOR DIFF CURVE
vp_avg_NO_SE<- df_avg_tp(data)
vp_avg_NO_SE<- vp_avg_NO_SE[vp_avg_NO_SE$Sample_Type != 'Diff',]

VPCL_AVG_NO_SE<- vipcal_avg(vp_avg_NO_SE)

###2.4.2 DIFF VURVE AND THEN VIPCAL no se
vp_avg_NO_SE<- df_avg_tp(data)
VPCL_AVG_NO_SE_Diff<- vipcal_avg(vp_avg_NO_SE)


#####2.5 VIPCAL_LMER
