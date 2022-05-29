source("vp_functions.R")
NJ1<- read.csv("NJ1.csv")


NJ1_sr<- df_sr_tp(NJ1)

#With the help of separate replicates, I can calculate viral production rates
#for all time point ranges, with three slopes, and then averaging them. 


#1. Linear regression slope for total viruses.
# VP and T0:T24

VP_T0_T24<- NJ1_sr[NJ1_sr$Sample_Type == 'VP',]
VP_T0_T24<- VP_T0_T24[VP_T0_T24$Time_Range == 'T0_T24',]
VP_T0_T24<- VP_T0_T24[VP_T0_T24$count == 'c_Viruses',]

a<- lm(data= VP_T0_T24, value ~ Timepoint)
summary(a)
plot(a)

slope<- a$coefficients[[2]]
intercept<- a$coefficients[[1]]


VPC_T0_T24<- NJ1_sr[NJ1_sr$Sample_Type == 'VPC',]
VPC_T0_T24<- VPC_T0_T24[VPC_T0_T24$Time_Range == 'T0_T24',]
VPC_T0_T24<- VPC_T0_T24[VPC_T0_T24$count == 'c_Viruses',]
b<- lm(data= VPC_T0_T24, value ~ Timepoint)
summary(b)
plot(b)



#for all the time ranges, calculate LM for all viruses, VP and VPC, per replicate



slope_lm_sr(NJ1_sr)
