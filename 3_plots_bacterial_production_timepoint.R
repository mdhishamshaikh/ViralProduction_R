source("vp_functions.R")

df<- read.csv("NJ1.csv")
pull<-df_lm_tp(NJ1)
plots_lm_tp(pull)


#WORK ON THIS



#function for dataframe 
df_lm_tp<- function(df){

df<- df[df$Sample_Type != '0.22',]
df<- gather(df, 'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value") %>%
  group_by(Location, Expt_No, Depth, Sample_Type, Timepoint, count ) %>%
  summarise(n =n(), mean=mean(value), sd=sd(value)) #calculating means and sd

if ('VPC' %in% df$Sample_Type){
  colnames_mean<- c("VP", "VPC")
  colnames_sd<- c("VP", "VPC")
  
} else {
  colnames_mean<- c("VP")
  colnames_sd<- c("VP")
  
}

df_mean<- df[,1:8] %>% #splitting the dataframe cause I haven't figure out how to spread teh table without adding NAs
  spread('Sample_Type', 'mean')
if (length(colnames_mean)==2){
  colnames(df_mean)[7:8]<- colnames_mean
} else if (length(colnames_mean)==1){
  colnames(df_mean)[7]<- colnames_mean
} 

df_sd<- df[,c(1:7,9)] %>%
  spread('Sample_Type', 'sd')
if (length(colnames_sd)==2){
  colnames(df_sd)[7:8]<- colnames_sd
} else if (length(colnames_sd)==1){
  colnames(df_sd)[7]<- colnames_sd
} 

if ('VPC' %in% df$Sample_Type){
  df_mean$Diff <- with(df_mean, VPC-VP) #calcualting Diff mean
  df_mean<- pivot_longer(df_mean, cols = c("VP", "VPC", "Diff"), names_to= 'Sample_Type', values_to='mean_value')
  df_sd$Diff <- with(df_sd, VPC+VP) #Calculating Diff sd, whcih si addition of the other sds
  df_sd<- pivot_longer(df_sd, cols = c("VP", "VPC", "Diff"), names_to='Sample_Type', values_to= 'sd_value')
}

df<- merge(df_mean, df_sd, by= c('Location', 'Expt_No', 'Depth',
                                    'Timepoint', 'count', 'n', 'Sample_Type')) %>%
  mutate(Microbe = if_else(count == 'c_Bacteria' | count == 'c_HNA' | count == 'c_LNA', "Bacteria", "Viruses"))%>%
  mutate(Subgroup = if_else(count == 'c_Bacteria' | count == 'c_Viruses', "Parent", "Subgroup"))

rm(df_mean)
rm(df_sd)



TP<- unique(df$Timepoint)
colnames<- c()
for (col in 2: length(TP)){
  a<- paste("T", TP[1], "_T", TP[col], sep = "")
  colnames[length(colnames)+1]<- a
}
df<- df%>%
  mutate("T0_T3" = case_when(Timepoint == '0' ~ "T0:T3",
                             Timepoint == '3' ~ "T0:T3"))%>%
  
  mutate("T0_T6" = case_when(Timepoint == '0' ~ "T0:T6",
                             Timepoint == '3' ~ "T0:T6", 
                             Timepoint == '6' ~ "T0:T6"))%>%
  
  mutate("T0_T17" = case_when(Timepoint == '0' ~ "T0:T17",
                              Timepoint == '3' ~ "T0:T17", 
                              Timepoint == '6' ~ "T0:T17",
                              Timepoint == '17' ~ "T0:T17"))%>%
  
  mutate("T0_T20" = case_when(Timepoint == '0' ~ "T0:T20",
                              Timepoint == '3' ~ "T0:T20",
                              Timepoint == '6' ~ "T0:T20",
                              Timepoint == '17' ~ "T0:T20",
                              Timepoint == '20' ~ "T0:T20"))%>%
  
  mutate("T0_T24" = case_when(Timepoint == '0' ~ "T0:T24",
                              Timepoint == '3' ~ "T0:T24",
                              Timepoint == '6' ~ "T0:T24",
                              Timepoint == '17' ~ "T0:T24",
                              Timepoint == '20' ~ "T0:T24",
                              Timepoint == '24' ~ "T0:T24")) %>%
  pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
  drop_na()


return(df)
}
