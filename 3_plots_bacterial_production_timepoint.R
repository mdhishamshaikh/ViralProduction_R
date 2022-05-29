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
  mutate(Subgroup = if_else(count == 'c_Bacteria' | count == 'c_Viruses', "Parent", "Subgroup"))%>%
  arrange(Location, Expt_No, Depth, Timepoint, count, n, Sample_Type)



rm(df_mean)
rm(df_sd)


TP<-c()
TP<- unique(df$Timepoint)
colnames<- c()
for (col in 2: length(TP)){
  a<- paste("T", TP[1], "_T", TP[col], sep = "")
  colnames[length(colnames)+1]<- a
}
colvalues<- c()
for (col in 2: length(TP)){
  a<- paste("T", TP[1], ":T", TP[col], sep = "")
  colvalues[length(colvalues)+1]<- a
}
ncol<-c()
ncol<- ncol(df2)
df[colnames]<- NA

df2<- df


#df<- df%>%

for (i in 1:5){
  for(j in 1:5){
  
    return(df2[, ncol+i] <- case_when(df2$Timepoint == TP[j] ~ colvalues[i]))
  
  if(j==1+i){
    next
  }
  
  }
  
}

for(i in 1:5){
  
  return(df2[, 12] <- case_when(df2$Timepoint == TP[i] ~ colvalues[i]))
  
  next
}




df2[,ncol+1]<- case_when(df2$Timepoint == TP[1] ~ colvalues[1],
                         df2$Timepoint == TP[2] ~ colvalues[1])


df2[,ncol+2]<- case_when(df2$Timepoint == TP[1] ~ colvalues[2],
                         df2$Timepoint == TP[2] ~ colvalues[2],
                         df2$Timepoint == TP[3] ~ colvalues[2])

df2[,ncol+3]<- case_when(df2$Timepoint == TP[1] ~ colvalues[3],
                         df2$Timepoint == TP[2] ~ colvalues[3],
                         df2$Timepoint == TP[3] ~ colvalues[3],
                         df2$Timepoint == TP[4] ~ colvalues[3])

df2[,ncol+4]<- case_when(df2$Timepoint == TP[1] ~ colvalues[4],
                         df2$Timepoint == TP[2] ~ colvalues[4],
                         df2$Timepoint == TP[3] ~ colvalues[4],
                         df2$Timepoint == TP[4] ~ colvalues[4],
                         df2$Timepoint == TP[5] ~ colvalues[4])

df2[,ncol+5]<- case_when(df2$Timepoint == TP[1] ~ colvalues[5],
                         df2$Timepoint == TP[2] ~ colvalues[5],
                         df2$Timepoint == TP[3] ~ colvalues[5],
                         df2$Timepoint == TP[4] ~ colvalues[5],
                         df2$Timepoint == TP[5] ~ colvalues[5],
                         df2$Timepoint == TP[6] ~ colvalues[5])


df2<- df2 %>%
  pivot_longer(cols = colnames, names_to = "Time_Range", values_to = "Time_Time")%>%
  drop_na()

rm('colnames', 'colvalues', 'TP', 'colnames_mean', 'colnames_sd', 'a', 'ncol')
 
  
  

return(df2)
}


print()
