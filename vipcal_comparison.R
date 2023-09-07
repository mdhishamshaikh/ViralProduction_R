
#1. VIPCAL Separate replicates
#Here we calculate VIPCAL on each replicate followed by avergaing VP and VPC. 
#Diff is calculated by subtracting average VPC and VP

i<- VPCL_SR_AVG

#2. VIPCAL AVERGAE replicates. VP and VPC are averaged first. VIPCAL is then calucalted.
#Diff is calcualted by subtracting VPC-VIPCAL and VP-VIPCAL
####SE was used for this method

j<- VPCL_AVG

#3. Differnce curve is calculated first by subtracting average VPC and VP per timepoint.
#VIPCAL is then calculated for VP, VPC and Diff
####SE was used for this method

k<- VPCL_AVG_Diff

#4. Just like #2 but VPCL is not calculted using SE

l<- VPCL_AVG_NO_SE 

#5. Just like #3 but VPCL is not calculated using SE

m<- VPCL_AVG_NO_SE_Diff

#6. LMER Extraction of DIFF followed by SE VIPCAL calc

n<- VIPCAL_SR_DIFF_LMER_SE

#7. LMER Extraction of Diff followed by No SE VIPCAL calc

o<- VIPCAL_SR_DIFF_LMER_No_SE


df_list<- list(i, j, k, l, m , n, o)
df_list<- list(VPCL_SR_avg, VPCL_AVG, VPCL_AVG_Diff, VPCL_AVG_NO_SE,
               VPCL_AVG_NO_SE_Diff, VIPCAL_SR_DIFF_LMER_SE, VIPCAL_SR_DIFF_LMER_No_SE)

df_list<- lapply(df_list, function(x) arrange( x, 'Location',
                                   'Expt_No',
                                   'Depth',
                                   factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
                                   factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
                                   factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
)


names(df_list)<- c('VPCL_SR_avg', 'VPCL_AVG', 'VPCL_AVG_Diff', 'VPCL_AVG_NO_SE',
             'VPCL_AVG_NO_SE_Diff', 'VIPCAL_SR_DIFF_LMER_SE', 'VIPCAL_SR_DIFF_LMER_No_SE')
                                                                                                                                                
vipcal_df<- data.frame()
for (i in df_names) {
  df<- i
  vipcal_df [, i]<- df$VP_Mean
}


vipcal_df<- lapply(df_list, function(x) x$VP_Mean)
vipcal_df<- data.frame(matrix(unlist(vipcal_df), ncol=length(vipcal_df)))
colnames(vipcal_df)<- names(df_list)

vipcal_df2<- pivot_longer(vipcal_df, cols = 1:2, values_to = 'VP', names_to = 'VP_Method')

ggplot(data = vipcal_df2, aes(x = VP_Method, y = VP)) +
  geom_point()


summary(lm(vipcal_df$VIPCAL_SR_DIFF_LMER_No_SE~ vipcal_df$VPCL_AVG_NO_SE_Diff))
summary(lm(vipcal_df$VIPCAL_SR_DIFF_LMER_SE~ vipcal_df$VPCL_AVG_Diff))


ggplot(data = vipcal_df, aes(x = VPCL_AVG_NO_SE_Diff, y = VIPCAL_SR_DIFF_LMER_SE))+
  geom_point()+
  geom_abline(intercept = 0)+
  
  geom_smooth(method = 'lm')


scale_x_continuous(limits = c(-100, 10e+4))+
  scale_y_continuous(limits = c(-100, 10e+4))+
summary(lm(vipcal_df$VPCL_AVG_NO_SE_Diff~ vipcal_df$VIPCAL_SR_DIFF_LMER_SE))


pivot_longer(df_list[1],cols = 'VP_Mean', names_to = 'VP_Type', values_to = 'VP')



list2env(df_list, .GlobalEnv)



na_list<- list()

na<- sum(is.na(i$Mean))
na_list[length(na_list)+1]<- na

na<- sum(is.na(j$VPCL_AVG))
na_list[length(na_list)+1]<- na

na<- sum(is.na(k$VPCL_AVG))
na_list[length(na_list)+1]<- na

na<- sum(is.na(l$VPCL_AVG))
na_list[length(na_list)+1]<- na

na<- sum(is.na(m$VPCL_AVG))
na_list[length(na_list)+1]<- na

na<- sum(is.na(n$VPCL_LMER_Diff_Slope))
na_list[length(na_list)+1]<- na

na<- sum(is.na(o$VPCL_LMER_Diff_Slope))
na_list[length(na_list)+1]<- na


na_list

na_df<- data.frame(df= c('i', 'j', 'k', 'l', 'm', 'n', 'o'),
                  na = unlist(na_list,))
na_df
summary(na_df)

barplot(na_df$na)


ggplot()+
  geom_point(aes(x = g$VP, y = f$VP, col = g$Sample_Type, shape = g$Population))+
  geom_abline(intercept = 0)+
  scale_x_continuous(limits = c(-100, 10e+4))+
  scale_y_continuous(limits = c(-100, 10e+4))



#m vs n

m1<- m %>% arrange('Location',
          'Expt_No',
          'Depth',
          factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
          factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
          factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
n1<- n%>% arrange('Location',
             'Expt_No',
             'Depth',
             factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
             factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
             factor(Time_Range, levels = c('T0_T3', 'T0_T6', 'T0_T9', 'T0_17', 'T0_T24')))
ggplot()+
  geom_point(aes(x = m1$VPCL_AVG, y = n1$VPCL_LMER_Diff_Slope, col = m$Sample_Type, shape = n$Population))+
  geom_abline(intercept = 0)

sum(is.na(m1)) 
sum(is.na(n1))
  
#If NAs were included as zeros
m2<- m1
(m2[is.na(m2)]<- 0)
m2

n2<- n1
n2[is.na(n2)]<-0
n2    

ggplot()+
  geom_point(aes(x = m2$VPCL_AVG, y = n2$VPCL_LMER_Diff_Slope, col = m2$Sample_Type, shape = n2$Population))+
  geom_abline(intercept = 0)    

ggplot()+
  geom_point(aes(x = m2$VPCL_AVG, y = n2$VPCL_LMER_Diff_Slope, col = m2$Sample_Type, shape = n2$Population))+
  geom_abline(intercept = 0)    +
scale_x_continuous(limits = c(-100, 10e+4))+
  scale_y_continuous(limits = c(-100, 10e+4))




