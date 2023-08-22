library(mnonr)
mnonr(n = 10000, p = 1, ms = 2, mk = 61, Sigma = matrix(c(1,0.5,0.5,1),2,2), initial = NULL)
plot(mnonr::mnonr(n=6,p=1,ms=3,mk=61,Sigma=matrix(c(1,0.5,0.5,1),1,1),initial=NULL))
plot(mnonr::mnonr(n=6,p=1,ms=3,mk=61,Sigma=matrix(c(1,0.5,0.5,1),1,1),initial=NULL))
plot(mnonr::mnonr(n=6,p=1,ms=3,mk=61,Sigma=matrix(c(1,0.5,0.5,1),1,1),initial=NULL))

set.seed(2023)


a<- runif(n = 18, min = -1, max = 1)
a

B<- as.data.frame(a)
B$rep = rep(1:3, 6)

B<- B%>% arrange(rep)
B
B$time = rep(1:6, 3)
B

plot(a)
A<- matrix(a, nrow = 6)
plot(A)

ggplot(B, aes(x = time, y = a  ))+
  geom_point()

plot_df<- function(){
  d<- runif(n=6, min = -1, max = 1)
  d
  D<- as.data.frame(d)
  D
  D$time<-rep(1:6,1)
  D
  D$SE<- runif(n = 6, min = 0, max = 200)
  D
  D$SE2<- (D$SE*D$d)/100
  D
  
  f<- ggplot(D, aes(x = time, y = d))+
    geom_point()+
    geom_errorbar(aes(ymin = d-SE2,
                      ymax = d+SE2))
  return(f)
}
plot_df()

set.seed(2023)
replicate(10,
{plot_df
})



rnorm2<- function(n, mean,se){
  mean+se*scale(rnorm(n))
}

mean_values<- rnorm2(3, -0.5126253, -0.91706747 )
mean_values
0.4829704, -1.3228406, -0.6980057

mean(c(0.4829704, -1.3228406, -0.6980057))
sd(c(0.4829704, -1.3228406, -0.6980057))
plotrix::std.error(mean_values)*sqrt(3)


-0.5126253
-0.5126253 + (-0.91706747) = -1.429693
-0.5126253 - (-0.91706747) = 0.4044422

-0.5126253, -1.429693, 0.4044422

mean(c(-0.5126253, -1.429693, 0.4044422))
sd(c(-0.5126253, -1.429693, 0.4044422))
plotrix::std.error(c(-0.5126253, -1.429693, 0.4044422))

0.9170676/1.732


#we can use to find linear regression values
mean(mean_values)
sd(mean_values)

0.91706747/(3^(1/2))

plotrix::std.error(mean_values)
0.5294692

rnorm2(3, VP, SD)



#Create a single dataframe to calculate viral prodcution rate
set.seed(2023)
A<-  data.frame(VP = runif(n=6, min = -1, max = 1))
A$Time<- 1:6
A$SD_percent<- runif(n = 6, min = 0, max = 100) #UP TO 100% ERROR
A$SD<- (A$VP*A$SD_percent)/100
A$SE<- A$SD/sqrt(3)
ggplot(A, aes(x = Time, y = VP))+
  geom_point()+
  geom_errorbar(aes(ymin = VP-(SD/sqrt(3)),
                    ymax = VP+(SD/sqrt(3))))


B<- apply(A, 1, function(x){x[1]+x[4]*scale(rnorm(3))} )
sd(B[,1])

B<- as.data.frame(B)
str(B)

colnames(B)<- 1:6

B<- B%>% pivot_longer(cols = everything(), names_to = 'Time',
                      values_to = 'VP') 
B<- B%>%
  arrange(Time) %>%
  mutate(Replicate = rep(1:3, 6))



############


#I can make a thousand dataframes using this and save it.
set.seed(2023)
simu_df<- data.frame(matrix(ncol = 4, nrow = 0))
colnames(simu_df)<- c("Timepoint", "Count", "Replicate", "Expt_No")
simu_df<- simu_df %>% mutate_all(as.double)


for (df_no in 1:1000){
  
  try(simu_df<- simu_df %>%full_join(B))
  
  A<- data.frame(Count = runif(n=6, min = -1, max = 1))
  A$Timepoint<- 1:6
  A$SD_percent<- runif(n = 6, min = 0, max = 100) #UP TO 100% ERROR
  A$SD<- (A$Count*A$SD_percent)/100
  #A$SE<- A$SD/sqrt(3)
  
  B<- apply(A, 1, function(x){x[1]+x[4]*scale(rnorm(3))} )
  B<- as.data.frame(B)
  colnames(B)<- 1:6
  B<- B%>% pivot_longer(cols = everything(), names_to = 'Timepoint',
                        values_to = 'Count') 
  B<- B%>%
    arrange(Timepoint) %>%
    mutate(Replicate = rep(1:3, 6))
  B$Expt_No = df_no
  
  
  B<- B%>% mutate_all(as.numeric)
  
  
}

# for (expt in 1:10){
#   plot<- ggplot(data = simu_df %>% filter(Expt_No == expt), aes(x = Time, y = VP))+
#     geom_point()+
#     geom_smooth(method = 'lm')
#   
#   return(plot)
# }

# lm(data = simu_df %>% filter(Expt_No == expt), Time ~ VP )



#Write functions for calcualting VP using the vp_Calc_functions
simu_df<- simu_df%>% arrange(Expt_No)
which(is.na(simu_df))


lm_vp<- list()
slope<- data.frame()

for (expt in 1:999) {
  print(expt)
  df2<- simu_df[simu_df$Expt_No == expt,]
  try(lm<- summary(lm(data = df2, Count ~ as.numeric(Timepoint))))
  slope<- c(expt, lm$coefficients[c(2,4)], lm$r.squared)
  lm_vp[[length(lm_vp)+1]] <- slope
 print('done')
}

slope<- data.frame(t(sapply(lm_vp, c)))
colnames(slope)<- c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared')
slope$VP_Type<- 'LM_AP'
ggplot(data = slope, aes(y = VP_Slope, x = Expt_No))+
  geom_point()+
  geom_smooth(method = 'lm')




#Done with linear. 
Try tomorrow for VIPCAL and VIPCAL_SE
After this , we could also try dividing this into VP and VPC. Randomly. 


#Making aDF_AVG for VIPCAL without SE


vipcal_vp<- list()
vpcl_vp_df<- data.frame()

for(expt in 1:999){
 
  df2<- simu_df[simu_df$Expt_No == expt,]
  df2<- df2 %>% group_by(Expt_No, Timepoint)%>%
    summarise(n = n(), mean= mean(Count), se = plotrix::std.error(Count))
  
  p<- peaks(c(+10e+10, df2$mean, -10e+10))-1
  v<- valleys(c(+10e+10, df2$mean, -10e+10))-1
  if(identical(length(p), length(v))){
    print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
  }else{
    print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
  }
  
  if(length(p)==0){
    vp<- 0
    
  } else if (length(p)==1) {
    vp<- (df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
    
  } else if (length(p)==2) {
    vp<- ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
            (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
    
  } else if (length(p)==3) {
    vp<-  ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
             (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
             (df2$mean[p[3]] - df2$mean[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
  }
  vipcal<- c(expt, vp)
  
  vipcal_vp[[length(vipcal_vp)+1]] <- vipcal
  
  
}


vpcl_vp_df<- data.frame(t(sapply(vipcal_vp, c)))
colnames(vpcl_vp_df)<- c( 'Expt_No',  'VP')
vpcl_vp_df$VP<- as.numeric(vpcl_vp_df$VP)

vpcl_vp_df
vpcl_vp_df$VP_Type <- 'VPCL_AR'


# VIPCAL woth SE


vipcal_vp2<- list()
vpcl_vp_df2<- data.frame()

for(expt in 1:999){
  
  df2<- simu_df[simu_df$Expt_No == expt,]
  df2<- df2 %>% group_by(Expt_No, Timepoint)%>%
    summarise(n = n(), mean= mean(Count), se = plotrix::std.error(Count))
  
  p<- peaks_se(c(+10e+10, df2$mean, -10e+10),
               c(0, df2$se, 0))-1
  v<- valleys_se(c(+10e+10, df2$mean, -10e+10),
                 c(0, df2$se, 0))-1
  
  if(identical(length(p), length(v))){
    print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
  }else{
    print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
  }
  
  if(length(p)==0){
    vp<- 0
    
  } else if (length(p)==1) {
    vp<- (df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
    
  } else if (length(p)==2) {
    vp<- ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
            (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
    
  } else if (length(p)==3) {
    vp<-  ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
             (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
             (df2$mean[p[3]] - df2$mean[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
  }
  
  if(length(p)==0){
    se<- 0
    
  } else if (length(p)==1) {
    se<- (df2$se[p[1]] + df2$se[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
    
  } else if (length(p)==2) {
    se<- ((df2$se[p[1]] + df2$se[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
            (df2$se[p[2]] + df2$se[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
    
  } else if (length(p)==3) {
    se<-  ((df2$se[p[1]] + df2$se[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
             (df2$se[p[2]] + df2$se[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
             (df2$se[p[3]] + df2$se[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
  }
  
  vipcal2<- c(expt, vp, se)
  
  vipcal_vp2[[length(vipcal_vp2)+1]] <- vipcal2
  
  
}


vpcl_vp_df2<- data.frame(t(sapply(vipcal_vp2, c)))
colnames(vpcl_vp_df2)<- c( 'Expt_No',  'VP', 'VP_SE')
vpcl_vp_df2[, c('VP', 'VP_SE')]<- lapply(vpcl_vp_df2[, c('VP', 'VP_SE')], as.numeric)

vpcl_vp_df2
vpcl_vp_df2$VP_Type <- 'VPCL_AR_SE'



plot(vpcl_vp_df$VP ~vpcl_vp_df2$VP)

slope
vpcl_vp_df
vpcl_vp_df2

simu_df2<- rbind(slope[,c(1,2,5)], vpcl_vp_df, vpcl_vp_df2[, c(1,2,4)])


ggbetweenstats(
  data = simu_df2,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

labs(title = "Kruskal-Wallis Test: VPCL Mean Diff")
simu_df3<- simu_df2 %>% filter(VP_Type != 'LM_AP')

t.test(VP ~ VP_Type, data = simu_df3)





#####I. Simulation####
#####1.0 Installing Packages####
source("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/sourcesourcebaby.R")
library(tidyverse)
library(ggstatsplot)
library(ggsci)
#library(moments)



####1.1 Creating a dataframe####
#Set number of datafraes you'd like to create

length_df<- 1000

set.seed(2023)
simu_df<- data.frame(matrix(ncol = 5, nrow = 0))
colnames(simu_df)<- c("Timepoint", "Count", "Replicate", "Expt_No", "Sample_Type")
simu_df<- simu_df %>%  mutate(across(!Sample_Type, as.numeric)) %>%
  mutate(across(Sample_Type, as.character))


for (df_no in 1:length_df){
  
  A<- data.frame(Count = runif(n=12, min = -1, max = 1))
  A$Timepoint<- rep(1:6,2)
  A$SD_percent<- runif(n = 12, min = 0, max = 100) #UP TO 100% ERROR
  A$SD<- (A$Count*A$SD_percent)/100
 
  
  B<- apply(A, 1, function(x){x[1]+x[4]*scale(rnorm(3))} )
  B<- as.data.frame(B)
  colnames(B)<- rep(1:6,2)
  B<- B%>% pivot_longer(cols = everything(), names_to = 'Timepoint',
                        values_to = 'Count') 
  B<- B%>%
    arrange(Timepoint) %>%
    mutate(Replicate = rep(1:6, 6))%>%
    mutate(Sample_Type = rep(c(rep('VP',3), rep('VPC',3)),6))
  
  B$Expt_No = df_no
  
  
  B<- B%>% mutate(across(!Sample_Type, as.numeric))
  
  try(simu_df<- simu_df %>%full_join(B))
  rm(A,B, df_no)
}



####1.2 LM_AP####

lm_vp<- list()
slope<- data.frame()

for (expt in 1:length_df) {
  print(expt)
  df2<- simu_df %>% filter(Expt_No == expt)
  
  for(type in c('VP', 'VPC')){
    
    df3<- df2 %>% filter(Sample_Type == type)
  
  try(lm<- summary(lm(data = df3, Count ~ as.numeric(Timepoint))))
  slope1<- c(expt, type, lm$coefficients[c(2,4)], lm$r.squared)
  lm_vp[[length(lm_vp)+1]] <- slope1
  }
  print('done')
  rm(df2, df3, expt, type, slope1, lm)
}

slope<- data.frame(t(sapply(lm_vp, c)))
colnames(slope)<- c('Expt_No', 'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared')
slope$VP_Type<- 'LM_AP'
slope<- slope%>% mutate(across(c('Expt_No', 'VP', 'VP_SE', 'VP_R_Squared'), as.numeric))


#calculate Lysogeny from here.

df<- slope


VP<- df[df$Sample_Type == "VP",]
VPC<- df[df$Sample_Type == "VPC",]
Diff<- as.data.frame(cbind(VP[,c('Expt_No', 'Sample_Type')], VPC[,'VP'] - VP[,'VP'], VPC[,'VP_SE'] + VP[,'VP_SE']))
Diff$Sample_Type<- "Diff"
colnames(Diff)[3]<- 'VP'
colnames(Diff)[4]<- 'VP_SE'
Diff$VP_R_Squared <- NA
Diff$VP_Type <- 'LM_AP'

LM_AP_output<- full_join(df, Diff, by = NULL)
rm(df, Diff, VP, VPC, slope, lm_vp)

# 
# ggplot(data = LM_AP_output, aes(x = Sample_Type, y = VP))+
#   geom_violin()+
#   geom_point()


####1.3 VIPCAL (VPCL_AR_No_SE)####



df<- simu_df
df<- df %>% group_by(Expt_No, Sample_Type, Timepoint )%>%
  summarise(n = n(), mean = mean(Count), se = plotrix::std.error(Count) )

VP<- df[df$Sample_Type == "VP",]
VPC<- df[df$Sample_Type == "VPC",]


Diff<- as.data.frame(cbind(VP[,c('Expt_No', 'Sample_Type', 'Timepoint')], VPC[,'mean'] - VP[,'mean'], VPC[,'se'] + VP[,'se']))
Diff$Sample_Type<- "Diff"
colnames(Diff)[4]<- 'mean'
colnames(Diff)[5]<- 'se'
Diff$VP_R_Squared <- NA


df<- full_join(df, Diff, by = NULL)


rm(Diff, VP, VPC)

vipcal_vp3<- list()
vpcl_vp_df3<- data.frame()



# %>% group_by(Expt_No, Sample_Type) %>%
#   arrange(Expt_No, Sample_Type, Timepoint)

for(expt in 1:length_df){
  
  for(type in c('VP', 'VPC', 'Diff')){
  
  df2<- df %>% filter(Expt_No == expt,
                       Sample_Type == type)
    
  
     # summarise(n = n(), mean= mean(Count), se = plotrix::std.error(Count))
  
  p<- peaks(c(+10e+10, df2$mean, -10e+10))-1
  v<- valleys(c(+10e+10, df2$mean, -10e+10))-1
  if(identical(length(p), length(v))){
    print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
  }else{
    print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
  }
  
  if(length(p)==0){
    vp<- 0
    
  } else if (length(p)==1) {
    vp<- (df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]])
    
  } else if (length(p)==2) {
    vp<- ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
            (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]))/2
    
  } else if (length(p)==3) {
    vp<-  ((df2$mean[p[1]] - df2$mean[v[1]])/(df2$Timepoint[p[1]] - df2$Timepoint[v[1]]) + 
             (df2$mean[p[2]] - df2$mean[v[2]])/(df2$Timepoint[p[2]] - df2$Timepoint[v[2]]) +
             (df2$mean[p[3]] - df2$mean[v[3]])/(df2$Timepoint[p[3]] - df2$Timepoint[v[3]]))/3
  }
  vipcal<- c(expt, type, vp)
  
  vipcal_vp3[[length(vipcal_vp3)+1]] <- vipcal
  
  
  }
}


vpcl_vp_df3<- data.frame(t(sapply(vipcal_vp3, c)))
colnames(vpcl_vp_df3)<- c( 'Expt_No', 'Sample_Type', 'VP')
vpcl_vp_df3$VP<- as.numeric(vpcl_vp_df3$VP)

vpcl_vp_df3
vpcl_vp_df3$VP_Type <- 'VPCL_AR_Diff_No_SE'
vpcl_vp_df3$VP_R_Squared<- NA
vpcl_vp_df3$VP_SE<- NA

VPCL_AR_Diff_No_SE_output_df<- vpcl_vp_df3

rm(df, df2, vipcal_vp3, vpcl_vp_df3, expt, type)

# ggplot(data = VPCL_AR_No_SE_output_df, aes(x = Sample_Type, y = VP))+
#   geom_violin()+
#   geom_point()


####1.4 VIPCAL SE (VPCL_AR_LMER_SE)####


vipcal_vp4<- list()
vpcl_vp_df4<- data.frame()



df<- simu_df
for(expt in 1:length_df){
  
  
  df3<- df %>% filter(Expt_No == expt)
  lmer_data<- data.frame()
  model<- lme4::lmer(data = df3, Count ~ Sample_Type*as.factor(Timepoint) + (1  | Replicate))
  emmeans<- emmeans::emmeans(model, ~ Sample_Type|as.factor(Timepoint))
  contrast<- pairs(emmeans) 
  dataf1<- data.frame(
    rep("Diff", length(unique(df$Timepoint))),
    summary(contrast)$Timepoint, 
    -(summary(contrast)$estimate),
    summary(contrast)$SE
  )
  colnames(dataf1)<- c( "Sample_Type", "Timepoint", "Mean", "SE")
  st<- summary(emmeans)[1]
  tp<- summary(emmeans)[2]
  mean<- summary(emmeans)[3]
  se<- summary(emmeans)[4]
  dataf2<- data.frame(st,tp,mean,se)
  colnames(dataf2)<- c("Sample_Type", "Timepoint", "Mean", "SE")
  
  
  df4<- rbind(lmer_data, dataf2, dataf1)
  rm(dataf1, dataf2, lmer_data, st, tp)
 
  
  for(type in c('VP', 'VPC', 'Diff')){
    
    df5<- df4 %>% filter(Sample_Type == type)
    
    
    # summarise(n = n(), mean= mean(Count), se = plotrix::std.error(Count))
    p<- peaks_se(c(+10e+10, df5$Mean, -10e+10),
                 c(0, df5$SE, 0))-1
    v<- valleys_se(c(+10e+10, df5$Mean, -10e+10),
                   c(0, df5$SE, 0))-1
    
     
   
    if(identical(length(p), length(v))){
      print(paste0("Number of peaks and valleys identified are the same: ", length(p)))
    }else{
      print("Number peaks and valleys are not the same. This will lead to erroneous viral production calculations")
    }
    
    if(length(p)==0){
      vp<- 0
      
    } else if (length(p)==1) {
      vp<- (df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]])
      
    } else if (length(p)==2) {
      vp<- ((df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
              (df5$Mean[p[2]] - df5$Mean[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]))/2
      
    } else if (length(p)==3) {
      vp<-  ((df5$Mean[p[1]] - df5$Mean[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
               (df5$Mean[p[2]] - df5$Mean[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]) +
               (df5$Mean[p[3]] - df5$Mean[v[3]])/(df5$Timepoint[p[3]] - df5$Timepoint[v[3]]))/3
    }
    
    
    if(length(p)==0){
      se<- 0
      
    } else if (length(p)==1) {
      se<- (df5$SE[p[1]] + df5$SE[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]])
      
    } else if (length(p)==2) {
      se<- ((df5$SE[p[1]] + df5$SE[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
              (df5$SE[p[2]] + df5$SE[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]))/2
      
    } else if (length(p)==3) {
      se<-  ((df5$SE[p[1]] + df5$SE[v[1]])/(df5$Timepoint[p[1]] - df5$Timepoint[v[1]]) + 
               (df5$SE[p[2]] + df5$SE[v[2]])/(df5$Timepoint[p[2]] - df5$Timepoint[v[2]]) +
               (df5$SE[p[3]] + df5$SE[v[3]])/(df5$Timepoint[p[3]] - df5$Timepoint[v[3]]))/3
    }
    vipcal<- c(expt, type, vp, se)
    
    vipcal_vp4[[length(vipcal_vp4)+1]] <- vipcal
    
    
    }
    
  rm(df5, vp, se, expt, type, vipcal, p, v)
}   



vpcl_vp_df4<- data.frame(t(sapply(vipcal_vp4, c)))
colnames(vpcl_vp_df4)<- c( 'Expt_No', 'Sample_Type', 'VP', 'VP_SE')
vpcl_vp_df4$VP<- as.numeric(vpcl_vp_df4$VP)

vpcl_vp_df4
vpcl_vp_df4$VP_Type <- 'VPCL_AR_Diff_LMER_SE'
vpcl_vp_df4$VP_R_Squared<- NA
VPCL_AR_Diff_LMER_SE_output_df <- vpcl_vp_df4
rm(vpcl_vp_df4, vipcal_vp4 )

# ggplot(data = VPCL_AR_LMER_SE_output_df, aes(x = Sample_Type, y = VP))+
#   geom_violin()+
#   geom_point()


simu_output_df<- rbind(LM_AP_output, VPCL_AR_Diff_No_SE_output_df, VPCL_AR_Diff_LMER_SE_output_df)

ggplot(data = simu_output_df, aes(x = Sample_Type, y = VP))+
  geom_violin()+
  geom_point()

ggbetweenstats(
  data = simu_output_df%>% filter(Sample_Type == 'VP'),
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test:  Simu VPC")

ggbetweenstats(
  data = simu_output_df%>% filter(Sample_Type == 'VP'),
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test:  Simu VPC")

ggbetweenstats(
  data = simu_output_df%>% filter(Sample_Type == 'Diff'),
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Simu Diff")



#https://github.com/GRousselet/rogme
###ROGME####

library(rogme)

simu_output_df$VP_Type<- factor(simu_output_df$VP_Type,
                                levels = c("LM_AP",
                                           "VPCL_AR_Diff_No_SE",
                                           "VPCL_AR_Diff_LMER_SE"))

simu_output_df2<- simu_output_df %>% filter(VP_Type == c("VPCL_AR_Diff_No_SE", "VPCL_AR_Diff_LMER_SE"))
simu_output_df2$VP_Type<- factor(simu_output_df2$VP_Type,
                                 levels = c("VPCL_AR_Diff_LMER_SE","VPCL_AR_Diff_No_SE"))
                          
                          
p<-  plot_scat2(simu_output_df2,
                formula = VP ~ VP_Type,
                xlabel = "VP_Type",
                ylabel = "VP",
                alpha = 1,
                shape = 21,
                colour = "grey10",
                fill = "grey90",
                size = 3) +
  scale_x_discrete(breaks=c("VPCL_AR_Diff_No_SE", "VPCL_AR_Diff_LMER_SE"),
                   labels=c("VPCL_AR_Diff_No_SE", "VPCL_AR_Diff_LMER_SE")) +
  theme(axis.text.y = element_text(angle = 90, hjust = .5))
p + coord_flip()


q<- ggplot(data = simu_output_df2,
           aes(fill = VP_Type, x = VP, col = VP_Type))+
  geom_density(alpha = 0)+
  facet_grid(rows = (simu_output_df2 ))+
  theme_bw()
q


plot_hd_bars(q,
             col = "black",
             q_size = 0.5,
             md_size = 1.5,
             alpha = 1)

p <- plot_hd_bars(p,
                  col = "black",
                  q_size = 0.5,
                  md_size = 1.5,
                  alpha = 1)
p <- p + coord_flip() #> flip axes
pscat <- p
pscat

sf<- shifthd_pbci(data = simu_output_df2 ,
                  formula = VP ~ VP_Type)


psf <- plot_sf(sf, plot_theme = 2)[[1]] 
psf


psf+xlim(c(-1, 1))

dasym <- asymhd(data = simu_output_df2, formula = VP ~ VP_Type, 
                q = seq(5,40,5)/100, alpha = .05, nboot = 100)

#> ggplot
diff_asym <- plot_diff_asym(data = dasym)[[1]]
diff_asym +
  xlim(c(-1, 0.5))



####Running combine dfunction on simu_df 

#As these are mandatory columns, I'll add them here
simu_df2<- simu_df
simu_df2[c("Location", "Depth", "Sample_Type", "Time_Range" )]<- NA

viral_production(data = simu_df2, bp_endpoint = F, output.dir = 'simulations')
