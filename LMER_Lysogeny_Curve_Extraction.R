
#I will use the LMER model to extrcat the difference curve from VP and VPC values.

source("0_vp_source.R")
source("vp_functions.R")
library(sjPlot)
#Importing Dataset

NJ1<- read.csv("NJ1.csv") %>%
  df_sr()%>%
  subset( Microbe == "Viruses" & Sample_Type != '0.22')
NJ1<- read.csv("NJ1.csv") %>%
  df_sr()%>%
  subset( count == "c_V3" & Sample_Type != '0.22')
#Changing the replicate for VPC to 4,5,6
for (i in c(1,2,3)){
NJ1$Replicate[NJ1$Replicate == i & NJ1$Sample_Type == 'VPC'] <- i+3
}

#Converting timepoint and sample type as factor

NJ1$Sample_Type<- as.factor(NJ1$Sample_Type)
NJ1$Timepoint2<- as.factor(NJ1$Timepoint)

NJ1

lmer_model<- function(df, value = count){
  
  lmer_data<- data.frame()
  model_plots<- list()
  n<-0
  for (rep in c(1,2,3)){
    df$Replicate[df$Replicate == rep & df$Sample_Type == 'VPC'] <- rep+3
  }
  
  for (virus in unique(df$count)){
  
  model<- lme4::lmer(data = df, value ~ Sample_Type*as.factor(Timepoint) + (1+ Sample_Type | Replicate))
  plot<- model_plot(model, df = df)
  emmeans<- emmeans::emmeans(model, ~ Sample_Type|as.factor(Timepoint))
  contrast<- pairs(emmeans) 
  df1<- data.frame(rep(virus, 6),
                   rep("Diff", 6),
                   summary(contrast)$Timepoint, 
                   -(summary(contrast)$estimate),
                   summary(contrast)$SE
  )
  colnames(df1)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
  st<- summary(emmeans)[1]
  tp<- summary(emmeans)[2]
  mean<- summary(emmeans)[3]
  se<- summary(emmeans)[4]
  pop<- rep(virus, 6)
  df2<- data.frame(pop,st,tp,mean,se)
  colnames(df2)<- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
  
  
  lmer_data<- rbind(lmer_data, df2, df1)
 # DF<- rbind(df2, df1) %>%
  #  arrange(Timepoint)
  #n<- n +1
  #lmer_data[[n]]<- DF
  #model_plots[[n]]<- plot
  
  }
  return(lmer_data)
  return(model_plots)
  
  
}
model_Data<- lmer_model(NJ1)
model_Data

d<- inner_join(model_Data[[1]], model_Data[[2]], by = c("Sample_Type", "Timepoint") )
d<-  gather(d[,c(1,6)], "Population" )


V<- NJ1[NJ1$count == vir,]

shapiro.test(NJ1[NJ1$count == "c_Viruses",]$value)
shapiro.test(NJ1[NJ1$count == "c_V1",]$value) #parametric
shapiro.test(NJ1[NJ1$count == "c_V2",]$value)
shapiro.test(NJ1[NJ1$count == "c_V3",]$value)



#NLME

nlme


lmer_model<- lme4::lmer(data = NJ1, value ~ Sample_Type*Timepoint2 + (1 + Sample_Type| Replicate))
summary(model)
model_plot(model, NJ1) #for quick visualization of model fit
emmeans<- emmeans::emmeans(model, ~ Sample_Type|Timepoint2)
summary(emmeans)
contrast<- pairs(emmeans)
summary(contrast)
summary(contrast)$estimate
summary(contrast)$SE
summary(contrast)$Timepoint2
Sample_Type2<- rep("Diff",6)
df1<- data.frame(Sample_Type2,
                 summary(contrast)$Timepoint2, 
                 -(summary(contrast)$estimate),
                 summary(contrast)$SE
                 )
colnames(df1)<- c("Sample_Type", "Timepoint", "Mean", "SE")
df1

summary(emmeans)
st<- summary(emmeans)[1]
tp<- summary(emmeans)[2]
mean<- summary(emmeans)[3]
se<- summary(emmeans)[4]
df<- data.frame(st,tp,mean,se)
colnames(df)<- c("Sample_Type", "Timepoint", "Mean", "SE")

DF<- rbind(df, df1) %>%
  arrange(Timepoint)

plot(x= df$Timepoint2, y = df$emmean, pch = 19)


ggplot(DF, aes(x =  Timepoint, y = Mean, group = Sample_Type, color = Sample_Type, shape = Sample_Type))+
  geom_point(size = 2.5)+
  geom_line(size = 0.5)+
  geom_hline(yintercept = 0, color= '#636363', size= 0.3, linetype = "dashed")+
  theme_bw()+
  geom_errorbar(aes(ymin=Mean + SE, ymax= Mean - SE), width = 0.5, size = 0.5)

DF
DF2<- DF[DF$Sample_Type== 'VP',]
peaks(DF2$Mean, DF2$SE)
valleys(DF2$Mean, DF2$SE)

DF3<- DF[DF$Sample_Type== 'VPC',]
peaks(DF3$Mean, DF3$SE)
valleys(DF3$Mean, DF3$SE)

DF4<- DF[DF$Sample_Type== 'Diff',]
peaks(DF4$Mean, DF4$SE)
valleys(DF4$Mean, DF4$SE)
#this works. 
#the question now is if we'll calculate viral production in VP using the SE from the model or individual SE?

for (ST in unique(DF$Sample_Type)) {
  print(ST)
  df<- DF[DF$Sample_Type == ST ,]
  print(df)
  p<- peaks(df$Mean, df$SE)
  v<- valleys(df$Mean, df$SE)
  print(paste(p,v))
  if (identical(length(p),length(q))) {
    print("Number of peaks and valleys are identical. Proceeding to calculating viral production ")
  } else {
    print("Number of peaks and valleys aren't identical")
  }
}



p<- peaks(data[,i],data[,i+3])
print(p)
v<- valleys(data[,i],data[,i+3])
print(v)

print (i)

if (identical(length(p),length(q))) {
  print("Number of peaks and valleys are identical. Proceeding to calculating viral production ")
} else {
  print("Number of peaks and valleys aren't identical")
}

if (length(p)==0) {
  print("No viral production")
} else if (length(p)==1) {
  
  vp<- (data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1])
  
} else if (length(p)==2) {
  vp<- ((data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1]) + 
          (data[p[2],i] - data[v[2],2])/(data[p[2],i] - data[v[2],1]))/2
  
} else if (length(p)==3) {
  vp<-  ((data[p[1],i] - data[v[1],i])/(data[p[1],1] - data[v[1],1]) + 
           (data[p[2],i] - data[v[2],i])/(data[p[2],1] - data[v[2],1]) +
           (data[p[3],i] - data[v[1],i])/(data[p[3],1] - data[v[3],1]))/3
}

print(vp)
viralp = c(viralp, vp)
}


MAKE A MASSIVE FLOW CHART FOR CALCULATIONS PER STATION. 
AND THEN RUN THAT OVER MADE UP DATAFRAMES