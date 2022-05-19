setwd("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R")

{
  library(readxl)
library(quantmod)
library(magrittr)
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(gridExtra)
}


{
data<- as.data.frame(read_excel("viralproduction_R.xlsx")) 
data<- data[, c(1:4,8:10)]

#Adding the first highest count (1E10) and the last lowest count (-1E10)
data[7,]<- c(-1e+2, rep(1e+10, times=6)) #high
data[8,]<- c(1e+2, rep(-1e+10, times=6)) #low
data<-as.data.frame(data[c(7, 1:6,8),])
}

####VIPCAL####
viralp<- c()
for ( i in 2:7){
  p<- quantmod::findPeaks(data[,i])-1
  print(p)
  v<- quantmod::findValleys(data[,i])-1
  print(v)

print (i)

if (identical(length(p),length(q))) {
  print("great")
}

if (length(p)==1) {
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

lmmodel <- lm(Timepoints ~ VP_a, lm_data)
summary(lmmodel)
lmmodel2 <- lm(VP_a ~ Timepoints, lm_data)
summary(lmmodel2)
lmmodel2$coefficients[1] #slope
lmmodel2$coefficients[2] #r2

check_model(lmmodel2)

lm_data<-data[2:7,]
#Linear Model
pkr::Slope(data[,1], data[,2])

ggplot(lm_data,aes(Timepoints, VP_a)) +
  geom_point() +
#  geom_line(method='lm')
  geom_smooth(method='lm', se = F) +
  stat_poly_eq(formula = VP_a ~ Timepoints,
               aes(label = paste(..eq.label.., ..rr.label..,
                                 sep = "~~~")),
               parse=TRUE)

par(mfrow=c(1,3))
mat <- matrix(c(1, 2, 3  # First, second
               ), # and third plot
             nrow = 1, ncol = 1,
             byrow = TRUE)

layout(mat = mat)
#Linear models 
lm_data<- data[1:6,]
lm_data$VP_mean <- apply(lm_data[,2:4], 1, mean)
lm_data$VP_sd <- apply(lm_data[,2:4], 1, sd)
lm_data$VPC_mean <- apply(lm_data[,5:7], 1, mean)
lm_data$VPC_sd <- apply(lm_data[,5:7], 1, sd)

for (i in 2:7) {
  model<- lm(lm_data[,i] ~ Timepoints, lm_data)
  slope<- model$coefficient[2]
  r2<- summary(model)$r.squared
  print(slope)
  print(r2)
}

check_model(lmmodel2)


a<-ggplot(lm_data,aes(Timepoints, VP_a)) +
  geom_point() +
  geom_smooth(method='lm', se=F, formula = y ~x) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..))  +
  theme_minimal()

b<-ggplot(lm_data,aes(Timepoints, VP_b)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..))  +
  theme_minimal()

c<-ggplot(lm_data,aes(Timepoints, VP_c)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..)) +
  theme_minimal()

d<-ggplot(lm_data,aes(Timepoints, VP_mean)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..)) +
  geom_errorbar(aes(ymin= VP_mean-VP_sd, ymax=VP_mean+VP_sd), width=.2,
                position=position_dodge(.9)) +
  theme_minimal()

e<-ggplot(lm_data,aes(Timepoints, VPC_a)) +
  geom_point() +
  geom_smooth(method='lm', se=F, formula = y ~x) +
  stat_regline_equation(label.y= 4.2e+06, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 4.18e+06, aes(label = ..rr.label..))  +
  theme_minimal()

f<-ggplot(lm_data,aes(Timepoints, VPC_b)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 4.2e+06, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 4.18e+06, aes(label = ..rr.label..))  +
  theme_minimal()

g<-ggplot(lm_data,aes(Timepoints, VPC_c)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 4.2e+06, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 4.18e+06, aes(label = ..rr.label..)) +
  theme_minimal()

h<-ggplot(lm_data,aes(Timepoints, VPC_mean)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 4.2e+06, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 4.18e+06, aes(label = ..rr.label..)) +
  geom_errorbar(aes(ymin= VPC_mean-VPC_sd, ymax=VPC_mean+VPC_sd), width=.2,
                position=position_dodge(.9)) +
  theme_minimal()

grob_list<- list(a,b,c,d,e,f,g,h)
#can arrange the above three plots in one plot
windows(width = 45, height = 20)
gridExtra::grid.arrange(grobs=grob_list, ncol=4)

colnames<- colnames(lm_data)[2:7]

for (i in colnames){
plot<- ggplot(lm_data,aes_string(x=lm_data$Timepoints, y=i)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 4.2e+06, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 4.18e+06, aes(label = ..rr.label..)) +
  theme_minimal() +
  labs(x = " Timepoints")
plot_list[[length(plot_list)+1]]<- plot
}
windows(width = 45, height = 20)
gridExtra::grid.arrange(grobs=plot_list, ncol=3)






lm_data<- tidyr::gather(lm_data, 'VP_a', 'VP_b', "VP_c", "Diff_a", "Diff_b", "Diff_c", 
                        key = "Replicate", value = "Counts")
{windows(width=25,height=16)                                                                                    
  par(mfrow=c(1,3), pty = 's')}

ggplot(lm_data,aes(x= Timepoints, y = Counts, color= as.factor(Replicate))) +
  geom_point() +
  geom_smooth(method='lm', se=F, formula = y ~x) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.1e+07, aes(label = ..rr.label..)) +
  facet_wrap(~Replicate)


par(mfrow=c(2,2))
