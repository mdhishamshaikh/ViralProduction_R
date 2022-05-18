setwd("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R")

{
  library(readxl)
library(quantmod)
library(magrittr)
  library(ggplot2)
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


ggplot(lm_data,aes(Timepoints, VP_a)) +
  geom_point() +
  geom_smooth(method='lm', se=F, formula = y ~x) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..))

ggplot(lm_data,aes(Timepoints, VP_b)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..))

ggplot(lm_data,aes(Timepoints, VP_c)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_regline_equation(label.y= 1.2e+07, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y= 1.18e+07, aes(label = ..rr.label..))


#Linear models 
lm_data<- data[2:7,]
for (i in 2:7) {
model<- lm(lm_data[,i] ~ Timepoints, lm_data)
slope<- model$coefficient[2]
r2<- summary(model)$r.squared
print(slope)
print(r2)
}

check_model(lmmodel2)


