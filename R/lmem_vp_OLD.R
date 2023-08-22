library("lme4")
library("emmeans")
library("ggplot2")
data<- read.csv("lme4_data.csv")

#creating factors

data$group<- as.factor(data$group)
data$time2<- as.factor(data$time)

data
data_m_sd<- read.csv("lme4_data_mean_sd.csv")
data_m_sd$group<- as.factor(data_m_sd$group)
data_m_sd$time2<- as.factor(data_m_sd$time)


ggplot(data = data_m_sd, aes(x= time2, y = y, color = group, group =group)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2)

## build basic model
modl<-lmer(y~group*time2 +(1|ID), data=data)

## see results
summary(modl)



warp.emm <- emmeans(modl, ~ group|time2)
warp.emm

## get contrasts for each time point
tw.emm <- pairs(warp.emm)
tw.emm

plot(summary(tw.emm)$estimate)
contrast(warp.emm, "consec", simple = "each", combine = TRUE, adjust = "mvt")

emm_group<- emmeans(modl, ~ time2|group)
emm_group

grp_pairs<- pairs(emm_group)
grp_pairs


#but i want to check if the contrast between the timepoint is different between lytic and lytic+lysogenic