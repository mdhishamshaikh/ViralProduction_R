## Calculating the difference curve for lysogeny using linear mixed effects models 
library(lme4)

lmer_data<- df_sr(NJ1)

#only calculating for total virus for now.

lmer_data<- lmer_data[lmer_data$count =="c_Viruses",]
lmer_data$Sample_Type <- as.factor(lmer_data$Sample_Type)
lmer_data$Timepoint2<- as.factor(lmer_data$Timepoint)


lmer_data[lmer_data$Sample_Type == "VPC",]$Replicate<- replace(lmer_data[lmer_data$Sample_Type == "VPC",]$Replicate, lmer_data[lmer_data$Sample_Type == "VPC",]$Replicate == c("1","2", "3"), c("4", "5", "6"))
lmer_data[order(lmer_data$Replicate),]



lyso_model<- lmer(value~Sample_Type*Timepoint2 + (1 | Replicate), data = lmer_data)
summary(lyso_model)

warm.emm_lyso<- emmeans(lyso_model, ~ Sample_Type|Timepoint2)
warm.emm_lyso

tw.emm_lyso<- pairs(warm.emm_lyso)
tw.emm_lyso



summary(tw.emm_lyso)$SE
as.numeric(unique(lmer_data$Timepoint))


summary(lmer_data)


plot(-(summary(tw.emm_lyso)$estimate))


peaks(-summary(tw.emm_lyso)$estimate, summary(tw.emm_lyso)$SE)
valleys(-summary(tw.emm_lyso)$estimate, summary(tw.emm_lyso)$SE)


#non equal variance model 
library(nlme)

lyso_model2<- nlme(value~Sample_Type*Timepoint2, weights = varIdent(form = ~ 1| lmer_data$Replicate), data = lmer_data)
