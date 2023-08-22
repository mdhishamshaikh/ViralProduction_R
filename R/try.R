library(tidyverse)
nj_vp<- read.csv("./VI3/vp_calc_all.csv")

nj_vp<- vp_calc
unique(nj_vp$VP_Type)

for (len in unique(nj_vp$VP_Type)){
  print(sum(nj_vp$VP_Type == len))
}
for (len in unique(nj_vp$Sample_Type)){
  print(sum(nj_vp$Sample_Type == len))
}
for (len in unique(nj_vp$Population)){
  print(sum(nj_vp$Population == len))
}
for (len in unique(nj_vp$Time_Range)){
  print(len)
  print(sum(nj_vp$Time_Range == len))
}


#MANOVA

hist(nj_vp$VP)

qqnorm(nj_vp$VP , main = "normal")
qqline(nj_vp$VP)

shapiro.test(nj_vp$VP[1:5000])

ks.test(nj_vp$VP, 'pnorm')




v1<- nj_vp%>%filter(Population == 'c_Viruses'& Sample_Type == 'Diff'& Time_Range == 'T0_T24'&
                    VP_Type == c('LM_AP', 'VPCL_AR_Diff_No_SE', 'VPCL_AR_Diff_LMER_SE'))
hist(v1$VP)
qqnorm(v1$VP, main = 'normal')
qqline(v1$VP)

shapiro.test(v1$VP)

model <- lm(VP ~ VP_Type + Sample_Type + Population + Time_Range + 
              VP_Type*Population + VP_Type*Sample_Type + VP_Type*Time_Range+
              Sample_Type*Population + Sample_Type*Time_Range+
              Population*Time_Range, data = nj_vp)
model <- lm(VP ~ VP_Type + Sample_Type  + Time_Range 
             + VP_Type*Sample_Type + VP_Type*Time_Range
               + Sample_Type*Time_Range, data = v1)
model <- lm(VP ~ VP_Type + Sample_Type  + Time_Range 
            , data = v1)
model <- lm(VP ~ VP_Type   
            , data = v1)
summary(model)
par(mfrow=c(2,2))
plot(model)

anova_model<- anova(model)
anova_model
library(emmeans)

marginal<- emmeans(model, ~ VP_Type)
pairs(marginal)
otp<- as.data.frame( pairs(marginal))

TukeyHSD(anova_model)

.















df<- nj_vp %>%
  filter(VP_Type == c('LM_AP', 'VPCL_AR_Diff_No_SE', 'VPCL_AR_Diff_LMER_SE') )
unique(df$VP_Type)


ggplot(data = df, aes(x = VP_Type, y = VP,fill = VP_Type ))+
  geom_violin(color= NA)+
  geom_jitter(size = 0.5, width = 0.1)+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_lancet()+
  geom_hline(yintercept = 0)+
  labs(title = 'Comparison of viral production calculation methods')+
  xlab('VP Methods')+
  ylab('Viral Production (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))
  
  ylim(-5e+04/2, +5e+04/2)

  
  df1<- df%>% filter(Sample_Type == 'VP') %>% select(-c(VP_SE, VP_R_Squared))%>%
    arrange(Location, Expt_No, Depth, Population, Sample_Type, Time_Range)%>%
    pivot_wider(names_from =  'VP_Type', values_from = 'VP')
  
  ggplot(data = df1)+
    geom_point(aes(x = LM_AP, y = VPCL_AR_Diff_LMER_SE))
  
                
  
  nj<- read.csv("./VI3/vp_calc_bp.csv")

  nj<- nj%>%filter(VP_Type == c('LM_AP', 'VPCL_AR_Diff_No_SE', 'VPCL_AR_Diff_LMER_SE') )
  ggplot(data = nj, aes(x = VP_Type, y = VP, col = Population, shape = VP_Type))+
    geom_point()
  
  