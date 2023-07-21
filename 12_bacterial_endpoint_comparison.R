library(tidyverse)
library(ggsci)
library(moments)
library(ggstatsplot)


nj_vp<- read.csv("./V5000/vp_calc_bp.csv")
nj_vp2<- read.csv("./V5000/vp_calc_all.csv")

select_df<- nj_vp %>%filter(Population == 'c_Viruses', VP_Type == 'LM_AP' | VP_Type == 'VPCL_AR_Diff_No_SE' | VP_Type == 'VPCL_AR_Diff_LMER_SE') %>%
  mutate(Kind = 'Endpoint')

t24_df<- nj_vp2 %>% filter(Population == 'c_Viruses', Time_Range == 'T0_T24', VP_Type == 'LM_AP' | VP_Type == 'VPCL_AR_Diff_No_SE' | VP_Type == 'VPCL_AR_Diff_LMER_SE') %>%
  mutate(Kind = 'T24')


df<- rbind(select_df, t24_df)
df<- df %>% arrange(Kind, VP_Type, Sample_Type,
                    Expt_No)

ggplot(df, aes(x = Kind, y = VP))+
  geom_violin()

ggbetweenstats(
  data = df,
  x = Kind,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Selected Mean - ALL Endpoint")

ggbetweenstats(
  data = df %>% filter(VP_Type == 'LM_AP'),
  x = Kind,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test:  Mean - LM_AP Endpoint")

ggbetweenstats(
  data = df %>% filter(VP_Type == 'VPCL_AR_Diff_No_SE'),
  x = Kind,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test:  Mean - VPCL_AR_Diff_No_SE Endpoint")

ggbetweenstats(
  data = df %>% filter(VP_Type == 'VPCL_AR_Diff_LMER_SE', Sample_Type == 'Diff'),
  x = Kind,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test:  Mean - VPCL_AR_Diff_LMER_SE Endpoint")

#no significance




####Line plots
df<- df %>% group_by(Location, Expt_No, Depth, Population, Sample_Type, VP_Type, Kind)
df2<- df%>%select(-c('VP_SE', 'VP_R_Squared', 'Time_Range')) %>% pivot_wider(names_from = 'Kind', values_from = 'VP')

plot(df2$Endpoint ~ df2$T24)

ggplot(data = df2, aes( x= Endpoint, y = T24, color = as.factor(Expt_No), shape = VP_Type))+
  geom_point()+
 # geom_smooth(method = 'lm')+
  geom_abline(slope =1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)


ggplot(data = df %>% filter(VP_Type=='VPCL_AR_Diff_LMER_SE', Sample_Type == 'VP'), aes(x= as.factor(Expt_No), y = VP, fill = Kind))+
  geom_bar(stat = 'identity', position = 'dodge')


df3<- df%>% separate(Time_Range, c(NA, "B") )


ggplot(data = df3, aes(y = VP, x = as.factor(Expt_No), fill = Kind ))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_grid(rows = vars(VP_Type),
             cols = vars(Sample_Type),
             scale = 'free')+
  geom_text(data = df3 %>% filter(Kind == 'Endpoint'),
            aes(x = as.factor(Expt_No), y = -2e+06, label = B
                                ),size = 3)+
  theme_bw()+
  scale_fill_locuszoom()

lms_df<- data.frame(
  Location = rep("NJ2020", 21),
  Expt_No = rep(1:7, 3),
  Depth = rep(1, 21),
  B = c('6_20', '0_24', '0_24', '0_9', '0_12', '0_12', '0_12',
        '0_12', '0_24', '0_9', '0_9', '0_12', '0_12', '0_9',
        rep(NA, 7)),
  Population = rep('c_Viruses', 21),
  Sample_Type = c(rep('VP', 7),
                  rep('VPC', 7),
                  rep('Diff', 7)),
  VP = c(138785.2381,
         2893346.085,
         359805.8512,
         2443283.582,
         88116.56006,
         -341937.2525,
         -68676.61692,
         36690.83156,
         932176.6718,
         148258.7065,
         1008813.077,
         -390405.1173,
         -417419.0273,
         42467.66169,
         -102094.4065,
         -1961169.413,
         -211547.1447,
         -1434470.505,
         -478521.6773,
         -75481.7748,
         111144.2786),
  VP_SE = rep(NA, 21),
  VP_R_Squared = rep(NA, 21),
  VP_Type = rep('LM-S', 21),
  Kind = rep('LM-S', 21)
  
)

lms_bep_df<- data.frame(
  Location = rep("NJ2020", 21),
  Expt_No = rep(1:7, 3),
  Depth = rep(1, 21),
  B = c('6_17', '0_6', '0_12', '0_3', '0_12', '0_12', '0_12',
        '0_17', '0_6', '0_6', '3_3', '0_12', '0_12', '0_9',
        rep(NA, 7)),
  Population = rep('c_Viruses', 21),
  Sample_Type = c(rep('VP', 7),
                  rep('VPC', 7),
                  rep('Diff', 7)),
  VP = c(162261.64,
         -741791.0448,
         -317142.1272,
         422459.1329,
         88116.56006,
         -341937.2525,
         -68676.61692,
         44891.38722,
         -129708.5999,
         148258.7065,
         30160798.83,
         -390405.1173,
         -417419.0273,
         42467.66169,
         -117370.2528,
         612082.4449,
         465400.8337,
         29738339.7,
         -478521.6773,
         -75481.7748,
         111144.2786),
  VP_SE = rep(NA, 21),
  VP_R_Squared = rep(NA, 21),
  VP_Type = rep('LM-S', 21),
  Kind = rep('Endpoint', 21)
  
)

df4<- rbind(lms_df,lms_bep_df, df3)

df4$Sample_Type <- factor(df4$Sample_Type, levels = c('VP', 'VPC','Diff'))
df4$VP_Type<- factor(df4$VP_Type, levels = c('LM-S', 'LM_AP', 'VPCL_AR_Diff_No_SE', 'VPCL_AR_Diff_LMER_SE'))

ggplot(data = df4, aes(y = VP, x = as.factor(Expt_No), fill = Kind ))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_grid(cols = vars(VP_Type),
             rows = vars(Sample_Type),
             scales = 'fixed')+
  geom_text(data = df4 %>% filter(Kind == 'Endpoint'| Kind == 'LM-S'),
            aes(x = as.factor(Expt_No), y = -1e+06, label = B,
              angle = 90  
            ),size = 3)+
  theme_bw()+
  scale_fill_lancet()

which(df4$VP > 1e+7)
df4$VP<- replace(df4$VP, c(32, 39), NA)



ggplot(data = df4 %>% filter(Kind != 'Endpoint'),
       aes(y = VP/1e+04, x = as.factor(Expt_No), fill = Kind, label = paste0(round(VP/1e+04, 2))))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_grid(cols = vars(VP_Type),
             rows = vars(Sample_Type),
             scales = 'fixed')+
  geom_text(data = df4 %>% filter(Kind != 'Endpoint'),
            aes(x = as.factor(Expt_No), y = -100, label = B,
                angle = 90
            ),size = 3)+
  theme_bw()+
  scale_fill_lancet()+
  ylab("VP X 10000")+
  xlab("Expt_No")+
  ylim(c(-500, 600))+
  geom_text(aes(y = ((4e+06 + 1e+3)/1e+04) * sign(VP)), angle = 90)


ggplot(data = df4 %>% filter(Kind == 'Endpoint'),
       aes(y = VP/1e+04, x = as.factor(Expt_No), fill = Kind, label = paste0(round(VP/1e+04, 2))))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_grid(cols = vars(VP_Type),
             rows = vars(Sample_Type),
             scales = 'fixed')+
  geom_text(data = df4 %>% filter(Kind == 'Endpoint'),
            aes(x = as.factor(Expt_No), y = -150, label = B,
                angle = 90
            ),size = 3)+
  theme_bw()+
  scale_fill_lancet()+
  ylab("VP X 10000")+
  xlab("Expt_No")+
  ylim(c(-500, 600))+
  geom_text(aes(y = (400 * sign(VP)), angle = 90))

ggplot(data = df4 ,
       aes(y = VP/1e+04, x = as.factor(Expt_No), fill = Kind))+
  geom_bar(stat = 'identity', position = 'dodge')+
  facet_grid(cols = vars(VP_Type),
             rows = vars(Sample_Type),
             scales = 'fixed')+
  geom_text(data = df4 %>% filter(Kind == 'Endpoint'),
            aes(x = as.factor(Expt_No), y = -150, label = B,
                angle = 90
            ),size = 3)+
  geom_text(data = df4 %>% filter(Kind != 'Endpoint'),
            aes(x = as.factor(Expt_No), y = -400, label = B,
                angle = 90
            ),size = 3)+
  theme_bw()+
  scale_fill_lancet()+
  ylab("VP X 10000")+
  xlab("Expt_No")+
  ylim(c(-500, 600))
  #geom_text(aes(y = (400 * sign(VP)), angle = 90))


