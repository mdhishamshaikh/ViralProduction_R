library(tidyverse)
library(ggsci)
library(moments)
library(ggstatsplot)


nj_vp<- read.csv("./V5000/vp_calc_all.csv")
nj_vp$VP_Type<- factor(nj_vp$VP_Type, levels = c("LM_AP",
                                                 "LM_SR_AVG",
                                                 "LM_AR",
                                                 "LM_AR_Diff",
                                                 "LM_AR_Diff_LMER",
                                                 "VPCL_SR_AVG",
                                                 "VPCL_AR_No_SE",
                                                 "VPCL_AR_SE",
                                                 "VPCL_AR_Diff_No_SE",
                                                 "VPCL_AR_Diff_SE",
                                                 "VPCL_AR_Diff_LMER_No_SE",
                                                 "VPCL_AR_Diff_LMER_SE"))




lm_df<- nj_vp[grep('LM_', nj_vp$VP_Type),] %>%
  filter(Population == 'c_Viruses')

ggplot(data = lm_df, aes(x = VP_Type, y = VP_SE,fill = VP_Type ))+
  geom_violin(color= 'black', linewidth = 1)+
  #geom_jitter(size = 0.5, width = 0.1)+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_lancet()+
  geom_hline(yintercept = 0)+
  labs(title ='Linear model viral production calculation methods - SE')+
  xlab('VP Methods')+
  ylab('Viral Production SE (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))


ggplot(data = lm_df, aes(x = VP_Type, y = VP,fill = VP_Type ))+
  geom_violin(color= 'black', linewidth = 1)+
  #geom_jitter(size = 0.5, width = 0.1)+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_lancet()+
  geom_hline(yintercept = 0)+
  labs(title ='Linear model viral production calculation methods - MEAN')+
  xlab('VP Methods')+
  ylab('Viral Production mean (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))


kruskal.test(VP~VP_Type, data = lm_df)
#https://statsandr.com/blog/kruskal-wallis-test-nonparametric-version-anova/
  
ggbetweenstats(
  data = lm_df,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: LM Mean")
  
ggbetweenstats(
  data = lm_df,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: LM SE")

hist(lm_df$VP)
qqnorm(lm_df$VP)
qqline(lm_df$VP)
shapiro.test(lm_df$VP) #NON-NORMAL DATA

skewness(lm_df$VP) #RIGHT SKEWED
kurtosis(lm_df$VP) #LEPTOKURTIC. More otliers than normal distribution. Normal distribution Kurtosis is 3.
jarque.test(lm_df$VP) #Cannot reject the hypothesis. Therefore, non-normal.

#Transformations

shapiro.test(log(lm_df$VP+abs(min(lm_df$VP))+1))
hist(1/(lm_df$VP))
skewness(1/(lm_df$VP))

vpcl_df<- nj_vp[grep('VPCL_', nj_vp$VP_Type),] %>%
  filter(Population == 'c_Viruses')
shapiro.test(vpcl_df$VP)

ggbetweenstats(
  data = vpcl_df,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean")

ggbetweenstats(
  data = vpcl_df,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE")


vpcl_df2<- vpcl_df[grep('_Diff_', vpcl_df$VP_Type),]


unique(vpcl_df2$VP_Type)


ggbetweenstats(
  data = vpcl_df2,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean")

ggbetweenstats(
  data = vpcl_df2,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE")



select_df<- nj_vp %>%filter(Population == 'c_Viruses', VP_Type == c('LM_AP',
                                                           'VPCL_AR_Diff_No_SE',
                                                           'VPCL_AR_Diff_LMER_SE'))

shapiro.test(select_df$VP)

ggbetweenstats(
  data = select_df,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Selected Mean")

ggbetweenstats(
  data = select_df,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Selected SE")


ggbetweenstats(
  data = select_df %>% filter(Sample_Type == 'Diff'),
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Selected Mean")

ggbetweenstats(
  data = select_df %>% filter(Sample_Type == 'Diff'),
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Selected SE")




#DO THIS WITH TRANSFORMATIONS
#WHERE YOU ADD THE ABSOLUTE VALUE OF MIN + 1

#(abs(min(select_df$VP))+1)

hist(log(select_df$VP + abs(min(select_df$VP)) +1))

select_df<- select_df %>% mutate(log_VP = log(select_df$VP + abs(min(select_df$VP)) +1) )
select_df<- select_df %>% mutate(log_SE = log(select_df$VP_SE + abs(min(select_df$VP_SE, na.rm = T)) +1) )
ggbetweenstats(
  data = select_df,
  x = VP_Type,
  y = log_VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean")

ggbetweenstats(
  data = select_df,
  x = VP_Type,
  y = log_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE")



#VPCL SELECT DF
vpcl_select_df<- nj_vp %>%filter(Population == 'c_Viruses', VP_Type == c('VPCL_AR_Diff_SE',
                                                                         'VPCL_AR_Diff_LMER_No_SE',
                                                                    'VPCL_AR_Diff_No_SE',
                                                                    'VPCL_AR_Diff_LMER_SE'))%>% mutate(VP_Type = as.factor(VP_Type))
shapiro.test(vpcl_select_df$VP) #not normal

ggbetweenstats(
  data = vpcl_select_df,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean")

ggbetweenstats(
  data = vpcl_select_df,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE")


vpcl_select_df_VP<- nj_vp %>%filter(Population == 'c_Viruses', 
                                 VP_Type == c('VPCL_AR_Diff_SE',
                                                                         'VPCL_AR_Diff_LMER_No_SE',
                                                                         'VPCL_AR_Diff_No_SE',
                                                                         'VPCL_AR_Diff_LMER_SE'),
                                 Sample_Type == 'VP')

ggbetweenstats(
  data = vpcl_select_df_VP,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean VP")

ggbetweenstats(
  data = vpcl_select_df_VP,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE VP")

vpcl_select_df_VPC<- nj_vp %>%filter(Population == 'c_Viruses', 
                                    VP_Type == c('VPCL_AR_Diff_SE',
                                                 'VPCL_AR_Diff_LMER_No_SE',
                                                 'VPCL_AR_Diff_No_SE',
                                                 'VPCL_AR_Diff_LMER_SE'),
                                    Sample_Type == 'VPC')

ggbetweenstats(
  data = vpcl_select_df_VPC,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean VPC")

ggbetweenstats(
  data = vpcl_select_df_VPC,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE VPC")

vpcl_select_df_Diff<- nj_vp %>%filter(Population == 'c_Viruses', 
                                     VP_Type == c('VPCL_AR_Diff_SE',
                                                  'VPCL_AR_Diff_LMER_No_SE',
                                                  'VPCL_AR_Diff_No_SE',
                                                  'VPCL_AR_Diff_LMER_SE'),
                                     Sample_Type == 'Diff')

ggbetweenstats(
  data = vpcl_select_df_Diff,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL Mean Diff")

ggbetweenstats(
  data = vpcl_select_df_Diff,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: VPCL SE Diff")


#do a simple t - test
tt_df_vpcl<- vpcl_select_df_Diff %>% select(-c(VP_SE, VP_R_Squared)) %>% pivot_wider(names_from = 'VP_Type', values_from = 'VP')%>%
  replace_na(list(VPCL_AR_Diff_No_SE = 0,
                  VPCL_AR_Diff_SE = 0,
                  VPCL_AR_Diff_LMER_No_SE = 0,
                  VPCL_AR_Diff_LMER_SE = 0))

t.test(tt_df_vpcl$VPCL_AR_Diff_SE , tt_df_vpcl$VPCL_AR_Diff_No_SE) #there is NO differenc ebtween the two
t.test(tt_df_vpcl$VPCL_AR_Diff_LMER_SE , tt_df_vpcl$VPCL_AR_Diff_LMER_No_SE)  #there is NO differenc ebtween the two
t.test(tt_df_vpcl$VPCL_AR_Diff_LMER_SE , tt_df_vpcl$VPCL_AR_Diff_SE) #there is NO differenc ebtween the two
t.test(tt_df_vpcl$VPCL_AR_Diff_No_SE , tt_df_vpcl$VPCL_AR_Diff_LMER_No_SE) #there is NO differenc ebtween the two

wilcox.test(tt_df_vpcl$VPCL_AR_Diff_SE , tt_df_vpcl$VPCL_AR_Diff_No_SE)
wilcox.test(tt_df_vpcl$VPCL_AR_Diff_LMER_SE , tt_df_vpcl$VPCL_AR_Diff_LMER_No_SE)  #there is NO differenc ebtween the two
wilcox.test(tt_df_vpcl$VPCL_AR_Diff_LMER_SE , tt_df_vpcl$VPCL_AR_Diff_SE) #there is NO differenc ebtween the two
wilcox.test(tt_df_vpcl$VPCL_AR_Diff_No_SE , tt_df_vpcl$VPCL_AR_Diff_LMER_No_SE)
wilcox.test(tt_df_vpcl$VPCL_AR_Diff_LMER_SE , tt_df_vpcl$VPCL_AR_Diff_No_SE)


#TRYING ROGME
library(devtools)
devtools::install_github("GRousselet/rogme")
library(rogme)


vpcl_select_df_Diff2 <- vpcl_select_df_Diff %>% filter(VP_Type == 'VPCL_AR_Diff_LMER_SE' | VP_Type == 'VPCL_AR_Diff_No_SE') %>% mutate(VP_Type = as.factor(VP_Type))


vp_df<- vpcl_select_df %>% filter(VP_Type == 'VPCL_AR_Diff_LMER_SE' | VP_Type == 'VPCL_AR_Diff_No_SE')%>% mutate(VP_Type = factor(VP_Type,
                                                                                                                 levels = c('VPCL_AR_Diff_No_SE')))
  
vp_df<- nj_vp %>% filter(Population == 'c_Viruses') %>% 
  filter(VP_Type == 'VPCL_AR_Diff_LMER_SE' | VP_Type == 'VPCL_AR_Diff_No_SE') %>%
  select(-c(VP_SE, VP_R_Squared)) %>%
  mutate(VP_Type = factor(VP_Type)) %>%
  pivot_wider(names_from = VP_Type, values_from = VP)%>%
  
  
                                                                                                                                                                                                                                                      'VPCL_AR_Diff_LMER_SE')))
ps <- plot_scat2(vp_df,
                 formula = VP ~ VP_Type,
                 alpha = 1,
                 shape = 21,
                 colour = "grey10",
                 fill = "grey90")
ps <- ps + coord_flip()
ps <- ps + labs(title = "No clear difference") +
  theme(plot.title = element_text(face = "bold", size = 20))
ps

sf<- shifthd(vp_df, formula =VP ~ VP_Type, nboot = 10)

sf <- plot_sf(sf, plot_theme = 2)
sf




p <- plot_scat2(vp_df,
               
                alpha = .3,
                shape = 21,
                colour = "grey10",
                fill = "grey90")

n = 300 #> number of observations per group
set.seed(6)
g1 <- rnorm(n)
g2 <- rnorm(n)

ks(g1,g2)

t.test(g1,g2)

#> make tibble
df <- mkt2(g1,g2)

#> -------------------------------------------------------
#> scatterplots alone
ps <- plot_scat2(df,
                 xlabel = "",
                 ylabel = "Scores (a.u.)",
                 alpha = 1,
                 shape = 21,
                 colour = "grey10",
                 fill = "grey90")
ps <- ps + coord_flip()
ps <- ps + labs(title = "No clear difference") +
  theme(plot.title = element_text(face = "bold", size = 20))
 ps

 sf <- shifthd(data = df, formula = obs ~ gr, nboot = 200)
 sf
 psf <- plot_sf(sf, plot_theme = 2) 
psf 
psf[[1]] <- psf[[1]] +
  labs(x = "Group 1 quantiles of scores (a.u.)",
       y = "Group 1 - group 2 \nquantile differences (a.u.)")

#> add labels for deciles 1 & 9
psf <- add_sf_lab(psf, sf, y_lab_nudge = .1, labres = 1, text_size = 4)
#> psf

#> -------------------------------------------------------
#> scatterplot + deciles + colour coded decile differences
p <- plot_scat2(df,
                xlabel = "",
                ylabel = "Scores (a.u.)",
                alpha = .3,
                shape = 21,
                colour = "grey10",
                fill = "grey90")
p
p <- plot_hd_links(p, sf[[1]],
                   q_size = 1,
                   md_size = 1.5,
                   add_rect = TRUE,
                   rect_alpha = 0.1,
                   rect_col = "grey50",
                   add_lab = TRUE) #> superimposed deciles + rectangle
p <- p + coord_flip() #> flip axes
 p
 
 cowplot::plot_grid(ps, p, psf[[1]], labels=c("A", "B", "C"), ncol = 1, nrow = 3,
                    rel_heights = c(1, 1, 1), 
                    label_size = 20, 
                    hjust = -0.5, 
                    scale=.95,
                    align = "v")
 
 
 
 
 #Independent groups
 vp_df<- nj_vp %>% filter(Population == 'c_Viruses') %>% 
   filter(VP_Type == 'VPCL_AR_Diff_LMER_SE' | VP_Type == 'VPCL_AR_Diff_No_SE') %>%
   select(-c(VP_SE, VP_R_Squared)) %>%
   mutate(VP_Type = factor(VP_Type)) %>%
   pivot_wider(names_from = VP_Type, values_from = VP)
 
 vp_df2<- vp_df %>% pivot_longer(cols = c('VPCL_AR_Diff_No_SE',
                                          'VPCL_AR_Diff_LMER_SE'),
                                 names_to = 'VP_Type',
                                 values_to = 'VP') %>%
   mutate(VP_Type = as.factor(VP_Type))
 
p<-  plot_scat2(vp_df2,
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
  

q<- ggplot(data = vp_df2, aes(fill = VP_Type, x = VP, col = VP_Type))+
  geom_density(alpha = 0)+
 facet_grid(rows = vp_df2$VP_Type)+
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

sf<- shifthd_pbci(data = vp_df2, formula = VP ~ VP_Type)
psf <- plot_sf(sf, plot_theme = 2)[[1]] 
psf
psf+xlim(c(-1, 1))

dasym <- asymhd(data = vp_df2, formula = VP ~ VP_Type, 
                q = seq(5,40,5)/100, alpha = .05, nboot = 100)

#> ggplot
diff_asym <- plot_diff_asym(data = dasym)[[1]]
diff_asym +
  xlim(c(-1, 0.5))


#https://garstats.wordpress.com/2019/02/21/hsf/

