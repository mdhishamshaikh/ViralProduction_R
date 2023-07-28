library(tidyverse)
library(ggstatsplot)
library(ggsci)
library(ggstatsplot)
library(cowplot)
#install.packages(parallel)
unique(output2$VP_Type)

make_a_csv <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/simulations/simu_output_df.csv")


#### 1. Linear Model Comparison ####
lm_vp_type <- c('LM_AP', 'LM_SR_AVG',
                'LM_AR', 'LM_AR_Diff',
                'LM_AR_Diff_LMER')

lm_output<- make_a_csv %>%
  filter(VP_Type  %in% c('LM_AP', 'LM_SR_AVG',
                         'LM_AR', 'LM_AR_Diff',
                         'LM_AR_Diff_LMER'))
unique(lm_output$VP_Type)

vpcl_output<- make_a_csv %>%
  filter(!(VP_Type  %in% c('LM_AP', 'LM_SR_AVG',
                           'LM_AR', 'LM_AR_Diff',
                           'LM_AR_Diff_LMER')))

unique(vpcl_output$VP_Type)
        
ggplot(data = lm_output,
       aes(x = VP_Type,
           y = VP,
           fill = VP_Type))+
  geom_violin()+
  geom_boxplot()

lm_comp<- ggbetweenstats(
    data = lm_output,
    x = VP_Type,
    y = VP,
    type = "nonparametric", # ANOVA or Kruskal-Wallis
    plot.type = "violin",
    pairwise.comparisons = TRUE,
    pairwise.display = "s",
    centrality.plotting = FALSE,
    bf.message = FALSE
  )+labs(x = 'VP Data Analyses Type',
         y = 'Viral Production Rate',
    title = "Kruskal-Wallis Test: LM VP")


 ggbetweenstats(
   data = lm_output,
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   ggsignif.args = list(textsize = 2, tip_length = 0.01, na.rm = TRUE, angle =90),
   pairwise.display = 'ns',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = lm_output,
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_lancet()+
   coord_flip()
 
 

a<-  ggbetweenstats(
   data = lm_output,
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = lm_output,
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 
 b<- ggbetweenstats(
   data = lm_output,
   x = VP_Type,
   y = VP_SE,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = lm_output,
               aes(x = VP_Type,
                   y = VP_SE,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 extract_stats(a)
 extract_stats(b)
 extract_subtitle(a)
 
 
 combine_plots(list(a,b))
 
 
 
 c<-  ggbetweenstats(
   data = vpcl_output,
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = vpcl_output,
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 
 d<- ggbetweenstats(
   data = vpcl_output,
   x = VP_Type,
   y = VP_SE,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 1, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 's',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data = vpcl_output,
               aes(x = VP_Type,
                   y = VP_SE,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.1)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 extract_stats(c)
 extract_stats(d)
 extract_subtitle(c)
 
 
 combine_plots(list(c,d))
 
 
 
 e<-  ggbetweenstats(
   data = vpcl_output %>%
     filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                           'VPCL_AR_Diff_LMER_SE')),
   x = VP_Type,
   y = VP,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 'ns',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data =  vpcl_output %>%
                 filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                       'VPCL_AR_Diff_LMER_SE')),
               aes(x = VP_Type,
                   y = VP,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.05)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 
 f<- ggbetweenstats(
   data =  vpcl_output %>%
     filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                           'VPCL_AR_Diff_LMER_SE')),
   x = VP_Type,
   y = VP_SE,
   violin.args = list(linewidth = 1,
                      fill = NA),
   type = 'nonparametric',
   plot.type = 'violin',
   point.args = list(alpha= 1, position = position_jitterdodge(jitter.width = 1)),
   pairwise.display = 'ns',
   centrality.type = 'parametric',
   centrality.plotting = F,
   centrality.label.args = list(alpha = 0),
   centrality.point.args = list(color = 'black',
                                size = 3)
   
 )+
   geom_violin(data =  vpcl_output %>%
                 filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                       'VPCL_AR_Diff_LMER_SE')),
               aes(x = VP_Type,
                   y = VP_SE,
                   fill = VP_Type))+
   geom_boxplot(size = 1,
                outlier.alpha = 0,
                width = 0.05)+
   scale_fill_npg()+
   theme_classic()+
   theme(legend.position = 'none')
 
 extract_stats(e)
 extract_stats(f)
 extract_subtitle(e)
 
 
 combine_plots(list(e,f))
 
 
 
 ggbetweenstats(data =  vpcl_output %>%
                  filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                        'VPCL_AR_Diff_LMER_SE')),
                x= VP_Type,
                y = VP,
                pairwise.display = 's',
                plot.type = 'violin',
                type = 'non-parametric',
                point.args = list(alpha= 0))+
   geom_point(data =  vpcl_output %>%
                filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                      'VPCL_AR_Diff_LMER_SE')),
             aes( x= VP_Type,
              y = VP,
              color = VP_Type),
              position = position_jitterdodge())
   
   
   ggplot(data =  vpcl_output %>%
            filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                  'VPCL_AR_Diff_LMER_SE')),
          x= VP_Type,
          y = VP)+
   geom_beeswarm()
 
 
   #VIPCAL vs VIPCAL-SE violin plot #####
   
   
   simu_output_df <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/simulations/simu_output_df.csv")
   
   vpcl_se_df<- simu_output_df %>%
     filter(VP_Type %in% c('VPCL_AR_Diff_LMER_SE', 'VPCL_AR_Diff_No_SE'))
   vpcl_se_df$VP_Type<- factor(vpcl_se_df$VP_Type, 
                               levels = c('VPCL_AR_Diff_No_SE', 'VPCL_AR_Diff_LMER_SE'),
                               labels = c('VIPCAL', 'VIPCAL-SE'))
   
   ggplot(data = vpcl_se_df,
          aes(x = VP_Type,
          y = VP,
          #fill = VP_Type,
          color = VP_Type))+
     geom_point(position = position_jitterdodge(dodge = 0.1),
                alpha = 0.2)+
    # geom_violin()+
     scale_color_manual(values = c("178b76ff", "#52495aff"))
 
 #
 
 simu_output_df <- read.csv("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/simulations/simu_output_df.csv")
 
 simu_output_df<- simu_output_df %>%
   filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                         'VPCL_AR_Diff_LMER_SE',
                         'LM_AP'))
simu_output_df$VP_Type<-   factor(simu_output_df$VP_Type, levels = c('VPCL_AR_Diff_No_SE',
                     'VPCL_AR_Diff_LMER_SE',
                     'LM_AP'))
 
 ggbetweenstats(data =  simu_output_df ,
                x= VP_Type,
                y = VP,
                pairwise.display = 's')
 
 
 
 ggbetweenstats(data =  vpcl_output %>%
                  filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                                        'VPCL_AR_Diff_LMER_SE')),
                x= VP_Type,
                y = VP,
                PA)
 
 
 
 
 geom_boxplot()+
   scale_fill_manual(values = NULL)
 
 ggplot(data = lm_output,
           aes(x = VP_Type,
               y = VP,
               color = Sample_Type,
               shape = Sample_Type))+
   ggbeeswarm::geom_beeswarm(cex = 0.3,
                             alpha = 0.2)
   
 
 geom_point(stat = 'identity',
              position = position_jitterdodge(),
              alpha = 0.1)+
   geom_violin(aes(fill = NULL),
               draw_quantiles = c(0.25, 0.5, 0.75))
 
 
 df6<- simu_output_df %>% 
   select(-c(VP_SE, VP_R_Squared, n))%>%group_by(Expt_No, Sample_Type)%>%
   filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                         'VPCL_AR_Diff_LMER_SE'))
 df6$Sample_Type<- factor(df6$Sample_Type,
                             levels = c("VP",
                                        "VPC",
                                        "Diff"))
 df7<- simu_output_df %>% 
   select(-c(VP_SE, VP_R_Squared, n))%>%group_by(Expt_No, Sample_Type)%>%
   filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                         'VPCL_AR_Diff_LMER_SE',
                         'LM_AP')) %>% pivot_wider(names_from = VP_Type, values_from = VP)
 df7$Sample_Type<- factor(df7$Sample_Type,
                          levels = c("VP",
                                     "VPC",
                                     "Diff"))
 lm_vp_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7 %>% filter(Sample_Type == 'VP')))
 lm_vp_fit$coefficients[[2]]
 
 lm_vpc_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7 %>% filter(Sample_Type == 'VPC')))
 lm_vpc_fit$coefficients[[2]]
 
 lm_diff_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7 %>% filter(Sample_Type == 'Diff')))
 lm_diff_fit$coefficients
 
 lm_fit<- summary(lm(VPCL_AR_Diff_No_SE ~ VPCL_AR_Diff_LMER_SE, data = df7))
 lm_fit$coefficients[[2]]
 
 
 
 
 vp_se_comp_legend<- ggplot(data = df7, aes(x = VPCL_AR_Diff_No_SE,
                                     y = VPCL_AR_Diff_LMER_SE,
                                     color = Sample_Type,
                                     fill = Sample_Type,
                                     shape = Sample_Type))+
   geom_point(alpha = 0.2)+
   theme_classic()+
   geom_abline(slope = 1, linewidth = 1)+
   scale_color_npg()+
   labs(x = 'VIPCAL',
        y = 'VIPCAL-SE')+
  geom_smooth(method = 'lm', linewidth =2, alpha = 0.2)+
   scale_color_manual(values = cols)+
   scale_fill_manual(values = cols)+
   scale_shape_manual(values = shape)+
   theme(axis.line = element_line(linewidth = 1),
         axis.text = element_text(face = 'bold', size = 15),
         axis.title = element_text(face = 'bold', size = 15),
         legend.box.background = element_rect(color = 'black'),
         legend.box.margin = margin(t = 1, l = 1),
         legend.title = element_text(face = 'bold'),
         legend.position = c(0.85, 0.25))+
   # xlim(c(-0.1, 2.2))+
   # ylim(c(-0.1, 2.2))+
   labs(fill = 'Treatment',
        color = 'Treatment',
        shape = 'Treatment')
 
 
 vp_se_comp_legend
 
 
vpcl_se_lm_plot<- vp_se_comp_legend+
   ggpubr::stat_regline_equation(show.legend = F)
 
vpcl_se_lm_plot

 
 #adding linear regression equations on the plot
 
 lm_eqn <- function(m){
   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                    list(a = format(unname(coef(m)[1]), digits = 3),
                         b = format(unname(coef(m)[2]), digits = 3),
                         r2 = format(summary(m)$r.squared, digits = 3)))
   eq<- as.expression(eq)
   
 }
 
 vp_se_comp_legend + geom_text(x = 25, y = 300, label = lm_eqn(vp_lm), parse = TRUE)
 

 
lm2<- function(df, x, z){
  m<- lm(data = df, z ~x)
  return(m)
}

vp_lm<- lm(VPCL_AR_Diff_LMER_SE ~ VPCL_AR_Diff_No_SE, subset(df7, Sample_Type == 'VP'))
vpc_lm<- lm(VPCL_AR_Diff_LMER_SE ~ VPCL_AR_Diff_No_SE, subset(df7, Sample_Type == 'VPC'))
diff_lm<- lm(VPCL_AR_Diff_LMER_SE ~ VPCL_AR_Diff_No_SE, subset(df7, Sample_Type == 'Diff'))



lm2(df7, x = 'VPCL_AR_Diff_LMER_SE', z = 'VPCL_AR_Diff_No_SE')





#split violins 
violin_df<- simu_output_df%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE')
violin_df$VP_Type2<- 'A'



ggplot(df6,
       aes(y = VP,
           x = Sample_Type,
           fill = VP_Type))+
  #geom_point()+
  geom_split_violin(nudge = 0.02,
                    color = 'transparent')+
  theme_classic()+
  scale_fill_manual(values = c("178b76ff", "#52495aff"),
                    labels = c('VIPCAL-SE', 'VIPCAL'))+
  scale_color_manual(values<- NA)+
  theme(axis.line = element_line(linewidth = 1),
        axis.text = element_text(face = 'bold', size = 15),
        axis.title = element_text(face = 'bold', size = 15),
        legend.box.background = element_rect(color = 'black'),
        legend.box.margin = margin(t = 1, l = 1),
        legend.title = element_text(face = 'bold'),
        legend.position = 'right')+
  labs(x = 'Treatment',
       y = 'Viral Production',
       fill = 'VP Type') #PLOT: VIPCAL VIPCAL SE split violin #####



#Perfeom Kruskal Wallis test between VIPCAL and VIPCAL-SE for VP, VPC and Diff 
v_df <- make_a_csv
v_df$VP_Type<- factor(v_df$VP_Type, levels = c('LM_AP', 'LM_SR_AVG',
                                             'LM_AR', 'LM_AR_Diff',
                                             'LM_AR_Diff_LMER', "VPCL_SR_AVG" ,
                                             "VPCL_AR_No_SE" , "VPCL_AR_SE" ,
                                             "VPCL_AR_Diff_No_SE" ,
                                             "VPCL_AR_Diff_SE"  ,
                                             "VPCL_AR_Diff_LMER_No_SE", "VPCL_AR_Diff_LMER_SE"
),
labels = c('LM-1', 'LM-2', 'LM-3', 'LM-4', 'LM-5',
           'VPCL-1', 'VPCL-2', 'VPCL-3',
           'VPCL-4', 'VPCL-5', 'VPCL-6',
           'VPCL-7'))


kruskal.test(VP ~ VP_Type, data = v_df)

library(FSA)
dunn_test<- dunnTest(VP ~ VP_Type, data = v_df,
         method = 'holm')
dunn_output<- as.data.frame(as.list(dunn_test)[[2]])

dunn2<- dunn_output

colnames(dunn2)



dunn2<- dunn2 %>% separate(Comparison, c("Comp1", "Comp2"), " - ")



dunn2<- dunn2 %>%
  arrange(Comp1, Comp2)

corr_p_plot<- ggplot(data = dunn2 , aes(x = Comp1, y = Comp2, fill = P.adj, color = P.adj))+
  geom_tile()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 10,
                                 face = 'bold'),
        axis.title = element_blank(),
        legend.position = c(0.8, 0.3),
        legend.title = element_text(size = 8,
                                    face = 'bold'),
        legend.text = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
        panel.grid = element_blank())+
  labs(fill = 'P(adjusted)',
       color = 'P(adjusted)')+
  scale_fill_gradient(high = '#ea7e82ff', low = '#52495aff')+
  scale_color_gradient(high = '#ea7e82ff', low = '#52495aff')

corr_p_plot



dunn2<- dunn2 %>%
  select(-c(Z, P.unadj)) %>%
  pivot_wider(names_from = Comp2,
              values_from = P.adj) %>%
  column_to_rownames(var = 'Comp1')%>%
  as.matrix()

ggcorrplot::ggcorrplot(dunn2)


#Complete set of 12 methods####




v_df <- make_a_csv
v_df$VP_Type<- factor(v_df$VP_Type, levels = c('LM_AP', 'LM_SR_AVG',
                          'LM_AR', 'LM_AR_Diff',
                          'LM_AR_Diff_LMER', "VPCL_SR_AVG" ,
                          "VPCL_AR_No_SE" , "VPCL_AR_SE" ,
                          "VPCL_AR_Diff_No_SE" ,
                          "VPCL_AR_Diff_SE"  ,
                          "VPCL_AR_Diff_LMER_No_SE", "VPCL_AR_Diff_LMER_SE"
                        ),
                      labels = c('LM-1', 'LM-2', 'LM-3', 'LM-4', 'LM-5',
                                 'VPCL-1', 'VPCL-2', 'VPCL-3',
                                 'VPCL-4', 'VPCL-5', 'VPCL-6',
                      'VPCL-7'))

make_a_csv2$VP_Type<- factor(make_a_csv2$VP_Type, levels = c('LM_AP', 'LM_SR_AVG',
                                               'LM_AR', 'LM_AR_Diff',
                                               'LM_AR_Diff_LMER', "VPCL_SR_AVG" ,
                                               "VPCL_AR_No_SE" , "VPCL_AR_SE" ,
                                               "VPCL_AR_Diff_No_SE" ,
                                               "VPCL_AR_Diff_SE"  ,
                                               "VPCL_AR_Diff_LMER_No_SE", "VPCL_AR_Diff_LMER_SE"
),
labels = c('LM-1', 'LM-2', 'LM-3', 'LM-4', 'LM-5',
           'VPCL-1', 'VPCL-2', 'VPCL-3',
           'VPCL-4', 'VPCL-5', 'VPCL-6',
           'VPCL-7'))

ggplot(data = v_df, aes(x = VP_Type, y = VP))+
  geom_violin()+
  geom_point(aes(color = Sample_Type),
             alpha = 0.1,
             position = position_jitter()) +
  scale_color_npg()


ggplot(data = v_df, aes(x = VP_Type, y = VP, fill = VP_Type, color = VP_Type))+
  geom_violin()+
  geom_boxplot(width = 0.1, fill = NA, color = 'white', linewidth = 1, outlier.alpha = 0)+
#  scale_fill_npg()+
#  coord_flip()+
  scale_x_discrete(limits = rev(levels(v_df$VP_Type)))+
  theme_classic()




# c('LM_AP', 'LM_SR_AVG',
#   'LM_AR', 'LM_AR_Diff',
#   'LM_AR_Diff_LMER', "VPCL_SR_AVG" ,
#   "VPCL_AR_No_SE" , "VPCL_AR_SE" ,
#   "VPCL_AR_Diff_No_SE" ,
#   "VPCL_AR_Diff_SE"  ,
#   "VPCL_AR_Diff_LMER_No_SE", "VPCL_AR_Diff_LMER_SE"
# )



#####All combined

  cols<- c(as.vector(paletteer::paletteer_d(palette = 'ggsci::nrc_npg')), as.vector(paletteer::paletteer_d(palette = 'ggsci::default_aaas')))
cols_wes<- c('#e19e34ff', '#e3782fff', '#bc4b2fff',
         '#4888a2ff','#433d67ff', '#52495aff', 
         '#405529ff', '#972029ff', '#86b0a5ff',
         '#178b76ff', '#ea7e82ff', '#df536bff')
  
  
  
  
vp_methods_simu<-   ggbetweenstats(
  data = v_df,
  x = VP_Type,
  y = VP,
  violin.args = list(linewidth = 0.5,
                     fill = NA, colour = NA),
  type = 'nonparametric',
  plot.type = 'violin',
  pairwise.comparisons = F,
  point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
  ggsignif.args = list(textsize = 2, vjust = 0.4, tip_length = 0.01, na.rm = TRUE, alpha = 0),
  pairwise.display = FALSE,
  centrality.type = 'parametric',
  centrality.plotting = F,
  centrality.label.args = list(alpha = 0),
  centrality.point.args = list(color = 'black',
                               size = 3),
  results.subtitle = F
  
)+
    
  geom_violin(data = v_df,
              aes(x = VP_Type,
                  y = VP,
                  fill = VP_Type),
              colour = NA,
              scale = 'count')+
  geom_boxplot(size = 0.8,
               outlier.alpha = 0,
               width = 0.05,
               fill = 'white',
               color = 'black')+
    geom_hline(yintercept = 0, color = 'black', linewidth  =1 , alpha = 0.2)+
  scale_fill_manual(values = cols_wes)+
  theme_classic()+
    labs(x = 'Data Analyses Method',
         y = 'Viral Production')+
  #  ylim(-2,2)+
  theme(legend.position = 'none',
        axis.text = element_text(size = 15,
                                 face = 'bold'),
        #axis.text.x = element_text(face = 'bold'),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.title = element_text(size = 12,
        #                           face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold',
                                    size = 15),
        axis.line = element_line(linewidth = 1),
        axis.line.y.right = element_blank())

  vp_methods_simu
  
  
vp_se_methods_simu<-  ggbetweenstats(
    data = make_a_csv,
    x = VP_Type,
    y = VP_SE,
    violin.args = list(linewidth = 1,
                       fill = NA, colour = NA),
    type = 'nonparametric',
    plot.type = 'violin',
    pairwise.comparisons = F,
    point.args = list(alpha= 0, position = position_jitterdodge(jitter.width = 1)),
    ggsignif.args = list(textsize = 2, vjust = 0.4, tip_length = 0.01, na.rm = TRUE, alpha = 0),
    pairwise.display = FALSE,
    centrality.type = 'parametric',
    centrality.plotting = F,
    centrality.label.args = list(alpha = 0),
    centrality.point.args = list(color = 'black',
                                 size = 3),
    results.subtitle = F
    
  )+
    geom_violin(data = make_a_csv,
                aes(x = VP_Type,
                    y = VP_SE,
                    fill = VP_Type),
                colour = NA,
                scale = 'width')+
    geom_boxplot(size = 1,
                 outlier.alpha = 0,
                 width = 0.05,
                 fill = 'white',
                 color = 'black')+
    geom_hline(yintercept = 0, color = 'black', linewidth  =1 , alpha = 0.2)+
    scale_fill_manual(values = cols_wes)+
    theme_classic()+
    labs(x = 'Data Analyses Method',
         y = 'Viral Production')+
   ylim(-2,2)+
    theme(legend.position = 'none',
          axis.text = element_text(size = 10),
          axis.text.x = element_text(face = 'bold'),
          axis.title = element_text(size = 12,
                                    face = 'bold'),
          axis.line = element_line(linewidth = 1),
          axis.line.y.right = element_blank())
  
vp_se_methods_simu
  
  
  #Making a table ####
  library(flextable)
  
#Lets make a datafrmae first
  
  table_df <- data.frame(name = c('A', 'B', 'C'),
                         Linear = c(1,0,0),
                         VIPCAL = c(0,1,1),
                         No_Rep = c(0,0,0),
                         Sep_rep = c(1,0,0),
                         AAvg_rep = c(0,1,1),
                         no_SE = c(1,0,0),
                         yes_SE = c(0,1,1),
                         no_diff = c(1,0,0),
                         diff_Sub = c(0,1,0),
                         diff_lmer = c(0,0,1))
  
 
ft  

ft<- add_header_row(ft, values = c('', 'Linear', 'VIPCAL',
                                   'No', 'Yes', 'No', 'Yes',
                                   'No', 'Yes'),
                    colwidths = c(1, 1, 1, 1, 2, 1, 1, 1, 2))
ft<- add_header_row(ft, values = c('Name', 'Type', 'Replicates',
                                   'SE', 'Diff_Curve'),
                    colwidths = c(1,2, 3,2,3))
ft

library(officer)
bigger_border = fp_border(color = 'black', width = 2)
big_border = fp_border(color = 'black', width = 1.5)
small_border = fp_border(color = 'gray', width = 1.5)
smaller_border = fp_border(color = 'gray', width = 1)

table_df <- readxl::read_excel("simulations/table/table_simulation.xlsx")
ft<- flextable(table_df) %>%
  #theme_box() %>%
  merge_v(j = ~ A1 + A2) %>%
  merge_h_range(i=c(1,2, 6, 7, 8), j1 = 2, j2 = 3)%>% 
  #rotate(j = 1,rotation = 'btlr', part = 'body', align = 'top') %>%
  #void(j = 1:3, part = 'header') %>%
  bg(j = 4:15, i =1, part = 'header', bg = cols) %>%
  color(j = 4:15, part = 'header', color = 'white') %>%
  align(part = 'header', align = 'center') %>%
  font(fontname = 'Arial')%>%
  fontsize(size = 10, part = 'all') %>%
  fontsize(size = 12, j = 4:15,  part = 'body') %>%
  #bold(j = 1:3)%>%
  bold(part = 'all') %>%
  height_all(height = 0.2, unit = 'cm',)%>%
  width(width = 2, unit = 'cm', j = 1:15) %>%
  # 
  # add_header_row(values = c("", ""),
  #                colwidths = c(3,12),
  #                top = F) %>%
  border_remove() %>%
  hline_top( part = 'all', border = bigger_border) %>%
  hline_bottom(part = 'all', j= 1:15, border = bigger_border) %>%
  # hline_top( part = 'header', border = big_border) %>%
  # hline(part = 'header', border = big_border) %>%
  hline(i = c(2,5,7, 10), j = 1:15, border = big_border) %>%
  hline(i = c(1,3,6, 8), j = 2:15, border = small_border) %>%
  hline(i = c(4, 9), j = 3:15, border = smaller_border) %>%
  padding(padding = 4, i = 1, part = 'header') %>%
  vline(j = 3, border = bigger_border) %>%
  fix_border_issues() %>%
  labelizor(part = 'all',
            labels = c("0" = "",
                       "1" = "✓",
                       "A1" = "Viral Production Analyses Methods",
                       "A2" = "Viral Production Analyses Methods",
                       "A3" = "Viral Production Analyses Methods"))%>%
  merge_h(i = 1, part = 'header')%>%
  align(j = 4:15, part = 'body', align = 'center') %>%
  valign(j = 4:15, valign = 'center', part = 'body') %>%
  line_spacing(j =1:3, part = 'all', space = 1.5)

ft
plot(ft, fit = 'fixed', just = 'center')

save_as_image(ft, path = 'ft_table.png', res = 800)

gr <- gen_grob(ft, fit = "fixed", just = "center")
dims <- dim(gr)
dims
svglite::svglite(filename = "ft_table.svg", width = dims$width + .1, height = dims$height + .1)
plot(gr)
dev.off()

ft_gg<- gen_grob(ft)
plot(ft_gg)



####Combining all the plots####

mock_df_plots
corr_p_plot
vp_methods_simu

plot_grid(NULL, corr_p_plot, NULL, )



plot_grid(mock_df_plots, plot_grid(corr_p_plot, vp_methods_simu, 
                                   nrow = 1, 
                                   rel_widths =  c(1,3),
                                   labels = c("B", "C"),
                                   label_y  = 1.1),
          ncol = 1,
          rel_heights = c(1,1),
          labels = c("A", ""))

plot_grid(plot_grid(mock_df_plots, plot_grid(corr_p_plot, vp_methods_simu, 
          nrow = 1, 
          rel_widths =  c(1,3)),
          ncol = 1,
          rel_heights = c(1,1)), ft_gg,
          nrow = 2, rel_heights = c(2,1))

  

library(patchwork)

mock_df_plots + plot_spacer()+  corr_p_plot + vp_methods_simu+
  plot_layout(widths = c(2,1), heights = c(1, 1),
              nrow =2)

