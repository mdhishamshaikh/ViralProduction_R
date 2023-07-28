library(tidyverse)
library(ggsci)
library(cowplot)
library(ggpubr)
library(patchwork)
library(colorspace)


bep_df <- read_csv("./V5000/vp_calc_bp.csv") #This dataframe only has the endpoint time ranges
summary
#We'll add a column to flag this data as the endpoint data
bep_df$Time_Range_Type <- "BEP"
#We will add T0-T24 time range to this dataframe

all_df <- read_csv("./V5000/vp_calc_all.csv")


bep_df<- rbind(bep_df, all_df %>% filter(Time_Range == 'T0_T24')
               %>% mutate(Time_Range_Type = 'T0_T24'))
unique(bep_df$Time_Range_Type)
unique(bep_df$VP_Type)


#Visualisation
vpcl_se_df <- bep_df %>% filter(VP_Type == 'VPCL_AR_Diff_LMER_SE',
                                Population == 'c_Viruses',
                                Sample_Type == 'VP')
ggplot(data = vpcl_se_df,
       aes(x = as.factor(Expt_No),
           y = VP,
           fill = Time_Range_Type,
           #group = Time_Range_Type,
           #color = Time_Range_Type,
           #shape = Time_Range_Type
           ))+
 #geom_point(alpha = 1, size = 3)+
  geom_bar(data = vpcl_se_df,
           aes(x = as.factor(Expt_No),
               y = VP,
               fill = Time_Range_Type,
               group = Time_Range_Type,
               #color = Time_Range_Type,
               #shape = Time_Range_Type
           ),
           stat = 'identity',
           width = 0.125,
           position = position_dodge(width = 0.5),
           color = NA,
           fill = 'black')+
  theme_bw()+
  geom_point(data = vpcl_se_df %>% filter(Time_Range_Type == 'BEP'),
             aes(x = as.factor(Expt_No),
                 y = VP,
                 fill = Time_Range_Type),
             position = position_nudge(x = -0.125),
             size = 3,
             shape = 19,
             color = 'coral')+
  geom_point(data = vpcl_se_df %>% filter(Time_Range_Type == 'T0_T24'),
             aes(x = as.factor(Expt_No),
                 y = VP,
                 fill = Time_Range_Type),
             position = position_nudge(x = 0.125),
             size = 3,
             shape = 17,
             color = 'red')
  
geom_segment( aes(x=as.factor(Expt_No), xend=as.factor(Expt_No), y=0, yend=VP))



#Geom_Segments ####


  

vpcl_se_df <- bep_df %>% filter(
                                Population == 'c_Viruses',
                                VP_Type %in% c("LM_AP",
                                             "VPCL_AR_Diff_No_SE",
                                             "VPCL_AR_Diff_LMER_SE")) %>%
  mutate(VP = replace_na(VP, 0))  %>%
  arrange(Location, Depth, Expt_No, Population, 
          Sample_Type, VP_Type, Time_Range_Type)
vpcl_se_df$Sample_Type <- factor(vpcl_se_df$Sample_Type,
                                 levels = c('VP', 'VPC', 'Diff'))

vpcl_se_df$VP_Type <- factor(vpcl_se_df$VP_Type,
                                 levels = c('LM_AP', 'VPCL_AR_Diff_No_SE', 'VPCL_AR_Diff_LMER_SE'),
                             label = c("LM", "VIPCAL", "VIPCAL-SE"))


bep_time_ranges<- vpcl_se_df %>% 
  filter(Time_Range_Type == 'BEP') %>%
  select(Expt_No, Time_Range) %>%
  unique()




gs_df1 <- vpcl_se_df %>%
  select(-c(VP_SE, VP_R_Squared, Location, Depth, Population)) %>%
  ungroup()%>%
  pivot_wider(id_cols = !Time_Range,
              names_from = Time_Range_Type,
              values_from = VP)%>%
  rowwise() %>%
  mutate(start = min(BEP, T0_T24),
         end = max(BEP, T0_T24),
         sign = sign(BEP-T0_T24)) 

  gs_df1<- merge(gs_df1, bep_time_ranges)

  

gs_df2<- vpcl_se_df %>%
  select(-c(VP_SE, VP_R_Squared)) %>%
  # mutate(Time_Range_Type2 = Time_Range_Type)%>%
  group_by(Location, 
           Depth, VP_Type,
           Population, Sample_Type) 




gs_df2$Time_Range<- factor(gs_df2$Time_Range, levels = c( "T0_T3" , "T0_T6",
                                                          "T0_T12", "T0_T17", "T0_T24"),
                           labels = c( "T0-T3" , "T0-T6",
                                       "T0-T12", "T0-T17", "T0-T24") )

gs_df1$Time_Range<- factor(gs_df1$Time_Range, levels = c( "T0_T3" , "T0_T6",
                                                          "T0_T12", "T0_T17", "T0_T24"),
                           labels = c( "T0-T3" , "T0-T6",
                                       "T0-T12", "T0-T17", "T0-T24") )


cols_points<- c('#4DBBD5FF', '#8491B4FF','#F39B7FFF', '#7E6148FF', '#3C5488FF')
cols_segments<- c('#E64B35FF', '#B09C85FF', '#00A087FF')
cols<-  c('-1' = '#E64B35FF', '0' = '#B09C85FF', '1'  ='#00A087FF',
          "T0_T3" = '#4DBBD5FF', "T0_T6" = '#8491B4FF', "T0_T12" = '#F39B7FFF', "T0_T17" = '#7E6148FF', "T0_T24" = '#3C5488FF')

bep_plot<- ggplot(data = gs_df2 ,
       aes(x = as.factor(Expt_No),
           y = VP))+
  geom_segment(data = gs_df1,
               aes(x = as.factor(Expt_No),
                   xend = as.factor(Expt_No),
                   y = start,
                   yend = end,
                   color = as.factor(sign)),
               linewidth = 1,
               show.legend = F)+
  geom_point(aes(color = Time_Range,
                 shape = Time_Range_Type),
             size = 3,
             show.legend = T) +
  scale_shape_manual(values = c(16, 17))+
  scale_color_manual( values = c('#F39B7FFF', '#7E6148FF', '#3C5488FF',
                                 '#00A087FF', '#B09C85FF', '#E64B35FF',
                                 '#4DBBD5FF', '#8491B4FF'))+
  theme_bw()+
  theme(#legend.position = 'none'
  )+
  facet_grid(Sample_Type ~ VP_Type,
             scales = 'free')+
  theme(strip.background = element_rect(fill = lighten("#323D5E",
                                                       amount = 0.0) ,
                                        color = NA),
        strip.text = element_text(face = 'bold',
                                  color = 'white',
                                  size = 10),
        #strip.background.x = element_rect(vjust = 3),
        legend.position = 'right',
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 7),
        panel.border = element_rect(linewidth =  2),
        panel.background = element_rect(fill = NA),
        #strip.text.x = element_text(vjust = 3),
        strip.clip = "on")+
  labs(x = 'Experiment',
       y = 'Viral Production') +
  guides(shape = guide_legend(title = 'Endpoints')) +
  theme(legend.position = 'none')



bep_plot


bep_plot1<- ggplot(data = gs_df2 ,
       aes(x = as.factor(Expt_No),
           y = VP))+
  geom_segment(data = gs_df1,
               aes(x = as.factor(Expt_No),
                   xend = as.factor(Expt_No),
                   y = start,
                   yend = end,
                   color = as.factor(sign)),
               linewidth = 1,
               show.legend = F)+
  geom_point(aes(color = Time_Range,
                 shape = Time_Range_Type),
             size = 3,
             show.legend = T) +
  scale_shape_manual(values = c(16, 17))+
  scale_color_manual( values = c('#F39B7FFF', '#7E6148FF', '#3C5488FF', '#E64B35FF', '#B09C85FF', '#00A087FF', '#4DBBD5FF', '#8491B4FF'))+
  theme_bw()+
  theme(#legend.position = 'none'
        )+
  facet_grid(Sample_Type ~ VP_Type,
             scales = 'free')+
  theme(strip.background = element_rect(fill = 'gray',
                                        color = NA),
        #strip.background.x = element_rect(vjust = 3),
        legend.position = 'right',
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        panel.border = element_rect(linewidth =  2),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 15,
                                  face = 'bold'),
      #strip.text.x = element_text(vjust = 3),
        strip.clip = "on")+
  labs(x = 'Experiment',
       y = 'Viral Production') +
  guides(shape = guide_legend(title = 'Endpoints')) +
  guides(color = 'none')


bep_plot1

bep_legend1<- get_legend(bep_plot1)
as_ggplot(get_legend(bep_plot1))






bep_plot2<- ggplot(data = gs_df2 ,
                   aes(x = as.factor(Expt_No),
                       y = VP))+
  # geom_segment(data = gs_df1,
  #              aes(x = as.factor(Expt_No),
  #                  xend = as.factor(Expt_No),
  #                  y = start,
  #                  yend = end,
  #                  color = as.factor(sign)),
  #              linewidth = 1,
  #              show.legend = F)+
  geom_point(aes(color = Time_Range,
                 shape = Time_Range_Type),
             size = 3,
             show.legend = T) +
  scale_shape_manual(values = c(19, 17))+
  scale_color_manual( values = cols_points)+
  theme_bw()+
  theme(#legend.position = 'none'
  )+
  facet_grid(Sample_Type ~ VP_Type,
             scales = 'free')+
  theme(strip.background = element_rect(fill = 'gray',
                                        color = NA),
        #strip.background.x = element_rect(vjust = 3),
        legend.position = 'right',
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 9),
        panel.border = element_rect(linewidth =  2),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 10,
                                  face = 'bold'),
        #strip.text.x = element_text(vjust = 3),
        strip.clip = "on")+
  labs(x = 'Experiment',
       y = 'Viral Production') +
  guides(color = guide_legend(title = 'Time Ranges')) +
  guides(
    shape = 'none'  
  )


bep_plot2

bep_legend2<- get_legend(bep_plot2)
as_ggplot(get_legend(bep_plot2))


gs_df1$sign <- factor( gs_df1$sign, levels = c(-1, 0 , 1),
                       label = c('T0-T24 > BEP',
                                 'T0-T24 = BEP',
                                 'T0-T24 < BEP'))

bep_plot3<- ggplot(data = gs_df2 ,
                   aes(x = as.factor(Expt_No),
                       y = VP))+
  geom_segment(data = gs_df1,
               aes(x = as.factor(Expt_No),
                   xend = as.factor(Expt_No),
                   y = start,
                   yend = end,
                   color = as.factor(sign)),
               linewidth = 1,
               show.legend =T)+
  # geom_point(aes(color = Time_Range,
  #                shape = Time_Range_Type),
  #            size = 3,
  #            show.legend = T) +
  # scale_shape_manual(values = c(19, 17))+
  scale_color_manual( values = cols_segments)+
 
  theme_bw()+
  facet_grid(Sample_Type ~ VP_Type,
             scales = 'free')+
  theme(strip.background = element_rect(fill = 'gray',
                                        color = NA),
        #strip.background.x = element_rect(vjust = 3),
        legend.position = 'right',
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15,
                                  face = 'bold'),
        legend.title = element_text(face = 'bold',
                                    size = 7),
      legend.text =  element_text(size = 9),
        panel.border = element_rect(linewidth =  2),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = 15,
                                  face = 'bold'),
        #strip.text.x = element_text(vjust = 3),
        strip.clip = "on")+
  labs(x = 'Experiment',
       y = 'Viral Production') +
  guides(color = guide_legend(title = '')) +
  guides(
    shape = 'none'  
  )


bep_plot3

bep_legend3<- get_legend(bep_plot3)
as_ggplot(get_legend(bep_plot3))


bep_legends<- cowplot::plot_grid(plotlist = list(bep_legend1, bep_legend2,  bep_legend3),
                             ncol= 1,
                             rel_widths = c(1, 1, 1),
                             axis = 'l')

bep_legends


bep_plot_final<- bep_plot + bep_legends +
  plot_layout(widths = c(8,2),
              ncol = 2)

bep_plot_final

bep_cr_plot<- bep_plot_final/cr_plot_final +
  plot_layout(ncol = 1, heights = c(4,2))+
  plot_annotation(tag_levels = c('A',
                                 'B'))

bep_cr_plot

#CR_PLOT_FINAL IS from 13_contatc_rates script
bep_cr_plot_legend<- plot_grid(NULL, bep_plot_final, NULL, cr_plot_final,
          nrow = 4,
          rel_heights = c(0.1, 2,0.1, 1),
          labels = c('A','', 'B', ''))


bep_cr_plot_legend

bep_cr_plot_legend<- plot_grid(bep_plot_final, NULL, cr_plot_final,
                               nrow = 1,
                               rel_widths = c(2,0.1, 1),
                               labels = c('A','', 'B'))


bep_cr_plot_legend


plot_grid(cr_plot, legends, bep_plot, bep_legends, 
          nrow = 1,
          rel_widths = c( 2, 1, 4, 1),
          labels = c( 'A', '', 'B', ''))









ggsave(file = 'bep_cr_NO_STRIP.png', width = 35,
       height = 20, unit = 'cm')


# 
# ggplot(data = gs_df1,
#        aes(x = BEP,
#            y = T0_T24,
#            color = as.factor(Expt_No),
#            shape = Time_Range))+
#   geom_point(size = 3,
#              alpha = 0.5)+
#   geom_abline(slope =1,
#               linewidth  =1)+
#   theme_bw()+
#   facet_grid(Sample_Type ~ VP_Type,
#              scales = 'fixed')+
#   theme(strip.background = element_rect(color = NA,
#                                         fill = NA),
#         legend.position = 'bottom')+
#   labs(x = 'Bacterial Endpoint',
#        y = 'T0-T24')+
#   guides(color = guide_legend('Experiment'),
#          shape = guide_legend('Endpoint Time Range'))
# 



