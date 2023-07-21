library(tidyverse)
library(ggsci)
library(cowplot)
library(grid)


df<- data.frame(Timepoint = rep(1:6,6),
                Sample_Type = c(rep('VP', 18),
                rep('VPC', 18)),
                Replicate = c(rep(1, 6), rep(2, 6), rep(3, 6),
                              rep(4, 6), rep(5, 6), rep(6, 6)),
                Count = c(1, 2, 5, 6, 7, 4,
                          1.1, 2.3, 5.3, 6.3, 7.5, 4.3,
                          0.9, 1.9, 5.8, 5.2, 6.9, 4.1,
                          1.2, 3.1, 4.0, 6.2, 5.4, 3.0,
                          1.3, 3.4, 3.9, 6.4, 5.9, 3.4,
                          0.8, 3.7, 3.6, 6.1, 5.8, 3.5)
                )


df <- df %>% group_by(Sample_Type, Replicate, Timepoint)
df$Sample_Type<- factor(df$Sample_Type, levels = c('VP', 'VPC', 'Diff'))
#Mean and SE
df2<- df%>% group_by(Sample_Type, Timepoint) %>%
  summarise(mean = mean(Count), se = plotrix::std.error(Count)) %>%
  as.data.frame()%>%
  mutate(Timepoint= as.factor(Timepoint))%>%
  mutate(Sample_Type= as.factor(Sample_Type))

df3_a<- df2 %>%
  select(-se)%>%
  pivot_wider(names_from = Sample_Type,
                            values_from = c(mean))%>%
  mutate(Diff = VPC - VP)%>%
  pivot_longer(cols = -Timepoint,
               names_to = 'Sample_Type',
               values_to = 'mean')
  

df3_b<- df2 %>%
  select(-mean)%>%
  pivot_wider(names_from = Sample_Type,
              values_from = c(se))%>%
  mutate(Diff = VPC - VP)%>%
  pivot_longer(cols = -Timepoint,
               names_to = 'Sample_Type',
               values_to = 'se')
 df3<- merge(df3_a, df3_b, by.x = c('Timepoint', 'Sample_Type'),
        by.y = c('Timepoint', 'Sample_Type'))
 
 df3$Sample_Type<- factor(df3$Sample_Type, levels = c('VP', 'VPC', 'Diff'))
  


bp<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')
  
bp1<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.1)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments',  shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')

bp2<- ggplot(data = df2, aes(x = as.factor(Timepoint), y = mean, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  geom_errorbar(data = df2, aes(ymin=mean-se-1, ymax=mean+se+1), width=.2, alpha = 0.5)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')

# bp2<- ggplot(data = df3, aes(x = as.factor(Timepoint), y = mean, color = Sample_Type, shape = Sample_Type))+
#   geom_point(alpha = 1)+
#   geom_errorbar(data = df3, aes(ymin=mean-se, ymax=mean+se), width=.2)+
#   theme_bw()+
#   xlab("Timepoints")+
#   ylab("Viral Counts")+
#   labs(color = 'Treatments', shape = 'Treatments')+
#   scale_color_lancet()+
#   theme(legend.position = 'bottom')

bp3<- ggplot(data = df3, aes(x = as.factor(Timepoint), y = mean, color = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5)+
  #geom_errorbar(data = df3, aes(ymin=mean-se, ymax=mean+se), width=.2)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_lancet()+
  theme(legend.position = 'bottom')


cols<- c("VP" = '#d32f27ff',
         "VPC" = '#4888a2ff',
         "Diff" = '#edb81dff')

shapes<- c("VP" = 15,
          "VPC" = 17,
          "Diff" = 19)

theme<- theme_bw()+
  theme(legend.position = 'none',
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))


#0_Data
lm0<- bp +
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('DATA', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
  

#1_ LM_AP
lm1<- bp +
  geom_line(stat = 'smooth', aes(group = Sample_Type) , method = 'lm', se = F, lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('LM-1', x = 0.2, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
lm1

#2_LM_SR_AVG
lm2<- bp+
  geom_line(stat = 'smooth',aes(group = Replicate), method = 'lm', se = F, lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme+ 
  annotation_custom(grobTree(textGrob('LM-2', x = 0.2, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
lm2

#3_LM_AR
lm3<- bp1 +
  geom_point(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type, alpha = 0.5))+
  geom_line(stat = 'smooth',data = df2, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), 
            method = 'lm', se = F, lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme+ 
  annotation_custom(grobTree(textGrob('LM-3', x = 0.2, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
lm3

#4_LM_AR_Diff

lm4 <- bp1+
  geom_point(data = df3, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type, alpha = 0.5))+
  geom_line(stat = 'smooth',data = df3, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), method = 'lm', se = F, lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme+ 
  annotation_custom(grobTree(textGrob('LM-4', x = 0.2, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
lm4

#5_LM_aR_Diff_LMER
lm5 <- bp1+
  geom_point(data = df3, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type, alpha = 0.5))+
  geom_line(stat = 'smooth',data = df3, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), method = 'lm', se = F, lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme+ 
  annotation_custom(grobTree(textGrob('LM-5', x = 0.2, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
lm5


#6_VPCL_SR_AVG

vpcl1<- bp+
  geom_line(aes(group = Replicate), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme+ 
  annotation_custom(grobTree(textGrob('VPCL-1', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl1


#7_VPCL_AR_No_SE

vpcl2<- bp1+
  geom_point(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type, alpha = 0.1))+
  geom_line(data = df2, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('VPCL-2', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl2
#8_VPCL_AR_SE

vpcl3<- bp1 +
  geom_errorbar(inherit.aes = F, data = df2, aes(x = as.factor(Timepoint), ymin= (mean-se-1), ymax= (mean+se+1), color = Sample_Type), width=.2, alpha = 0.5)+
  geom_point(data = df2, aes(x = as.factor(Timepoint), y = mean, shape = Sample_Type, alpha = 0.1))+
  geom_line(data = df2, aes(x = as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('VPCL-3', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl3  

#9_VPCL_AR_Diff_No_SE

vpcl4<- bp3 +
  geom_point( data = df, aes(x = as.factor(Timepoint), y = Count, shape = Sample_Type ), alpha = 0.1)+ 
  geom_line(aes(x =  as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('VPCL-4', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl4
#10_VPCL_AR_Diff_SE

vpcl5<- bp3 +
  geom_point(data = df, aes(x = as.factor(Timepoint), y = Count, shape = Sample_Type), alpha = 0.1)+
  geom_errorbar(data = df3, aes(x = as.factor(Timepoint),ymin=mean-se-1, ymax=mean+se+1), width=.2, alpha = 0.5)+
  geom_line(aes(x =  as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('VPCL-5', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl5
#11_VPCL_AR_Diff_LMER_No_SE

vpcl6<- bp3 +
  geom_point(data = df, aes(x = as.factor(Timepoint), y = Count, shape = Sample_Type), alpha = 0.1)+ 
  geom_line(aes(x =  as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('VPCL-6', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl6
#12_VPCL_AR_Diff_LMER_SE

vpcl7<- bp3 +
  geom_point(data = df, aes(x = as.factor(Timepoint), y = Count, shape = Sample_Type), alpha = 0.1)+
  geom_errorbar(data = df3, aes(x = as.factor(Timepoint), ymin=mean-se-1, ymax=mean+se+1), width=.2, alpha = 0.5)+
  geom_line(aes(x =  as.factor(Timepoint), y = mean, group = Sample_Type), lwd = 1, alpha = 0.5)+
  scale_colour_manual(breaks = c('VP','VPC','Diff'),
                      values = cols)+
  scale_shape_manual(breaks = c('VP','VPC','Diff'),
                     values = shapes)+
  theme + 
  annotation_custom(grobTree(textGrob('VPCL-7', x = 0.25, y = 0.9,
                                      gp = gpar(fontzie = 12,
                                                fontface = 'bold'))))
vpcl7

vpcl7 +theme(legend.position = 'right',
             legend.text = element_text(size = 10,
                                        face = 'bold'),
             legend.title = element_text(size = 10,
                                        face = 'bold'))


library(cowplot)
lm0_grid<- plot_grid(NULL, lm0, NULL,
                     nrow = 3,
                     rel_heights = c(0.5,1.0, 0.5))
  



lm_grid<- plot_grid(lm1, lm2, lm3, lm4, lm5, NULL, get_legend(vpcl7 +theme(legend.position = 'right',
                                                                           legend.text = element_text(size = 10,
                                                                                                      face = 'bold'),
                                                                           legend.title = element_text(size = 10,
                                                                                                       face = 'bold'))),
          nrow = 1)


vpcl_grid<- plot_grid(vpcl1, vpcl2, vpcl3, vpcl4, vpcl5, vpcl6, vpcl7,
                      nrow = 1)

mock_df_plots<- plot_grid(lm0_grid, plot_grid(lm_grid, vpcl_grid, nrow = 2),
          ncol=2,
          rel_widths = c(1,7))


mock_df_plots
# 
# ggsave("mock_df_plots.svg", width = 800, height =200,  unit = 'px', dpi = 800)
# #I used the xport option to save 1400 x 350


















ggplot(data = df , aes(x = as.factor(Timepoint), y = Count, color = Sample_Type))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  
  theme(legend.position = 'bottom')+
    geom_line(stat = 'smooth',method = 'lm', aes(group = Sample_Type), se = F)+
  scale_color_lancet()





bp2

bp3

bp  +
  geom_line(aes(group = Replicate))
