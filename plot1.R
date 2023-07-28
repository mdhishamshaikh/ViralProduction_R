library(tidyverse)
library(cowplot)


make_a_csv<- read_csv("./simulations/simu_output_df.csv")


#violinplots
####Geom Split Violin function####
{
  GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self,
                          data,
                          ...,
                          # add the nudge here
                          nudge = 0,
                          draw_quantiles = NULL) {
      data <- transform(data,
                        xminv = x - violinwidth * (x - xmin),
                        xmaxv = x + violinwidth * (xmax - x))
      grp <- data[1, "group"]
      newdata <- plyr::arrange(transform(data,
                                         x = if (grp %% 2 == 1) xminv else xmaxv),
                               if (grp %% 2 == 1) y else -y)
      newdata <- rbind(newdata[1, ],
                       newdata,
                       newdata[nrow(newdata), ],
                       newdata[1, ])
      newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
      
      # now nudge them apart
      newdata$x <- ifelse(newdata$group %% 2 == 1,
                          newdata$x - nudge,
                          newdata$x + nudge)
      
      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        
        quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                             draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)),
                           setdiff(names(data), c("x", "y")),
                           drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
                         grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                        quantile_grob))
      }
      else {
        ggplot2:::ggname("geom_split_violin",
                         ggplot2::GeomPolygon$draw_panel(newdata, ...))
      }
    }
  )
  geom_split_violin <- function(mapping = NULL,
                                data = NULL,
                                stat = "ydensity",
                                position = "identity",
                                # nudge param here
                                nudge = 0,
                                ...,
                                draw_quantiles = NULL,
                                trim = TRUE,
                                scale = "area",
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
    
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = stat,
                   geom = GeomSplitViolin,
                   position = position,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(trim = trim,
                                 scale = scale,
                                 # don't forget the nudge
                                 nudge = nudge,
                                 draw_quantiles = draw_quantiles,
                                 na.rm = na.rm,
                                 ...))
  }
  
}

simu_output_df<- read.csv('./simulations/simu_output_df.csv')
violin_df<- simu_output_df%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE')



violin_df$VP_Type <- factor(violin_df$VP_Type, levels = c("LM_SR_AVG", "VPCL_AR_Diff_No_SE"),
                            labels = c('Linear Regression', 'VIPCAL'))
violin_df$Sample_Type <- factor(violin_df$Sample_Type, levels = c("VP", "VPC", "Diff"))

violin_df$VP_Type2<- 'A'
lm_vs_vpcl_violin<- ggplot(violin_df,
                           aes(y = VP,
                               x = VP_Type2,
                               fill = VP_Type))+
  #geom_point()+
  geom_split_violin(nudge = 0.02,
                    color = 'transparent')+
  theme_bw()+
  scale_fill_manual(values = c("#ff9c00ff", '#52495aff'))+
  scale_color_manual(values<- NA)+
  theme(#legend.position =  c(0.5, 0.8),
    legend.text = element_text(size = 8),
    legend.title = element_text(face = 'bold',
                                size = 8),
    legend.position = 'bottom' ,
    legend.direction = 'vertical',
    legend.title.align = 0.5,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10,
                              face = 'bold'),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2))+
  labs(y = 'Viral Production Rate',
       fill = 'Data Analyses Methods')+
  ylim(c(-0.5, 2.2))



lm_vs_vpcl_violin<- (lm_vs_vpcl_violin + theme(legend.position = 'none')) +
  inset_element(get_legend(lm_vs_vpcl_violin),
                0.2, 0.8, 0.8, 0.95) #PLOT: LM VS VIPCAL split violin ####
lm_vs_vpcl_violin

plot_grid(get_legend(lm_vs_vpcl_violin),
          lm_vs_vpcl_violin + theme(legend.position = 'none'),
          nrow = 2,
          rel_heights = c(1,5)
)

lm_vs_vpcl_violin_sample_type<- ggplot(violin_df,
                           aes(y = VP,
                               x = Sample_Type,
                               fill = VP_Type))+
  #geom_point()+
  geom_split_violin(nudge = 0.02,
                    color = 'transparent')+
  theme_bw()+
  scale_fill_manual(values = c("#ff9c00ff", '#52495aff'))+
  scale_color_manual(values<- NA)+
  theme(#legend.position =  c(0.5, 0.8),
    legend.text = element_text(size = 8),
    legend.title = element_text(face = 'bold',
                                size = 8),
    legend.position = 'bottom' ,
    legend.direction = 'vertical',
    legend.title.align = 0.5,
    #axis.title.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.text = element_text(size = 10,
                             face = 'bold'),
    axis.title = element_text(size = 10,
                              face = 'bold'),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=2))+
  labs(y = 'Viral Production Rate',
       x = 'Treatment',
       fill = 'Data Analyses Methods')+
  ylim(c(-0.5, 2.2))

lm_vs_vpcl_violin_sample_type

lm_vs_vpcl_violin_sample_type<- (lm_vs_vpcl_violin_sample_type + theme(legend.position = 'none')) +
  inset_element(get_legend(lm_vs_vpcl_violin_sample_type),
                0.1, 0.8, 0.4, 0.95) #PLOT: LM VS VIPCAL split violin ####
lm_vs_vpcl_violin_sample_type
#violin plot - difference between LM and VIPCAL(VP, VPC, Diff)



#Combined Split Violin####

violin_df2<- violin_df
violin_df2<- violin_df2 %>%
  group_by(Expt_No, Sample_Type) %>%
  summarise()


#checking if LM_sR_AVG and VPCL_AR_Diff_No_SE

vs_df<- simu_output_df%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE')

vs_df<- vs_df %>% select(-VP_SE, -VP_R_Squared, -n)%>%
  group_by(Expt_No, Sample_Type) %>%
  pivot_wider(names_from = VP_Type,
              values_from = VP)

ggplot(data = vs_df,
       aes(y = VPCL_AR_Diff_No_SE,
           x = LM_SR_AVG))+
  geom_point()+
  geom_abline(slope = 1)


for(i in 1:length(vs_df$LM_SR_AVG)){
  
  if(vs_df[i,]$LM_SR_AVG > vs_df[i,]$VPCL_AR_Diff_No_SE){
    print(i)
    print(vs_df[i,])
  }
  
}#only two of the expts have a higher LM_SR_AVG. Expt No: 291, 394
cols<- c("VP" = '#d32f27ff',
         "VPC" = '#4888a2ff',
         "Diff" = '#edb81dff')
shape<- c("VP" = 15,
          "VPC" = 17,
          "Diff" = 19)

vs_df$Sample_Type <- factor(vs_df$Sample_Type, levels = c('VP', 'VPC', 'Diff'))


lm_vs_vpcl_lm<- ggplot(data = vs_df,
                       aes(y = VPCL_AR_Diff_No_SE,
                           x = LM_SR_AVG,
                           group = Sample_Type,
                           col = Sample_Type,
                           shape = Sample_Type))+
  geom_point(alpha = 0.2)+
  geom_abline(slope = 1, lwd = 1)+
  geom_smooth(method = 'lm', linewidth = 2, alpha = 0.3)+
  scale_color_manual(values = cols)+
  theme_bw()+
  labs(x = 'Linear Regression Method\n(LM-SR-AVG)',
       y = 'VICAL Method\n(VPCL-AR-Diff-No-SE)',
       col = 'Treatment',
       shape = 'Treatment')+
  theme( panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
         panel.background = element_rect(fill = NA),
         axis.title = element_text(face = 'bold',
                                   size = 10),


#Mock Plots for LM and VIPCAL #####
         axis.text = element_text(size = 10),
         legend.title = element_text(face = 'bold',
                                     size = 10),
         legend.text = element_text(size = 10),
         legend.position = c(0.15, 0.7))

lm_vs_vpcl_lm #PLOT: lm vs vipcal Linear Regression####

library(tidyverse)

#### 1.0 Make mock plots to explain LM and VIPCAL methods.

library(tidyverse)
library(ggsci)
{
  
  df<- data.frame(Timepoint = rep(c('T1', 'T2', 'T3',
                                    'T4', 'T5', 'T6'),6),
                  Sample_Type = c(rep('VP', 18),
                                  rep('VPC', 18)),
                  Replicate = c(rep(1, 6), rep(2, 6), rep(3, 6),
                                rep(4, 6), rep(5, 6), rep(6, 6)),
                  Count = c(2.0, 0.7, 5.0, 6.9, 8.7, 5.0,
                            2.1, 0.3, 5.3, 6.3, 8.3, 5.3,
                            2.9, 1.9, 5.8, 7.2, 7.9, 5.1,
                            1.2, 2.1, 2.0, 6.2, 8.4, 3.0,
                            1.3, 2.4, 1.9, 6.4, 8.2, 3.4,
                            0.8, 2.7, 1.6, 6.1, 8.8, 3.5)
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
  
  
}




#linear background plot with VP and VPC
bp<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_classic()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme(legend.position = 'none',
        axis.text = element_text(size = 20, 
                                 face = 'bold'),
        axis.title = element_blank())+
  ylim(0,10)

bp




#####MAIN PLOTS####
library(ggsci)
{
  
  df<- data.frame(Timepoint = rep(c('T1', 'T2', 'T3',
                                    'T4', 'T5', 'T6'),6),
                  Sample_Type = c(rep('VP', 18),
                                  rep('VPC', 18)),
                  Replicate = c(rep(1, 6), rep(2, 6), rep(3, 6),
                                rep(4, 6), rep(5, 6), rep(6, 6)),
                  Count = c(2.0, 0.7, 5.0, 6.9, 8.7, 5.0,
                            2.1, 0.3, 5.3, 6.3, 8.3, 5.3,
                            2.9, 1.9, 5.8, 7.2, 7.9, 5.1,
                            1.2, 2.1, 2.0, 6.2, 8.4, 3.0,
                            1.3, 2.4, 1.9, 6.4, 8.2, 3.4,
                            0.8, 2.7, 1.6, 6.1, 8.8, 3.5)
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
  
  
}




cols<- c("VP" = '#d32f27ff',
         "VPC" = '#4888a2ff',
         "Diff" = '#edb81dff')
shape<- c("VP" = 15,
          "VPC" = 17,
          "Diff" = 19)

#linear background plot with VP and VPC
bp<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.5,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_classic()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme(legend.position = 'none',
        axis.text = element_text(size = 20, 
                                 face = 'bold'),
        axis.title = element_blank())+
  ylim(0,10)

bp



theme_mockplots<- theme(legend.position = 'none',
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=2),
                        panel.background = element_rect(fill = NA),
                        axis.text = element_text(size = 10, 
                                                 face = 'bold'),
                        axis.title = element_blank())

#lm - vp

lm_1<- ggplot(data = df %>% filter(Sample_Type == 'VP'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = .8,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1, size = 0.5)

#LM - VPC
lm_2<- ggplot(data = df %>% filter(Sample_Type == 'VPC'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.8,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1, size = 0.5)

#LM -VP VPC
lm_3<- ggplot(data = df , aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.8,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = c(15, 17))+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_smooth(aes(group = Replicate), method = 'lm', se = F, alpha = 0.1, size = 0.5)



#vpcl - vp

vpcl_1<-ggplot(data = df %>% filter(Sample_Type == 'VP'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.2,
             size =2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = shape)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_point(data = df3 %>% filter(Sample_Type == 'VP') ,aes(x = as.factor(Timepoint), y = mean, 
                                                             color = Sample_Type, fill = Sample_Type, 
                                                             shape = Sample_Type), alpha = 1, size =2.0) +
  geom_line(data = df3 %>% filter(Sample_Type == 'VP'),
            aes(y = mean, group = Sample_Type),alpha = 1, size = 0.5)

#vpcl -vpc

vpcl_2<-ggplot(data = df %>% filter(Sample_Type == 'VPC'), aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.2,
             size = 2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = shape)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_point(data = df3 %>% filter(Sample_Type == 'VPC') ,aes(x = as.factor(Timepoint), y = mean, 
                                                              color = Sample_Type, fill = Sample_Type, 
                                                              shape = Sample_Type), alpha = 1, size = 2.0) +
  geom_line(data = df3 %>% filter(Sample_Type == 'VPC'),
            aes(y = mean, group = Sample_Type),alpha = 1, size = 0.5)

#vpcl - vp vpc diff

vpcl_3<- ggplot(data = df, aes(x = as.factor(Timepoint), y = Count, color = Sample_Type, fill = Sample_Type, shape = Sample_Type))+
  geom_point(alpha = 0.2,
             size = 2.0,
             position = position_jitter(width = 0.1))+
  scale_shape_manual(values = shape)+
  theme_bw()+
  xlab("Timepoints")+
  ylab("Viral Counts")+
  labs(color = 'Treatments', shape = 'Treatments')+
  scale_color_manual(values = cols)+
  theme_mockplots+
  ylim(-4,10)+
  geom_point(data = df3 ,aes(x = as.factor(Timepoint), y = mean, 
                             color = Sample_Type, fill = Sample_Type, 
                             shape = Sample_Type), alpha = 1, size = 2.0) +
  geom_line(data = df3 %>% filter(Sample_Type != 'Diff'),
            aes(y = mean, group = Sample_Type),alpha = 0.2, size = 0.5)+
  geom_line(data = df3 %>% filter(Sample_Type == 'Diff'),
            aes(y = mean, group = Sample_Type),alpha = 1, size = 0.5)


mock_plots_lm_vpcl<- plot_grid(lm_1, lm_2, lm_3,
          vpcl_1, vpcl_2, vpcl_3,
          nrow = 2) #PLOT: LM VS VIPCAL MOCKPLOTS####

mock_plots_lm_vpcl

og_plots_list<- list(lm_1, vpcl_1, lm_2, vpcl_2, lm_3, vpcl_3)
og_plots_list_name<- c('LM_1', 'VPCL_1',
                       'LM_2', 'VPCL_2',
                       'LM_3', 'VPCL_3')

for (i in 1:6){
  print(i)
  
  ggsave(filename = paste0(og_plots_list_name[i], ".svg"),
         plot = og_plots_list[[i]],
         device = 'svg',
         width = 4,
         height = 4,
         dpi = 800,
         units = 'in')
  
}

#### Literature Study ####
library(tidyverse)
library(ggsci)
library(cowplot)


data<- readxl::read_excel("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/2023_VP_Assay_MethodComparison/2023_Literature_Survey/LM_VIPCAL_Literature_Survey/../2022_VP_Assay_Method_Comparison_Literature_Survey_HMS.xlsx")
data<- data[1:38,]
data<- data%>%arrange(Year_Published)



df1<- data%>% filter(Method == 'LM')%>%
  group_by(Year_Published) %>%  arrange(Year_Published)%>%
  summarise(LM_Lytic = sum(Lytic), LM_Lysogenic = sum(Lysogenic))
df2<- data%>%filter(Method == 'VIPCAL')%>%
  group_by(Year_Published)%>% arrange(Year_Published)%>%
  summarise(VPCL_Lytic = sum(Lytic), VPCL_Lysogenic = sum(Lysogenic))

cumulative_df<- full_join(df1,df2)%>%
  arrange(Year_Published)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  mutate(across(!Year_Published, cumsum))%>%
  pivot_longer(cols = 2:5,
               names_to = 'Type',
               values_to = 'Study_No')
cumulative_df$Type<- factor(cumulative_df$Type, levels = c("LM_Lytic", "LM_Lysogenic",  "VPCL_Lytic", "VPCL_Lysogenic"),
                            labels =  c("LM-Lytic", "LM-Lysogenic",   "VPCL-Lytic", "VPCL-Lysogenic"))

ggplot(data = cumulative_df, aes(x = Year_Published, y = Study_No, fill = Type ))+
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_line(lwd = 0.5 )+
  scale_fill_npg()+
  theme_classic()+
  scale_x_continuous(breaks = c(2005:2020))+
  xlab("Year of Publishing")+
  ylab("No of Experiments Performed")+
  labs(title ="Number of Virus Reduction Assays performed over the years")


d1<- data%>%filter(Method == 'LM')%>%
  filter(Lytic != 0)%>%
  arrange(Year_Published)%>%
  group_by(Year_Published)%>%
  summarise(LM_Lytic = n())%>%
  add_row(Year_Published = 2015, LM_Lytic = 0)%>%
  add_row(Year_Published = 2016, LM_Lytic = 0)


d2<- data%>%filter(Method == 'LM')%>%
  filter(Lysogenic != 0)%>%
  arrange(Year_Published)%>%
  group_by(Year_Published)%>%
  summarise(LM_Lysogenic = n())%>%
  add_row(Year_Published = 2015, LM_Lysogenic = 0)%>%
  add_row(Year_Published = 2016, LM_Lysogenic = 0)

d3<- data%>%filter(Method == 'VIPCAL')%>%
  filter(Lytic != 0)%>%
  arrange(Year_Published)%>%
  group_by(Year_Published)%>%
  summarise(VPCL_Lytic = n())%>%
  add_row(Year_Published = 2015, VPCL_Lytic = 0)%>%
  add_row(Year_Published = 2016, VPCL_Lytic = 0)

d4<- data%>%filter(Method == 'VIPCAL')%>%
  filter(Lysogenic != 0)%>%
  arrange(Year_Published)%>%
  group_by(Year_Published)%>%
  summarise(VPCL_Lysogenic = n())%>%
  add_row(Year_Published = 2015, VPCL_Lysogenic = 0)%>%
  add_row(Year_Published = 2016, VPCL_Lysogenic = 0)


study_Df<- full_join(d1,d2)%>%
  full_join(d3)%>%
  full_join(d4)%>%
  mutate(across(everything(), ~replace_na(.x,0)))%>%
  arrange(Year_Published)%>%
  mutate(across(!Year_Published, cumsum))%>%
  pivot_longer(cols = 2:5,
               names_to = 'Type',
               values_to = 'Study_No')


study_Df$Type<- factor(study_Df$Type, levels = c("LM_Lytic", "LM_Lysogenic",  "VPCL_Lytic", "VPCL_Lysogenic"),
                       labels =  c("LM-Lytic", "LM-Lysogenic",   "VPCL-Lytic", "VPCL-Lysogenic"))


cols_black_pink<-  c('#474747',
                     '#CC527A',
                     '#363636',
                     '#E8175D'
)
cols_red_blue_pastel<- c('#355C7D',
                         '#6C5B7B',
                         '#C06C84',
                         '#F67280')


cols_orange_black<- c('#ff9c00ff', '#eac996ff', '#363636ff', '#605d5dff')


studies<- ggplot(data = study_Df, aes(x = Year_Published, y = Study_No, fill = Type ))+
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label = Study_No), position = position_dodge(width =2))+
  #geom_line(lwd = 0.5 )+
  scale_fill_manual(values =cols_orange_black)+
  theme_classic()+
  scale_x_continuous(breaks = c(2005, 2010, 2015, 2020))+
  xlab("Year")+
  ylab("No of Studies\nPublished")+
  #labs(title ="Number of published studies using Virus Reduction Approach over the years")+
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = 'bold',
                                 size = 8),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(face = 'bold',
                                  size = 8),
        axis.line = element_line(linewidth = 1))


assays_legend<- ggplot(data = cumulative_df, aes(x = Year_Published, y = Study_No, fill = Type ))+
  geom_bar(stat = 'identity', position = 'dodge')+
  #geom_text(aes(label = Study_No), position = position_dodge(width =2))+
  #geom_line(lwd = 0.5 )+
  scale_fill_manual(values =cols_orange_black)+
  theme_classic()+
  scale_x_continuous(breaks = c(2005, 2010, 2015, 2020))+
  xlab("Year")+
  ylab("No of Assays\nPerformed")+
  # labs(title ="Number of assays performed using Virus Reduction Approach over the years")+
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = 'bold',
                                 size = 8),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(face = 'bold',
                                  size = 8),
        axis.line = element_line(linewidth = 1))
legend <- get_legend(assays_legend) 

assays<- assays_legend +
  theme(legend.position = 'none')


plot_grid(NULL, legend, NULL,
          nrow = 3
)


vp_ls_plot<- plot_grid(studies, NULL, plot_grid(NULL, legend, NULL,
                                                nrow = 3
), NULL, assays, 
ncol = 5,
rel_widths =  c(1,0.25, 0.25, 0.15,1))




ls_lm_vpcl_plot<- plot_grid(plot_grid(studies, assays,
                                      nrow = 2), legend,
                            nrow = 2,
                            rel_heights = c(6,1))


ls_lm_vpcl_plot #PLOT: Literature Study ####





####COMBINE EVERYTHING####

lm_vpcl_violin_plots<- plot_grid(lm_vs_vpcl_violin, lm_vs_vpcl_violin_sample_type)
lm_vpcl_violin_plots

row3<- plot_grid(lm_vpcl_violin_plots, ls_lm_vpcl_plot)

row2<- plot_grid(mock_plots_lm_vpcl, lm_vs_vpcl_lm,
          rel_widths = c(1.5,1))
plot_grid(row2, row3,
          nrow =2,
          rel_heights = c(1,1))
