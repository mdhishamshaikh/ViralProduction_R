library(tidyverse)

#Let's do ROGME

library(rogme)


df<- read_csv("./simulations/simu_output_df.csv")
df<- df %>% filter(VP_Type == 'VPCL_AR_Diff_LMER_SE' | VP_Type == 'VPCL_AR_Diff_No_SE') %>% 
  mutate(VP_Type = factor(VP_Type, levels = c('VPCL_AR_Diff_LMER_SE', 'VPCL_AR_Diff_No_SE'),
                          labels = c('VIPCALSE', 'VIPCAL')))
unique(df$VP_Type)

kruskal.test(VP ~ VP_Type, data = df)

df<- df[c(1:100,5901:6000) ,]# %>% filter(Sample_Type == 'Diff')

ps <- plot_scat2(df,
                 formula = VP ~ VP_Type,
                 alpha = 0.2,
                 shape = 21,
                 colour = 'pink',
                 fill = 'orange')
ps <- ps + coord_flip()+

  theme(plot.title = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 2))
ps

hd_bars<- plot_hd_bars(ps,
             col = "black",
             q_size = 0.5,
             md_size = 1.5,
             alpha = 1)

hd_bars




#sf <- shifthd_pbci(data = df[sample(nrow(df), size=100), ], formula = VP ~ VP_Type)

#> plot shift function
psf <- plot_sf(sf, plot_theme = 1)[[1]] +
  scale_x_continuous(breaks = seq(4, 6, 0.5)) +
  scale_y_continuous(breaks = seq(-6, 6, 2), limits = c(-6, 6))
psf


cowplot::plot_grid(hd_bars, psf)



hd_bars<- plot_hd_bars(ps,
                       col = "black",
                       q_size = 0.5,
                       md_size = 1.5,
                       alpha = 1)+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 1))
hd_bars+
  labs(x = 'VP Type')





sf <- shifthd(data = df, formula = VP ~ VP_Type)

sf
psf <- plot_sf(sf, plot_theme = 2)[[1]] +
  scale_x_continuous(breaks = seq(4, 6, 0.5)) +
  labs(x = 'Deciles for VIPCAL-SE',
       y = 'Decile Differences: VIPCAL-SE - VIPCAL')

psf











color1<- c('#d32f27ff','#4888a2ff', '#edb81dff')
color1<- color1[color]

#Plotting functions

plot_scat2_hd_bars_vp<- function(df2 = df, x_axis_title = 'VP', y_axis_title = 'VP_Type', color = 1){
  
  color1<- c('#d32f27ff','#4888a2ff', '#edb81dff')
  color1<- color1[color]
  
  
  df_rogme<- df2 
  
  ps <- plot_scat2(df_rogme,
                   formula = VP ~ VP_Type,
                   alpha = 0.3,
                   shape = 21,
                  colour = color1,
                   fill = color1)+ 
    coord_flip()+
    theme(plot.title = element_text(face = "bold", size = 20),
          axis.line = element_line(linewidth = 2))
  ps
  
  hd_bars<- plot_hd_bars(ps,
                         col = "black",
                         q_size = 0.5,
                         md_size = 1.5,
                         alpha = 1)+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1),
          axis.text.y = element_text(angle = 90),
          axis.text = element_text(size = 6))+
    labs(x = paste(y_axis_title),
         y = paste(x_axis_title))
  
  sf <- shifthd(data = df_rogme, formula = VP ~ VP_Type)
  
  
  hd_links<- plot_hd_links(hd_bars, sf[[1]],
                q_size = 0.3,
                md_size = 1.0,
                add_rect = TRUE,
                rect_alpha = 0.05,
                rect_col = "grey50",
                add_lab = TRUE,
                text_size = 2,
                link_alpha = c(0.2, 0.5)) 
  
  return(hd_links)
}

psf_vp<- function(df2 = df, y_axis_title = 'Decile Differences: VIPCAL-SE - VIPCAL', x_axis_title = "" ){
  
  df_rogme<- df2 
  
  sf <- shifthd(data = df_rogme, formula = VP ~ VP_Type)
  
  sf
  psf <- plot_sf(sf, plot_theme = 2, symb_size = 2)[[1]] 
  
  psf<- psf +
    scale_x_continuous(breaks = seq(4, 6, 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1),
          axis.text.y = element_text(angle = 90),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 10))+
    labs(x = paste('Deciles for VIPCAL-SE', x_axis_title),
         y = y_axis_title)
  
  
  
  psf
  
  return(psf)
}



df_vp<- df %>% filter(Sample_Type == 'VP')
df_vpc<- df %>% filter(Sample_Type == 'VPC')
df_diff<- df %>% filter(Sample_Type == 'Diff')


links_vp_plot<- cowplot::plot_grid(plot_scat2_hd_bars_vp(df_vp, 'VP', "", 1 ),
plot_scat2_hd_bars_vp(df_vpc, 'VPC', "VP Type", 2 ),
plot_scat2_hd_bars_vp(df_diff, 'Diff', "", 3 ),
nrow =3,
labels = c(1,3,5),
label_fontface = 'bold')
links_vp_plot


shift_vp_plot<- cowplot::plot_grid(psf_vp(df_vp, "", x_axis_title = "(VP)" ),
                   psf_vp(df_vpc, x_axis_title = "(VPC)" ),
                   psf_vp(df_diff, "", x_axis_title = "(Diff)" ),
                   nrow =3,
                   labels = c(2,4,6),
                   label_fontface = 'bold')

shift_vp_plot

cowplot::plot_grid(links_vp_plot, shift_vp_plot,
                   labels = 1:6)
#saved by exporting at 800 x 600



plot_list<- list()

for(type in unique(df$Sample_Type)){

  df_rogme<- df %>% filter(Sample_Type == type)
  
  ps <- plot_scat2(df_rogme,
                   formula = VP ~ VP_Type,
                   alpha = 0.2,
                   shape = 21,
                   colour = 'pink',
                   fill = 'orange')+ 
    coord_flip()+
    theme(plot.title = element_text(face = "bold", size = 20),
          axis.line = element_line(linewidth = 2))
  ps
  
  hd_bars<- plot_hd_bars(ps,
                         col = "black",
                         q_size = 0.5,
                         md_size = 1.5,
                         alpha = 1)+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 1))+
    labs(x = 'VP Type')
  hd_bars
  
  sf <- shifthd(data = df_rogme, formula = VP ~ VP_Type)
  
  sf
  psf <- plot_sf(sf, plot_theme = 2)[[1]] +
    scale_x_continuous(breaks = seq(4, 6, 0.5)) +
    labs(x = 'Deciles for VIPCAL-SE',
         y = 'Decile Differences: VIPCAL-SE - VIPCAL')
  
  psf
  
  

plot_list[[length(plot_list)+1]]<- cowplot::plot_grid(hd_bars, psf)
rm(hd_bars, psf, df_rogme, ps)
}

cowplot::plot_grid(plot_list[[1]],
                   plot_list[[2]],
                   plot_list[[3]],
                   nrow = 3)



##VIPCAL vs VIPCAL-SE ####

simu_output_df<- read_csv("./simulations/simu_output_df.csv")

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


vp_lm<- lm(VPCL_AR_Diff_LMER_SE ~ VPCL_AR_Diff_No_SE, subset(df7, Sample_Type == 'VP'))
vpc_lm<- lm(VPCL_AR_Diff_LMER_SE ~ VPCL_AR_Diff_No_SE, subset(df7, Sample_Type == 'VPC'))
diff_lm<- lm(VPCL_AR_Diff_LMER_SE ~ VPCL_AR_Diff_No_SE, subset(df7, Sample_Type == 'Diff'))

vp_se_comp_legend + geom_text(x = 25, y = 300, label = lm_eqn(vp_lm), parse = TRUE)



lm2<- function(df, x, z){
  m<- lm(data = df, z ~x)
  return(m)
}




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

