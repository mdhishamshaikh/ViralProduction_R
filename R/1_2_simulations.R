#####I. Simulation####
#####1.0 Installing Packages####y
source("C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/1_1_simulations_source.R")


####1.0 Creating a dataframe####
# #Set number of dataframes you'd like to create
simu_length<- 1000
{


  set.seed(2023)
  simu_df<- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(simu_df)<- c("Timepoint", "Count", "Replicate", "Expt_No", "Sample_Type")
  simu_df<- simu_df %>%  mutate(across(!Sample_Type, as.numeric)) %>%
    mutate(across(Sample_Type, as.character))


  for (df_no in 1:simu_length){

    A<- data.frame(Count = runif(n=12, min = 1, max = 2))
    A$Timepoint<- rep(1:6,2)
    A$SD_percent<- runif(n = 12, min = 0, max = 100) #UP TO 100% ERROR
    A$SD<- (A$Count*A$SD_percent)/100


    B<- apply(A, 1, function(x){x[1]+x[4]*scale(rnorm(3))} )
    B<- as.data.frame(B)
    colnames(B)<- rep(1:6,2)
    B<- B%>% pivot_longer(cols = everything(), names_to = 'Timepoint',
                          values_to = 'Count')
    B<- B%>%
      arrange(Timepoint) %>%
      mutate(Replicate = rep(1:6, 6))%>%
      mutate(Sample_Type = rep(c(rep('VP',3), rep('VPC',3)),6))

    B$Expt_No = df_no


    B<- B%>% mutate(across(!Sample_Type, as.numeric))

    try(simu_df<- simu_df %>%full_join(B))
    rm(A,B, df_no)
  }
}

####2.0 Simulation####

#Define the length of sim_df you'd like to analyse
length_df<- unique(simu_df$Expt_No)

#Define variables
var<- c('Expt_No', 'Sample_Type')
var2<- c('Expt_No', 'Sample_Type', 'Timepoint')



#Running the simulation


make_a_csv<- simulation_function(input_df = simu_df)



make_a_csv<- make_a_csv %>%
  filter(VP_Type != 'VPCL_AR_Diff_LMER_SE')%>%
  full_join(VPCL_AR_Diff_LMER_SE_output)

make_a_csv$VP_Type<- factor(make_a_csv$VP_Type, levels = c('LM_AP', 'LM_SR_AVG',
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
# make_a_csv2<- make_a_csv2 %>%
#   filter(VP_Type != 'VPCL-7')

unique(make_a_csv$VP_Type)

write.csv(make_a_csv,
          file = './simulations/simu_output_df.csv',
          row.names = F)




####VISUALIZATION####
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
           fill = Sample_Type))+
  geom_violin(draw_quantiles = 1)+
  scale_fill_npg()

ggbetweenstats(
  data = lm_output,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: LM VP")

ggbetweenstats(
  data = lm_output,
  x = VP_Type,
  y = VP_SE,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: LM VP_SE")



####vISUAL;IZATION#####

#Comparison
ggbetweenstats(
  data = simu_output_df,
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Simu")


ggbetweenstats(
  data = simu_output_df%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE'),
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "s",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Simu Diff") #significantly different


ggbetweenstats(
  data = simu_output_df%>% filter(Sample_Type == 'Diff'),
  x = VP_Type,
  y = VP,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "ns",
  centrality.plotting = FALSE,
  bf.message = FALSE
)+labs(title = "Kruskal-Wallis Test: Simu Diff")


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
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2))+
  labs(y = 'Viral Production Rate',
       fill = 'Data Analyses Methods')

plot_grid(get_legend(lm_vs_vpcl_violin),
          lm_vs_vpcl_violin + theme(legend.position = 'none'),
          nrow = 2,
          rel_heights = c(1,5)
          )

#checking if LM_sR_AVG and VPCL_AR_Diff_No_SE

vs_df<- simu_output_df%>% filter(VP_Type == 'LM_SR_AVG'|  VP_Type == 'VPCL_AR_Diff_No_SE')

vs_df<- vs_df%>% arrange() %>% select(-VP_SE, -VP_R_Squared)%>%
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

}#only two of the expts have a higher LM_SR_AVG
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
        axis.text = element_text(size = 10),
        legend.title = element_text(face = 'bold',
                                    size = 10),
        legend.text = element_text(size = 10),
        legend.position = c(0.15, 0.8))

lm_vs_vpcl_lm



#Adding density

dens1<-ggplot(vs_df, aes(x = LM_SR_AVG))+
  geom_density(fill = '#ff9c00ff', color =NA,
    alpha = 1)+
  theme_void()

dens2<- ggplot(vs_df, aes(x = VPCL_AR_Diff_No_SE))+
  geom_density(fill = '#52495aff', color =NA,
               alpha = 1)+
  theme_void()+
  coord_flip()

dens1 + plot_spacer() + lm_vs_vpcl_lm + dens2 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(6, 1), heights = c(1, 6))




ggplot(simu_output_df,
       aes(x = VP,
           fill = VP_Type))+
geom_histogram(bins = 100)

ggplot(data = simu_output_df,
       aes(x = VP,
           y = VP_Type,
           color = VP_Type, fill = VP_Type)) +
  geom_density_ridges(alpha = 0.8, scale = 5) +
  scale_fill_viridis(option = "A", discrete = TRUE) +
  scale_color_viridis(option = "A", discrete = TRUE) +
  theme_bw()
# 
