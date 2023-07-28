library(cowplot)
library(ggplot2)

#vipcal VS vipcal-se LM #####

cols<- c("VP" = '#d32f27ff',
         "VPC" = '#4888a2ff',
         "Diff" = '#edb81dff')
shapes<- c("VP" = 15,
          "VPC" = 17,
          "Diff" = 19)

vp_se_comp<- ggplot(data = df7, aes(x = VPCL_AR_Diff_No_SE,
                                           y = VPCL_AR_Diff_LMER_SE,
                                           color = Sample_Type,
                                           fill = Sample_Type,
                                           shape = Sample_Type))+
  geom_point(alpha = 0.2)+
  theme_classic()+
  geom_abline(slope = 1, linewidth = 1)+
  geom_hline(yintercept = 0, col = 'black', linetype = 5)+
  labs(x = 'VIPCAL',
       y = 'VIPCAL-SE')+
  geom_smooth(method = 'lm', linewidth =2, alpha = 0.2)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  scale_shape_manual(values = shapes)+
  theme(axis.line = element_line(linewidth = 1),
        axis.text = element_text(face = 'bold', size = 15),
        axis.title = element_text(face = 'bold', size = 15),
        legend.box.background = element_blank(),
        legend.box.margin = margin(t = 1, l = 1),
        legend.title = element_text(face = 'bold'),
        legend.position = c(0.85, 0.25),
        legend.text = element_text(size = 10))+
  # xlim(c(-0.1, 2.2))+
  # ylim(c(-0.1, 2.2))+
  labs(fill = 'Treatment',
       color = 'Treatment',
       shape = 'Treatment')+
  stat_poly_eq(use_label(c("eq", "R2")))


#vp_se_comp_no_legend +geom_point(alpha = 0.8, size = 3)+geom_smooth(method = 'lm', linewidth =1.5, alpha = 0.1)

#get_legend(vp_se_comp_no_legend +geom_point(alpha = 0.8, size = 3)+geom_smooth(method = 'lm', linewidth =1.5, alpha = 0.1))

vp_se_comp_plot<- plot_grid(get_legend(vp_se_comp +geom_point(alpha = 0.8, size = 3)+
                         geom_smooth(method = 'lm', linewidth =1.5, alpha = 0.1)+
                         theme(legend.position = 'top',
                               legend.text = element_text(size = 10))),
          vp_se_comp + theme(legend.position = 'none'),
          ncol = 1,
          rel_heights = c(0.5,6))


vp_se_comp_top_legend<- vp_se_comp +theme(legend.position = 'top')








####2.0 Split Violin plots simulations - vipca vs vipcal-se ####
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

df6<- simu_output_df %>% 
  select(-c(VP_SE, VP_R_Squared, n))%>%group_by(Expt_No, Sample_Type)%>%
  filter(VP_Type %in% c('VPCL_AR_Diff_No_SE',
                        'VPCL_AR_Diff_LMER_SE'))
df6$Sample_Type<- factor(df6$Sample_Type,
                         levels = c("VP",
                                    "VPC",
                                    "Diff"))



split_violin_vp_se_comp <- ggplot(df6,
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
        legend.box.background = element_blank(),
        legend.box.margin = margin(t = 1, l = 1),
        legend.title = element_text(face = 'bold'),
        legend.position = 'top',
        legend.text = element_text(size = 10))+
  labs(x = 'Treatment',
       y = 'Viral Production',
       fill = 'VP Type')

split_violin_vp_se_comp

split_violin_vp_se_comp<- plot_grid(get_legend(split_violin_vp_se_comp +
                                         theme(legend.text = element_text(size = 10))),
                                    split_violin_vp_se_comp + theme(legend.position = 'none'),
                            ncol = 1,
                            rel_heights = c(0.5,6))
                            
                            
cowplot::plot_grid(split_violin_vp_se_comp, NULL, vp_se_comp_plot,
                   labels = c("A", "", "B"),
                   rel_widths = c(1.5, 0.1, 2),
                   nrow = 1)
