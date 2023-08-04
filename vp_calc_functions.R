###---Additional source file for viral_production_step2.R---###
# Different methods for calculating viral production

#1. Linear Regression - All Points (LM_AP) LM1

vp_lm_ap <- function(data){
  
  data_sr<- df_sr_tp(data)
  
  slopes_LM_AP<- slope_all_points(data_sr) %>%
    calc_diff_lm_AP()%>%
    group_by(Location, Expt_No, Depth, Time_Range, Population, Sample_Type) %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')),
            factor(Time_Range, levels = c("T0_T3", "T0_T6", "T0_T17", "T0_T20", "T0_T24")))
  
  slopes_LM_AP$VP_Type<- 'LM_AP'
  return(slopes_LM_AP)
}

#2. Linear Regression - Separate Replicates (LM_SR) LM2

vp_lm_sr<- function(data){
  
  data_sr<- df_sr_tp(data)
  
  slopes_LM_SR<- slope_lm_sr(data_sr)%>%
    group_by(Location, Expt_No, Depth, Sample_Type, Population, Time_Range )%>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  slopes_LM_SR$VP_Type<- 'LM_SR'
  return(slopes_LM_SR)
}

vp_lm_sr_avg<- function(data){
  
  data_sr<- df_sr_tp(data)
  
  slopes_LM_SR<- slope_lm_sr(data_sr)
  
  slopes_LM_SR_avg<- slopes_LM_SR %>%
    group_by(Location, Expt_No, Depth, Time_Range, Population, Sample_Type) %>%
    summarise(VP_Mean=mean(VP_Slope), VP_SE=plotrix::std.error((VP_Slope)), VP_R_Squared = mean(VP_R_Squared)) %>%
    rename(VP = VP_Mean)%>%
    calc_diff_lm_AR() %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            ) 
  
  
  colnames(slopes_LM_SR_avg)[colnames(slopes_LM_SR_avg) == 'VP_Mean']<- 'VP'
  slopes_LM_SR_avg$VP_Type<- 'LM_SR_AVG'
  return(slopes_LM_SR_avg)
}

#3. Linear Regression - Averaged Replicates (LM_AR) LM3

vp_lm_ar<- function(data){
  
  data_avg<- df_avg_tp(data) %>%
    subset(Sample_Type != 'Diff')
  
  slopes_LM_AR<- slope_lm_avg(data_avg) %>%
    calc_diff_lm_AR() %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3')))
  
  slopes_LM_AR$VP_Type<- 'LM_AR'
  
  return(slopes_LM_AR)
}

#4. Linear Regression - Averaged Replicate - Difference Curve (LM_AR_Diff) LM4 (DIFF = Subtraction)

vp_lm_ar_diff<- function(data){
  
  data_avg<- df_avg_tp(data) 
  
  slopes_LM_AR_Diff<- slope_lm_avg(data_avg) %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  
  slopes_LM_AR_Diff$VP_Type <- 'LM_AR_Diff'
  colnames(slopes_LM_AR_Diff)[colnames(slopes_LM_AR_Diff) == 'VP_Slope'] <- 'VP'
  return(slopes_LM_AR_Diff)
}

#5. Linear Regression - Averaged Replicate - Difference Curve - LMER (LM_AR_Diff_LMER) LM5 (DIFF = LMER)

vp_lm_ar_diff_lmer<- function(data){
  
  data_sr<- df_sr_tp(data)
  
  slopes_LM_AR_Diff_LMER<- slope_lm_ar_diff_lmer(data_sr)%>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  
  slopes_LM_AR_Diff_LMER$VP_Type <- 'LM_AR_Diff_LMER'
  colnames(slopes_LM_AR_Diff_LMER)[colnames(slopes_LM_AR_Diff_LMER) == 'VP_Slope'] <- 'VP'
  return(slopes_LM_AR_Diff_LMER)
}

#6. VIPCAL - Separate Replicates (VPCL_SR) VPCL1

vp_vpcl_sr<- function(data){
  data_sr<- df_sr_tp(data)
  
  VPCL_SR<- vipcal_sr(data_sr) %>%
    arrange( 'Location',
             'Expt_No',
             'Depth',
             factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
             factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
             )
  
  VPCL_SR$VP_Type<- 'VPCL_SR'
  
  return(VPCL_SR)
}

vp_vpcl_sr_avg<- function(data){
  data_sr<- df_sr_tp(data)
  
  VPCL_SR_AVG<- vipcal_sr(data_sr)  %>%
    group_by(Location, Expt_No, Depth, Sample_Type, Population, Time_Range ) %>%
    summarise(VP_Mean=mean(VP), VP_SE=plotrix::std.error((VP))) 
  
  colnames(VPCL_SR_AVG)[colnames(VPCL_SR_AVG) == 'VP_Mean']<- 'VP'
  
  
    VPCL_SR_AVG<- VPCL_SR_AVG %>% 
      calc_diff_lm_AR() %>%
      arrange( 'Location',
             'Expt_No',
             'Depth',
             factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
             factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
             )
  
  VPCL_SR_AVG$VP_Type<- 'VPCL_SR_AVG'
  
  return(VPCL_SR_AVG)
}

#7. VIPCAL - Averaged Replicate - No SE (VPCL_AR_No_SE) VPCL2

vp_vpcl_ar_no_se<- function(data){
  
  data_avg<- df_avg_tp(data) %>%
    subset(Sample_Type != 'Diff')
  
  VPCL_AR_No_SE<- vipcal_avg(data_avg) %>%
    calc_diff_vpcl_AR() %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  
  VPCL_AR_No_SE$VP_Type<- 'VPCL_AR_No_SE'
  
  return(VPCL_AR_No_SE)
  
}

#8. VIPCAL - Averaged Replicates - Difference Curve - No SE (VPCL_AR_Diff_No_SE) VPCL4 (DIFF = Subtraction)

vp_vpcl_ar_diff_no_se<- function(data){
  
  data_avg<- df_avg_tp(data) 
  
  VPCL_AR_Diff_No_SE<- vipcal_avg(data_avg) %>%
        arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  
  VPCL_AR_Diff_No_SE$VP_Type<- 'VPCL_AR_Diff_No_SE'
  
  return(VPCL_AR_Diff_No_SE)
  
}

#9. VIPCAL - Averaged Replicate - SE (VPCL_AR_SE) VPCL3

vp_vpcl_ar_se<- function(data){
  
  data_avg<- df_avg_tp(data) %>%
    subset(Sample_Type != 'Diff')
  
  VPCL_AR_SE<- vipcal_avg_se(data_avg) %>%
    calc_diff_vpcl_AR() %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  
  VPCL_AR_SE$VP_Type<- 'VPCL_AR_SE'
  
  return(VPCL_AR_SE)
  
}

#10. VIPCAL - Averaged Replicates - Difference Curve - SE (VPCL_AR_Diff_SE) VPCL5 (DIFF = Substraction)

vp_vpcl_ar_diff_se<- function(data){
  
  data_avg<- df_avg_tp(data) 
  
  VPCL_AR_Diff_SE<- vipcal_avg_se(data_avg) %>%
    arrange('Location',
            'Expt_No',
            'Depth',
            factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
            factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
            )
  
  VPCL_AR_Diff_SE$VP_Type<- 'VPCL_AR_Diff_SE'
  
  return(VPCL_AR_Diff_SE)
  
}

#11. VIPCAL - Averaged Replicates - Difference Curve - LMER - No SE (VPCL_AR_Diff_LMER_No_SE) VPCL6 (DIFF = LMER)

vp_vpcl_ar_diff_lmer_no_se<- function(data){
  
  data_sr<- df_sr_tp(data)
  
  VIPCAL_SR_DIFF_LMER_No_SE<- vipcal_sr_diff_no_SE(data_sr) %>%
    arrange( 'Location',
             'Expt_No',
             'Depth',
             factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
             factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
             )
  
  VIPCAL_SR_DIFF_LMER_No_SE$VP_Type<- 'VPCL_AR_Diff_LMER_No_SE'
  
  return(VIPCAL_SR_DIFF_LMER_No_SE)
}

#12. VIPCAL - Averaged Replicates - Difference Curve - LMER - SE (VPCL_AR_Diff_LMER_SE) VPCL7 (DIFF = LMER)

vp_vpcl_ar_diff_lmer_se<- function(data){
  
  data_sr<- df_sr_tp(data)
  
  VIPCAL_SR_DIFF_LMER_SE<- vipcal_sr_diff_SE(data_sr) %>%
    arrange( 'Location',
             'Expt_No',
             'Depth',
             factor(Sample_Type, levels = c('VP', 'VPC', 'Diff')),
             factor(Population, levels = c('c_Viruses', 'c_V1', 'c_V2', 'c_V3'))
             )
  
  VIPCAL_SR_DIFF_LMER_SE$VP_Type<- 'VPCL_AR_Diff_LMER_SE'
  
  return(VIPCAL_SR_DIFF_LMER_SE)
}

   
#List of all these functions
vp_calc_funct_list<- list(vp_lm_ap, vp_lm_sr_avg, vp_lm_ar,
                          vp_lm_ar_diff, vp_lm_ar_diff_lmer,
                          vp_vpcl_sr_avg, vp_vpcl_ar_no_se, 
                          vp_vpcl_ar_diff_no_se, vp_vpcl_ar_se,
                          vp_vpcl_ar_diff_se, vp_vpcl_ar_diff_lmer_no_se,
                          vp_vpcl_ar_diff_lmer_se)

vp_calc_funct_list_sr<- list(vp_lm_sr, vp_vpcl_sr)
