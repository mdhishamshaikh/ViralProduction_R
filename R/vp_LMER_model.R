#' Linear Mixed-Effects Model
#' 
#' @description
#' The viral reduction assay has two types of samples: VP and VPC. In VP samples, count of bacteriophages in the lytic phase
#' can be measured (`lytic viral production`). On the other hand in VPC samples, count of bacteriophages in both the lytic and
#' lysogenic phase can be measured (`lytic + lysogenic viral production`) since treatment with antibiotic `mitomycin-C`
#' forces lysogenic bacteriophages to go into the lytic phase. To retrieve the lysogenic viral production, a difference
#' curve is used. The difference curve can be estimated by subtracting the VP samples count from the VPC samples count, `difference curve estimation by subtraction`. 
#' But also, a \code{LMER model} can be used which will incorporate both fixed- and random effect terms. The influence 
#' of the sample type and time point on the viral count will be considered plus the variability between the different replicates (random-effect), 
#' `difference curve estimation by LMER model`. 
#' 
#' See [lme4::lmer] for more details on the LMER model.
#' 
#' @param DF Data frame with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe] for more details.
#' 
#' @return Data frame with the mean viral count for VP samples, VPC samples and the difference (DIFF samples) for each population at the different time points of the assay.
#' 
#' @name vp_LMER_model
#' @rdname vp_LMER
#' 
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_all)
#' 
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
#' 
#' vp_LMER_model(DF_SR)
#' }
vp_LMER_model <- function(DF){ 
  result_DF <- data.frame() 
  
  # VP and VPC have same replicate reference, needs to be unique for model
  for (replicate in unique(DF$Replicate)){
    DF$Replicate[DF$Replicate == replicate & DF$Sample_Type == 'VPC'] <- replicate + 3
  }
  
  for (type in unique(DF$Population)){
    lmer_model <- lme4::lmer(data = DF, Count ~ Sample_Type * as.factor(Timepoint) + (1 | Replicate))
    
    estimate_marginal_means <- emmeans::emmeans(lmer_model, ~ Sample_Type | as.factor(Timepoint))
    results_VP_VPC_samples <- data.frame(rep(type, length(unique(DF$Timepoint))),
                                         summary(estimate_marginal_means)$Sample_Type,
                                         summary(estimate_marginal_means)$Timepoint,
                                         summary(estimate_marginal_means)$emmean,
                                         summary(estimate_marginal_means)$SE)
    colnames(results_VP_VPC_samples) <- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    
    contrasts_emmeans <- graphics::pairs(estimate_marginal_means)
    results_DIFF_samples <- data.frame(rep(type, length(unique(DF$Timepoint))),
                                       rep("Diff", length(unique(DF$Timepoint))),
                                       summary(contrasts_emmeans)$Timepoint,
                                       -(summary(contrasts_emmeans)$estimate),
                                       summary(contrasts_emmeans)$SE)
    colnames(results_DIFF_samples) <- c("Population", "Sample_Type", "Timepoint", "Mean", "SE")
    
    result_DF <- rbind(result_DF, results_VP_VPC_samples, results_DIFF_samples)
  }
  return(result_DF)
}