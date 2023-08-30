#' Linear Mixed-Effects Model
#' 
#' @description
#' To calculate lysogenic production, a difference curve is used. The average of increments of the difference curve 
#' represent the lysogenic production. The difference curve can be estimated by subtracting the slope of VP samples,
#' represent lytic viral production, from the slope of VPC samples, representing lytic + lysogenic viral production.
#' On the other hand, a \code{LMER model} can be used which will incorporate both fixed- and random effects terms. The
#' influence of the sample type and time point on the viral count will be considered plus the variability between the different replicates (random-effect).
#' See [lme4::lmer] for more details on the LMER model.
#' 
#' @param DF Data frame with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe].
#' 
#' @return Data frame with the mean viral count for VP samples, VPC samples and the difference (DIFF samples) for each population at each time point.
#' 
#' @name vp_LMER_model
#' @rdname vp_LMER
#' 
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
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