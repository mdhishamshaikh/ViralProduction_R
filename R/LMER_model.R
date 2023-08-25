LMER_model <- function(DF){ 
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