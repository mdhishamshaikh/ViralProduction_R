#' Calculate viral production with linear regression
#' 
#' @description
#' `Linear regression` employs first-order regressions (slopes) of viral abundance over time to calculate viral production.
#' Lytic viral production can be derived as the slope of the VP samples. For the lysogenic viral production,
#' subtraction of the VP slope from the VPC slope is needed. A difference curve estimation by subtraction or
#' LMER model is also an option to determine the lysogenic viral production rate. 
#' 
#' See [stats::lm] for more details on the linear model used to determine the slope. 
#' 
#' `determine_vp_linear_allpoints` considers all points, no distinguishing between replicates.
#' 
#' `determine_vp_linear_separate_replicates` uses a separate replicate treatment, distinguishing between replicates. 
#' 
#' `determine_vp_linear_average_replicates` uses an average replicate treatment, average over the replicates,
#' difference curve estimation by subtraction. 
#' 
#' `determine_vp_linear_LMER_model` uses an average replicate treatment, average over the replicates is included in LMER model, 
#' difference curve estimation by LMER model.
#' 
#' See [viralprod::vp_LMER_model] for more details about the LMER model.
#'
#' @param SR_dataframe Data frame with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe] for more details.
#' @param AVG_dataframe Data frame with the viral counts and time ranges, see [viralprod::vp_average_replicate_dataframe] for more details.
#' 
#' @return Data frame with the viral production rate and the absolute viral production for each population at the different time points of the assay.
#' 
#' @name determine_vp_linear_regression
#' @rdname vp_LM
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_all)
#' 
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020_all)
#' 
#' determine_vp_linear_allpoints(DF_SR)
#' 
#' determine_vp_linear_separate_replicates(DF_SR)
#' 
#' determine_vp_linear_average_replicates(DF_AVG)
#' 
#' determine_vp_linear_LMER_model(DF_SR)
#' }
determine_vp_linear_allpoints <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      for (sample in unique(SR_dataframe$Sample_Type)){
        DF <- SR_dataframe %>%
          dplyr::filter(.data$tag == combi_tag, .data$Population == virus, .data$Sample_Type == sample)
        
        for (time in unique(DF$Time_Range)){
          DF2 <- DF %>%
            dplyr::filter(.data$Time_Range == time)
          
          linear_model <- summary(stats::lm(data = DF2, Count ~ as.numeric(Timepoint)))
          absolute_vp_values <- linear_model$coefficients[2] * (DF2$Timepoint[length(DF2$Timepoint)])
          
          result <- c(combi_tag, time, virus, sample, linear_model$coefficients[2], absolute_vp_values, linear_model$coefficients[4], linear_model$r.squared)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_linear <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_linear) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 
                                         'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')
  viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')] <- lapply(viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  viral_production_linear <- viral_production_linear %>%
    dplyr::left_join(SR_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_linear)
}


#' @rdname vp_LM
determine_vp_linear_separate_replicates <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      for (sample in unique(SR_dataframe$Sample_Type)){
        for (rep in unique(SR_dataframe$Replicate)){
          DF <- SR_dataframe %>%
            dplyr::filter(.data$tag == combi_tag, .data$Population == virus, 
                          .data$Sample_Type == sample, .data$Replicate == rep)
          
          for (time in unique(DF$Time_Range)){
            DF2 <- DF %>%
              dplyr::filter(.data$Time_Range == time)
            
            linear_model <- summary(stats::lm(data = DF2, Count ~ as.numeric(Timepoint)))
            absolute_vp_values <- linear_model$coefficients[2] * (DF2$Timepoint[length(DF2$Timepoint)])
            
            result <- c(combi_tag, time, virus, sample, rep, linear_model$coefficients[2], absolute_vp_values, linear_model$coefficients[4], linear_model$r.squared)
            result_list[[length(result_list) + 1]] <- result
          }
        }
      }
    }
  }
  viral_production_linear <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_linear) <- c('tag', 'Time_Range', 'Population', 'Sample_Type',
                                         'Replicate', 'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')
  viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')] <- lapply(viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  viral_production_linear <- viral_production_linear %>%
    dplyr::left_join(SR_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_linear)
}


#' @rdname vp_LM
determine_vp_linear_average_replicates <- function(AVG_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(AVG_dataframe$tag)){
    for (virus in unique(AVG_dataframe[AVG_dataframe$Microbe == 'Viruses', ]$Population)){
      for (sample in unique(AVG_dataframe$Sample_Type)){
        DF <- AVG_dataframe %>%
          dplyr::filter(.data$tag == combi_tag, .data$Population == virus, .data$Sample_Type == sample)
        
        for (time in unique(DF$Time_Range)){
          DF2 <- DF %>%
            dplyr::filter(.data$Time_Range == time)
          
          linear_model <- summary(stats::lm(data = DF2, Mean ~ as.numeric(Timepoint)))
          absolute_vp_values <- linear_model$coefficients[2] * (DF2$Timepoint[length(DF2$Timepoint)])
          
          result <- c(combi_tag, time, virus, sample, linear_model$coefficients[2], absolute_vp_values, linear_model$coefficients[4], linear_model$r.squared)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_linear <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_linear) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 
                                         'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')
  viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')] <- lapply(viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  viral_production_linear <- viral_production_linear %>%
    dplyr::left_join(AVG_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_linear)
}


#' @rdname vp_LM
determine_vp_linear_LMER_model <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      DF <- SR_dataframe %>%
        dplyr::filter(.data$tag == combi_tag, .data$Population == virus)
      
      for (time in unique(DF$Time_Range)){
        DF2 <- DF %>%
          dplyr::filter(.data$Time_Range == time)
        
        DF_with_LMER_model <- vp_LMER_model(DF2)
        
        for (sample in unique(DF_with_LMER_model$Sample_Type)){
          linear_model <- summary(stats::lm(data = DF_with_LMER_model[DF_with_LMER_model$Sample_Type == sample, ], 
                                            Mean ~ as.numeric(Timepoint)))
          absolute_vp_values <- linear_model$coefficients[2] * (DF_with_LMER_model$Timepoint[length(DF_with_LMER_model$Timepoint)])
          
          result <- c(combi_tag, time, virus, sample, linear_model$coefficients[2], absolute_vp_values, linear_model$coefficients[4], linear_model$r.squared)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_linear <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_linear) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 
                                         'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')
  viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')] <- lapply(viral_production_linear[, c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared')], as.numeric)
  
  viral_production_linear <- viral_production_linear %>%
    dplyr::left_join(SR_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_linear)
}
