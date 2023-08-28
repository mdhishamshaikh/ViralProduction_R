#' Determine viral production by VIPCAL with separate replicate treatment
#' 
#' Viral production rate is calculated by using a separate replicate treatment, distinguishing between
#' replicates. VIPCAL determines viral production rate as the average of increments between the viral counts
#' on different time points of the assay. First, peaks and valleys in the count data are calculated followed by
#' averaging the increases over time to get the viral production rate.
#'
#' @param SR_dataframe Dataframe with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe] for more details.
#'
#' @return Dataframe with the viral production rate and the absolute viral production for each population and replicate at given time range of the assay.
#' 
#' @name determine_vp_VIPCAL_separate_replicates 
#' @rdname vp_VPCL_SR
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
#' 
#' viral_production_VIPCAL_SR <- determine_vp_VIPCAL_separate_replicates(DF_SR)
#' }
determine_vp_VIPCAL_separate_replicates <- function(SR_dataframe){
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
            
            index_peaks <- vp_determine_peaks(c(+10e+10, DF2$Count, -10e+10))
            index_valleys <- vp_determine_valleys(c(+10e+10, DF2$Count, -10e+10))
            
            if (length(index_peaks) == 0){
              viral_production <- 0
              abs_vp <- 0
            }else {
              total_vp <- 0
              total_abs_vp <- 0
              
              for (index in 1:length(index_peaks)){
                viral_production_index <- (DF2$Count[index_peaks[index]] - DF2$Count[index_valleys[index]]) / (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
                total_vp <- total_vp + viral_production_index
                
                abs_vp_index <- DF2$Count[index_peaks[index]] - DF2$Count[index_valleys[index]]
                total_abs_vp <- total_abs_vp + abs_vp_index
              }
              viral_production <- total_vp / length(index_peaks)
              abs_vp <- total_abs_vp
            }
            result <- c(combi_tag, time, virus, sample, rep, viral_production, abs_vp)
            result_list[[length(result_list) + 1]] <- result
          }
        }
      }
    }
  }
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('tag', 'Time_Range', 'Population', 'Sample_Type',
                                         'Replicate', 'VP', 'abs_VP')
  viral_production_VIPCAL[, c('VP', 'abs_VP')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP')], as.numeric)
  
  viral_production_VIPCAL <- viral_production_VIPCAL %>%
    dplyr::left_join(SR_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_VIPCAL)
}


#' Determine viral production by VIPCAL with average replicate treatment
#' 
#' Viral production rate is calculated by using average replicate treatment, average over the replicates
#' is taken. VIPCAL determines viral production rate as the average of increments between the viral counts
#' on different time points of the assay. First, peaks and valleys in the count data are calculated followed by
#' averaging the increases over time to get the viral production rate.
#'
#' @param AVG_dataframe Dataframe with the viral counts and time ranges, see [viralprod::vp_average_replicate_dataframe] for more details.
#'
#' @return Dataframe with the viral production rate and the absolute viral production for each population averaged over the replicates at given time range of the assay.
#' 
#' @name determine_vp_VIPCAL_average_replicates
#' @rdname vp_VPCL_AR
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020)
#' 
#' viral_production_VIPCAL_AR <- determine_vp_VIPCAL_average_replicates(DF_AVG)
#' }
determine_vp_VIPCAL_average_replicates <- function(AVG_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(AVG_dataframe$tag)){
    for (virus in unique(AVG_dataframe[AVG_dataframe$Microbe == 'Viruses', ]$Population)){
      for (sample in unique(AVG_dataframe$Sample_Type)){
        DF <- AVG_dataframe %>%
          dplyr::filter(.data$tag == combi_tag, .data$Population == virus, .data$Sample_Type == sample)
        
        for (time in unique(DF$Time_Range)){
          DF2 <- DF %>%
            dplyr::filter(.data$Time_Range == time)
          
          index_peaks <- vp_determine_peaks(c(+10e+10, DF2$Mean, -10e+10))
          index_valleys <- vp_determine_valleys(c(+10e+10, DF2$Mean, -10e+10))
          
          if (length(index_peaks) == 0){
            viral_production <- 0
            abs_vp <- 0
          }else {
            total_vp <- 0
            total_abs_vp <- 0
            
            for (index in 1:length(index_peaks)){
              viral_production_index <- (DF2$Mean[index_peaks[index]] - DF2$Mean[index_valleys[index]]) / (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
              total_vp <- total_vp + viral_production_index
              
              abs_vp_index <- DF2$Mean[index_peaks[index]] - DF2$Mean[index_valleys[index]]
              total_abs_vp <- total_abs_vp + abs_vp_index
            }
            viral_production <- total_vp / length(index_peaks)
            abs_vp <- total_abs_vp
          }
          result <- c(combi_tag, time, virus, sample, viral_production, abs_vp)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP')
  viral_production_VIPCAL[, c('VP', 'abs_VP')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP')], as.numeric)
  
  viral_production_VIPCAL <- viral_production_VIPCAL %>%
    dplyr::left_join(AVG_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_VIPCAL)
}


#' Determine viral production by VIPCAL-SE with average replicate treatment
#' 
#' Viral production rate is calculated by using average replicate treatment, average over the replicates 
#' is taken. VIPCAL-SE determines viral production rate as the average of increments between the viral counts
#' on different time points of the assay. VIPCAL-SE works the same as VIPCAL but takes the standard error into
#' account. When determining peaks and valleys, only true increments are kept. A peak/valley is only defined if
#' there is no overlap between the standard errors. By doing this, sufficient differences in count values, 
#' representing possible increments, are preserved.
#'
#' @param AVG_dataframe Dataframe with the viral counts and time ranges, see [viralprod::vp_average_replicate_dataframe] for more details.
#'
#' @return Dataframe with the viral production rate and the absolute viral production for each population averaged over the replicates at given time range of the assay. Standard error is taken into account.
#' 
#' @name determine_vp_VIPCAL_average_replicates_SE 
#' @rdname vp_VPCL_AR_SE
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020)
#' 
#' viral_production_VIPCAL_AR_SE <- determine_vp_VIPCAL_average_replicates_SE(DF_AVG)
#' }
determine_vp_VIPCAL_average_replicates_SE <- function(AVG_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(AVG_dataframe$tag)){
    for (virus in unique(AVG_dataframe[AVG_dataframe$Microbe == 'Viruses', ]$Population)){
      for (sample in unique(AVG_dataframe$Sample_Type)){
        DF <- AVG_dataframe %>%
          dplyr::filter(.data$tag == combi_tag, .data$Population == virus, .data$Sample_Type == sample)
        
        for (time in unique(DF$Time_Range)){
          DF2 <- DF %>%
            dplyr::filter(.data$Time_Range == time)
          
          index_peaks <- vp_determine_peaks_with_se(c(+10e+10, DF2$Mean, -10e+10),
                                                    c(0, DF2$SE, 0))
          index_valleys <- vp_determine_valleys_with_se(c(+10e+10, DF2$Mean, -10e+10),
                                                        c(0, DF2$SE, 0))
          
          if (length(index_peaks) == 0){
            viral_production <- 0
            abs_vp <- 0
            se <- 0
          }else {
            total_vp <- 0
            total_abs_vp <- 0
            total_se <- 0
            
            for (index in 1:length(index_peaks)){
              viral_production_index <- (DF2$Mean[index_peaks[index]] - DF2$Mean[index_valleys[index]]) / (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
              total_vp <- total_vp + viral_production_index
              
              abs_vp_index <- DF2$Mean[index_peaks[index]] - DF2$Mean[index_valleys[index]]
              total_abs_vp <- total_abs_vp + abs_vp_index
              
              se_index <- (DF2$SE[index_peaks[index]] + DF2$SE[index_valleys[index]]) / (DF2$Timepoint[index_peaks[index]] - DF2$Timepoint[index_valleys[index]])
              total_se <- total_se + se_index
            }
            viral_production <- total_vp / length(index_peaks)
            abs_vp <- total_abs_vp
            se <- total_se / length(index_peaks)
          }
          result <- c(combi_tag, time, virus, sample, viral_production, abs_vp, se)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP', 'VP_SE')
  viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')], as.numeric)
  
  viral_production_VIPCAL <- viral_production_VIPCAL %>%
    dplyr::left_join(AVG_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_VIPCAL)
}


#' Determine viral production by VIPCAL with difference curve estimation by LMER model
#' 
#' Viral production rate is calculated by using average replicate treatment, average over the replicates 
#' is taken and the difference curve is estimated by the LMER model instead by subtraction. See [viralprod::LMER_model] for 
#' more details about LMER model. VIPCAL determines viral production rate as the average of increments between the viral counts
#' on different time points of the assay. First, peaks and valleys in the count data are calculated followed by
#' averaging the increases over time to get the viral production rate.
#'
#' @param SR_dataframe Dataframe with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe] for more details.
#'
#' @return Dataframe with the viral production rate and the absolute viral production for each population averaged over the replicates at given time range of the assay. Difference curve estimated by LMER model.
#' 
#' @name determine_vp_VIPCAL_LMER_model
#' @rdname vp_VPCL_LMER
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
#' 
#' viral_production_VIPCAL_LMER <- determine_vp_VIPCAL_LMER_model(DF_SR)
#' }
determine_vp_VIPCAL_LMER_model <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      DF <- SR_dataframe %>%
        dplyr::filter(.data$tag == combi_tag, .data$Population == virus)
      
      for (time in unique(DF$Time_Range)){
        DF2 <- DF %>%
          dplyr::filter(.data$Time_Range == time)
        
        DF_with_LMER_model <- LMER_model(DF2)
        
        for(sample in unique(DF_with_LMER_model$Sample_Type)){
          DF3 <- DF_with_LMER_model %>%
            dplyr::filter(.data$Sample_Type == sample)
          
          index_peaks <- vp_determine_peaks(c(+10e+10, DF3$Mean, -10e+10))
          index_valleys <- vp_determine_valleys(c(+10e+10, DF3$Mean, -10e+10))
          
          if (length(index_peaks) == 0){
            viral_production <- 0
            abs_vp <- 0
          }else {
            total_vp <- 0
            total_abs_vp <- 0
            
            for (index in 1:length(index_peaks)){
              viral_production_index <- (DF3$Mean[index_peaks[index]] - DF3$Mean[index_valleys[index]]) / (DF3$Timepoint[index_peaks[index]] - DF3$Timepoint[index_valleys[index]])
              total_vp <- total_vp + viral_production_index
              
              abs_vp_index <- DF3$Mean[index_peaks[index]] - DF3$Mean[index_valleys[index]]
              total_abs_vp <- total_abs_vp + abs_vp_index
            }
            viral_production <- total_vp / length(index_peaks)
            abs_vp <- total_abs_vp
          }
          result <- c(combi_tag, time, virus, sample, viral_production, abs_vp)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP')
  viral_production_VIPCAL[, c('VP', 'abs_VP')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP')], as.numeric)
  
  viral_production_VIPCAL <- viral_production_VIPCAL %>%
    dplyr::left_join(SR_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_VIPCAL)
}


#' Determine viral production by VIPCAL-SE with difference curve estimation by LMER model
#' 
#' VIral production rate is calculated by using average replicate treatment, average over the replicates 
#' is taken and the difference curve is estimated by the LMER model instead by subtraction. See [viralprod::LMER_model] for 
#' more details about LMER model. VIPCAL-SE determines viral production rate as the average of increments between the viral counts
#' on different time points of the assay. VIPCAL-SE works the same as VIPCAL but takes the standard error into
#' account. When determining peaks and valleys, only true increments are kept. A peak/valley is only defined if
#' there is no overlap between the standard errors. By doing this, sufficient differences in count values, 
#' representing possible increments, are preserved.
#'
#' @param SR_dataframe Dataframe with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe] for more details.
#'
#' @return Dataframe with the viral production rate and the absolute viral production for each population averaged over the replicates at given time range of the assay. Difference curve estimated by LMER model and standard error is taken into account.
#' 
#' @name determine_vp_VIPCAL_LMER_model_SE
#' @rdname vp_VPCL_LMER_SE
#'
#' @examples \dontrun{
#' data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
#' 
#' viral_production_VIPCAL_LMER_SE <- determine_vp_VIPCAL_LMER_model_SE(DF_SR)
#' }
determine_vp_VIPCAL_LMER_model_SE <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      DF <- SR_dataframe %>%
        dplyr::filter(.data$tag == combi_tag, .data$Population == virus)
      
      for (time in unique(DF$Time_Range)){
        DF2 <- DF %>%
          dplyr::filter(.data$Time_Range == time)
        
        DF_with_LMER_model <- LMER_model(DF2)
        
        for(sample in unique(DF_with_LMER_model$Sample_Type)){
          DF3 <- DF_with_LMER_model %>%
            dplyr::filter(.data$Sample_Type == sample)
          
          index_peaks <- vp_determine_peaks_with_se(c(+10e+10, DF3$Mean, -10e+10),
                                                    c(0, DF3$SE, 0))
          index_valleys <- vp_determine_valleys_with_se(c(+10e+10, DF3$Mean, -10e+10),
                                                        c(0, DF3$SE, 0))
          
          if (length(index_peaks) == 0){
            viral_production <- 0
            abs_vp <- 0
            se <- 0
          }else {
            total_vp <- 0
            total_abs_vp <- 0
            total_se <- 0
            
            for (index in 1:length(index_peaks)){
              viral_production_index <- (DF3$Mean[index_peaks[index]] - DF3$Mean[index_valleys[index]]) / (DF3$Timepoint[index_peaks[index]] - DF3$Timepoint[index_valleys[index]])
              total_vp <- total_vp + viral_production_index
              
              abs_vp_index <- DF3$Mean[index_peaks[index]] - DF3$Mean[index_valleys[index]]
              total_abs_vp <- total_abs_vp + abs_vp_index
              
              se_index <- (DF3$SE[index_peaks[index]] + DF3$SE[index_valleys[index]]) / (DF3$Timepoint[index_peaks[index]] - DF3$Timepoint[index_valleys[index]])
              total_se <- total_se + se_index
            }
            viral_production <- total_vp / length(index_peaks)
            abs_vp <- total_abs_vp
            se <- total_se / length(index_peaks)
          }
          result <- c(combi_tag, time, virus, sample, viral_production, abs_vp, se)
          result_list[[length(result_list) + 1]] <- result
        }
      }
    }
  }
  viral_production_VIPCAL <- data.frame(t(sapply(result_list, c)))
  colnames(viral_production_VIPCAL) <- c('tag', 'Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP', 'VP_SE')
  viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')] <- lapply(viral_production_VIPCAL[, c('VP', 'abs_VP', 'VP_SE')], as.numeric)
  
  viral_production_VIPCAL <- viral_production_VIPCAL %>%
    dplyr::left_join(SR_dataframe %>%
                       dplyr::select(dplyr::all_of(c('tag', 'Location', 'Station_Number', 'Depth'))), 
                     by = 'tag', relationship = 'many-to-many') %>%
    dplyr::distinct()
  
  return(viral_production_VIPCAL)
}
