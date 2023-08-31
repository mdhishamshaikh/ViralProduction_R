#' Calculate viral production with VIPCAL
#' 
#' @description
#' `VIPCAL` uses the average of increments to determine the viral production rate. After determining peaks and valleys
#' in the count data, averaging the increases will result in the viral production rate. Again,
#' lytic viral production can be derived from the VP samples, the average of increments of the difference curve is needed
#' for lysogenic viral production. `VIPCAL-SE` goes one step further and takes the standard error into account, since this
#' has a major influence on the results. With the determination of peaks and valleys, only true increments are returned with
#' VIPCAL-SE. A peak/valley is only defined if there is no overlap between the standard errors. By adding this requirement,
#' sufficient differences in count values are needed so that there is a true increase in the viral count.
#' 
#' `determine_vp_VIPCAL_separate_replicates` uses a separate replicate treatment, distinguishing between replicates.  
#' 
#' `determine_vp_VIPCAL_average_replicates` uses an average replicate treatment, average over the replicates.
#' 
#' `determine_vp_VIPCAL_average_replicates_SE` uses VIPCAL-SE and an average replicate treatment, average over the replicates.
#' 
#' `determine_vp_VIPCAL_LMER_model` uses an average replicate treatment, average over the replicates is included in LMER model. 
#' Difference curve estimation by LMER model, See [viralprod::vp_LMER_model] for more details about LMER model.
#' 
#' `determine_vp_VIPCAL_LMER_model_SE` Uses VIPCAL-SE and an average replicate treatment, average over the replicates is included in LMER model. 
#' Difference curve estimation by LMER model, See [viralprod::vp_LMER_model] for more details about LMER model.
#'
#' @param SR_dataframe Data frame with the viral counts and time ranges, see [viralprod::vp_separate_replicate_dataframe] for more details.
#' @param AVG_dataframe Data frame with the viral counts and time ranges, see [viralprod::vp_average_replicate_dataframe] for more details.
#' 
#' @return Data frame with the viral production rate and the absolute viral production for each population at the different time points of the assay.
#' 
#' @name determine_vp_VIPCAL
#' @rdname vp_VIPCAL
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_all)
#' 
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020_all)
#' 
#' determine_vp_VIPCAL_separate_replicates(DF_SR)
#' 
#' determine_vp_VIPCAL_average_replicates(DF_AVG)
#' 
#' determine_vp_VIPCAL_average_replicates_SE(DF_AVG)
#' 
#' determine_vp_VIPCAL_LMER_model(DF_SR)
#' 
#' determine_vp_VIPCAL_LMER_model_SE(DF_SR)
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


#' @rdname vp_VIPCAL
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


#' @rdname vp_VIPCAL
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


#' @rdname vp_VIPCAL
determine_vp_VIPCAL_LMER_model <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      DF <- SR_dataframe %>%
        dplyr::filter(.data$tag == combi_tag, .data$Population == virus)
      
      for (time in unique(DF$Time_Range)){
        DF2 <- DF %>%
          dplyr::filter(.data$Time_Range == time)
        
        DF_with_LMER_model <- vp_LMER_model(DF2)
        
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


#' @rdname vp_VIPCAL
determine_vp_VIPCAL_LMER_model_SE <- function(SR_dataframe){
  result_list <- list()
  
  for (combi_tag in unique(SR_dataframe$tag)){
    for (virus in unique(SR_dataframe[SR_dataframe$Microbe == 'Viruses', ]$Population)){
      DF <- SR_dataframe %>%
        dplyr::filter(.data$tag == combi_tag, .data$Population == virus)
      
      for (time in unique(DF$Time_Range)){
        DF2 <- DF %>%
          dplyr::filter(.data$Time_Range == time)
        
        DF_with_LMER_model <- vp_LMER_model(DF2)
        
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
