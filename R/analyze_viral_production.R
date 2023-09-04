#' Analyze viral production results
#' 
#' @description
#' Luef et al. (2009) created an online tool, `Viral Production Calculator`, that assesses and analyses viral production.
#' Next to calculating viral production and proportion of lysogenic bacteria in the sample, other parameters for estimating
#' virus-induced mortality are determined including lysis rate of bacteria, viral turnover time, organic carbon release and others.
#' 
#' Based on the original data and the output of [viralprod::calculate_viral_production], the same parameters are
#' calculated and added to the previous output. Another data frame with the description and units of each of the variables
#' is also generated. 
#' 
#' Reference to paper:
#' 
#' Luef, B., Luef, F., and Peduzzi, P. (2009). Online program ‘VIPCAL’ for calculating lytic viral production
#' and lysogenic cells based on a viral reduction approach. Environmental microbiology reports, 1(1):78–85.
#' doi:10.1111/j.1758-2229.2008.00008.x.
#'
#' @param vp_results Data frame with the viral production calculation results. 
#' @param data Data frame with the output of the flow cytometry.
#' @param original_abundances Data frame with the abundances of bacterial and virus population in the original sample. 
#' @param burst_sizes Vector with three different burst sizes. The burst size refers to the number of new viral particles released from an infected bacterial cell. 
#' @param bacterial_secondary_production Value for the bacterial secondary production, how much new bacterial biomass is produced as a result of bacterial growth and reproduction. 
#' @param nutrient_content_bacteria List with the amount of organic carbon, nitrogen and phosphor released by a bacteria, preferred a aquatic, North Sea bacteria.
#' @param nutrient_content_viruses List with the amount of organic carbon, nitrogen and phosphor released by a marine virus (bacteriophage).
#' @param write_csv If \code{TRUE}, the output data frames will be saved as csv files in the same folder of the viral production calculation.
#' If no csv files are wanted, set to \code{FALSE}. (Default = \code{TRUE})
#'
#' @return Data frame with the analyzed viral production results, different parameters for estimating virus-induced mortality are calculated. A second data frame with some more information about the variables will also be available in the global environment. 
#' @export
#' 
#' @name analyze_viral_production
#' @rdname analyze_viral_production
#'
#' @examples \dontrun{
#' # Setup
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' original_abundances_NJ2020 <- read.csv(system.file('extdata',
#' 'NJ2020_original_abundances.csv', package = "viralprod"))
#' 
#' calculate_viral_production(data_NJ2020_all, write_csv = F, SR_calc = F, BP_endpoint = F)
#' 
#' # Perform
#' analyze_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, 
#' original_abundances_2020) # Output files
#' 
#' analyze_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, 
#' original_abundances_2020, write_csv = F)
#' 
#' analyze_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, original_abundances_2020,
#' burst_sizes = c(15,30,50), bacterial_secondary_producition = 1000, 
#' nutrient_content_bacteria = list(C = 20, N = 15, P = 5),
#' nutrient_content_virus = list(C = 5, N = 3, P = 1)) # Set parameters, no default values
#' }
analyze_viral_production <- function(vp_results = data.frame(),
                                     data = data.frame(),
                                     original_abundances = data.frame(),
                                     burst_sizes = c(),
                                     bacterial_secondary_production = NULL,
                                     nutrient_content_bacteria = list(),
                                     nutrient_content_viruses = list(),
                                     write_csv = TRUE){
  ## 1. Checks
  # Check for valid output directory if csv file needs to be written
  # Output directory of calculation viral production is saved in global environment, save analyze csv file in same folder
  if (write_csv == T){
    if(!exists("output_dir", where = globalenv(), inherits = FALSE)){
      stop('There exists no output_dir value in the global environment, please define output folder to save analyzing results!')
      
    } else if(!file.exists(.GlobalEnv$output_dir)){
      warning('The output folder does not exists, results will be saved in given folder!')
      analyze_vp_results_path <- paste0(.GlobalEnv$output_dir, '/')
      
    } else {
      analyzed_vp_results_path <- paste0(.GlobalEnv$output_dir, '/')
    }
  }
  
  # Check if all three data frames aren't empty
  if (any(sapply(list(vp_results, data, original_abundances), function(df) rlang::is_empty(df)))){
    stop('Not able to proceed analyzing since one of the input data frames is empty!')
  }
  
  ## 2. Setup
  # Add bacterial abundance at T0 and the bacterial and viral abundances of the original sample to vp_results data frame
  B_T0_df <- vp_average_replicate_dataframe(data) %>%
    dplyr::filter(.data$Timepoint == 0, .data$Population == 'c_Bacteria', .data$Sample_Type == 'VP') %>%
    dplyr::select(dplyr::all_of(c('Station_Number', 'Timepoint', 'Population', 'Sample_Type', 'Mean'))) %>%
    dplyr::distinct()
  
  original_abundances_df <- original_abundances %>%
    dplyr::select(dplyr::all_of(c('Station_Number', 'Total_Bacteria', 'Total_Viruses')))
  
  .GlobalEnv$analyzed_vp_results_df <- vp_results %>%
    dplyr::left_join(dplyr::select(B_T0_df, 'Station_Number', 'Mean'), by = c('Station_Number')) %>%
    dplyr::rename(B_0 = "Mean") %>%
    dplyr::left_join(dplyr::select(original_abundances_df, 'Station_Number', 'Total_Bacteria', 'Total_Viruses'), by = c('Station_Number')) %>%
    dplyr::rename(B_OS = "Total_Bacteria",
                  V_OS = "Total_Viruses")
  
  # If no input for hyperparameters is given, default values are used
  if (rlang::is_empty(burst_sizes) | !rlang::is_bare_numeric(burst_sizes)){
    burst_sizes <- c(10,25,40)
    print('Default values used for burst size!')
  }
  
  if (!rlang::is_bare_numeric(bacterial_secondary_production)){
    bacterial_secondary_production <- 0.0027e6
    print('Default value used for bacterial secondary production!')
  }
  
  if (rlang::is_empty(nutrient_content_bacteria) | !all(sapply(nutrient_content_bacteria, rlang::is_bare_numeric))){
    # Article: Fagerbakke, KM & Heldal, Mikal & Norland, S. (1996). 
    # Content of Carbon, Nitrogen, Oxygen, Sulfur and Phosphorus in Native Aquatic and Cultured Bacteria. 
    # Aquatic Microbial Ecology - AQUAT MICROB ECOL. 10. 15-27. 10.3354/ame010015
    nutrient_content_bacteria <- list(C = 19e-15, N = 5e-15, P = 0.8e-15)
    print('Default values used for nutrient content of bacteria!')
  }
  
  if (rlang::is_empty(nutrient_content_viruses) | !all(sapply(nutrient_content_viruses, rlang::is_bare_numeric))){
    # Article: Jover, L., Effler, T., Buchan, A. et al. The elemental composition of virus particles: 
    # implications for marine biogeochemical cycles. Nat Rev Microbiol 12, 519–528 (2014). 10.1038/nrmicro3289
    C_V <- (1664612 / 6.022e23) * 12.01
    N_V <- (563058 / 6.022e23) * 14.01
    P_V <- (92428 / 6.022e23) * 30.97
    nutrient_content_viruses <- list(C = C_V, N = N_V, P = P_V) 
    print('Default values used for nutrient condent of viruses!')
  }
  
  ## 3. Analyze viral production results
  # 3.1 Correct values to represent original sample abundances
  .GlobalEnv$analyzed_vp_results_df <- .GlobalEnv$analyzed_vp_results_df %>%
    dplyr::mutate(c_VP = .data$VP * (.data$B_OS / .data$B_0),
                  c_abs_VP = .data$abs_VP * (.data$B_OS / .data$B_0),
                  c_VP_SE = abs(.data$c_VP) * (.data$VP_SE / abs(.data$VP)))
  
  # 3.2 Lytic and lysogenic viral production
  # 3.3 Lysis and lysogenic rate of bacteria
  # 3.4 Percentage of bacterial production lysed and bacterial loss per day
  for (bs in burst_sizes){
    col_name_3_2 <- paste0('P_Cells_BS_', bs)
    .GlobalEnv$analyzed_vp_results_df[[col_name_3_2]] <- .GlobalEnv$analyzed_vp_results_df$abs_VP * (100 / (.GlobalEnv$analyzed_vp_results_df$B_0 * bs))
    
    col_name_3_3 <- paste0('Rate_BS_', bs)
    .GlobalEnv$analyzed_vp_results_df[[col_name_3_3]] <- .GlobalEnv$analyzed_vp_results_df$c_VP / bs
    
    col_name_3_4_1 <- paste0('P_BP_Lysed_BS_', bs)
    .GlobalEnv$analyzed_vp_results_df[[col_name_3_4_1]] <- .GlobalEnv$analyzed_vp_results_df[[col_name_3_3]] / bacterial_secondary_production
    
    col_name_3_4_2 <- paste0('P_B_Loss_BS_', bs)
    .GlobalEnv$analyzed_vp_results_df[[col_name_3_4_2]] <- ((.GlobalEnv$analyzed_vp_results_df[[col_name_3_3]] * 100) / .GlobalEnv$analyzed_vp_results_df$B_OS) * 24
  }
  
  # 3.5 Viral turnover time
  .GlobalEnv$analyzed_vp_results_df <- .GlobalEnv$analyzed_vp_results_df %>%
    dplyr::mutate(V_TT = .data$c_VP / .data$V_OS)
  
  # 3.6 Nutrient release
  for (nutrient in 1:length(nutrient_content_viruses)){
    col_name_nutrient_virus <- paste0('DO', names(nutrient_content_viruses[nutrient]), '_V')
    .GlobalEnv$analyzed_vp_results_df[[col_name_nutrient_virus]] <- .GlobalEnv$analyzed_vp_results_df$VP * nutrient_content_viruses[[nutrient]]
  }
  
  for (bs in burst_sizes){
    for (nutrient in 1:length(nutrient_content_bacteria)){
      current_bs_column <- paste0('Rate_BS_', bs)
      col_name_nutrient_bacteria <- paste0('DO', names(nutrient_content_bacteria[nutrient]), '_B_BS_', bs)
      .GlobalEnv$analyzed_vp_results_df[[col_name_nutrient_bacteria]] <- .GlobalEnv$analyzed_vp_results_df[[current_bs_column]] * nutrient_content_bacteria[[nutrient]]
      
      # Add column with total nutrient release (bacteria + virus)
      col_name_nutrient_virus <- paste0('DO', names(nutrient_content_viruses[nutrient]), '_V')
      col_name_nutrient_total <- paste0('Total_DO', names(nutrient_content_bacteria[nutrient]), '_BS_', bs)
      .GlobalEnv$analyzed_vp_results_df[[col_name_nutrient_total]] <- .GlobalEnv$analyzed_vp_results_df[[col_name_nutrient_bacteria]] + .GlobalEnv$analyzed_vp_results_df[[col_name_nutrient_virus]]
    }
  }
  
  ## 4. Write output and legend if wanted
  .GlobalEnv$analyzed_vp_results_dictionary <- data.frame(
    Variable = colnames(.GlobalEnv$analyzed_vp_results_df),
    Unit = c('/', '/', 'm', 'h', '/', '/', '#VLP (virus-like particles)/mLh', '#VLP/mL', '/', '/', '/',
             '#VLP/mL', '#VLP/mL', '#VLP/mL', '#VLP/mLh', '#VLP/mL', '/', '%', '#VLP/mLh', '%', '%',
             '%', '#VLP/mLh', '%', '%', '%', '#VLP/mLh', '%', '%', '1/h', 'g C/mLh', 'g N/mLh', 'g P/mLh',
             'g C/mLh', 'g C/mLh', 'g N/mLh', 'g N/mLh', 'g P/mLh', 'g P/mLh', 'g C/mLh', 'g C/mLh', 'g N/mLh',
             'g N/mLh', 'g P/mLh', 'g P/mLh', 'g C/mLh', 'g C/mLh', 'g N/mLh', 'g N/mLh', 'g P/mLh', 'g P/mLh'),
    Description = c(
      'Location of the experiment',
      'Number of station where experiment is conducted',
      'Depth at which experiment is performed',
      'Timepoints of sampling data, expressed in a time range: starting at T = 0 until time of measuring X => T0_TX',
      'Population Types: c_Viruses covers the entire virus population, while c_V1, c_V2, and c_V3 represent subpopulations',
      'Sample Types: VP rfor lytic viral production, Diff for lysogenic viral production, VPC for both',
      'Viral production rate: the mean viral production rate during the current Time_Range',
      'Absolute viral production at current timepoint',
      'The standard error on the viral production rate',
      'R-squared value: goodness of fit of the linear regression model',
      'Calculation method of viral production',
      'Bacterial abundance at the beginning of the experiment (T0)',
      'Bacterial abundance in the original sample',
      'Viral abundance in the original sample',
      'Corrected viral production rate: viral production rate in the original sample',
      'Corrected absolute viral production: absolute viral production in the original sample',
      'Corrected standard error on the viral production rate',
      'Percentage of cells for given burst size: % lytically infected cells for VP samples, % lysogenic cells for Diff samples',
      'Rate of bacteria for given burst size: lysis rate of bacteria for VP samples, lysogenic rate of bacteria for Diff samples',
      'Percentage of bacterial production lysed: the quantity of bacterial biomass that undergoes lysis',
      'Percentage of bacterial loss per day: the rate at which bacteria are removed due to viral lysis',
      'Percentage of cells for given burst size: % lytically infected cells for VP samples, % lysogenic cells for Diff samples',
      'Rate of bacteria for given burst size: lysis rate of bacteria for VP samples, lysogenic rate of bacteria for Diff samples',
      'Percentage of bacterial production lysed: the quantity of bacterial biomass that undergoes lysis',
      'Percentage of bacterial loss per day: the rate at which bacteria are removed due to viral lysis',
      'Percentage of cells for given burst size: % lytically infected cells for VP samples, % lysogenic cells for Diff samples',
      'Rate of bacteria for given burst size: lysis rate of bacteria for VP samples, lysogenic rate of bacteria for Diff samples',
      'Percentage of bacterial production lysed: the quantity of bacterial biomass that undergoes lysis',
      'Percentage of bacterial loss per day: the rate at which bacteria are removed due to viral lysis',
      'Viral turnover time: time to replacte the current virus population by new viruses',
      'Dissolved organic carbon release of viruses',
      'Dissolved organic nitrogen release of viruses',
      'Dissolved organic phosphorous release of viruses',
      'Dissolved organic carbon release of bacteria for given burst size',
      'Total dissolved organic carbon release for given burst size',
      'Dissolved organic nitrogen release of bacteria for given burst size',
      'Total dissolved organic nitrogen release for given burst size', 
      'Dissolved organic phosphorous release of bacteria for given burst size',
      'Total dissolved organic phosphorous release of bacteria for given burst size',
      'Dissolved organic carbon release of bacteria for given burst size',
      'Total dissolved organic carbon release for given burst size',
      'Dissolved organic nitrogen release of bacteria for given burst size',
      'Total dissolved organic nitrogen release for given burst size', 
      'Dissolved organic phosphorous release of bacteria for given burst size',
      'Total dissolved organic phosphorous release of bacteria for given burst size',
      'Dissolved organic carbon release of bacteria for given burst size',
      'Total dissolved organic carbon release for given burst size',
      'Dissolved organic nitrogen release of bacteria for given burst size',
      'Total dissolved organic nitrogen release for given burst size', 
      'Dissolved organic phosphorous release of bacteria for given burst size',
      'Total dissolved organic phosphorous release of bacteria for given burst size'
    )
  )
  
  if (write_csv == T){
    utils::write.csv(.GlobalEnv$analyzed_vp_results_df, file.path(analyze_vp_results_path, 'vp_results_ANALYZED.csv'), row.names = F)
    utils::write.csv(.GlobalEnv$analyzed_vp_results_dictionary, file.path(analyze_vp_results_path, 'ANALYZED_Legend.csv'), row.names = F)
  }
}
