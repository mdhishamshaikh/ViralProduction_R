# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate') %in% names(x)))
  expect_true(any(startsWith(names(x), 'c_')))
}

expect_original_abundances_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Station_Number', 'Total_Bacteria', 'Total_Viruses') %in% names(x)))
}

expect_input_vp_results <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared', 'VP_Method') %in% names(x)))
}

expect_output_analyzed_vp <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('B_0', 'B_OS', 'V_OS', 'V_TT') %in% names(x)))
  expect_true(any(startsWith(names(x), 'c_')))
  expect_true(any(startsWith(names(x), 'P_')))
  expect_true(any(startsWith(names(x), 'Rate_')))
  expect_true(any(startsWith(names(x), 'DO')))
  expect_true(any(startsWith(names(x), 'Total_DO')))
}

expect_output_analyzed_vp_dictionary <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Variable', 'Unit', 'Description') %in% names(x)))
}


# Perform
test_that("Viral production can be analyzed", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', 
                                          package = "viralprod"))
  expect_data_input(data_NJ2020_all)
  
  original_abundances_NJ2020 <- read.csv(system.file('extdata','NJ2020_original_abundances.csv', 
                                                     package = "viralprod"))
  expect_original_abundances_input(original_abundances_NJ2020)
  
  calculate_viral_production(data_NJ2020_all, SR_calc = F, BP_endpoint = F, write_output = F)
  expect_input_vp_results(.GlobalEnv$vp_results_output_df)
  
  analyze_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, original_abundances_NJ2020, 
                           write_output = F)
  expect_output_analyzed_vp(.GlobalEnv$analyzed_vp_results_df)
  expect_output_analyzed_vp_dictionary(.GlobalEnv$analyzed_vp_results_dictionary)
  
  analyze_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, original_abundances_NJ2020,
                           write_output = T, output_dir = '') %>% expect_error()
})
