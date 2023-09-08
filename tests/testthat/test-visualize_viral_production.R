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

# Perform
test_that("Viral production can be visualized", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', 
                                          package = "viralprod"))
  expect_data_input(data_NJ2020_all)
  
  original_abundances_NJ2020 <- read.csv(system.file('extdata','NJ2020_original_abundances.csv', 
                                                     package = "viralprod"))
  expect_original_abundances_input(original_abundances_NJ2020)
  
  calculate_viral_production(data_NJ2020_all, write_output = F)
  expect_input_vp_results(.GlobalEnv$vp_results_output_df)
  
  visualize_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, original_abundances_NJ2020, 
                             write_output = T, output_dir = '') %>% expect_error()
  visualize_viral_production(.GlobalEnv$vp_results_output_df, data_NJ2020_all, original_abundances_NJ2020, write_output = F)
  expect_equal(length(.GlobalEnv$plot_list), 11)
})
