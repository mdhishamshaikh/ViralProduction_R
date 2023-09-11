# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('viralprod','data.frame'))
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate') %in% names(x)))
  expect_true(any(startsWith(names(x), 'c_')))
}

expect_original_abundances_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('viralprod_analyze','data.frame'))
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
  x <- new_viralprod_class(data_NJ2020_all)
  expect_data_input(x)
  
  original_abundances_NJ2020 <- read.csv(system.file('extdata','NJ2020_original_abundances.csv', 
                                                     package = "viralprod"))
  y <- new_viralprod_class_2(original_abundances_NJ2020)
  expect_original_abundances_input(y)
  
  calculate_viral_production(x, SR_calc = F, BP_endpoint = F, write_output = F)
  expect_input_vp_results(.GlobalEnv$vp_results_output_df)
  
  visualize_viral_production(x = data_NJ2020_all, vp_results = vp_results_output_df,
                             original_abundances = original_abundances_NJ2020, write_output = F) %>% expect_message()
  
  visualize_viral_production(x, vp_results = .GlobalEnv$vp_results_output_df,
                             original_abundances = y, write_output = F)
  expect_equal(length(.GlobalEnv$plot_list), 11)
  
  visualize_viral_production(x, vp_results = .GlobalEnv$vp_results_output_df,
                             original_abundances = y, write_output = T, output_dir = '') %>% expect_error()
})
