# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', !!!populations_to_analyze) %in% names(x)))
}

expect_df_output_separate_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('tag', 'Population', 'Count', 'Microbe') %in% names(x)))
}

expect_df_output_average_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('tag', 'Population', 'n', 'Mean', 'SE', 'Microbe') %in% names(x)))
}

# Perform
test_that("Correct dataframe can be made", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
  vp_check_populations(data_NJ2020_all)
  expect_data_input(data_NJ2020_all)
  
  vp_separate_replicate_dataframe(data_NJ2020_all) %>% expect_df_output_separate_replicates()
  vp_separate_replicate_dataframe(data_NJ2020_all) %>% ncol() %>% expect_equal(12)
  SR_dataframe <- vp_separate_replicate_dataframe(data_NJ2020_all)
  expect_true(all(unique(SR_dataframe$Population) %in% populations_to_analyze))
  
  vp_separate_replicate_dataframe(data_NJ2020_all, add_timepoints = F) %>% expect_df_output_separate_replicates()
  vp_separate_replicate_dataframe(data_NJ2020_all, add_timepoints = F) %>% ncol() %>% expect_equal(10)
  
  vp_separate_replicate_dataframe(data_NJ2020_all, keep_0.22_samples = T) %>% expect_df_output_separate_replicates()
  test <- vp_separate_replicate_dataframe(data_NJ2020_all, keep_0.22_samples = T)
  expect_true(any(test$Sample_Type == "0.22"))
  
  vp_average_replicate_dataframe(data_NJ2020_all) %>% expect_df_output_average_replicates()
  vp_average_replicate_dataframe(data_NJ2020_all) %>% ncol() %>% expect_equal(13)
  AVG_dataframe <- vp_average_replicate_dataframe(data_NJ2020_all)
  expect_true(all(unique(AVG_dataframe$Population) %in% populations_to_analyze))
  
  vp_average_replicate_dataframe(data_NJ2020_all, add_timepoints = F) %>% expect_df_output_average_replicates()
  vp_average_replicate_dataframe(data_NJ2020_all, add_timepoints = F) %>% ncol() %>% expect_equal(11)
})
