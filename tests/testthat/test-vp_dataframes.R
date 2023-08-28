# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate',
                    'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3') %in% names(x)))
}

expect_df_output_separate_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('tag', 'Population', 'Count', 'Microbe') %in% names(x)))
}

expect_df_output_average_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('tag', 'Population', 'n', 'Mean', 'SE', 'Microbe', 'Subgroup') %in% names(x)))
}

# Perform
test_that("vp_separate_replicate_dataframe", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  expect_data_input(data_NJ2020)
  
  vp_separate_replicate_dataframe(data_NJ2020) %>% expect_df_output_separate_replicates() 
  vp_separate_replicate_dataframe(data_NJ2020) %>% ncol() %>% expect_equal(12)
  vp_separate_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% expect_df_output_separate_replicates() 
  vp_separate_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% ncol() %>% expect_equal(10)
  vp_separate_replicate_dataframe(data_NJ2020, keep_0.22_samples = T) %>% expect_df_output_separate_replicates() 
  
  test <- vp_separate_replicate_dataframe(data_NJ2020, keep_0.22_samples = T)
  expect_true(any(test$Sample_Type == "0.22")) 
})

test_that("vp_average_replicate_dataframe",{
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  
  vp_average_replicate_dataframe(data_NJ2020) %>% expect_df_output_average_replicates()
  vp_average_replicate_dataframe(data_NJ2020) %>% ncol() %>% expect_equal(14)
  vp_average_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% expect_df_output_average_replicates()
  vp_average_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% ncol() %>% expect_equal(12)
})
