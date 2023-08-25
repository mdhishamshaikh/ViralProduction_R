# Config
expect_df_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate',
                    'c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3') %in% names(x)))
}

expect_SRdf_output <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('Population', 'Count', 'Microbe') %in% names(x)))
}

expect_AVGdf_output <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('Population', 'n', 'Mean', 'SE', 'Microbe', 'Subgroup') %in% names(x)))
}

# Perform

test_that("vp_separate_replicate_dataframe", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  expect_df_input(data_NJ2020)
  
  vp_separate_replicate_dataframe(data_NJ2020) %>% expect_SRdf_output() 
  vp_separate_replicate_dataframe(data_NJ2020) %>% ncol() %>% expect_equal(11)
  vp_separate_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% expect_SRdf_output() 
  vp_separate_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% ncol() %>% expect_equal(9)
  vp_separate_replicate_dataframe(data_NJ2020, keep_0.22_samples = T) %>% expect_SRdf_output() 
  
  test <- vp_separate_replicate_dataframe(data_NJ2020, keep_0.22_samples = T)
  expect_true(any(test$Sample_Type == "0.22")) 
})

test_that("vp_average_replicate_dataframe",{
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  expect_df_input(data_NJ2020)
  
  vp_average_replicate_dataframe(data_NJ2020) %>% expect_AVGdf_output()
  vp_average_replicate_dataframe(data_NJ2020) %>% ncol() %>% expect_equal(13)
  vp_average_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% expect_AVGdf_output()
  vp_average_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% ncol() %>% expect_equal(11)
})
