# Config
expect_df <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true(all(c('Population', 'Count', 'Microbe') %in% names(x)))
}

# Perform
#.GlobalEnv$data_NJ2020 <- read.csv('./inst/NJ2020_subset.csv')

test_that("vp_separate_replicate_dataframe", {
  data_NJ2020 <- read.csv(paste0(path.package('viralprod'), "/inst/NJ2020_subset.csv"))
  
  vp_separate_replicate_dataframe(data_NJ2020) %>% expect_df() 
  vp_separate_replicate_dataframe(data_NJ2020) %>% ncol() %>% expect_equal(11)
  vp_separate_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% expect_df() 
  vp_separate_replicate_dataframe(data_NJ2020, add_timepoints = F) %>% ncol() %>% expect_equal(9)
  vp_separate_replicate_dataframe(data_NJ2020, keep_0.22_samples = T) %>% expect_df() 
  
  test <- vp_separate_replicate_dataframe(data_NJ2020, keep_0.22_samples = T)
  expect_true(any(test$Sample_Type == "0.22")) 
})
