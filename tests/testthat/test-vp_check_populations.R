# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate') %in% names(x)))
  expect_true(any(startsWith(names(x), 'c_')))
}

expect_output_check_populations <- function(x){
  expect_type(x, 'character')
}

# Perform
test_that("Right populations are analyzed", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
  expect_data_input(data_NJ2020_all)
  
  vp_check_populations(data_NJ2020_all)
  expect_output_check_populations(.GlobalEnv$populations_to_analyze)
  length(.GlobalEnv$populations_to_analyze) %>% expect_equal(length(which(startsWith(names(data_NJ2020_all),'c_'))))
  
  data_NJ2020_less <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_less_populations.csv', package = "viralprod"))
  expect_data_input(data_NJ2020_less)
  
  vp_check_populations(data_NJ2020_less)
  expect_output_check_populations(.GlobalEnv$populations_to_analyze)
  length(.GlobalEnv$populations_to_analyze) %>% expect_equal(length(which(startsWith(names(data_NJ2020_less),'c_'))))
  
  data_NJ2020_more <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_more_populations.csv', package = "viralprod"))
  expect_data_input(data_NJ2020_more)
  
  vp_check_populations(data_NJ2020_more)
  expect_output_check_populations(.GlobalEnv$populations_to_analyze)
  length(.GlobalEnv$populations_to_analyze) %>% expect_equal(length(which(startsWith(names(data_NJ2020_more),'c_'))))
  
  data_NJ2020_without_cViruses <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_without_cViruses.csv', package = "viralprod"))
  expect_data_input(data_NJ2020_without_cViruses)
  
  vp_check_populations(data_NJ2020_without_cViruses) %>% expect_error()
})
