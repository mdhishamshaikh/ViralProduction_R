# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Sample_Type', 'Timepoint') %in% names(x)))
}

expect_output_bacterial_endpoint <- function(x){
  expect_type(x, 'character')
}

expect_output_bacterial_endpoint_for_visual <- function(x){
  expect_type(x, 'integer')
}

# Perform
test_that("Bacterial endpoint can be determinated", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
  vp_check_populations(data_NJ2020_all)
  
  # Bacterial endpoint is determined per station/experiment
  subset_data <- subset(data_NJ2020_all, data_NJ2020_all$Station_Number == 2)
  expect_data_input(subset_data)
  
  vp_bacterial_endpoint(subset_data) %>% expect_output_bacterial_endpoint()
  vp_bacterial_endpoint(subset_data, visual = T) %>% expect_output_bacterial_endpoint_for_visual()
})
