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
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  subset_data <- subset(data_NJ2020, data_NJ2020$Station_Number == 2) # Bacterial endpoint is determined per station/experiment
  expect_data_input(subset_data)
  
  vp_bacterial_endpoint(subset_data) %>% expect_output_bacterial_endpoint()
  vp_bacterial_endpoint(subset_data, visual = T) %>% expect_output_bacterial_endpoint_for_visual()
})
