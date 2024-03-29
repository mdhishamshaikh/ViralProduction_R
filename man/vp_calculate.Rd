% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vp_calculate.R
\name{vp_calculate}
\alias{vp_calculate}
\alias{vp_calculate.default}
\alias{vp_calculate.viralprod}
\title{Wrapper function to calculate viral production}
\usage{
vp_calculate(x, ...)

\method{vp_calculate}{default}(x, ...)

\method{vp_calculate}{viralprod}(
  x,
  ...,
  methods = c(1:12),
  SR_calc = TRUE,
  BP_endpoint = TRUE,
  write_output = TRUE,
  output_dir = ""
)
}
\arguments{
\item{x}{Data frame with the output of the flow cytometry, has to have the \code{viralprod} class.}

\item{...}{Arguments passed on to the next function.}

\item{methods}{Integer vector with the indexes of \code{list_of_methods}. Indexes determine which methods will be performed within the viral production calculation.}

\item{SR_calc}{If \code{TRUE}, separate replicate treatment results will be stored and saved as an separate data frame.
Set to \code{FALSE}, if separate replicate treatment results are not wanted. (Default = \code{TRUE})}

\item{BP_endpoint}{If \code{TRUE}, the bacterial endpoint will be taken into account and only those results will be saved in a new
data frame. If not wanted, set to \code{FALSE}. (Default = \code{TRUE})}

\item{write_output}{If \code{TRUE}, the output data frames will be saved as csv files in a folder specified by \code{output_dir}.
If no csv files are wanted, set to \code{FALSE}. (Default = \code{TRUE})}

\item{output_dir}{String that refers to the location of folder to save the data frames as csv files.}
}
\value{
Depending on setting of parameters, different data frames with the viral production calculation will be available in the global environment.
}
\description{
Wrapper function that performs viral production calculation. Given the output of the flow cytometry step, viral production
is calculated for all methods of linear regression and VIPCAL. Based on the parameters, different data frames will be available
in the global environment and written as csv files.
\enumerate{
\item \code{vp_results_ALL.csv}: Contains the viral production results for all samples.
\item \code{vp_results_T24.csv}: Contains the viral production results for all samples at the end of the assay (T0_T24).
\item \code{vp_results_SR.csv}: Contains the viral production results for all samples of the separate replicate treatment, VP and VPC samples only.
\item \code{vp_results_BP.csv}: Contains the viral production results for all samples with the bacterial endpoint into account. Results for that time range will be stored in data frame.
}

More details about the calculation methods:
\itemize{
\item Linear regression variants: \link{vp_methods_linear}
\item VIPCAL variants: \link{vp_methods_VIPCAL}
}
}
\examples{
\dontrun{
# Setup
data_NJ2020_all <- read.csv(system.file('extdata', 
'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))

# Perform
# Add S3 class
x <- vp_class_count_data(data_NJ2020_all)

# Default method
vp_calculate(x = data_NJ2020_all, write_output = F)

# S3 class, viralprod, method
# Write output files
vp_calculate(x, 
output_dir = paste0(system.file(“extdata”, package = “viralprod”), 
“/NJ2020_vp_results”))

# No output files
vp_calculate(x, write_output = F)

# No bacterial endpoint and separate replicate treatment results
vp_calculate(x, write_output = F, SR_calc = F, BP_endpoint = F)

# Sub selection of the methods
vp_calculate(x, write_output = F, methods = c(2,3,8,12))
}
}
