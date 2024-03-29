% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vp_list_of_methods.R
\name{vp_list_of_methods}
\alias{vp_list_of_methods}
\title{List of methods used for viral production calculation}
\usage{
vp_list_of_methods()
}
\value{
A list, with the different variants of linear regression and VIPCAL to calculate viral production.
}
\description{
Generate a comprehensive list of available variants for viral production calculation

More details on the different methods:
\itemize{
\item Variants of linear regression: \link{vp_methods_linear}
\item Variants of VIPCAL: \link{vp_methods_VIPCAL}
}
}
\examples{
\dontrun{
viralprod::vp_list_of_methods()

# If you want to run for example the first method
data_NJ2020_all <- read.csv(system.file('extdata', 
'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))

vp_check_populations(data_NJ2020_all)

list_of_methods[[1]](data_NJ2020_all)
}
}
