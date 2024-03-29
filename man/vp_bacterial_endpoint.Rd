% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vp_bacterial_endpoint.R
\name{vp_bacterial_endpoint}
\alias{vp_bacterial_endpoint}
\title{Determine bacterial endpoint of assay}
\usage{
vp_bacterial_endpoint(data, visual = FALSE)
}
\arguments{
\item{data}{Data frame with the output of the flow cytometry.}

\item{visual}{If \code{FALSE}, a character with the time range to stop the assay is returned. The character value is from the form: T0_TX. If \code{TRUE}, an integer with the index of the time point to stop the assay is returned. (Default = \code{FALSE})}
}
\value{
A character or integer value defining the time range to stop the assay to retrieve comparable results.
}
\description{
In the VP samples, an increase of collision rates between the bacteriophages (viruses) and bacteria is noticed.
This is probably due to the net increase in bacterial growth that is established during the assay. In VPC samples,
this phenomenon isn't presented since treatment with antibiotic \code{Mitomycin-C} inhibits the
growth of bacteria. Such increased collision rates can result in an overestimated lytic viral production
and a negative lysogenic viral production.

To mitigate this, we define the bacterial endpoint as the point at which the bacterial generation time becomes
less than the assay duration, indication the appropriate moment to stop the assay and obtain less biased results.
The generation time is calculated based of the net increase of the total bacteria population in
the VP samples of the assay.
}
\examples{
\dontrun{
data_NJ2020_all <- read.csv(system.file('extdata', 
'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))

vp_check_populations(data_NJ2020_all)

# Bacterial endpoint is determined per station/experiment
subset_data <- subset(data_NJ2020_all, data_NJ2020_all$Station_Number == 2)

vp_bacterial_endpoint(subset_data)
vp_bacterial_endpoint(subset_data, visual = T)
}
}
