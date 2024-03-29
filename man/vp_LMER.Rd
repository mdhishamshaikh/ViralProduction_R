% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vp_LMER_model.R
\name{vp_LMER_model}
\alias{vp_LMER_model}
\title{Linear Mixed-Effects Model}
\usage{
vp_LMER_model(DF)
}
\arguments{
\item{DF}{Data frame with the viral counts and time ranges, see \link{vp_separate_replicate_dataframe} for more details.}
}
\value{
Data frame with the mean viral count for VP samples, VPC samples and the difference (DIFF samples) for each population at the different time points of the assay.
}
\description{
The viral reduction assay has two types of samples: VP and VPC samples.
In VP samples, count of bacteriophages in the lytic phase
can be measured (\verb{lytic viral production}), while in VPC samples, count of bacteriophages in both the lytic and
lysogenic phase can be measured (\verb{lytic + lysogenic viral production}) since treatment with antibiotic \code{Mitomycin-C}
forces lysogenic bacteriophages to go into the lytic phase. To retrieve the lysogenic viral production, a difference
curve is used. This difference curve can be derived in two ways: \code{subtraction} or the application of a
\code{Linear Mixed-Effects Model (LMER)}. In the subtraction approach, viral abundance in VP samples is subtracted from
that in VPC samples. On the other hand, the LMER model incorporates both fixed and random effect terms.
It performs a maximum likelihood estimation by considering the fixed effect terms, such as sample type and timepoint,
while also accounting for the random variability among different replicates.

See \link[lme4:lmer]{lme4::lmer} for more details on the LMER model.
}
\examples{
\dontrun{
data_NJ2020_all <- read.csv(system.file('extdata', 
'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))

vp_check_populations(data_NJ2020_all)

DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)

vp_LMER_model(DF_SR)
}
}
