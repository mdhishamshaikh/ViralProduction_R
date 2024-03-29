% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vp_visuals.R
\name{vp_visuals}
\alias{vp_visuals}
\alias{plot_overview_counts_over_time}
\alias{plot_collision_rates}
\alias{plot_comparison_methods}
\alias{plot_VIPCAL_vs_VIPCAL_SE}
\alias{plot_percentage_cells}
\alias{plot_nutrient_release}
\title{Visualizations of viral production data}
\usage{
plot_overview_counts_over_time(data)

plot_collision_rates(data, original_abundances)

plot_comparison_methods(vp_results)

plot_VIPCAL_vs_VIPCAL_SE(vp_results)

plot_percentage_cells(analyzed_vp_results_bacterial_endpoint)

plot_nutrient_release(analyzed_vp_results_T0_T24)
}
\arguments{
\item{data}{Data frame with the output of the flow cytometry.}

\item{original_abundances}{Data frame with the abundances of bacterial and virus population in the original sample.}

\item{vp_results}{Data frame with the viral production calculation results, available in global environment
after running \code{calculate_viral_production}.}

\item{analyzed_vp_results_bacterial_endpoint}{Data frame with the analyzed viral production results, available in global environment
after running \code{analyze_viral_production}. Important is that the analyzing step is ran on the \code{vp_results_ouput_BP_df}
data frame from \code{calculate_viral_production} instead of the \code{vp_results_ouput_BP_df} so that the bacterial
endpoint is taken into account.}

\item{analyzed_vp_results_T0_T24}{Data frame with the analyzed viral production results, available in global environment
after running \code{analyze_viral_production}. Important is that the analyzing step is ran on the \code{vp_results_ouput_T24_df}
data frame from \code{calculate_viral_production} instead of the \code{vp_results_ouput_BP_df} so that the nutrient release
at the end of the assay is calculated.}
}
\value{
Plot objects will be stored in variable \code{plot_list} in the global environment.
}
\description{
A major step in data analyses is \verb{data visualization}. Some suggestions to visualize the original viral count data,
calculated viral production data, or analyzed viral production data.

\code{plot_overview_counts_over_time}: plots the bacterial and viral counts, retrieved from the flow cytometry, in
function of the different time ranges of the assay. The bacterial endpoint is also highlighted in pink.

\code{plot_collision_rates}: plots the difference in collision rates over time between VP and VPC samples. The
collision rate is defined as the frequency of physical encounter between a bacteriophage and a bacterial cell in the sample.

\code{plot_comparison_methods}: compares the different methods to calculate viral production. The calculated viral
production is compared between the linear regression variants and the VIPCAL variants, also an overview of
all the methods is made. Next to that, a plot that compares linear regression with VIPCAL and VIPCAL-SE is also produced.

\code{plot_VIPCAL_vs_VIPCAL_SE}: a separate visualization for comparison between VIPCAL and VIPCAL-SE. The calculated viral
production values of both methods are compared plus a \verb{Robust Graphical Method (ROGME)} is used for the comparison between the two.

\code{plot_percentage_cells}: plots the percentage of lytically infected and lysogenic cells in function of the
different burst sizes. Important to note, that this is the percentage of cells at the bacterial endpoint of the assay,
so the point where ideally the assay is stopped to retrieve less biased results between samples.

\code{plot_nutrient_release}: plots the nutrient release for organic carbon, nitrogen and phosphor in function of the burst
size at the end of the assay (T0_T24).
}
\examples{
\dontrun{
# Setup
data_NJ2020_all <- read.csv(system.file('extdata', 
'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))

original_abundances_NJ2020 <- read.csv(system.file('extdata',
'NJ2020_original_abundances.csv', package = "viralprod"))

calculate_viral_production(data_NJ2020_all, write_csv = F)

.GlobalEnv$plot_list <- list()

# Visuals
plot_overview_counts_over_time(data_NJ2020_all)
plot_collision_rates(data_NJ2020_all, original_abundances_NJ2020)
plot_comparison_methods(vp_results_output_df) 
plot_VIPCAL_vs_VIPCAL_SE(vp_results_output_df)

analyze_viral_production(vp_results_output_BP_df, data_NJ2020_all, 
original_abundances_NJ2020, write_csv = F)

plot_percentage_cells(analyzed_vp_results_df)

analyze_viral_production(vp_results_output_T24_df, data_NJ2020_all, 
original_abundances_NJ2020, write_csv = F)

plot_nutrient_release(analyzed_vp_results_df)
}
}
