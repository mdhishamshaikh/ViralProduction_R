#' List of methods used for viral production calculation
#' 
#' @description
#' A list of all the different variants of linear regression and VIPCAL available to calculate the viral
#' production rate.
#' 
#' More details on the different variants:
#' 
#' - Variants of linear regression: [viralprod::vp_methods_linear]
#' - Variants of VIPCAL: [viralprod::vp_methods_VIPCAL]
#' 
#' @return A list with the different variants of linear regression and VIPCAL available to calculate viral production.
#' 
#' @name vp_list_of_methods
#' @rdname vp_list_of_methods
#' 
#' @export
vp_list_of_methods <- function(){
  methods_to_calculate_viral_production <- list(
    "vp_linear_allpoints",
    "vp_linear_separate_replicates",
    "vp_linear_average_replicates",
    "vp_linear_average_replicates_diff",
    "vp_linear_average_replicates_diff_LMER",
    "vp_VIPCAL_separate_replicates",
    "vp_VIPCAL_average_replicates",
    "vp_VIPCAL_average_replicates_SE",
    "vp_VIPCAL_average_replicates_diff",
    "vp_VIPCAL_average_replicates_diff_SE",
    "vp_VIPCAL_average_replicates_diff_LMER",
    "vp_VIPCAL_average_replicates_diff_LMER_SE"
  )
  list_of_methods <- list()
  
  for (method in methods_to_calculate_viral_production){
    method_function <- get(method)
    list_of_methods[[length(list_of_methods) + 1]] <- method_function
  }
  names(list_of_methods) <- methods_to_calculate_viral_production
  
  .GlobalEnv$list_of_methods <- list_of_methods
}
