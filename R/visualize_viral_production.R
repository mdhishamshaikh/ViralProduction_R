visualize_viral_production <- function(){
  if(!exists("output_dir", where = globalenv(), inherits = FALSE)){
    stop('There exists no output_dir value in the global environment, please define output folder to save figures!')
    
  } else if(!file.exists(.GlobalEnv$output_dir)){
    warning('The output folder does not exists, figures will be saved in given folder!')
    .GlobalEnv$visuals_path <- file.path(.GlobalEnv$output_dir, "Figures")
    
    if(!dir.exists(.GlobalEnv$visuals_path)){
      dir.create(.GlobalEnv$visuals_path)
    }
  } else {
    .GlobalEnv$visuals_path <- file.path(.GlobalEnv$output_dir, "Figures")
    
    if(!dir.exists(.GlobalEnv$visuals_path)){
      dir.create(.GlobalEnv$visuals_path)
    }
  }
  
  
}
