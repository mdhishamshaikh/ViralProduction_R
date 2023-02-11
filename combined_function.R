data<- read.csv("NJ1.csv")

viral_production<- function(data, method = c(1:12), SR_calc = T, bp_endpoint = T){
  
  if(file.exists('viral_production_calculations')){
    print("The viral_production_calculations folder already exists")
    
    stop("Please delete the folder before proceeding")
  } else{
  dir.create('viral_production_calculations')
  }
  
  
  
  #create an empty dataframe
 {
  output_df<- data.frame(matrix(ncol = 10, nrow =0))
  colnames(output_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                          'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
  
  output_df[, c('Location', 'Expt_No', 'Depth', 
                'Time_Range', 'Population',
                'Sample_Type',  'VP_Type') ] <- lapply(c('Location', 'Expt_No', 
                                                         'Depth', 'Time_Range', 
                                                         'Population',
                                                         'Sample_Type', 'VP_Type'),
                                                       as.character)
  
  output_df[, c('VP', 'VP_SE', 'VP_R_Squared') ] <- lapply(c('VP', 'VP_SE', 
                                                             'VP_R_Squared'), 
                                                           as.numeric)
  
  }
  
  for(mtd in method){
    
    tryCatch(
      expr = {
        print(paste("Processing using method", mtd))
        
        output_df<- output_df %>% 
          full_join(vp_calc_funct_list[[mtd]](data))
        
      },
      
      error = function(e){
        err<- paste(Sys.time(),
                    paste("Error in Analysis using Method:", mtd),
                    e)
        print(err)
        vp_error_list[[length(vp_error_list)+1]]<<- err
        
        
      },
      
      warning = function(w){
        warn<- paste(Sys.time(),
                     paste("Warning in Analysis using Method:", mtd),
                     w)
        print(warn)
        vp_warn_list[[length(vp_warn_list)+1]]<<- warn
      },
      
      finally = {
        
        print(paste0("Analyses done for method ", mtd, ". Please check vp_error_list and vp_warn_list for any error and warnings" ))
      }
      
    )
    
  }
  
  results_path<- "./viral_production_calculations"
  write.csv(output_df, file.path(results_path, 'vp_calc_all.csv'), row.names = F)
  
  try(if(SR_calc == T){
    
    {
      output_df_sr<- data.frame(matrix(ncol = 11, nrow =0))
      colnames(output_df_sr)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                              'Sample_Type', 'Replicate', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
      
      output_df_sr[, c('Location', 'Expt_No', 'Depth', 
                    'Time_Range', 'Population',
                    'Sample_Type', 'Replicate', 'VP_Type') ] <- lapply(c('Location', 'Expt_No', 
                                                             'Depth', 'Time_Range', 
                                                             'Population',
                                                             'Sample_Type', 'Replicate',
                                                             'VP_Type'),
                                                           as.character)
      
      output_df_sr[, c('VP', 'VP_SE', 'VP_R_Squared') ] <- lapply(c('VP', 'VP_SE', 
                                                                 'VP_R_Squared'), 
                                                               as.numeric)
      
    }
    
 
      print("Started analyses using method 2.1")
     try( m <- vp_lm_sr(data))
        output_df_sr<- output_df_sr%>%
        full_join(m)
   
    

      print("Started analyses using method 6.1")
     try( n<- vp_vpcl_sr(data))
      output_df_sr<- output_df_sr%>%
        full_join(n)

  
    
  })
  
  write.csv(output_df_sr, file.path(results_path, 'vp_calc_sr.csv'), row.names = F)
  
  try(if(bp_endpoint == T){
    
    vp_calc_bp<- output_df %>%
      filter(Time_Range == bacterial_endpoint_range(data))
  }) #This is different for different dataframes. Therefore, I need to add a for loop here
  
  
 
  
  write.csv(vp_calc_bp, file.path(results_path, 'vp_calc_bp.csv'), row.names = F)
  return(output_df)
  
  
  
}

NJ1_vp<- viral_production(data)

NJ2020_vp<- viral_production(NJ2020)

ggplot(data = NJ2020_vp, aes(x = VP_Type, y = VP,fill = VP_Type ))+
  geom_violin(color= NA)+
#  geom_jitter(size = 0.5, width = 0.1)+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_rickandmorty()+
  geom_hline(yintercept = 0)+
  labs(title = 'Comparison of viral production calculation methods')+
  xlab('VP Methods')+
  ylab('Viral Production (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))+
  ylim(-2.5e+06, +2.55e+06)



write.csv(NJ1_vp, file.path(results_path, 'NJ1_vp.csv'), row.names = F)




NJ1_vp_bp<- NJ1_vp%>%
  filter(Time_Range == bacterial_endpoint_range(data))


NJ1_vp_24<- NJ1_vp%>%
  filter(Time_Range == 'T0_T24')

NJ1_vp_20<- NJ1_vp%>%
  filter(Time_Range == 'T0_T20')

comp<- full_join(NJ1_vp_bp, NJ1_vp_20)%>%
  full_join(NJ1_vp_24)%>%
  filter(Population == 'c_Viruses')%>%
   filter(VP_Type == 'LM_AP' |VP_Type == 'VPCL_AR_Diff_LMER_SE' | VP_Type == 'VPCL_AR_Diff_No_SE')%>%
  select(-c('VP_SE', 'VP_R_Squared'))%>%
  pivot_wider(values_from = VP, names_from = Time_Range)



ggplot(data = comp)+
  geom_jitter(aes(x= T0_T6 , y = T0_T24, 
                 shape = Sample_Type, col = VP_Type))+
  geom_abline(intercept = 0)


abline(0,1)
abline(h = 0, v = 0)

ggplot(data = NJ1_vp, aes(x = VP_Type, y = VP,fill = VP_Type, col = Sample_Type, shape = Population ))+
  geom_point()+
 
  theme_bw()+
  scale_color_lancet(alpha = 0.75)+
  geom_hline(yintercept = 0)

ggplot(data = NJ1_vp, aes(x = VP_Type, y = VP,fill = VP_Type ))+
  geom_violin(color= NA)+
  geom_jitter(size = 0.5, width = 0.1)+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_rickandmorty()+
  geom_hline(yintercept = 0)+
  labs(title = 'Comparison of viral production calculation methods')+
  xlab('VP Methods')+
  ylab('Viral Production (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))

ggplot(data = NJ1_vp[NJ1_vp$Sample_Type == 'Diff',], aes(x = VP_Type, y = VP,fill = VP_Type ))+
  geom_violin()+
  geom_jitter(size = 0.5, width = 0.1)+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_rickandmorty()+
  geom_hline(yintercept = 0)+
  labs(title = 'Comparison of viral production calculation methods')+
  xlab('VP Methods')+
  ylab('Viral Production (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))





ggplot(data = NJ1_vp[NJ1_vp$Time_Range == 'T0_T17' & NJ1_vp$Sample_Type == 'VPC',], aes(x = VP_Type, y = VP, col = Sample_Type, shape = Population ))+
  geom_point()+
  theme_classic()+
  scale_color_lancet(alpha = 0.75)+
  geom_hline(yintercept = 0)


TRUE %in% NJ1_vp[NJ1_vp$Time_Range == 'T0_T17',]$VP == 0
NJ1_vp$VP == 0

NJ1_vp_long<- NJ1_vp[,c(1:7,10)] %>% pivot_wider(names_from = 'VP_Type', values_from = 'VP')

unique(NJ1_vp$VP_Type)

ggplot(NJ1_vp_long[NJ1_vp_long$Time_Range == 'T0_T17',], aes(x = VPCL_AR_Diff_No_SE, y = VPCL_AR_Diff_LMER_SE ))+
  geom_point(aes(col = Sample_Type, shape= Population))+
  geom_abline(intercept = 0)+
  geom_smooth(method = 'lm')

ggplot(NJ1_vp_long, aes(x = VPCL_AR_Diff_No_SE, y = VPCL_AR_Diff_LMER_SE ))+
  geom_point(aes(shape = Sample_Type, col= Time_Range))+
  geom_abline(intercept = 0)+
  geom_smooth(method = 'lm')

ggplot(NJ1_vp_long[NJ1_vp_long$Time_Range == 'T0_T17',], aes(x = VPCL_AR_Diff_No_SE, y = VPCL_AR_Diff_LMER_SE ))+
  geom_point(aes(shape = Sample_Type, col= Population))+
  geom_abline(intercept = 0)+
  geom_smooth(method = 'lm')


ab<- full_join(a,b)
full_join(output_df,a
          ) %>%
  full_join(b)







output_df<- data.frame(matrix(ncol = 10, nrow =0))
colnames(output_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                        'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')

output_df[, c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
              'Sample_Type',  'VP_Type') ] <- lapply(c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                                                                                     'Sample_Type', 'VP_Type'), as.character)
output_df[, c('VP', 'VP_SE', 'VP_R_Squared') ] <- lapply(c('VP', 'VP_SE', 'VP_R_Squared'), as.numeric)
summary(output_df)




empty<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
          'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
