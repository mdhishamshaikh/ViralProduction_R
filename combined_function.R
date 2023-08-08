source("sourcesourcebaby.R")
source("vp_calc_functions.R")

data<- read.csv("NJ1.csv")

viral_production<- function(data, method = c(1:12), SR_calc = T, method_sr = c(1:2), 
                            bp_endpoint = T, output.dir = "viral_production_calculations",
                            write_csv = T){
  
  if(file.exists(output.dir)){
    print(paste0("The ", output.dir, " folder already exists"))
    
    stop("Please define another output directory before proceeding")
  } else{
  dir.create(output.dir)
  }
  
  results_path<- paste0("./", output.dir, "/")
  
  #create an empty dataframe
 {
  output_df<- data.frame(matrix(ncol = 10, nrow =0))
  colnames(output_df)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                          'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
  
  output_df <- output_df %>%
    mutate_at(c('Location', 'Expt_No',
                'Depth', 'Time_Range',
                'Population','Sample_Type', 
                'VP_Type'), as.character) %>%
    mutate_at(c('VP', 'VP_SE', 
                'VP_R_Squared'), as.numeric)
  
  }
  
  .GlobalEnv$vp_error_list<- list()
  .GlobalEnv$vp_warn_list<- list()
  
  for(mtd in method){
    
    tryCatch(
      expr = {
        print(paste0("Processing using method: ", names(vp_calc_funct_list)[mtd]))
        
        try(output_df<- output_df %>% 
          full_join(vp_calc_funct_list[[mtd]](data)))
        
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
  
  if(write_csv == T){
  write.csv(output_df, file.path(results_path, 'vp_calc_all.csv'), row.names = F)
  }
  
  
  
  try(if(SR_calc == T){
    
    {
      output_df_sr<- data.frame(matrix(ncol = 11, nrow =0))
      colnames(output_df_sr)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                              'Sample_Type', 'Replicate', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
      
      output_df_sr <- output_df_sr %>%
        mutate_at(c('Location', 'Expt_No', 'Depth', 
                    'Time_Range', 'Population',
                    'Sample_Type', 'Replicate', 
                    'VP_Type'), as.character) %>%
        mutate_at(c('VP', 'VP_SE', 
                    'VP_R_Squared'), as.numeric)
      
      
    }
    
    for(mtd_sr in method_sr){
      
      tryCatch(
        expr = {
          print(paste("Processing using SR method", mtd_sr))
          
          output_df_sr<- output_df_sr %>% 
            full_join(vp_calc_funct_list_sr[[mtd_sr]](data))
          
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
   
     #  print("Started analyses using method 2.1")
     # try( m <- vp_lm_sr(data))
     #    output_df_sr<- output_df_sr%>%
     #    full_join(m)
     # 
     # 
     # 
     #  print("Started analyses using method 6.1")
     # try( n<- vp_vpcl_sr(data))
     #  output_df_sr<- output_df_sr%>%
     #    full_join(n)
    if(write_csv == T){
    write.csv(output_df_sr, file.path(results_path, 'vp_calc_sr.csv'), row.names = F)
    }
    
  })
  
  
  try(if(bp_endpoint == T){
    
    {
      vp_calc_bp<- data.frame(matrix(ncol = 10, nrow =0))
      colnames(vp_calc_bp)<- c('Location', 'Expt_No', 'Depth', 'Time_Range', 'Population',
                               'Sample_Type', 'VP', 'VP_SE', 'VP_R_Squared', 'VP_Type')
      
      vp_calc_bp <- vp_calc_bp %>%
        mutate_at(c('Location', 'Expt_No',
                    'Depth', 'Time_Range',
                    'Population','Sample_Type', 
                    'VP_Type'), as.character) %>%
        mutate_at(c('VP', 'VP_SE', 
                    'VP_R_Squared'), as.numeric)
      
      bep_df<- data %>% unite(c('Location', 'Expt_No', 'Depth'), 
                              col = "tag", remove = F)
      
      endpoint_list<- list()
      
    }
    
    for (combi_tag in unique(bep_df$tag)){
      
      bep_df1<- bep_df%>% filter(tag == combi_tag)
      endpoint_list[[length(endpoint_list)+1]]<- c(combi_tag , bacterial_endpoint_range(bep_df1))
      
    }
    
    for (i in 1: length(unique(bep_df$tag))){
      
      
      
      vp_calc_bp<- vp_calc_bp %>% 
        full_join(
          output_df%>%
            unite("tag", c('Location', 'Expt_No', 'Depth'), remove = F) %>%
            filter(tag == endpoint_list[[i]][1] & Time_Range == endpoint_list[[i]][2])%>%
            select(-tag)%>%
            mutate_at(c('Location', 'Expt_No', 
                        'Depth'), as.character)
        )
      
      
    }
    
    if(write_csv == T){
    write.csv(vp_calc_bp, file.path(results_path, 'vp_calc_bp.csv'), row.names = F)
    }
  }) #This is different for different dataframes. Therefore, I need to add a for loop here
  
  
 
  
  
  return(output_df)
  
  
  
}




vp_calc_H <- viral_production(data, output.dir ="vp_calc_H", write_csv = T)


time_range_factor<- as.vector(c())

TP<- unique(df2$Timepoint)

for (tpr in 2:length(TP)){
  
  a<- paste0("T",TP[1],"_T",TP[tpr])
 
  time_range_factor<- append(time_range_factor, a)
}


time_range_factor

levels<- c(T1 = time_range_factor[1],
           T2 = time_range_factor[2],
           T3 = time_range_factor[3],
           T4 = time_range_factor[4],
           T5 = time_range_factor[5])
levels(time_range_factor)

fct_recode(time_range_factor, !!!levels)
 

time_range_factor<- 
  
  fct_recode(time_range_factor, !!!levels)  %>% fct_relevel(sort)

levels(time_range_factor)


df2<- df2%>%arrange(factor(Time_Range, levels = levels(time_range_factor) ))


tprange<- paste0("'", "ll")


NJ1_vp<- viral_production(data)

NJ2020_vp<- viral_production(NJ2020)
nj_vp<- read.csv("./V5000/vp_calc_all.csv")

ggplot(data = nj_vp, aes(x = VP_Type, y = VP,fill = VP_Type ))+
  geom_violin(color= NA)+
geom_jitter(size = 0.5, width = 0.1, aes(col = Expt_No))+
  #geom_boxplot(aes(fill = Sample_Type))+
  #geom_point( size = 0.5)+
  theme_bw()+
  scale_fill_viridis_d()+
  geom_hline(yintercept = 0)+
  labs(title = 'Comparison of viral production calculation methods')+
  xlab('VP Methods')+
  ylab('Viral Production (per hour)')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(face = 'bold' ),
        axis.title = element_text(face = 'bold'),
        title = element_text(face = 'bold'))

  ylim(-5e+04/2, +5e+04/2)



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
