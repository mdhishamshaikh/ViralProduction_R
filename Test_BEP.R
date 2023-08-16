# Test bep with different tags
# H_code
source('sourcesourcebaby.R')
bep_df<- data_all_t %>% unite(c('Location', 'Expt_No', 'Depth'), 
                        col = "tag", remove = F)

endpoint_list<- list()

for (combi_tag in unique(bep_df$tag)){
  
  bep_df1<- bep_df%>% filter(tag == combi_tag)
  endpoint_list[[length(endpoint_list)+1]]<- c(combi_tag , bacterial_endpoint_range(bep_df1))
  
}
print(endpoint_list)

# My code
source('viral_production_step2_source.R')
bep_df2<- data_all %>% unite(c('Location', 'Station_Number', 'Depth'), 
                        col = "tag", remove = F)

endpoint_list2<- list()

for (combi_tag in unique(bep_df2$tag)){
  
  bep_df3<- bep_df2%>% filter(tag == combi_tag)
  endpoint_list2[[length(endpoint_list2)+1]]<- c(combi_tag , bacterial_endpoint(bep_df3, visual = F))
  
}
print(endpoint_list2)




