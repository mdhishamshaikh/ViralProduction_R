# Test bep with different tags
data2 <- read.csv("NJ2020.csv") # For H_code
names(data2)[names(data2) == 'Expt_No'] <- 'Station_Number'  # For My code


# H_code
bep_df<- data2 %>% unite(c('Location', 'Expt_No', 'Depth'), 
                        col = "tag", remove = F)

endpoint_list<- list()

for (combi_tag in unique(bep_df$tag)){
  
  bep_df1<- bep_df%>% filter(tag == combi_tag)
  endpoint_list[[length(endpoint_list)+1]]<- c(combi_tag , bacterial_endpoint_range(bep_df1))
  
}
print(endpoint_list)

# My code
bep_df2<- data2 %>% unite(c('Location', 'Station_Number', 'Depth'), 
                        col = "tag", remove = F)

endpoint_list2<- list()

for (combi_tag in unique(bep_df2$tag)){
  
  bep_df3<- bep_df2%>% filter(tag == combi_tag)
  endpoint_list2[[length(endpoint_list2)+1]]<- c(combi_tag , bacterial_endpoint(bep_df3))
  
}
print(endpoint_list2)




