
{
   library(tidyverse)
  library(ggsci)
  library(scales)
}


NJ1 <- as.data.frame(read_csv("NJ1.csv"))


#For VP, plot total bacteria, viruses, HNA, LNA, V1, V2, V3

VP<- NJ1[NJ1$Sample_Type=='VP',] %>%
  gather('c_Bacteria', 'c_HNA', 'c_LNA', 'c_Viruses', 'c_V1', 'c_V2', 'c_V3', key="count", value="value")%>%
  mutate(Microbe = if_else(VP$count == 'c_Bacteria' | VP$count == 'c_HNA' | VP$count == 'c_LNA', "Bacteria", "Viruses"))



count_factor<- factor(VP$count,
                      levels = c("c_Bacteria", "c_HNA", "c_LNA",
                                 "c_Viruses", "c_V1", "c_V2", "c_V3"))
labels<- c("Total Bacteria", "HNA Bacteria", "LNA Bacteria",
           "Total Viruses", "V1 Viruses", "V2 Viruses", "V3 Viruses")
ggplot(VP, aes(x=Timepoint, y=value, color=count_factor, shape= count_factor)) +
  theme_minimal()+
  theme(plot.title = element_text(face = 'bold'),
        legend.title = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold'))+
  geom_point()+
  scale_color_lancet(name= 'Populations',
                     labels=labels) +
  scale_shape_manual(name = 'Populations', 
                     values = c(16,16,16,17,17,17,17),
                       labels=labels)+
  labs(title = 'NJ1 - Viral Production', subtitle = 'Overview - Bacterial and Viral counts',
       x= 'Sampling Timepoints\n (in hours)', y='FCM Counts \n (in millions)')+
  scale_y_continuous(labels = unit_format(unit = NULL, scale = 1e-6))






