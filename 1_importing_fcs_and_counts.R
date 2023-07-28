library(tidyverse)
library(readxl)



####SetUp work directory####
work_dir<- "C:/Users/hisham.shaikh/OneDrive - UGent/Projects/FCM_R/ViralProduction_R/"
project_title<- 'NJ2020'



#### Import Metadata####
metadata <- read_excel("Metadata/Metadata_NJ2020.xlsx")
files_to_transfer<- metadata$Sample_Name
origin<- "C:/Users/hisham.shaikh/OneDrive - UGent/Fieldwork/202008_Biweekly_Sampling_NIOZ/NJ2020_RawData/NJ2020_AllData/"
destination<- "Data/"

metadata_processing("Metadata/Metadata_NJ2020.xlsx") #provide the location of the files

import_fcs(origin) #provide path where all the fcs files are stored

ref_fcs_create("vi201104.052")

bacterial_gate(bacterial_gate = 'same', fcs_file = metadata$Sample_Name)

get_bv_stats()

counts <- as.data.frame(read_csv("./results/NJ2020_counts.csv"))



{
  NJ2020<- counts_metadata[counts_metadata$Sample_Type != 'TE',]
  NJ2020<- NJ2020[, c('Sample_Name', 'Staining_Protocol', 'Expt_Date', 'Date_Measurement',
                'Location', 'Expt_No', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', 'c_Bacteria', 'c_HNA', 'c_LNA', 
                'c_Viruses', 'c_V1', 'c_V2', 'c_V3', 'VBR', 'HNAperLNA')]
  NJ2020<- NJ2020[
    with(NJ2020,
         order(
           NJ2020[, 'Staining_Protocol'],
           NJ2020[, 'Location'],
           NJ2020[, 'Expt_No'],
           NJ2020[, 'Depth'],
           NJ2020[, 'Timepoint'],
           NJ2020[, 'Replicate'],
           NJ2020[, 'Sample_Type']
           
         )),
  ]
}
write.csv(NJ2020, "./results/NJ2020.csv", row.names=F)


viral_count_overview_plots()
bacterial_count_overview_plots()

ggplot(data = NJ2020, aes(x = Sample_Type, y = c_Viruses, col = as.factor(Expt_No), group = interaction(Sample_Type,Expt_No)))+
  geom_violin(aes(group = Sample_Type, fill = Sample_Type), alpha = 0.15, color = 'NA')+
  labs(fill = "Sample Type")+
  scale_fill_lancet()+
  geom_boxplot( fill = 'white', outlier.alpha =0)+
  guides(color = guide_legend(title = "Expt No"))+
  geom_jitter(alpha = 1.0, aes(shape = as.factor(Timepoint)), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.9))+
  labs(shape = "Timepoint")+
  scale_color_lancet()+
  theme_classic()+
  labs(title = "Overview of Viral Counts per Sample Type")+
  xlab("Sample Type")+
  ylab("Total Virus Count")+
  scale_shape_manual(values = c(19, 0, 1, 2, 3, 4, 7, 8, 10, 15))

