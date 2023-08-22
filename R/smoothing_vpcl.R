library(tidyverse)

length(unique(simu_df$Expt_No))

df1<- simu_df %>% filter(Expt_No == 1,
                         Sample_Type == 'VP')
plot(df1$Count ~ df1$Timepoint)


lowess(df1$Timepoint, df1$Count)

plot(df1$Timepoint, df1$Count)

lines(lowess(df1$Timepoint, df1$Count), col = 'red')
lines(lowess(df1$Timepoint, df1$Count, f = 0.3), col = 'purple')
lines(lowess(df1$Timepoint, df1$Count, f = 3), col = 'orange')
lines(lowess(df1$Timepoint, df1$Count, f = 5), col = 'green')
lines(lowess(df1$Timepoint, df1$Count, f = 1), col = 'pink')


vp.lo<- loess(Count ~ Timepoint, df1)
vp.lo.pred<- predict(vp.lo, data.frame(Timepoint = 1:6), se = T)
class(vp.lo.pred)
vp.lo.pred[[1]]
df2<- data.frame(
  vp.lo.mean = vp.lo.pred[[1]],
  vp.lo.se = vp.lo.pred[[2]],
  Sample_Type = 'VP',
  Timepoint = 1:6
)

df1<- df1 %>% group_by(Sample_Type, Timepoint)%>%
  summarise(vpcl_mean = mean(Count), vpcl_se = plotrix::std.error(Count))


df<- full_join(df1, df2)

df<- df %>% pivot_longer(cols = contains('mean'), names_to = 'mean_method',
                    values_to = 'mean')  %>% 
  pivot_longer(cols = contains('se'), names_to = 'se_method',
                                                          values_to = 'se')
df


ggplot(data = df,
       aes(x = Timepoint,
           y = mean,
           color = mean_method))+
  geom_point(alpha = 0.1)
