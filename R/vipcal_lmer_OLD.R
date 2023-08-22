## Calculating the difference curve for lysogeny using linear mixed effects models 
library(lme4)
NJ1<- read.csv("NJ1.csv")
lmer_data<- df_sr(NJ1)

#only calculating for total virus for now.

lmer_data<- lmer_data[lmer_data$count =="c_Viruses",]
lmer_data$Sample_Type <- as.factor(lmer_data$Sample_Type)
lmer_data$Timepoint2<- as.factor(lmer_data$Timepoint)


lmer_data[lmer_data$Sample_Type == "VPC",]$Replicate<- replace(lmer_data[lmer_data$Sample_Type == "VPC",]$Replicate, lmer_data[lmer_data$Sample_Type == "VPC",]$Replicate == c("1","2", "3"), c("4", "5", "6"))
lmer_data[order(lmer_data$Replicate),]



lyso_model<- lmer(value~Sample_Type*Timepoint2 + (1 | Replicate), data = lmer_data)
summary(lyso_model)

warm.emm_lyso<- emmeans::emmeans(lyso_model, ~ Sample_Type|Timepoint2)
warm.emm_lyso

tw.emm_lyso<- pairs(warm.emm_lyso)
tw.emm_lyso



summary(tw.emm_lyso)$SE
as.numeric(unique(lmer_data$Timepoint))


summary(lmer_data)


plot(-(summary(tw.emm_lyso)$estimate))


peaks(-summary(tw.emm_lyso)$estimate, summary(tw.emm_lyso)$SE)
valleys(-summary(tw.emm_lyso)$estimate, summary(tw.emm_lyso)$SE)


#non equal variance model 
library(nlme)

lyso_model2<- nlme(value~Sample_Type*Timepoint2, weights = varIdent(form = ~ 1| lmer_data$Replicate), data = lmer_data)


####Attempting generalized linear mixed models####
library(lme4)

m<- glmer(data = lmer_data, formula = value~Sample_Type*Timepoint2 + (1 | Replicate), family = "poisson")
summary(m)
warm.m<- emmeans(m, ~ Sample_Type|Timepoint2)
warm.m

n<- lmer(data = lmer_data, formula = value~ Sample_Type + Timepoint2+ Sample_Type*Timepoint2 + ( 1 + Sample_Type  | Replicate))
summary(n)

n.res<- rstandard(n)
warm.m<- emmeans(n, ~ Sample_Type|Timepoint2)
warm.m




nlme<- nlme(value ~ Sample_Type + Sample_Type*Timepoint2,
            data = lmer_data,
            fixed = Sample_Type*Timepoint2 ~ 1,
            random = Replicate ~ 1,
            groups = ~Sample_Type/Replicate,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))


m1<- lm(data = lmer_data, formula = value~ Sample_Type + Timepoint2 + Replicate)
summary(m1)
m2<- lm(data = lmer_data, formula = value~ Sample_Type + Timepoint2)
summary(m2)
m3<- lm(data = lmer_data, formula = value~ Sample_Type + Timepoint2 + Replicate + Sample_Type*Timepoint2)
summary(m3)
m4<- lm(data = lmer_data, formula = value~ Sample_Type + Timepoint2 +  Sample_Type*Timepoint2)
summary(m4)
m5<- lmer(data = lmer_data, formula = value~ Sample_Type + Timepoint2 + Sample_Type*Timepoint2 + (1 | Replicate))
summary(m5)
m6<- lmer(data = lmer_data, formula = value~  Sample_Type*Timepoint2 + (1 | Replicate))
summary(m6)
m7<- lmer(data = lmer_data, formula = value~  Sample_Type*Timepoint2 + (1 + value| Replicate))
summary(m7)
m8<- lmer(data = lmer_data, formula = value~  Sample_Type*Timepoint2 + (1 + Sample_Type| Replicate))
summary(m8)

AIC.table <- MuMIn::model.sel(m1, m2, m3, m4, m5, m6, m7,
                              m8)
(AIC.table <- AIC.table[, c("df", "logLik", "AICc", "delta")])

#plotting residuals  for m5, m6, m7, m8
par(mfrow = c(2,2))
for( model in 5:8){
  a<- paste0("m", model)
 
hist(resid(get(a))) 
  
}

#plotting residuals vs predeicted for m5, m6, m7, m8
par(mfrow = c(2,2))
for( model in 5:8){
  a<- paste0("m", model)
   print(a)
   plot(resid(get(a)) ~ fitted(get(a)), 
        xlab = "Predicted values", 
        ylab = "Normalized residuals", 
        main = paste("LMEM", a , sep = "_"))
   abline(h = 0, lty = 2)
  
}

#plotting residuals vs replicates for m5, m6, m7, m8
par(mfrow = c(2,2))
for( model in 5:8){
  a<- paste0("m", model)
  print(a)
  boxplot(resid(get(a)) ~ Replicate, 
          data = lmer_data, xlab = "Replicate",
          ylab = "Normalized residuals",
          main = paste("LMEM", a , sep = "_"))
 
  abline(h = 0, lty = 2)
  
}

#plotting predicted vs observed for m5, m6, m7, m8
par(mfrow = c(2,2))
for( model in 5:8){
  a<- paste0("m", model)
  plot(predict(get(a)) ~ lmer_data$value, 
       xlab = "Predicted values", 
       ylab = "Normalized residuals", 
       main = paste("LMEM", a , sep = "_"))
  abline(a= 0, b = 1)
  
}

#m7
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(m7) ~ fitted(m7), xlab = "Predicted values", ylab = "Normalized residuals", main = "")
abline(h = 0, lty = 2)

par(mfrow = c(1, 3), mar = c(4, 4, 0.5, 0.5))
boxplot(resid(m7) ~ Replicate, data = lmer_data, xlab = "Replicate",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
hist(resid(m7))

#m6
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(m6) ~ fitted(m6), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)


boxplot(resid(m6) ~ Replicate, data = lmer_data, xlab = "Replicate",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
hist(resid(m6))
hist(predict(m6))
hist(lmer_data$value)
plot(predict(m6),
     lmer_data$value)
abline(a= 0, b = 1)
#m8
par(mar = c(4, 4, 0.5, 0.5))
plot(resid(m8) ~ fitted(m8), xlab = "Predicted values", ylab = "Normalized residuals")
abline(h = 0, lty = 2)


boxplot(resid(m8) ~ Replicate, data = lmer_data, xlab = "Replicate",
        ylab = "Normalized residuals")
abline(h = 0, lty = 2)
hist(resid(m8))

#https://r.qcbs.ca/workshop07/book-en/step-3.-model-validation.html

plot(lmer_data$value ~  lmer_data$Timepoint)


m6_means<- emmeans(m6, ~ Sample_Type|Timepoint2)
m6_means
plot(summary(m6_means)$SE ~ summary(m6_means)$Timepoint2)
ggplot2::ggplot(data = m6_means, aes(summary(m6_means)$emmean, summary(m6_means)$Timepoint2)) 




#Trying Pavpop for nlme

devtools::install_github('larslau/pavpop')
library(pavpop)
data<- lmer_data[, c(6,4,5,8)]
tp<- unique(data$Timepoint)
value<- array(dim = c(6,6))
data[data$Replicate == "1",]$value
for (i in unique(data$Replicate)){
  value[,as.numeric(i)]<- data[data$Replicate == i ,]$value
}

identical(data[data$Replicate == "6",]$value, value[,6])
n<- ncol(value)
tp_range<- c(0, 25)
t<-replicate(ncol(value),tp, simplify = F)
y<- lapply(1:n, function (x) value[, x])

kts <- seq(tp_range[1], tp_range[2], length = 15)
basis_fct <- make_basis_fct(kts = kts, type = 'increasing', intercept = TRUE,
                            control = list(boundary = tp_range))
tw <- seq(tp_range[1], tp_range[2], length = 6)
warp_fct <- make_warp_fct('smooth', tw, control = list(wright = 'extrapolate'))
mw <- attr(warp_fct, 'mw')

# Set up covariance functions
warp_cov_par <- c(tau = 10)
warp_cov <- make_cov_fct(Brownian, noise = FALSE, param = warp_cov_par, type = 'motion',
                         range = tp_range)

amp_cov_par <- c(scale = 200, range = 10, smoothness = 2)
amp_cov <- make_cov_fct(Matern, noise = TRUE, param = amp_cov_par)
# Estimate in the model

# Bounds of parameters
# NOTE: Prediction of velocities is only meaningful
#       when the smoothness parameter is > 0.5
lower <- c(1e-2, 1e-2, 0.5001, 1e-2)
upper <- c(1000, Inf, Inf, Inf)

res <- pavpop(y, t, basis_fct, warp_fct, amp_cov, warp_cov, homeomorphisms = 'soft',
              like_optim_control = list(lower = lower, upper = upper))
# Plot results
#

t_p <- seq(range(t)[1], range(t)[2], length = 100)

# Functional fixed effect
theta <- basis_fct(t_p) %*% res$c

# Display data with predictions
plot(t_p, theta, ylim = range(y), type = 'n', main = 'Original counts and predicted',
     xlab = 'Time (hours)', ylab = 'count')
for (i in 1:n) {
  points(t[[i]], y[[i]], pch = 19, cex = 0.3, col = rainbow(n)[i])
  lines(t_p, predict_curve(t_p, t[[i]], y[[i]], basis_fct, res$c, warp_fct, res$w[, i],
                           amp_cov, res$amp_cov_par),
        lwd = 0.5, col = rainbow(n)[i])
}
lines(t_p, theta, ylim = range(y), lwd = 2, lty = 2)
# Display predicted warping functions
plot(t_p, t_p, type = 'l', lwd = 2, lty = 2, main = 'Warping functions',
     xlab = 'Age (years)', ylab = 'Biological age (years)')
for (i in 1:n) lines(t[[i]], warp_fct(res$w[,i], t[[i]]), lwd = 0.4, col = rainbow(n)[i])


