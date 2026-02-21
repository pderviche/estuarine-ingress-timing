#################################################################################### 
# Tracking estuarine ingress timing of juvenile dog snappers using otolith chemistry
#################################################################################### 

#Authors: Patrick Derviche, ,*, Jean R. S. Vitule, Michael A. Dance, Mario V. Condini,  
#Eduardo S. Costa, Fabian Sá, Felippe A. Daros, Maurício Hostim–Silva, Marcelo Soeth


###################
# R script made by Patrick Derviche
###################

#################################################
#Summary
#################################################


# Intro.........................................................................line 43
# 1 Fish variables..............................................................line 63
# 1.1 General results - fish....................................................line 129
# 1.2 General results - otolith signatures......................................line 167
# 1.3 Estimating the size of the individual based on the laser transect.........line 178
# 1.4 Add the estimated TL......................................................line 244

# 2 Ba - GAMM otolith transect..................................................line 259
# 2.1 Ba - Plot - GAMM otolith transect.........................................line 308
# 2.2 Predict fish TL based on transect.........................................line 386

# 3 Ba - GAMM estimated TL......................................................line 420

# 4 Sr - GAMM...................................................................line 475
# 4.1 Sr - Plot.................................................................line 552

# 5 Mg - GAMM...................................................................line 561
# 5.1 Mg - Plot.................................................................line 598

# 6 All elements................................................................line 638

# 7 Water chemistry.............................................................line 717


####
# Intro
####

# Set working directory
setwd("C:/Users/patri/OneDrive/Documentos/UFES/Tese/estuarine-ingress-timing/Estuarine, Coastal and Shelf Science/github")

# Load libraries
if (!require(pacman)) install.packages("pacman")

pacman::p_load(
  devtools, ggplot2, ggpubr, vegan, dplyr, corrplot, ggbiplot,
  scatterplot3d, ggpmisc, FactoMineR, factoextra, GGally,
  MVN, cluster, ggdendro, gplots, NbClust,
  ggeffects, mgcv, gratia, DHARMa)

# Clean R environment 
rm(list = ls())


##############################################################
# 1 Fish variables
##############################################################

####
# Read datasets
####


otolith <- read.csv2("otolith_chemistry.csv") 

str(otolith)
length(table(otolith$ID))

otolith <- otolith[, c("ID","Month","Year","Season","Ontogeny","Site","TL","Weight","Elapsed_Time","Mg24","Sr87","Ba138")]

otolith <- otolith %>%
  mutate(
    TL = as.numeric(TL),
    Weight = as.numeric(Weight),
    Elapsed_Time = as.numeric(Elapsed_Time),
    Mg24 = as.numeric(Mg24),
    Sr87 = as.numeric(Sr87),
    Ba138 = as.numeric(Ba138),
    Year = as.numeric(Year),
    Ontogeny = as.factor(Ontogeny),
    Month = as.factor(Month),
    Site = as.factor(Site),
    Season = as.factor(Season),
    Month = factor(Month, levels = c('JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN')),
  )

otolith <- otolith %>% filter(Ontogeny == "Juvenile")
otolith <- otolith %>% filter(Site == "Sao Mateus")

str(otolith)
length(table(otolith$ID))
unique(otolith$ID)

otolith$Elapsed_Time <- round(otolith$Elapsed_Time, 1)
otolith$Ba138 <- round(otolith$Ba138, 1)
otolith$Sr87 <- round(otolith$Sr87, 1)
otolith$Mg24 <- round(otolith$Mg24, 1)

otolith %>%
  group_by(Month) %>%
  summarise(n_ID = n_distinct(ID))

otolith %>%
  group_by(Month) %>%
  summarise(mean = mean(TL),
            sd = sd(TL))

otolith$transect <- otolith$Elapsed_Time*10

otolith %>%
  group_by(Month) %>%
  summarise(n_ID = n_distinct(ID))

otolith %>%
  count(ID)

otolith$Ontogeny <- NULL
otolith$Site <- NULL


####
# 1.1 General results - fish
####

#Boxplot TL
mean_TL <- otolith %>%
  group_by(ID, Month) %>%
  summarise(mean_TL = mean(TL, na.rm = TRUE), .groups = "drop")

ggplot(mean_TL, aes(x = Month, y = mean_TL, fill = Month)) +
  geom_boxplot(alpha = 0.2, color = "black") +
  labs(
    x = "Sampled month",
    y = "Mean TL (mm)",
    title = " "
  ) +
  theme_minimal()+
  theme(legend.position = "none")

mean(otolith$TL)
sd (otolith$TL)
summary (otolith)
summary (otolith$TL)

#The TL of the juvenile dog snappers was 177.2 ± 30.9 mm (mean ± standard deviation), ranging from 72.0 to 260.0 mm

# ANOVA
anova <- aov(TL ~ Month, data = otolith)
summary(anova) #Significant differencea, p value = <2e-16 ***

shapiro.test(anova$residuals) #p-value < 2.2e-16, we cannot assume the normality
leveneTest(anova) #p-value < 2.2e-16 ***, we cannot assume the normality

kruskal.test(TL ~ Month, data = otolith)
#Kruskal-Wallis chi-squared = 1031.5, df = 11, p-value < 2.2e-16



####
# 1.2 General results - otolith signatures
####

summary (otolith)
summary (otolith)
summary (otolith)
summary (otolith)



####
# 1.3 Estimating the size of the individual based on the laser transect
####


str(otolith)
unique(otolith$ID)

otolith_tail <- otolith %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%
  ungroup()

# Fit GAM
m_gam <- gam(
  TL ~ s(transect, k = 10),
  family = Gamma(link = "log"),
  method = "REML",
  data = otolith_tail)

# Model summary / diagnostics
summary(m_gam)
res_m_gam <- simulateResiduals(m_gam, plot = TRUE)

# Create prediction grid
newdat <- data.frame(
  transect = seq(min(otolith_tail$transect),
                 max(otolith_tail$transect),
                 length.out = 200))

# Predict on response scale (TL scale) with 95% CI
pred <- predict(m_gam, newdata = newdat, type = "link", se.fit = TRUE)
newdat$fit <- exp(pred$fit)
newdat$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
newdat$upr <- exp(pred$fit + 1.96 * pred$se.fit)

# Plot: points + fitted GAM curve + 95% CI ribbon
fig_TL_transect <- ggplot(otolith_tail, aes(transect, TL)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_ribbon(
    data = newdat,
    aes(x = transect, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,     fill = "#00BFFF",
    alpha = 0.2
  ) +
  geom_line(
    data = newdat,
    aes(x = transect, y = fit),
    inherit.aes = FALSE,color = "#00BFFF",
    linewidth = 1
  ) +
  labs(
    y = expression(Total~length~(TL)),
    x = expression(Distance~from~core~(mu*m))) +
  theme_bw()

fig_TL_transect

ggsave(
  filename = "fig_TL_transect.png",
  plot = fig_TL_transect,
  width = 9,        
  height = 7,    
  units = "cm",
  dpi = 900)

######
# 1.4 Add the estimated TL
######

otolith$estimated_TL <- predict(
  m_gam,
  newdata = otolith,
  type = "response")

otolith$estimated_TL <- round(otolith$estimated_TL, 1)





##############################################################
# 2 Ba - GAMM otolith transect
##############################################################

library(qgam)


summary(otolith$Ba138)
ggplot(otolith, aes(x = transect, y = Ba138, group = ID)) +
  geom_line(alpha = 0.2) +
  theme_bw() 

# Ensure factors
otolith <- otolith %>%
  mutate(Month = factor(Month),
         ID    = factor(ID))

# ---- Candidate GAMMs (Gaussian; random intercept for Month)
# Full model: smooths for all continuous predictors + random effect for Month
gamm_Ba <- gam(Ba138 ~ s(transect) +
                 s(Month, bs = "re") + s(ID, bs = "re"),
               data = otolith, method = "REML", family = Gamma(link = "log"))


# ---- Summary / ANOVA-like tests
summary(gamm_Ba)          # EDFs, F-tests for smooths, R^2, deviance explained
anova(gamm_Ba)            # compares parametric terms (none here), OK to view

# ---- Validation
gam.check(gamm_Ba)        # qq, residuals vs fits, k-index, etc.

# DHARMa residual simulation (works with mgcv::gam)
res_Ba138 <- simulateResiduals(gamm_Ba, plot = TRUE)
testDispersion(res_Ba138)         # dispersion
testUniformity(res_Ba138)         # KS test
testOutliers(res_Ba138)


# ---- R2-like metrics
# From mgcv summary:
#   best_gamm$gcv.ubre ; summary(best_gamm)$r.sq ; summary(best_gamm)$dev.expl
gam_sum_Ba138 <- summary(gamm_Ba)
cat("R^2 (mgcv):", round(gam_sum_Ba138$r.sq, 3), 
    "  Deviance explained:", round(100*gam_sum_Ba138$dev.expl, 1), "%\n")





##############################################################
# 2.1 Ba - Plot
##############################################################

plot(gamm_Ba, pages = 1, shade = TRUE)

draw(gamm_Ba, select = "s(transect)")+
  theme_bw()

###

newdat_Ba <- data.frame(
  transect = seq(min(otolith$transect, na.rm = TRUE),
                 max(otolith$transect, na.rm = TRUE),
                 length.out = 4088),
  Month = otolith$Month[1],  # valid levels needed, but excluded
  ID    = otolith$ID[1])

pred_Ba <- predict(
  gamm_Ba,
  newdata = newdat_Ba,
  se.fit  = TRUE, 
  type    = "response",
  unconditional = TRUE)


newdat_Ba$fit <- pred_Ba$fit
newdat_Ba$se  <- pred_Ba$se.fit


fig3 <- ggplot(newdat_Ba, aes(x = transect, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = fit - 2*se, ymax = fit + 2*se), alpha = 0.10, fill = "blue") +
  labs(y = expression("Predicted Ba:Ca ("*mu*"mol "*mol^{-1}*")"), x = expression(Distance~from~core~(mu*m))) + 
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 1.77, ymax = 2.77,
           fill = "grey70", alpha = 0.3) +
  annotate("text",
           x = Inf, y = (1.77 + 2.77)/2,
           label = "Transition zone",
           color = "black",
           hjust = 1.1,
           size = 3,
           fontface = "bold") +
  geom_vline(xintercept = 337, linetype = "dashed", color = "darkgreen", linewidth = 0.5) +
  annotate("text",
           x = 337,
           y = 9.5,  # ajuste conforme o topo do seu eixo
           label = "Estuarine ingress", color = "darkgreen",
           vjust = 0.7,
           hjust = -0.2,
           size = 3,
           fontface = "bold")+
  theme_bw()

fig3

ggsave(
  filename = "Fig3.png",
  plot = fig3,
  width = 12,        
  height = 9,    
  units = "cm",
  dpi = 900)


first_cross <- newdat_Ba %>%
  dplyr::arrange(transect) %>%
  dplyr::filter(fit >= 2.77) %>%
  dplyr::slice(1)

first_cross
#337 um transect
#TL 90.7



######
# 2.2 Predict fish TL based on transect = 337
######

transect_337 <- data.frame(transect = 337)


# Prediction on the link scale
pred_337 <- predict(
  m_gam,
  newdata = transect_337,
  type = "link",
  se.fit = TRUE)

pred_337

# Estimate on the response scale
fit_resp <- exp(pred_337$fit)
SE_resp <- fit_resp * pred_337$se.fit

# 95% Confidence Interval
lwr_resp <- exp(pred_337$fit - 1.96 * pred_337$se.fit)
upr_resp <- exp(pred_337$fit + 1.96 * pred_337$se.fit)

# SD from Gamma distribution
phi <- summary(m_gam)$scale   # scale parameter
SD_resp <- sqrt(phi * fit_resp^2)

fit_resp #90.7
SD_resp #12.9
SE_resp #7.8
lwr_resp #76.6
upr_resp #107.4

##############################################################
# 3 Ba - GAMM estimated TL
##############################################################

library(qgam)


# ---- Candidate GAMMs (Gaussian; random intercept for Month)
# Full model: smooths for all continuous predictors + random effect for Month
gamm_Ba_estimated_TL <- gam(Ba138 ~ s(estimated_TL) +
                              s(Month, bs = "re") + s(ID, bs = "re"),
                            data = otolith, method = "REML", family = Gamma(link = "log"))

# ---- Summary / ANOVA-like tests
summary(gamm_Ba_estimated_TL)          # EDFs, F-tests for smooths, R^2, deviance explained
anova(gamm_Ba_estimated_TL)            # compares parametric terms (none here), OK to view

# ---- Validation
gam.check(gamm_Ba_estimated_TL)        # qq, residuals vs fits, k-index, etc.

# DHARMa residual simulation (works with mgcv::gam)
res_Ba138_estimated_TL <- simulateResiduals(gamm_Ba_estimated_TL, plot = TRUE)
testDispersion(res_Ba138_estimated_TL)         # dispersion
testUniformity(res_Ba138_estimated_TL)         # KS test
testOutliers(res_Ba138_estimated_TL)


# ---- R2-like metrics
# From mgcv summary:
#   best_gamm$gcv.ubre ; summary(best_gamm)$r.sq ; summary(best_gamm)$dev.expl
gam_sum_Ba138_estimated_TL <- summary(gamm_Ba_estimated_TL)
cat("R^2 (mgcv):", round(gam_sum_Ba138_estimated_TL$r.sq, 3), 
    "  Deviance explained:", round(100*gam_sum_Ba138_estimated_TL$dev.expl, 1), "%\n")


###

newdat_Ba_gamm_Ba_estimated_TL <- data.frame(
  estimated_TL = seq(min(otolith$estimated_TL, na.rm = TRUE),
                     max(otolith$estimated_TL, na.rm = TRUE),
                     length.out = 4088),
  Month = otolith$Month[1],  # valid levels needed, but excluded
  ID    = otolith$ID[1])

pred_Ba <- predict(
  gamm_Ba_estimated_TL,
  newdata = newdat_Ba_gamm_Ba_estimated_TL,
  se.fit  = TRUE, 
  type    = "response",
  unconditional = TRUE)


newdat_Ba$fit <- pred_Ba$fit
newdat_Ba$se  <- pred_Ba$se.fit

##############################################################
# 4 Sr - GAMM
##############################################################


summary(otolith$Sr87)
str(otolith)


ggplot(otolith, aes(x = estimated_TL, y = Sr87, group = ID)) +
  geom_line(alpha = 0.2) +
  theme_bw() 


# ---- Candidate GAMMs (Gaussian; random intercept for Month)
# Full model: smooths for all continuous predictors + random effect for Month
gamm_Sr <- gam(Sr87 ~ s(estimated_TL) +
                 s(Month, bs = "re") + s(ID, bs = "re"),
               data = otolith, method = "REML", family = Gamma(link = "log"))


# ---- Summary / ANOVA-like tests
summary(gamm_Sr)          # EDFs, F-tests for smooths, R^2, deviance explained
anova(gamm_Sr)            # compares parametric terms (none here), OK to view

# ---- Validation
gam.check(gamm_Sr)        # qq, residuals vs fits, k-index, etc.

# DHARMa residual simulation (works with mgcv::gam)
res_Sr87 <- simulateResiduals(gamm_Sr, plot = TRUE)
testDispersion(res_Sr87)         # dispersion
testUniformity(res_Sr87)         # KS test
testOutliers(res_Sr87)


# ---- R2-like metrics
# From mgcv summary:
#   best_gamm$gcv.ubre ; summary(best_gamm)$r.sq ; summary(best_gamm)$dev.expl
gam_sum_Sr87 <- summary(gamm_Sr)
cat("R^2 (mgcv):", round(gam_sum_Sr87$r.sq, 3), 
    "  Deviance explained:", round(100*gam_sum_Sr87$dev.expl, 1), "%\n")

library(ggeffects)
library(mgcv)
library(gratia)


##############################################################
# 4.1 Sr - Plot
##############################################################

plot(gamm_Sr, pages = 1, shade = TRUE)

draw(gamm_Sr, select = "s(estimated_TL)")+
  theme_bw()

###

newdat_Sr <- data.frame(
  estimated_TL = seq(min(otolith$estimated_TL, na.rm = TRUE),
                     max(otolith$estimated_TL, na.rm = TRUE),
                     length.out = 4088),
  Month = otolith$Month[1],  # valid levels needed, but excluded
  ID    = otolith$ID[1]
)

pred_Sr <- predict(
  gamm_Sr,
  newdata = newdat_Sr,
  se.fit  = TRUE, 
  type    = "response",
  unconditional = TRUE
)


newdat_Sr$fit <- pred_Sr$fit
newdat_Sr$se  <- pred_Sr$se.fit


ggplot(newdat_Sr, aes(x = estimated_TL, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = fit - 2*se, ymax = fit + 2*se), alpha = 0.15, fill = "blue") +
  labs(y = expression("Predicted Sr:Ca ("*mu*"mol:"*mol*")"), x = "Elapsed Time") +
  theme_bw()


##############################################################
# 5 Mg - GAMM
##############################################################

summary(otolith$Mg24)
ggplot(otolith, aes(x = estimated_TL, y = Mg24, group = ID)) +
  geom_line(alpha = 0.2) +
  theme_bw() 


# ---- Candidate GAMMs (Gaussian; random intercept for Month)
# Full model: smooths for all continuous predictors + random effect for Month
gamm_Mg <- gam(Mg24 ~ s(estimated_TL) +
                 s(Month, bs = "re") + s(ID, bs = "re"),
               data = otolith, method = "REML", family = Gamma(link = "log"))


# ---- Summary / ANOVA-like tests
summary(gamm_Mg)          # EDFs, F-tests for smooths, R^2, deviance explained
anova(gamm_Mg)            # compares parametric terms (none here), OK to view

# ---- Validation
gam.check(gamm_Mg)        # qq, residuals vs fits, k-index, etc.

# DHARMa residual simulation (works with mgcv::gam)
res_Mg24 <- simulateResiduals(gamm_Mg, plot = TRUE)
testDispersion(res_Mg24)         # dispersion
testUniformity(res_Mg24)         # KS test
testOutliers(res_Mg24)


# ---- R2-like metrics
# From mgcv summary:
#   best_gamm$gcv.ubre ; summary(best_gamm)$r.sq ; summary(best_gamm)$dev.expl
gam_sum_Mg24 <- summary(gamm_Mg)
cat("R^2 (mgcv):", round(gam_sum_Mg24$r.sq, 3), 
    "  Deviance explained:", round(100*gam_sum_Mg24$dev.expl, 1), "%\n")
##############################################################
# 5.1 Mg - Plot
##############################################################

plot(gamm_Mg, pages = 1, shade = TRUE)

draw(gamm_Mg, select = "s(estimated_TL)")+
  theme_bw()

###

newdat_Mg <- data.frame(
  estimated_TL = seq(min(otolith$estimated_TL, na.rm = TRUE),
                     max(otolith$estimated_TL, na.rm = TRUE),
                     length.out = 4088),
  Month = otolith$Month[1],  # valid levels needed, but excluded
  ID    = otolith$ID[1]
)

pred_Mg <- predict(
  gamm_Mg,
  newdata = newdat_Mg,
  se.fit  = TRUE, 
  type    = "response",
  unconditional = TRUE
)


newdat_Mg$fit <- pred_Mg$fit
newdat_Mg$se  <- pred_Mg$se.fit


ggplot(newdat_Mg, aes(x = estimated_TL, y = fit)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = fit - 2*se, ymax = fit + 2*se), alpha = 0.15, fill = "blue") +
  labs(y = expression("Predicted Mg:Ca ("*mu*"mol:"*mol*")"), x = "Elapsed Time") +
  theme_bw()



##############################################################
# 6 All elements
##############################################################

# ---- Transformações ----

# Mg: 150–500  →  0–10
newdat_Mg$fit_scaled <- newdat_Mg$fit / 50
newdat_Mg$se_scaled  <- newdat_Mg$se  / 50

# Sr: 1800–3000  →  0–10
newdat_Sr$fit_scaled <- (newdat_Sr$fit - 1800) / (3000 - 1800) * 10
newdat_Sr$se_scaled  <- newdat_Sr$se / (3000 - 1800) * 10


fig4 <- ggplot() +
  geom_vline(xintercept = 90.7, linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Ba
  geom_ribbon(data = newdat_Ba,
              aes(x = estimated_TL, ymin = fit - 2*se, ymax = fit + 2*se, fill = "Ba:Ca"),
              alpha = 0.15) +
  geom_line(data = newdat_Ba,
            aes(x = estimated_TL, y = fit, color = "Ba:Ca"),
            linewidth = 0.8) +
  
  # Mg
  geom_ribbon(data = newdat_Mg,
              aes(x = estimated_TL, ymin = fit_scaled - 2*se_scaled, ymax = fit_scaled + 2*se_scaled,
                  fill = "Mg:Ca"),
              alpha = 0.15) +
  geom_line(data = newdat_Mg,
            aes(x = estimated_TL, y = fit_scaled, color = "Mg:Ca"),
            linewidth = 0.8) +
  
  # Sr
  geom_ribbon(data = newdat_Sr,
              aes(x = estimated_TL, ymin = fit_scaled - 2*se_scaled, ymax = fit_scaled + 2*se_scaled,
                  fill = "Sr:Ca"),
              alpha = 0.15) +
  geom_line(data = newdat_Sr,
            aes(x = estimated_TL, y = fit_scaled, color = "Sr:Ca"),
            linewidth = 0.8) +
  
  scale_color_manual(
    name = " ",
    breaks = c("Ba:Ca", "Mg:Ca", "Sr:Ca"),
    values = c("Ba:Ca" = "#E41A1C",
               "Mg:Ca" = "#4DAF4A",
               "Sr:Ca" = "#40BBD2")
  ) +
  scale_fill_manual(
    name = " ",
    breaks = c("Ba:Ca", "Mg:Ca", "Sr:Ca"),
    values = c("Ba:Ca" = "#E41A1C",
               "Mg:Ca" = "#4DAF4A",
               "Sr:Ca" = "#40BBD2")
  ) +
  
  scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, 1),
    name = expression(Mg:Ca[otolith]/50 ~~ "and" ~~ Ba:Ca[otolith] ~~ "(" * mu * mol~mol^{-1} * ")"),
    sec.axis = sec_axis(~ . * (3000 - 1800)/10 + 1800,
                        name = expression(Sr:Ca[otolith]~(mu*mol~mol^{-1})))
  ) +
  labs(x = expression(Estimated~fish~TL~(mm))) +
  theme_bw()

fig4

ggsave(
  filename = "Fig4.png",
  plot = fig4,
  width = 17,        
  height = 12,    
  units = "cm",
  dpi = 900)

##############################################################
# 7 Water chemistry
##############################################################

#Set working directory
water <- read.csv2("water_chemistry.csv")
water <- as.data.frame(water)
str(water)

gradient <- water %>%
  filter(site == "Sao Mateus") %>%
  filter(tide == "gradient")  

gradient_Ba <- gradient[,c('Sample', 'Month', 'Year', 'Season', 'BaCa', 'Temp', 'Sal')]

#Set as numeric
gradient_Ba <- mutate_at(gradient_Ba, vars(c('BaCa','Temp', 'Sal')), as.numeric)

gradient_Ba$BaCam <- gradient_Ba$BaCa/1000 #transform umol mol-1 to mmol mol-1

Ba <- ggplot(data = gradient_Ba, aes(y = BaCam, x = Sal)) +
  geom_point(size=2,color='#00BFFF', alpha = 0.3)  +
  geom_smooth(method = "gam",method.args = list(family = Gamma(link = "log")),
              se = TRUE,color='#00BFFF',
              fill = "#00BFFF",
              alpha = 0.2)+
  labs(title = " ",x = "Salinity", y = expression(Ba:Ca[water]~ (mmol~mol^-1)))  +
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 35, by = 5))
Ba

ggsave(
  filename = "Fig2.png",
  plot = Ba,
  width = 9,        
  height = 7,    
  units = "cm",
  dpi = 900)

summary(gradient_Ba)

##############################################################
#END
##############################################################