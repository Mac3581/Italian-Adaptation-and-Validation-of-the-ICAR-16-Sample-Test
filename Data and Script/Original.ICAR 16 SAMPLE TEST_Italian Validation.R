#Packages
library(lattice) 
library(DescTools)
library(RcmdrMisc)
library(sjPlot)
library(pwr)
library(effsize)
library(psych)
library(MVN)
library(car)
library(performance)
library(lmtest)
library(ggplot2)
library(lavaan)
library(semPlot)
library(semTools)
library(polycor)
library(REdaS)
library(psych)
library(readr)

###############################################################################

#load_Data.frame
icar <- read_csv("icar_p.csv")
#View the DF
View(icar)

#############################################################################
#Data Cleaning

#Exclude all the non-Italians Participants from the dataset
icar<-subset(icar, icar$Nationality=="Italiana")
#########################################################################

#Demographic variables


#1. Sex Distribution
icar$sex2<-as.factor(icar$Sex)
str(icar$sex2)#factor, with 3 levels. we consider intersex as NA
icar$sex2<- ifelse(icar$sex2=="Intersex", "NA", icar$sex2)
str(icar$sex2)#chr, with 1=female; 3= male
icar$sex2 <- ifelse(icar$sex2 == "1", "F", 
                    ifelse(icar$sex2 == "3", "M", NA))
icar$sex2<-as.factor(icar$sex2)
#Describing sex distribution
Desc(icar$sex2, plotit = F)



#2. Age Distribution
Desc(icar$Age, main = "Age", plotit = F)
Desc(icar$Age~icar$sex2, plotit = F)



#3. Education Level Distribution
icar$Education_level2<-as.factor(icar$Education_level)
str(icar$Education_level2)
Desc(icar$Education_level2, plotit = F)


#Translate education levels from Italian to English

# "Licenza media" with "Diploma"
icar$Education_level2 <- factor(icar$Education_level2)
levels(icar$Education_level2)[levels(icar$Education_level2) == "Licenza media"] <- "Diploma"

# Rename categories
icar$Education_level2 <- factor(icar$Education_level2,
                                levels = c("Diploma", 
                                           "Laurea di primo livello",
                                           "Laurea specialistica/magistrale",
                                           "Dottorato/Specializzazioni post laurea"),
                                labels = c("Undergraduate",
                                           "Bachelor",
                                           "Master degree",
                                           "Ph.D/Post graduate specialization"))

# Check
Desc(icar$Education_level2, plotit = F)
Desc(icar$Education_level2~icar$sex2, plotit = F)


#4. Geographical Distribution
Desc(icar$Region, plotit = F, maxrows = Inf)

df_Macro_Regions<- data.frame(
  Region = c("Calabria", "Campania", "Toscana", "Emilia-Romagna", "Lazio", 
             "Abruzzo", "Marche", "Veneto", "Piemonte", "Liguria", 
             "Trentino-Alto Adige", "Sicilia", "Lombardia", "Sardegna", "Umbria", "Friuli-Venezia Giulia"),
  Latitude = factor(c("South", "South", "Center", "North", "Center",
                      "South", "Center", "North", "North", "North",
                      "North", "South", "North", "South", "Center", "North"), 
                    levels = c("South", "Center", "North")), 
  stringsAsFactors = FALSE
)

icar$latitude<-df_Macro_Regions$Latitude[match(icar$Region, df_Macro_Regions$Region)]


str(icar$latitude)


Desc(icar$latitude, plotit = F)

Desc(icar$latitude~icar$sex2, plotit = F)

################################################################################
################################################################################

#DATA ANALYSES


#5. Preliminary Analysis: icar 16 total score and items

Desc(icar$ICAR, plotit = F)#Total score distribution
shapiro.test(icar$ICAR)#ICAR-16 distribution normality test
qqPlot(icar$ICAR)

#ICAR 16 Mean and standard deviation: CTT

mean<-sapply(icar[, 13:28], mean)
sd<-sapply(icar[, 13:28], sd)


# table
table <- data.frame(
  Item = names(icar[, 13:28]),
  mean = mean,
  sd= sd,
  row.names = NULL
)

# view the table
print(table)


################################################################################
################################################################################

#6. Condon & Revelle ICAR16_NORM (Condon & Revelle, 2014)

Icar_norm<- sum(.67, .69, .73, .61,
                .62, .59, .62, .42,
                .52, .60, .62, .36, 
                .17, .21, .29, .17)


cor_int <- .81/(16 - .81*(16-1))
sd_items <- c(.49, .49, .48, .49, .5, .49, .48, .48, .37, .41, .46, .37, .47, .46, .44, .49)

sum_var <- sum(sd_items^2)
sum_sd <- sum(sd_items)

var_total <- sum_var + cor_int * (sum_sd^2 - sum_var)
sd_total <- sqrt(var_total)

################################################################################
################################################################################
#7. Descriptive Statistics of CRT-long

Desc(icar$CRT, plotit = F)


meancrt<-sapply(icar[, c(25:30)], mean)
sdcrt<-sapply(icar[, c(25:30)], sd)


# table
table2 <- data.frame(
  Item = names(icar[, c(25:30)]),
  mean = meancrt,
  sd= sdcrt,
  row.names = NULL
)

# view the table
print(table2)

################################################################################
################################################################################

#8. Polychoric Matrix correlation of ICAR 16 items

coricar<-c('VR4' , 'VR17' , 'VR16' , 'VR19' , 'R3D3' , 'R3D4' , 'R3D6' , 'R3D8' , 'LN7' , 'LN33' , 'LN34' , 'LN58' , 'MR45' , 'MR46' , 'MR47' , 'MR55')
teticar<- psych::polychoric(icar[, 8:23][, coricar])

psych::cor.plot(teticar$rho, main = "ICAR 16", upper = F, diag = F)


################################################################################
#9. Structural Validity: 

#CFA 1
library(lavaan)
library(semTools)

dev.off()          #only if necessary
MVN::mvn(icar[, 8:23], mvnTest = "mardia", multivariatePlot = "qq")# Multivariate Normality not respected
#DWLS has to be used as estimator

#Model Estimation


#Unidimensional

Unidimensional<- 'g=~ VR4 + VR17 + VR16 + VR19 + R3D3 + R3D4 + R3D6 + R3D8 + LN7 + LN33 + LN34 + LN58 + MR45 + MR46 + MR47 + MR55'
            
Uni_CFA<-lavaan::cfa(Unidimensional, data = icar[, 8:23], ordered = T)
summary(Uni_CFA, standardize= TRUE, ci= TRUE)
Fit_uni<-fitmeasures(Uni_CFA, fit.measures=c ("df", "chisq.scaled" , "pvalue.scaled", "cfi.scaled", "tli.scaled","rmsea.ci.lower.scaled", "rmsea.scaled","rmsea.ci.upper.scaled" ,"srmr")) 
print(Fit_uni)
REL_uni<-semTools::reliability(Uni_CFA, return.total = T)#Reliability indexes of the uni_factor model 
print(REL_uni)
lavaangui::lavaangui(Uni_CFA)


#Bi_factor model

Bi<- '
          VR=~ VR4 + VR17 + VR16 + VR19 
          R3D=~ R3D3 + R3D4 + R3D6 + R3D8
          LN=~ LN7 + LN33 + LN34 + LN58
          MR=~ MR45 + MR46 + MR47 + MR55 
          g =~ VR4 + VR17 + VR16 + VR19 + R3D3 + R3D4 + R3D6 + R3D8 + LN7 + LN33 + LN34 + LN58 + MR45 + MR46 + MR47 + MR55
          '

Bi_CFA<-lavaan::cfa(Bi, data = icar[, 8:23], ordered = T, orthogonal= T)
summary(Bi_CFA, standardize= TRUE, ci= TRUE)
Fit_bi<-fitmeasures(Bi_CFA, fit.measures=c ("df", "chisq.scaled" , "pvalue.scaled", "cfi.scaled", "tli.scaled","rmsea.ci.lower.scaled", "rmsea.scaled","rmsea.ci.upper.scaled" ,"srmr")) 
print(Fit_bi)

REL_bi<-semTools::reliability(Bi_CFA, return.total = T)
print(REL_bi)
lavaangui::lavaangui(Bi_CFA)


#Second Order Latent Variable Model (CONDON & REVELLE STRUCTURE MODEL)

SOLV<- '
          VR=~ VR4 + VR17 + VR16 + VR19 
          R3D=~ R3D3 + R3D4 + R3D6 + R3D8
          LN=~ LN7 + LN33 + LN34 + LN58
          MR=~ MR45 + MR46 + MR47 + MR55 
          g =~ VR + R3D + LN + MR
          g~~g
          '

SOLV_CFA<-lavaan::cfa(SOLV, data = icar[, 8:23], ordered = T, orthogonal= F)
summary(SOLV_CFA, standardize= TRUE, ci= TRUE)

Fit_solv<-fitmeasures(SOLV_CFA, fit.measures=c ("df", "chisq.scaled" , "pvalue.scaled", "cfi.scaled", "tli.scaled","rmsea.ci.lower.scaled", "rmsea.scaled","rmsea.ci.upper.scaled" ,"srmr")) 
print(Fit_solv)

REL_SOLV<-semTools::reliability(SOLV_CFA, return.total = T) 
print(REL_SOLV)

lavaangui::lavaangui(SOLV_CFA)
#Compare the two CFA models: 1) Model fit indices; 2) Chi-Square Difference Test (for nested models). 
#if models are nested do not use AIC and BIC

lavaan::lavTestLRT(SOLV_CFA, Uni_CFA)



################################################################################
#####################################################################

#10. Reliability: internal consistency

#Hierarchical Omega, Total Omega and Cronbach Alpha

omegaicar<- psych::omega(icar[, 8:23], nfactors = 4, fm = "minres")

print(omegaicar)

#What about CRT-long Reliability?
omegacrt<- psych::omega(icar[, 25:30], nfactors = 1, fm= "minres")
print(omegacrt)

#Condon & Revelle 2014, icar16: hierarchical omega : 0.6; total Omega : 0.83

###############################################################################

#11. CFA2 and ICAR + CRT correlation: convergent validity (ICAR 16 + CRT-Long)

#tetracoric matrix correlation icar + crt

crticar<-c('VR4' , 'VR17' , 'VR16' , 'VR19' , 'R3D3' , 'R3D4' , 'R3D6' , 'R3D8' , 'LN7' , 'LN33' , 'LN34' , 'LN58' , 'MR45' , 'MR46' , 'MR47' , 'MR55', 'CRT1', 'CRT2', 'CRT3', 'CRT4', 'CRT5', 'CRT6')
tet<- psych::polychoric(icar[, c(8:23, 25:30)][, crticar])
psych::cor.plot(tet$rho, main = "icar & crt", upper = F, diag = F)


#Convergent Validity by correlation method
#correlation between icar scores and crt scores

cor.test(icar$ICAR, icar$CRT, method = "spearman") #positive correlation 0.616 (moderate/strong)


#Unidimensional model: g loaded tests


CV_uni<- '  g=~ R3D3 + R3D4 + R3D6 + R3D8 + LN7 + LN33 + LN34 + LN58 + MR45 + MR46 + MR47 + MR55 + VR4 + VR17 + VR16 + VR19 + CRT1 + CRT2 + CRT3 + CRT4 + CRT5+ CRT6
          '

Unidimensionale<-cfa(CV_uni, data = icar[, c(8:23, 25:30)], ordered = T)
summary(Unidimensionale, standardize= TRUE, ci= TRUE)
fitmeasures(Unidimensionale, fit.measures=c ("df.scaled", "chisq.scaled" , "pvalue.scaled", "cfi.scaled", "tli.scaled","rmsea.ci.lower.scaled", "rmsea.scaled","rmsea.ci.upper.scaled" ,"srmr"))
REL_unicrticr<-semTools::reliability(Unidimensionale, return.total = T)
print(REL_unicrticr)

#Orthogonal bidimensional

CV3<- '   CRT=~ CRT1 + CRT2 + CRT3 + CRT4 + CRT5+ CRT6 
          ICAR=~ R3D3 + R3D4 + R3D6 + R3D8 + LN7 + LN33 + LN34 + LN58 + MR45 + MR46 + MR47 + MR55 + VR4 + VR17 + VR16 + VR19
          ICAR~~0*CRT
          
'

Distinti_e_non_correlati<-cfa(CV3, data = icar[, c(8:23, 25:30)], ordered = T)
summary(Distinti_e_non_correlati, standardize= TRUE, ci= TRUE)
fitmeasures(Distinti_e_non_correlati, fit.measures=c ("df.scaled", "chisq.scaled" , "pvalue.scaled", "cfi.scaled", "tli.scaled","rmsea.ci.lower.scaled", "rmsea.scaled","rmsea.ci.upper.scaled" ,"srmr"))

REL_CRTICR2<- semTools::reliability(Distinti_e_non_correlati, return.total = T)
print(REL_CRTICR2)


#CV bidimensional oblique


CV1<- '   CRT=~ CRT1 + CRT2 + CRT3 + CRT4 + CRT5+ CRT6 
          ICAR=~ R3D3 + R3D4 + R3D6 + R3D8 + LN7 + LN33 + LN34 + LN58 + MR45 + MR46 + MR47 + MR55 + VR4 + VR17 + VR16 + VR19
          ICAR~~CRT
'

Distinti_ma_liberi_di_Correlare<-cfa(CV1, data = icar[, c(8:23, 25:30)] , ordered = T)
summary(Distinti_ma_liberi_di_Correlare, standardize= TRUE, ci= TRUE)
fitmeasures(Distinti_ma_liberi_di_Correlare, fit.measures=c ("df.scaled", "chisq.scaled" , "pvalue.scaled", "cfi.scaled", "tli.scaled","rmsea.ci.lower.scaled", "rmsea.scaled","rmsea.ci.upper.scaled" ,"srmr"))
REL_CRTICR<- semTools::reliability(Distinti_ma_liberi_di_Correlare, return.total = T)
print(REL_CRTICR)

# CI 95%
std_sol <- standardizedSolution(Distinti_ma_liberi_di_Correlare)
# CI 95%
std_sol[std_sol$op == "~~" & std_sol$lhs == "CRT" & std_sol$rhs == "ICAR", ]

#upper limits: <0.80= two different scales; <= 0.80 e <= 0.90 =  two different scales, but correlated; 
# >= 0.90 = two highly correlated scales measuring the same psychological construct. 
# >=1 = indistinguishable scales

lavaangui::lavaangui(Distinti_ma_liberi_di_Correlare)


# compare two nested models
anova(Unidimensionale, Distinti_e_non_correlati, Distinti_ma_liberi_di_Correlare)  
################################################################################
################################################################################

#12. Unidimensional IRT 

library(mirt)

#Unidimensionality

psych::unidim(icar[, 8:23]) 

#Monotonicity
library(mokken)
#Compute Monotonicity of Questionnaire Items
monotonicity_result <- check.monotonicity(icar[, 8:23])
monotonicity_result$violations #It detects the name of items that violates monotonicity. It returns NULL when no items violate the assumption.


#compute IRT models
mod_1pl <- mirt(icar[, c(8:23)], model= 1, itemtype = "1PL")

mod_2pl <- mirt(icar[, c(8:23)], model= 1, itemtype = "2PL")

mod_3pl <- mirt(icar[, c(8:23)], model= 1, itemtype = "3PL")

#model comparison: SABIC and BIC to prioritize parsimony
anova(mod_1pl, mod_2pl, mod_3pl)
M2(mod_1pl); M2(mod_2pl); M2(mod_3pl)


#Local Dependency
residuals2 <- residuals(mod_1pl, type="Q3") 
residuals_matrix<-cor.plot(residuals2, upper = FALSE)
# find the residual correlations |Q3| > 0.2 

cor_matrix <- residuals2 


which(abs(cor_matrix) > 0.2 & lower.tri(cor_matrix), arr.ind = TRUE)

coef(mod_1pl, IRTpars = TRUE, simplify= T)


plot(mod_1pl,type= "trace")# plot icc

#model fit
M2(mod_1pl)



#TIF
theta_vals <- seq(-4, 4, by = 0.01)
test_info <- testinfo(mod_1pl, theta_vals)
tif_df <- data.frame(theta = theta_vals, info = test_info)
theta<-seq(-4, 4, by=0.1)
tif_values<-mirt::testinfo(mod_1pl, Theta = theta)
tif_values

#tif values:  Peak at 40
which.max(tif_values) 
theta[40] <- (-4 + 0.1*(40-1)) 
print(theta[40])# theta= -0.1 more precision


#TIF
ggplot(tif_df, aes(x = theta, y = info)) +
  geom_line(color = "lightblue", size = 1.2) +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "grey40", linewidth = 1)+
  labs(
    title = "Test Information Function (TIF)",
    x = "Ability (θ)",
    y = "Information"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+ 
  scale_x_continuous(breaks = c(-4: +4))



#Test information reliability
theta <- seq(-4, 4, by = 0.01)  
tif_mirt <- mirt::testinfo(mod_1pl, Theta = theta)
total_info_mirt <- sum(tif_mirt) * 0.01  
print(total_info_mirt)  
percentuale_info <- (total_info_mirt / 40) * 100 #40 = (k*2.5)
print(percentuale_info)


theta_large <- seq(-10, 10, by = 0.01) 
tif_large <- mirt::testinfo(mod_1pl, Theta = theta_large)
total_info_full <- sum(tif_large) * 0.01  
print(total_info_full)


total_info_full <- 15.99
R <- 1 - 1/(1 + total_info_full)
R #94% of estimated variance is true variance


################################################################################

#13. uniform and not uniform DIF 
#DIF 
library(difR)

icar$sex3<- ifelse(icar$sex2=="F", 0, 1)


log2 <- difLogistic(
  Data = as.matrix(icar[, c(8:23)]), 
  group = icar$sex3, 
  focal.name = "0", 
  match = icar$ICAR, 
  type = "both", 
  p.adjust.method = "bonferroni"
)
summary(log2)

log2$names
log2$adjusted.p
log2$DIFitems

# 
dif_table2 <- data.frame(
  Item = log2$names,
  Adjusted_p = log2$adjusted.p
)

#
print(dif_table2)

dif_table2$DIF_detected <- ifelse(dif_table2$Adjusted_p < 0.05, "Yes", "No")

dif_table2$DIF_detected

dif_table12 <- data.frame(
  Item = dif_table2$Item,
  Adjusted_p = dif_table2$Adjusted_p, DIF= dif_table2$DIF_detected
)

print(dif_table12)


# Uniform DIF
log_udif <- difLogistic(
  Data = as.matrix(icar[, c(8:23)]),
  group = icar$sex3,
  focal.name = "0",
  match = icar$ICAR,
  type = "udif",
  p.adjust.method = "bonferroni"
)

# Non Uniform DIF
log_nudif <- difLogistic(
  Data = as.matrix(icar[, c(8:23)]),
  group = icar$sex3,
  focal.name = "0",
  match = icar$ICAR,
  type = "nudif",
  p.adjust.method = "bonferroni"
)

log_udif$adjusted.p[log_udif$names == "MR47"]
log_nudif$adjusted.p[log_nudif$names == "MR47"]

log_udif$logitPar[log_udif$names == "MR47", ]


plot(
  log_udif,
  plot = "itemCurve",
  item = "MR47",
  group.names = c("Male", "Female")
)

################################################################################
###############################################################################

#14. ICAR sex measurement invariance:


mgcfa1<- '  g=~ R3D3 + R3D4 + R3D6 + R3D8 + LN7 + LN33 + LN34 + LN58 + MR45 + MR46 + MR47 + MR55 + VR4 + VR17 + VR16 + VR19 
          '


#configural factorial invariance
mgcfa.configural<-cfa(mgcfa1, data=icar[, c(8:23, 32)], group = "sex2", ordered = T)
summary(mgcfa.configural, standardize= TRUE, ci= TRUE)

#weak factorial (metric) invariance
mgcfa.Weak <- cfa(mgcfa1, data=icar[, c(8:23, 32)], group = "sex2", ordered = T, group.equal="loadings")
summary(mgcfa.Weak, standardize= TRUE, ci= TRUE)

#strong factorial (scalar) invariance
mgcfa.strong <- cfa(mgcfa1, data=icar[, c(8:23, 32)], group = "sex2", ordered = T, group.equal= c("loadings", "intercepts"))
summary(mgcfa.strong, standardize= TRUE, ci= TRUE)

# Modification indices of scalar model
mi_strong <- modindices(mgcfa.strong, sort. = TRUE, maximum.number = 30)
modindices(mgcfa.strong, sort. = T, maximum.number = 10)

#Free intercept of MR47
partial_thresholds <- c("MR47|t1")

mgcfa.partial <- cfa(
  mgcfa1,
  data          = icar[, c(8:23, 32)],
  group         = "sex2",
  ordered       = TRUE,
  group.equal   = c("loadings", "intercepts"),
  group.partial = partial_thresholds
)

summary(mgcfa.partial, standardize = TRUE, ci = TRUE)
fitMeasures(mgcfa.partial, c("cfi.scaled", "tli.scaled", 
                             "rmsea.scaled", "srmr",
                             "chisq.scaled", "df.scaled"))

# Confronto con scalare completo
lavTestLRT(mgcfa.partial, mgcfa.strong)

#strict factorial invariance
mgcfa.strict <- cfa(mgcfa1, data=icar[, c(8:23, 32)], group = "sex2", ordered = T, group.equal= c("loadings", "intercepts", "residuals"), group.partial = partial_thresholds)
summary(mgcfa.strict, standardize= TRUE, ci= TRUE)

# configural, metric, scalar and strict fit indexes

fitmeasures(mgcfa.configural, c("cfi.scaled", "tli.scaled", "rmsea.scaled","rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled" ,"srmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))

fitmeasures(mgcfa.Weak, c("cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled" ,"srmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))

fitmeasures(mgcfa.strong, c("cfi.scaled", "tli.scaled", "rmsea.scaled","rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled" , "srmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))

fitMeasures(mgcfa.partial, c("cfi.scaled", "tli.scaled", 
                             "rmsea.scaled", "srmr", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "pvalue.scaled",
                             "chisq.scaled", "df.scaled"))
fitmeasures(mgcfa.strict, c("cfi.scaled", "tli.scaled", "rmsea.scaled","rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled" , "srmr", "chisq.scaled", "df.scaled", "pvalue.scaled"))

#cut point delta scaled cfi
delta_cfi1<- (0.896 - 0.885) 

delta_cfi2<-(0.875 - 0.896)

delta_cfi2p <-(0.880- 0.896)

delta_cfi3<- (0.880 -0.880) 

delta_cfi1#from configural to weak
delta_cfi2# from weak to strong
delta_cfi2p #from strong to partial strong
delta_cfi3# from strong to strict

#Dmacs

# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("ddueber/dmacs")

library(dmacs)
lavaan_dmacs(mgcfa.partial, RefGroup = "F")
________________________________________________________________________________
