
# install packages

#install.packages("survminer")

require(data.table)
require(survival)
require(ggplot2)
require(magrittr)
require(survminer)

# read in data
g = read.csv(file="C:/Users/hs_90/Dropbox/Coursera/Survival analysis in R/data/simulated data.csv",header=TRUE, sep=',')

setwd("C:/Users/hs_90/Dropbox/Coursera/Survival analysis in R/R script")

# To run these packages, we of course need some variables to put into them. My preferred way to do this is to turn each column of the data set, which we've called "g", into a variable and tell R what kind of variable it is. 

gender <- as.factor(g[,"gender"]) # R calls categorical variables factors
fu_time <- g[,"fu_time"] # continuous variable (numeric) - follow up time
death <- g[,"death"] # binary variable (numeric) 
age <- g[,"age"]
copd <- g[, "copd"]
prior_dnas <- g[,"prior_dnas"]

# Kaplan-Meier plot: 
km_fit <- survfit(Surv(fu_time, death) ~ 1)

plot(km_fit)

summary(km_fit, times = c(1:7,30,60,90*(1:10))) 


# lets extend this by splitting the curve by gender
km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 

plot(km_gender_fit)

summary(km_gender_fit, times = c(1:7,30,60,90*(1:10))) 

## to compare survival by gender we run a log-rank test
survdiff(Surv(fu_time, death) ~ gender, rho=0) # With rho = 0, which is the default so we don't need to write this bit, it yields the log-rank or Mantel-Haenszel test.

# compare by age grouping 65 and over
g$age_grp = ifelse(g$age >= 65, 1, 0)
age_grp = g[, "age_grp"]
table(g[,"age"], age_grp, exclude = NULL) # check the numbers 

km_age_fit = survfit(Surv(fu_time, death) ~ age_grp)
plot(km_age_fit)

survdiff(Surv(fu_time, death) ~ age_grp, rho=0)


# Simple Cox model 
cox = coxph(Surv(fu_time, death) ~ age, data=g)
summary(cox)

ethnicgroup = factor(g[,"ethnicgroup"])
cox2 = coxph(Surv(fu_time, death) ~ ethnicgroup)
summary(cox2)

# make a category for the missing ethnic group 
levels(ethnicgroup) = c(levels(ethnicgroup), '8')  # add level 8 to the factor
ethnicgroup[is.na(ethnicgroup)] = '8' # Change NA to "None"

# rerun Cox model
cox3 = coxph(Surv(fu_time, death) ~ ethnicgroup)
summary(cox3)

# explore the variable a bit
summary(g$age)
addmargins(table(g$ethnicgroup, exclude = NULL)) # addmargins to get total 
addmargins(table(ethnicgroup, exclude = NULL)) # where we've categorised missing as 8
prop.table(table(g$ethnicgroup, exclude=NULL)) # % of col total


cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + prior_dnas + ethnicgroup)
summary(cox)

# Non-convergence

quintile = factor(g[, "quintile"])
cox2 = coxph(Surv(fu_time,death) ~ age + gender + copd + quintile + ethnicgroup)
summary(cox2)

t = table(quintile, death)
t
prop.table(t,1) # % of row total


## change reference 
quintile <- relevel(quintile, ref = 2) # quintile 1 as the ref cat again

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile + ethnicgroup)

summary(cox)

## combine categories
quintile_5groups <- g[,'quintile'] # best start with the original data set, not from "quintile" 

quintile_5groups[quintile_5groups==0] <- 5 

quintile_5groups <- factor(quintile_5groups) 

table(quintile_5groups, exclude=NULL) 

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile_5groups + ethnicgroup) 

summary(cox) 


## drop quintile 0 patients
quintile_5groups <- g[,'quintile'] 

quintile_5groups[quintile_5groups==0] <- NA # set the zeroes to missing 

quintile_5groups <- factor(quintile_5groups) 

table(quintile_5groups, exclude=NULL) 

cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + quintile_5groups + ethnicgroup) 

summary(cox) 

# Test for proportionality assumption

cox3 = coxph(Surv(fu_time, death) ~ gender)
summary(cox3)

proportionality_assumption_test = cox.zph(cox3)
print(proportionality_assumption_test)

plot(proportionality_assumption_test)

# KM plot for gender
km_fit <- survfit(Surv(fu_time, death) ~ gender) 
plot(km_fit, xlab = "time", ylab = "Survival probability")


# Other types of residuals 

res.cox <- coxph(Surv(fu_time, death) ~ age) 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 


res.cox <- coxph(Surv(fu_time, death) ~ age) 
ggcoxdiagnostics(res.cox, type = "deviance", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 


# Martingale residual
require(survminer)
ggcoxfunctional(Surv(fu_time, death) ~ age + log(age) + sqrt(age)) 


# Test - testing the proportionality assumption with another variable

cox_test = coxph(Surv(fu_time, death) ~ copd)
summary(cox_test)

proportionality_assumption_test2 = cox.zph(cox_test)
print(proportionality_assumption_test2)

plot(proportionality_assumption_test2)


fit <- coxph(Surv(fu_time, death) ~ gender + tt(gender)) # "tt" is the time-transform function 
summary(fit) 


# Running a multiple Cox model 
ihd <- factor(g[,'ihd']) 

valvular <- factor(g[,'valvular_disease']) 

pvd <- factor(g[,'pvd']) 

stroke <- factor(g[,'stroke']) 

copd<- factor(g[,'copd'])

pneumonia <- factor(g[,'pneumonia']) 

ht <- factor(g[,'hypertension'])

renal <- factor(g[,'renal_disease']) 

ca <- factor(g[,'cancer']) 

mets <- factor(g[,'metastatic_cancer']) 

mental_health <- factor(g[,'mental_health']) 

los <- g[,'los']

prior_dna <- g[,'prior_dnas']

# generate cognitive impairment variable (senility and dementia combined)

cog_imp <- as.factor(ifelse(g$dementia == 1 | g$senile == 1, 1, 0))

# run the full model 

cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + 
               
               valvular + pvd + stroke + copd + pneumonia + ht + renal + 
               
               ca + mets + mental_health + cog_imp + los + prior_dna) 

summary(cox) 

## Keeping only the significant variables - backward elimination 
cox <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + mets + cog_imp) 
summary(cox) 

table(cog_imp) 
t <- table(cog_imp,death)
t
round(100*prop.table(t,1),digits=1) 


# Testing the proportionality assumption on the remaining variables
fit <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + 
               mets + cog_imp) # test them all in the same model 

temp <- cox.zph(fit)  

print(temp) 



#########################################################################################
# Final code 

# Install and load relevant packages 

install.packages("survival") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'survival' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  install.packages("ggplot2") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'ggplot2' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##  install.packages("survminer") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'survminer' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 
##   

install.packages("ggfortify") 

## Installing package into  
## (as 'lib' is unspecified) 

## package 'ggfortify' successfully unpacked and MD5 sums checked 
##  
## The downloaded binary packages are in 


require(survival)  

## Loading required package: survival 

## Warning: package 'survival' was built under R version 3.5.1 

require(ggplot2)  

## Loading required package: ggplot2 

## Warning: package 'ggplot2' was built under R version 3.5.1 

require(survminer) 

## Loading required package: survminer 

## Warning: package 'survminer' was built under R version 3.5.1 

## Loading required package: ggpubr 

## Warning: package 'ggpubr' was built under R version 3.5.1 

## Loading required package: magrittr 

require(ggfortify) 

## Loading required package: ggfortify 

## Warning: package 'ggfortify' was built under R version 3.5.1 

# Load dataset 
g <- read.csv(file = paste0(getwd(),"/simulated HF mort data for GMPH (1K) final.csv"), header=TRUE, sep=',') 


# Define variables 
gender <- factor(g[,"gender"]) 
fu_time <- g[,"fu_time"]  
death <-  g[,"death"] 
age <- g[,"age"] 
copd <- factor(g[,"copd"]) 
ethnicgroup <- factor(g[,"ethnicgroup"]) 
quintile <- factor(g[,"quintile"]) 
ihd <- factor(g[,'ihd']) 
valvular <- factor(g[,'valvular_disease']) 
pvd <- factor(g[,'pvd']) 
stroke <- factor(g[,'stroke']) 
pneumonia <- factor(g[,'pneumonia']) 
renal <- factor(g[,'renal_disease']) 
ca <- factor(g[,'cancer']) 
mets <- factor(g[,'metastatic_cancer']) 
mental_health <- factor(g[,'mental_health']) 
ht <- factor(g[,"hypertension"]) 
cog_imp <- factor(g[,"senile"]) 
prior_dnas <- g[,"prior_dnas"] 



# Plotting a Kaplan-Meier curve 
###################### 

# 1. Generate the survival curve 
km_fit <- survfit(Surv(fu_time, death) ~ 1) 

# 2b. Alternative plot with ggplot2 
autoplot(km_fit) + theme_bw() # theme_bw() is a predesigned "theme" which makes the plot prettier 

# Output the probability of survival at certain times after hospital admission 
summary(km_fit, times = c(1:7,30,60,90*(1:10))) 

# Plotting a Kaplan-Meier curve by gender 
###################### 

# 1. Generate the survival curve 
km_gender_fit <- survfit(Surv(fu_time, death) ~ gender) 

# 2. Plot the curve 
plot(km_gender_fit)

# 2b. Alternative plot with ggplot2 
autoplot(km_gender_fit) + theme_bw() 

# Perform log rank test to see whether survival varies by gender 
survdiff(Surv(fu_time, death) ~ gender, rho = 0) 

## Call: 
## survdiff(formula = Surv(fu_time, death) ~ gender, rho = 0) 
##  
##            N Observed Expected (O-E)^2/E (O-E)^2/V 
## gender=1 548      268      271    0.0365     0.082 
## gender=2 452      224      221    0.0448     0.082 
##  
##  Chisq= 0.1  on 1 degrees of freedom, p= 0.8 

# Testing whether those over the age of 65 have different survival to those under it 
###################### 

# 1. Dichotomise age into categorical (binary in this case) variable 
age_65plus <- ifelse(g[,'age']>=65, 1, 0) 

# 2. Perform log rank test 
survdiff(Surv(fu_time, death) ~ age_65plus, rho = 0) 

## Call: 
## survdiff(formula = Surv(fu_time, death) ~ age_65plus, rho = 0) 
##  
##                N Observed Expected (O-E)^2/E (O-E)^2/V 
## age_65plus=0 115       18       67     35.85      41.7 
## age_65plus=1 885      474      425      5.65      41.7 
##  
##  Chisq= 41.7  on 1 degrees of freedom, p= 1e-10 

###################### 

# Plot survival curve by age above or below 65 
###################### 

# 1. Generate survival curve 
km_old_fit <- survfit(Surv(fu_time, death) ~ age_65plus) 

# 2. Plot 
plot(km_old_fit) 

# 2b. Alternative plot in ggplot2 
autoplot(km_old_fit) + theme_bw()

###################### 

# Run Cox regression model with age as predictor (continuous variable) 
###################### 

# 1. Generate model 
cox <- coxph(Surv(fu_time, death) ~ age, data = g) 

# 2. Summarise model 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age, data = g) 
##  
##   n= 1000, number of events= 492  
##  
##         coef exp(coef) se(coef)     z Pr(>|z|)     
## age 0.056005  1.057602 0.005193 10.78   <2e-16 *** 
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##     exp(coef) exp(-coef) lower .95 upper .95 
## age     1.058     0.9455     1.047     1.068 
##  
## Concordance= 0.651  (se = 0.015 ) 
## Rsquare= 0.129   (max possible= 0.997 ) 
## Likelihood ratio test= 138  on 1 df,   p=<2e-16 
## Wald test            = 116.3  on 1 df,   p=<2e-16 
## Score (logrank) test = 115.7  on 1 df,   p=<2e-16 

###################### 



# Run Cox regression model with quintile as predictor (categorical variable) 
# Changing the reference group to first quintile 
# Removing the zero quintile altogether 
###################### 


# 1. Summarise the variable 
table(quintile, exclude = NULL) 

## quintile 
##    0    1    2    3    4    5 <NA>  
##    4  138  205  211  220  216    6 

# 2. Check levels 
levels(quintile) 

## [1] "0" "1" "2" "3" "4" "5" 

# 3. Generate model 
cox <- coxph(Surv(fu_time, death) ~ quintile) # warning 

## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, : 
## Loglik converged before variable 1,2,3,4,5 ; beta may be infinite. 

# 4. Summarise model 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ quintile) 
##  
##   n= 994, number of events= 489  
##    (6 observations deleted due to missingness) 
##  
##                coef exp(coef)  se(coef)     z Pr(>|z|) 
## quintile1 1.511e+01 3.665e+06 1.202e+03 0.013     0.99 
## quintile2 1.479e+01 2.645e+06 1.202e+03 0.012     0.99 
## quintile3 1.518e+01 3.912e+06 1.202e+03 0.013     0.99 
## quintile4 1.500e+01 3.261e+06 1.202e+03 0.012     0.99 
## quintile5 1.503e+01 3.383e+06 1.202e+03 0.013     0.99 
##  
##           exp(coef) exp(-coef) lower .95 upper .95 
## quintile1   3665107  2.728e-07         0       Inf 
## quintile2   2645345  3.780e-07         0       Inf 
## quintile3   3911559  2.557e-07         0       Inf 
## quintile4   3261463  3.066e-07         0       Inf 
## quintile5   3382505  2.956e-07         0       Inf 
##  
## Concordance= 0.544  (se = 0.015 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 13.22  on 5 df,   p=0.02 
## Wald test            = 8.39  on 5 df,   p=0.1 
## Score (logrank) test = 10.79  on 5 df,   p=0.06 

# 5. Make the first quintile the reference group 
quintile <- relevel(quintile, ref = "1") 

# 6. Regenerate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ quintile) # warning 

## Warning in fitter(X, Y, strats, offset, init, control, weights = weights, : 
## Loglik converged before variable 1 ; beta may be infinite. 

summary(cox) # still an issue where quintile = 0 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ quintile) 
##  
##   n= 994, number of events= 489  
##    (6 observations deleted due to missingness) 
##  
##                 coef  exp(coef)   se(coef)      z Pr(>|z|)   
## quintile0 -1.511e+01  2.728e-07  1.202e+03 -0.013   0.9900   
## quintile2 -3.261e-01  7.218e-01  1.537e-01 -2.121   0.0339 * 
## quintile3  6.508e-02  1.067e+00  1.495e-01  0.435   0.6633   
## quintile4 -1.167e-01  8.899e-01  1.501e-01 -0.777   0.4369   
## quintile5 -8.024e-02  9.229e-01  1.490e-01 -0.538   0.5903   
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##           exp(coef) exp(-coef) lower .95 upper .95 
## quintile0 2.728e-07  3.665e+06    0.0000       Inf 
## quintile2 7.218e-01  1.385e+00    0.5340    0.9756 
## quintile3 1.067e+00  9.370e-01    0.7962    1.4306 
## quintile4 8.899e-01  1.124e+00    0.6631    1.1942 
## quintile5 9.229e-01  1.084e+00    0.6891    1.2360 
##  
## Concordance= 0.544  (se = 0.015 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 13.22  on 5 df,   p=0.02 
## Wald test            = 8.39  on 5 df,   p=0.1 
## Score (logrank) test = 10.79  on 5 df,   p=0.06 

# 7. Inspecting quintile variable 
table(quintile, g$death) # Only 4 entries for quintile = 0 and 100% didn't die 

##          
## quintile   0   1 
##        1  60  78 
##        0   4   0 
##        2 111  94 
##        3 105 106 
##        4 116 104 
##        5 109 107 

# 8. Removing quintile = 0 entries as there are only 4 of them 
quintile_5groups <- quintile 
quintile_5groups[quintile_5groups == 0] <- NA # set the zeroes to missing 
quintile_5groups <- factor(quintile_5groups) # this removes 0 as a level as it is an empty category 

# 9. Regenerating the model and summarising 
cox <- coxph(Surv(fu_time, death) ~ quintile_5groups) 
summary(cox) # still an issue where quintile = 0 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ quintile_5groups) 
##  
##   n= 990, number of events= 489  
##    (10 observations deleted due to missingness) 
##  
##                       coef exp(coef) se(coef)      z Pr(>|z|)   
## quintile_5groups2 -0.32606   0.72176  0.15374 -2.121   0.0339 * 
## quintile_5groups3  0.06508   1.06724  0.14949  0.435   0.6633   
## quintile_5groups4 -0.11668   0.88987  0.15010 -0.777   0.4369   
## quintile_5groups5 -0.08024   0.92289  0.14903 -0.538   0.5903   
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##                   exp(coef) exp(-coef) lower .95 upper .95 
## quintile_5groups2    0.7218      1.385    0.5340    0.9756 
## quintile_5groups3    1.0672      0.937    0.7962    1.4306 
## quintile_5groups4    0.8899      1.124    0.6631    1.1942 
## quintile_5groups5    0.9229      1.084    0.6891    1.2360 
##  
## Concordance= 0.542  (se = 0.015 ) 
## Rsquare= 0.009   (max possible= 0.997 ) 
## Likelihood ratio test= 8.62  on 4 df,   p=0.07 
## Wald test            = 8.39  on 4 df,   p=0.08 
## Score (logrank) test = 8.45  on 4 df,   p=0.08 

###################### 

# Run Cox regression model with ethnic group as predictor (categorical variable) 
# Including missing values as another category 
###################### 

# 1. Summarise variable 
table(ethnicgroup, exclude = NULL) 

## ethnicgroup 
##    1    2    3    9 <NA>  
##  889   17   34   17   43 

# 2. Generate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ ethnicgroup) 
##  
##   n= 957, number of events= 471  
##    (43 observations deleted due to missingness) 
##  
##                  coef exp(coef) se(coef)      z Pr(>|z|)    
## ethnicgroup2 -0.06428   0.93774  0.32000 -0.201  0.84078    
## ethnicgroup3 -1.19586   0.30244  0.41108 -2.909  0.00362 ** 
## ethnicgroup9  0.07394   1.07674  0.35706  0.207  0.83596    
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##              exp(coef) exp(-coef) lower .95 upper .95 
## ethnicgroup2    0.9377     1.0664    0.5008    1.7558 
## ethnicgroup3    0.3024     3.3064    0.1351    0.6769 
## ethnicgroup9    1.0767     0.9287    0.5348    2.1679 
##  
## Concordance= 0.516  (se = 0.007 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 12.99  on 3 df,   p=0.005 
## Wald test            = 8.55  on 3 df,   p=0.04 
## Score (logrank) test = 9.61  on 3 df,   p=0.02 

# 3. Add another category (8) 
levels(ethnicgroup) <- c(levels(ethnicgroup),"8")  

# 4. Redefine NA as another group, 8 
ethnicgroup[is.na(ethnicgroup)] <- "8" 

# 5. Regenerate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ ethnicgroup) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ ethnicgroup) 
##  
##   n= 1000, number of events= 492  
##  
##                  coef exp(coef) se(coef)      z Pr(>|z|)    
## ethnicgroup2 -0.06573   0.93638  0.31999 -0.205  0.83725    
## ethnicgroup3 -1.19368   0.30310  0.41107 -2.904  0.00369 ** 
## ethnicgroup9  0.08160   1.08502  0.35706  0.229  0.81923    
## ethnicgroup8 -0.02353   0.97675  0.22363 -0.105  0.91621    
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##              exp(coef) exp(-coef) lower .95 upper .95 
## ethnicgroup2    0.9364     1.0679    0.5001    1.7532 
## ethnicgroup3    0.3031     3.2992    0.1354    0.6784 
## ethnicgroup9    1.0850     0.9216    0.5389    2.1846 
## ethnicgroup8    0.9767     1.0238    0.6301    1.5140 
##  
## Concordance= 0.518  (se = 0.008 ) 
## Rsquare= 0.013   (max possible= 0.997 ) 
## Likelihood ratio test= 12.95  on 4 df,   p=0.01 
## Wald test            = 8.53  on 4 df,   p=0.07 
## Score (logrank) test = 9.58  on 4 df,   p=0.05 

###################### 



# Investigating our variables in order to best perform a Cox model with multiple predictors 
# Checking for missing values 
# Running a multiple Cox regression 
###################### 

# 1. Summarising age 
summary(g$age) # no NAs 

##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  
##   29.00   73.00   80.00   78.73   87.00  102.00 

# 2. Gender 
gender_table <- table(gender, exclude = NULL) 
addmargins(gender_table) # no NAs 

## gender 
##    1    2  Sum  
##  548  452 1000 

round(100 * prop.table(gender_table), digits = 1) # Percentages rounded to 1 decimal place 

## gender 
##    1    2  
## 54.8 45.2 

# 3. Chronic Obstructive Pulmonary Disease (COPD) 
copd_table <- table(copd, exclude = NULL)  
addmargins(copd_table) # no NAs 

## copd 
##    0    1  Sum  
##  758  242 1000 

round(100 * prop.table(copd_table), digits = 1) # Percentages rounded to 1 decimal place 

## copd 
##    0    1  
## 75.8 24.2 

# 4. Prior OPD appointments missed  
prior_dnas_table <- table(prior_dnas, exclude = NULL)  
addmargins(prior_dnas_table) # no NAs 

## prior_dnas 
##    0    1    2    3    4    5    6    7    8   10  Sum  
##  732  156   50   34   17    3    3    2    1    2 1000 

round(100 * prop.table(prior_dnas_table), digits = 1) # Percentages rounded to 1 decimal place 

## prior_dnas 
##    0    1    2    3    4    5    6    7    8   10  
## 73.2 15.6  5.0  3.4  1.7  0.3  0.3  0.2  0.1  0.2 

# 5. Ethnic group 
ethnicgroup_table <- table(ethnicgroup, exclude = NULL)  
addmargins(ethnicgroup_table) # 4.3% NA 

## ethnicgroup 
##    1    2    3    9    8  Sum  
##  889   17   34   17   43 1000 

round(100 * prop.table(ethnicgroup_table), digits = 1) # Percentages rounded to 1 decimal place 

## ethnicgroup 
##    1    2    3    9    8  
## 88.9  1.7  3.4  1.7  4.3 

# 6. Generate and summarise model 
cox <- coxph(Surv(fu_time, death) ~ age + gender + copd + prior_dnas + ethnicgroup) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age + gender + copd +  
##     prior_dnas + ethnicgroup) 
##  
##   n= 1000, number of events= 492  
##  
##                   coef exp(coef)  se(coef)      z Pr(>|z|)     
## age           0.061999  1.063961  0.005516 11.241  < 2e-16 *** 
## gender2      -0.253460  0.776111  0.094349 -2.686  0.00722 **  
## copd1         0.136649  1.146425  0.103880  1.315  0.18836     
## prior_dnas    0.163461  1.177579  0.039832  4.104 4.07e-05 *** 
## ethnicgroup2 -0.307915  0.734978  0.353009 -0.872  0.38307     
## ethnicgroup3 -0.823643  0.438830  0.414301 -1.988  0.04681 *   
## ethnicgroup9  0.408255  1.504190  0.360737  1.132  0.25775     
## ethnicgroup8 -0.045372  0.955642  0.225204 -0.201  0.84033     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##              exp(coef) exp(-coef) lower .95 upper .95 
## age             1.0640     0.9399    1.0525    1.0755 
## gender2         0.7761     1.2885    0.6451    0.9338 
## copd1           1.1464     0.8723    0.9352    1.4053 
## prior_dnas      1.1776     0.8492    1.0891    1.2732 
## ethnicgroup2    0.7350     1.3606    0.3680    1.4681 
## ethnicgroup3    0.4388     2.2788    0.1948    0.9884 
## ethnicgroup9    1.5042     0.6648    0.7417    3.0504 
## ethnicgroup8    0.9556     1.0464    0.6146    1.4859 
##  
## Concordance= 0.667  (se = 0.015 ) 
## Rsquare= 0.155   (max possible= 0.997 ) 
## Likelihood ratio test= 168.4  on 8 df,   p=<2e-16 
## Wald test            = 141.7  on 8 df,   p=<2e-16 
## Score (logrank) test = 140  on 8 df,   p=<2e-16 

###################### 


# Investigating whether the assumptions of the Cox model are being broken 
# Testing for proportional hazards assumption (with gender as predictor variable) 
###################### 

# 1. Generate model fit 
fit <- coxph(Surv(fu_time, death) ~ gender) 

# 2. Apply the test to the model 
temp <- cox.zph(fit)     

# 3. Display results 
print(temp) 

##            rho chisq     p 
## gender2 0.0493  1.19 0.275 

# 4. Plot the curves 
plot(temp)     
hist(g$age) 

# 4b. Alternative plot in ggplot 
ggcoxzph(temp) 

###################### 

# Generating other diagnostic plots for Cox Proportional Hazards model 
###################### 

# 1. Define model 
res.cox <- coxph(Surv(fu_time, death) ~ age) 

# Generate diagnostic plots 

# 2. Plotting the estimated changes in the regression coefficients on deleting each patient 
ggcoxdiagnostics(res.cox, type = "dfbeta", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

# 3. Plotting deviance residuals 
ggcoxdiagnostics(res.cox, type = "deviance", 
                 linear.predictions = FALSE, ggtheme = theme_bw()) 

# 4. Plotting Martingale residuals 
fit <- coxph(Surv(fu_time, death) ~ age + log(age) + sqrt(age)) 
ggcoxfunctional(fit, data = g) # note we must specify original dataframe 

###################### 

# Testing proportionality assumption 
# Testing for a statistical relationship between gender and time 
###################### 

# 1. Generate model with time-transform function (tt) 
fit <- coxph(Surv(fu_time, death) ~ gender + tt(gender))  

# 2. Summarise 
summary(fit) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ gender + tt(gender)) 
##  
##   n= 1000, number of events= 492  
##  
##               coef exp(coef) se(coef)      z Pr(>|z|) 
## gender2     0.6405    1.8974   0.9800  0.654    0.513 
## tt(gender) -0.2003    0.8185   0.3182 -0.629    0.529 
##  
##            exp(coef) exp(-coef) lower .95 upper .95 
## gender2       1.8974      0.527    0.2779    12.953 
## tt(gender)    0.8185      1.222    0.4387     1.527 
##  
## Concordance= 0.497  (se = 0.197 ) 
## Rsquare= 0   (max possible= 0.997 ) 
## Likelihood ratio test= 0.49  on 2 df,   p=0.8 
## Wald test            = 0.48  on 2 df,   p=0.8 
## Score (logrank) test = 0.48  on 2 df,   p=0.8 

###################### 

# Backwards elimination to choose predictors for Cox regression 
###################### 



# 1. Run the full model with all of your predictors 
cox <- coxph(Surv(fu_time, death) ~ age + gender + ethnicgroup + ihd + 
               valvular + pvd + stroke + copd + pneumonia + ht + renal + 
               ca + mets + mental_health + cog_imp + los + prior_dna) 
summary(cox) 

## Call:
## coxph(formula = Surv(fu_time, death) ~ age + gender + ethnicgroup + 
##     ihd + valvular + pvd + stroke + copd + pneumonia + ht + renal + 
##     ca + mets + mental_health + cog_imp + los + prior_dna)
##
##   n= 1000, number of events= 492 
##
##                     coef exp(coef)  se(coef)      z Pr(>|z|)    
## age             0.058348  1.060083  0.005801 10.058  < 2e-16 ***
## gender2        -0.215988  0.805745  0.097592 -2.213  0.02688 *  
## ethnicgroup2   -0.159816  0.852301  0.351947 -0.454  0.64976    
## ethnicgroup3   -0.724770  0.484436  0.416059 -1.742  0.08151 .  
## ethnicgroup9    0.437472  1.548787  0.363859  1.202  0.22924    
## ethnicgroup8   -0.086505  0.917131  0.231729 -0.373  0.70892    
## ihd1            0.174550  1.190710  0.096181  1.815  0.06955 .  
## valvular1       0.194524  1.214732  0.108805  1.788  0.07380 .  
## pvd1            0.035442  1.036077  0.161851  0.219  0.82667    
## stroke1        -0.007433  0.992594  0.301432 -0.025  0.98033    
## copd1           0.103113  1.108616  0.106823  0.965  0.33441    
## pneumonia1      0.302242  1.352889  0.141232  2.140  0.03235 *  
## ht1            -0.057556  0.944069  0.096372 -0.597  0.55036    
## renal1          0.153977  1.166464  0.109547  1.406  0.15985    
## ca1             0.280847  1.324251  0.206595  1.359  0.17402    
## mets1           2.194757  8.977824  0.399639  5.492 3.98e-08 ***
## mental_health1 -0.066662  0.935511  0.181811 -0.367  0.71387    
## cog_imp1        0.327423  1.387388  0.149461  2.191  0.02847 *  
## los             0.011553  1.011620  0.003257  3.547  0.00039 ***
## prior_dna       0.107472  1.113460  0.041807  2.571  0.01015 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
##
##                exp(coef) exp(-coef) lower .95 upper .95
## age               1.0601     0.9433    1.0481    1.0722
## gender2           0.8057     1.2411    0.6655    0.9756
## ethnicgroup2      0.8523     1.1733    0.4276    1.6989
## ethnicgroup3      0.4844     2.0643    0.2143    1.0949
## ethnicgroup9      1.5488     0.6457    0.7591    3.1602
## ethnicgroup8      0.9171     1.0904    0.5823    1.4444
## ihd1              1.1907     0.8398    0.9861    1.4377
## valvular1         1.2147     0.8232    0.9814    1.5035
## pvd1              1.0361     0.9652    0.7544    1.4229
## stroke1           0.9926     1.0075    0.5498    1.7921
## copd1             1.1086     0.9020    0.8992    1.3668
## pneumonia1        1.3529     0.7392    1.0258    1.7843
## ht1               0.9441     1.0592    0.7816    1.1403
## renal1            1.1665     0.8573    0.9411    1.4458
## ca1               1.3243     0.7551    0.8833    1.9853
## mets1             8.9778     0.1114    4.1020   19.6492
## mental_health1    0.9355     1.0689    0.6551    1.3360
## cog_imp1          1.3874     0.7208    1.0351    1.8596
## los               1.0116     0.9885    1.0052    1.0181
## prior_dna         1.1135     0.8981    1.0259    1.2085
##
## Concordance= 0.707  (se = 0.012 )
## Likelihood ratio test= 225.9  on 20 df,   p=<2e-16
## Wald test            = 216.9  on 20 df,   p=<2e-16
## Score (logrank) test = 235.1  on 20 df,   p=<2e-16


# 2. Run the model with only significant predictors 
cox <- coxph(Surv(fu_time, death) ~ age + gender + valvular + pneumonia + mets + cog_imp) 
summary(cox) 

## Call: 
## coxph(formula = Surv(fu_time, death) ~ age + gender + valvular +  
##     pneumonia + mets + cog_imp) 
##  
##   n= 1000, number of events= 492  
##  
##                 coef exp(coef)  se(coef)      z Pr(>|z|)     
## age         0.059147  1.060931  0.005441 10.871  < 2e-16 *** 
## gender2    -0.260853  0.770394  0.093497 -2.790  0.00527 **  
## valvular1   0.229659  1.258170  0.107118  2.144  0.03204 *   
## pneumonia1  0.481375  1.618299  0.135192  3.561  0.00037 *** 
## mets1       2.494299 12.113233  0.364373  6.845 7.62e-12 *** 
## cog_imp1    0.241623  1.273314  0.184232  1.312  0.18968     
## --- 
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
##  
##            exp(coef) exp(-coef) lower .95 upper .95 
## age           1.0609    0.94257    1.0497    1.0723 
## gender2       0.7704    1.29804    0.6414    0.9253 
## valvular1     1.2582    0.79480    1.0199    1.5521 
## pneumonia1    1.6183    0.61793    1.2416    2.1093 
## mets1        12.1132    0.08255    5.9307   24.7409 
## cog_imp1      1.2733    0.78535    0.8874    1.8271 
##  
## Concordance= 0.686  (se = 0.015 ) 
## Rsquare= 0.17   (max possible= 0.997 ) 
## Likelihood ratio test= 186.9  on 6 df,   p=<2e-16 
## Wald test            = 179.5  on 6 df,   p=<2e-16 
## Score (logrank) test = 193  on 6 df,   p=<2e-16 

# 3. Test proportionality assumption on these predictors 
cox.zph(cox)  

##                 rho  chisq     p 
## age        -0.04409 0.9098 0.340 
## gender2     0.05474 1.4133 0.235 
## valvular1  -0.04565 1.0202 0.312 
## pneumonia1 -0.04397 0.9540 0.329 
## mets1       0.00615 0.0184 0.892 
## cog_imp1    0.05985 1.8140 0.178 
## GLOBAL           NA 5.4044 0.493 