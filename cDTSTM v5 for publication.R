###########################################################################################################
# Example analysis for paper ASSESSMENT OF VALUE OF NEW INTERVENTIONS SHOULD INCLUDE HEALTH EQUITY IMPACT #
#                                                                                                         #
# An aducanumab-informed depiction of the concept of DCEA with health equity impact evaluated across      #
# social subgroups defined according to race/ethnicity. Estimates were obtained by using                  #
# aducanumab-specific relative treatment effects, and age and race-specific mild cognitive impairment     #
# prevalence and background progression rates in the open-source health economic model for AD by          #
# Green et al 2019                                                                                        #  
#                                                                                                         #
# Time-inhomogenous Markov state-transition model implemented in HESIM                                    #                                                           #  
###########################################################################################################

rm(list=ls())


#devtools::install_github("hesim-dev/hesim")
library(hesim)
library(survival)
library(flexsurv)
library(icenReg)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)







ad_model <- function(
sexes = c("Female", "Male"),
ages = seq(60, 80, 10),
races = c("White", "African American", "Asian", "Hispanic"),
weights = c(0.154, 0.105, 0.066, 0.053, 0.029, 0.016, 0.022, 0.013, 0.007, 0.046, 0.025, 0.014, 
            0.143, 0.091, 0.043, 0.043, 0.02, 0.008, 0.018, 0.01, 0.005, 0.041, 0.019, 0.009),   # Distribution of MCI target population across age and race
tx_RR_input = c(0.776, 0.828, 0.391, 0.828,
                0.622, 0.541, 0.132, 0.541,
                0.967, 1.267, 1.157, 1.267),
tx_cost = 10000,
tx_duration = 30,
follow_up = 30,
PSA = 125,
wtp_dcea = 100000,
atkinson_dcea = 11,
kolm_dcea = 0.15
)
{

 
  

################  
# MODEL SET-UP #
################  
  
strategies <- data.frame(
  strategy_id = 1:2,
  strategy_name = c("NoTx", "Tx")
)
patient_id <- c(1:(length(sexes)*length(races)*length(ages)))
grp_id <- c(1:(length(sexes)*length(races)*length(ages)))
sex <- rep(sexes, each = length(ages)*length(races))
race <- rep(races, each = length(ages))
age <- rep(ages, times = length(sexes)*length(races))
patient_wt <- weights[1:length(patient_id)]/sum(weights[1:length(patient_id)])
patient_name <- str_c(sex, " ",race, " ", age)

patients <- data.table(patient_id,grp_id,sex,race,age,patient_name,patient_wt)


hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients
)



##############
# PARAMETERS #
##############

# Pre-dementia disease progression; from MCI to AD dementia
## Data from Green et al 2019 (Supplement 1) NACC data 2005-2017
tdata <- data.frame(time  =c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,9),
                    time2 =c(1,Inf,2,Inf,3,Inf,4,Inf,5,Inf,6,Inf,7,Inf,8,Inf,Inf,Inf),
                    n     =c(648,646,364,452,191,303,97,274,55,182,44,123,13,90,6,54,10,1))

fit <- ic_np(cbind(tdata$time, tdata$time2), data = tdata, weights = tdata$n)

## Fit parametric survival function (Weibull)
icbfit <- ic_bayes(Surv(time, time2, type = 'interval2') ~ 1, dist = "weibull",
                   data = tdata, weights = tdata$n)
summary(icbfit)

## Create matrix for coefficients of survival function parameters
icbfit_log_shape <- as.matrix(unlist(icbfit$mcmcList[1:1000,1,]))

### transform scale from icenREG (i.e  S(t)=exp(-(t/scale)^shape) ) to WeibullPH notation (i.e. S(t)=exp(-scale*t^shape) )
icbfit_log_scale <- as.matrix(unlist(icbfit$mcmcList[1:1000,2,]))
icbfit_log_scale <- log((1/exp(icbfit_log_scale))^(exp(icbfit_log_shape)))

log_shape <- mean(icbfit_log_shape)
log_scale <- mean(icbfit_log_scale)

## add effect of age and race to MCI to AD conversion
### effect of age with log HRs (adjustment of scale parameter; ignoring uncertainty in HRs based on limited information available in Green et al 2019)
age_HR <- matrix(c(1, 1.18, 1.39, 1.45, 1.68, 1.66),
                 ncol = 6, nrow = 1, byrow = TRUE)
age_logHR <- log(age_HR)

age_dist <- c(252, 515, 694, 843, 799, 450) ### distribution over age cohorts according to Green et al, 2019
average_HR <- sum(age_dist * age_HR)/sum(age_dist)
log_scale <- log_scale-log(average_HR) ### revised log scale, now corresponding to 60-64 age group

### Effect of race with log HRs (no adjustment of scale parameter because unclear race distribution in Green et al 2019. ~3/4 is white in NACC, assume log-scale relfects white)
### Based on NACC: Michaud TL, Su D, Siahpush M, Murman DL. The Risk of Incident Mild Cognitive Impairment and Progression to Dementia Considering Mild Cognitive Impairment Subtypes. Dement Geriatr Cogn Dis Extra. 2017 Feb 2;7(1):15-29. 
### Hispanic assumed the same as White due to no info in Michaud
race_HR <-   matrix(c(1, 1.39, 1.19, 1,
                      1, 0.60, 0.14, 1,
                      1, 3.22, 9.97, 1),
                    ncol = 4, nrow = 3, byrow = TRUE)
colnames(race_HR) <- c("White", "African American", "Asian", "Hispanic")
rownames(race_HR) <- c("HR", "95%CI-low", "95%CI-high")
race_logHR <- log(race_HR[1,])

### combine estimates
MCItoADdem_coeff <- matrix(c(log_shape, log_scale, age_logHR[2], age_logHR[3], age_logHR[4], age_logHR[5], age_logHR[6], race_logHR[2], race_logHR[3], race_logHR[4]),
                           ncol = 10,nrow = 1, byrow = TRUE)
colnames(MCItoADdem_coeff) <- c("log_shape", "log_scale", "age6569", "age7074", "age7579", "age8084", "age8589", "African American", "Asian", "Hispanic")
rownames(MCItoADdem_coeff) <- NULL

### Create variance-covariance matrix 
sd_log_shape <- sd(icbfit_log_shape)
sd_log_scale <- sd(icbfit_log_scale)
cor <- cor(icbfit_log_shape,icbfit_log_scale)

race_logHR_var <- ((log(race_HR[3,])-log(race_HR[2,]))/3.92)^2

MCItoADdem_vcov <- matrix(c(
                            sd_log_shape^2, sd_log_shape * sd_log_scale * cor, rep(0,8),
                            sd_log_shape * sd_log_scale * cor, sd_log_scale^2, rep(0,8),
                            rep(0,50),
                            rep(0,7),race_logHR_var[2],0,0,
                            rep(0,7),0,race_logHR_var[3],0,
                            rep(0,7),0,0,race_logHR_var[4]
                            ),
                            ncol = 10, nrow = 10, byrow = TRUE)
colnames(MCItoADdem_vcov) <- rownames(MCItoADdem_vcov) <-  c("log_shape", "log_scale", "age6569", "age7074", "age7579", "age8084", "age8589", "African American", "Asian", "Hispanic") 



# Proportions of landing states at conversion to dementia 
## Data from Green et al 2019 (Supplement 1)
landingstates_ADdem <- matrix(
  c(163, 34, 14, 229, 125, 42, 16, 15, 7, 16, 6, 2, 38, 26, 8, 9, 4, 4, 0, 0, 0, 1, 0, 0, 0, 1, 0),
  nrow = 1, byrow = TRUE)

state_names <- c(
  "ADdem111", "ADdem121", "ADdem131", "ADdem112", "ADdem122", "ADdem132", "ADdem113", "ADdem123", "ADdem133", "ADdem211",
  "ADdem221", "ADdem231", "ADdem212", "ADdem222", "ADdem232", "ADdem213", "ADdem223", "ADdem233", "ADdem311", "ADdem321",  
  "ADdem331","ADdem312", "ADdem322", "ADdem332", "ADdem313", "ADdem323", "ADdem333")
colnames(landingstates_ADdem) <- state_names
rownames(landingstates_ADdem) <- NULL



# AD-Dementia transitions
## Data from Green et al 2019 (Supplement 2)
transitions_ADdem <- matrix(
  c(
    247, 37, 3, 203, 73, 20, 11, 11, 0, 31, 4, 1, 86, 25, 8, 8, 2, 4, 1, 0, 0, 2, 0, 0, 0, 0, 1,
    32, 16, 2, 43, 39, 17, 3, 4, 0, 5, 2, 2, 9, 10, 7, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 5, 2, 6, 9, 8, 2, 0, 1, 3, 1, 1, 1, 5, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    48, 12, 0, 687, 255, 44, 145, 67, 25, 5, 2, 0, 261, 84, 24, 136, 98, 37, 0, 0, 0, 3, 0, 0, 5, 8, 5,
    15, 9, 5, 228, 228, 63, 45, 83, 30, 0, 0, 0, 74, 80, 33, 43, 75, 49, 0, 0, 0, 0, 3, 0, 2, 7, 7,
    6, 5, 1, 44, 64, 49, 17, 21, 31, 0, 0, 2, 4, 21, 16, 11, 26, 20, 0, 0, 0, 1, 0, 0, 2, 2, 0,
    0, 0, 0, 25, 5, 4, 63, 40, 23, 0, 0, 0, 15, 5, 1, 95, 36, 19, 0, 0, 0, 0, 0, 0, 6, 3, 2,
    0, 0, 0, 12, 9, 3, 45, 65, 24, 0, 1, 0, 12, 10, 4, 44, 70, 39, 0, 0, 0, 0, 0, 0, 2, 7, 3,
    0, 0, 0, 2, 5, 8, 18, 25, 42, 0, 1, 0, 3, 1, 1, 11, 27, 34, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    10, 0, 0, 11, 3, 2, 3, 0, 0, 11, 3, 0, 27, 8, 3, 7, 3, 2, 0, 0, 0, 6, 1, 0, 3, 0, 1,
    0, 0, 0, 1, 2, 1, 0, 1, 0, 2, 0, 0, 5, 2, 2, 1, 2, 3, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 3, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0,
    2, 0, 1, 57, 15, 6, 12, 11, 4, 5, 0, 0, 230, 76, 15, 134, 80, 32, 0, 0, 0, 21, 19, 2, 41, 35, 13,
    1, 0, 1, 17, 18, 6, 4, 13, 4, 3, 1, 0, 62, 56, 17, 46, 83, 41, 1, 0, 0, 4, 5, 3, 18, 20, 18,
    1, 0, 0, 5, 7, 6, 4, 4, 1, 0, 1, 1, 11, 20, 17, 17, 17, 39, 0, 0, 0, 1, 1, 2, 5, 5, 6,
    1, 0, 0, 7, 1, 0, 25, 12, 4, 1, 0, 0, 28, 8, 0, 228, 118, 43, 0, 0, 0, 3, 2, 0, 63, 53, 21,
    0, 0, 0, 2, 4, 1, 8, 17, 9, 1, 1, 0, 9, 13, 5, 84, 149, 90, 1, 0, 0, 1, 1, 2, 44, 59, 39,
    0, 0, 0, 3, 0, 0, 4, 4, 8, 0, 0, 0, 4, 4, 1, 24, 57, 101, 0, 0, 0, 1, 0, 2, 16, 26, 52,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 3, 2, 0, 11, 8, 8,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 1, 3, 1, 3, 8, 4,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 2, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 2, 0, 0, 78, 43, 29,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 5, 1, 0, 0, 0, 1, 2, 0, 39, 76, 44,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 1, 0, 0, 20, 28, 73
),
nrow = 27, byrow = TRUE)
colnames(transitions_ADdem) <- rownames(transitions_ADdem) <- state_names



# Mortality
## U.S. Census Bureau. Table NC-EST2019-ASR5H Annual Estimates of the Resident Population by Sex, Age, Race Alone or in Combination, and Hispanic Origin for the United States: 
## April 1, 2010 to July 1, 2019. Retrieved from https://www.census.gov/newsroom/press-kits/2020/population-estimates-detailed.html
mort_tbl <- matrix(
  c(
  50, 0.003159, 0.004967, 0.004628, 0.007450, 0.001895602, 0.002533271, 0.002076, 0.003693201,
  51, 0.003422, 0.005381, 0.005018, 0.008062, 0.002053318, 0.002744360, 0.002289, 0.004024453,
  52, 0.003734, 0.005893, 0.005474, 0.008786, 0.002240421, 0.003005408, 0.002513, 0.004401453,
  53, 0.004092, 0.006504, 0.005997, 0.009630, 0.002455462, 0.003317027, 0.002741, 0.004832045,
  54, 0.004473, 0.007169, 0.006565, 0.010559, 0.002684041, 0.003655973, 0.002975, 0.005311445,
  55, 0.004857, 0.007832, 0.007138, 0.011504, 0.002962632, 0.004150851, 0.003218, 0.005830034,
  56, 0.005234, 0.008482, 0.007718, 0.012466, 0.003192734, 0.004495342, 0.003479, 0.006377555,
  57, 0.005612, 0.009151, 0.008339, 0.013520, 0.003423602, 0.004850247, 0.003768, 0.006954584,
  58, 0.006002, 0.009863, 0.009023, 0.014710, 0.003661426, 0.005227289, 0.004093, 0.007561236,
  59, 0.006416, 0.010626, 0.009772, 0.016039, 0.003913770, 0.005631660, 0.004456, 0.008204376,
  60, 0.006871, 0.011452, 0.010579, 0.017471, 0.003710293, 0.006184011, 0.004860, 0.008917011,
  61, 0.007359, 0.012308, 0.011409, 0.018951, 0.003973871, 0.006646269, 0.005291, 0.009685997,
  62, 0.007867, 0.013163, 0.012240, 0.020484, 0.004248167, 0.007107782, 0.005728, 0.010457966,
  63, 0.008393, 0.013999, 0.013055, 0.022056, 0.004532271, 0.007559210, 0.006157, 0.011199611,
  64, 0.008958, 0.014846, 0.013873, 0.023682, 0.004837227, 0.008016683, 0.006594, 0.011928449,
  65, 0.009568, 0.015747, 0.014752, 0.025507, 0.005357879, 0.009133199, 0.007068, 0.012687281,
  66, 0.010302, 0.016844, 0.015721, 0.027435, 0.005769204, 0.009769591, 0.007613, 0.013541680,
  67, 0.011165, 0.018017, 0.016723, 0.029266, 0.006252196, 0.010450070, 0.008247, 0.014522870,
  68, 0.012216, 0.019337, 0.017771, 0.030803, 0.006841052, 0.011215634, 0.008991, 0.015671542,
  69, 0.013432, 0.020871, 0.018982, 0.032082, 0.007521696, 0.012105328, 0.009845, 0.016979223,
  70, 0.014833, 0.022427, 0.020418, 0.033425, 0.008306585, 0.013007673, 0.010806, 0.018436566,
  71, 0.016219, 0.024018, 0.021918, 0.035003, 0.009082412, 0.013930655, 0.011873, 0.020005114,
  72, 0.018152, 0.026452, 0.024027, 0.037203, 0.010164849, 0.015342326, 0.013063, 0.021660171,
  73, 0.019904, 0.028594, 0.025480, 0.040182, 0.011146140, 0.016584260, 0.014395, 0.023394244,
  74, 0.022051, 0.031490, 0.027795, 0.042821, 0.012348742, 0.018264089, 0.015907, 0.025269354,
  75, 0.024299, 0.034416, 0.030169, 0.046725, 0.015308252, 0.022026193, 0.017619, 0.027323818,
  76, 0.027123, 0.038177, 0.032943, 0.050461, 0.017087400, 0.024433126, 0.019637, 0.029788317,
  77, 0.030175, 0.042325, 0.035305, 0.054404, 0.019010467, 0.027088192, 0.021995, 0.032897148,
  78, 0.033527, 0.046248, 0.038829, 0.058845, 0.021122215, 0.029598677, 0.024664, 0.036370687,
  79, 0.037195, 0.051102, 0.042354, 0.063064, 0.023432550, 0.032704971, 0.027593, 0.040438816,
  80, 0.041367, 0.056288, 0.046095, 0.068031, 0.026060998, 0.036024177, 0.030848, 0.044811063,
  81, 0.046474, 0.062478, 0.051296, 0.073965, 0.029278743, 0.039985776, 0.034871, 0.050022066,
  82, 0.052049, 0.069313, 0.057979, 0.080399, 0.032790961, 0.044360514, 0.039259, 0.055627927,
  83, 0.058430, 0.076720, 0.062654, 0.088236, 0.036810905, 0.049100823, 0.044315, 0.061725955,
  84, 0.066267, 0.086706, 0.068704, 0.094411, 0.041748280, 0.055491819, 0.050572, 0.069985680,
  85, 0.074403, 0.096170, 0.074764, 0.102986, 0.046874061, 0.065395575, 0.057118, 0.077853374,
  86, 0.081415, 0.105892, 0.083635, 0.112188, 0.051291725, 0.072006597, 0.062800, 0.085975155,
  87, 0.092232, 0.119068, 0.092281, 0.122034, 0.058106402, 0.080965969, 0.071634, 0.097044937,
  88, 0.104231, 0.133492, 0.101665, 0.132537, 0.065665230, 0.090774738, 0.081526, 0.109246172,
  89, 0.117470, 0.149188, 0.111818, 0.143706, 0.074006276, 0.101448061, 0.092549, 0.122619100,
  90, 0.131997, 0.166155, 0.122764, 0.155542, 0.083158178, 0.112985265, 0.104769, 0.137186244,
  91, 0.147837, 0.184363, 0.134522, 0.168040, 0.093137177, 0.125367088, 0.118235, 0.152948216,
  92, 0.164991, 0.203755, 0.147101, 0.181185, 0.103944146, 0.138553424, 0.132981, 0.169879869,
  93, 0.183431, 0.224238, 0.160502, 0.194957, 0.115561746, 0.152481886, 0.149014, 0.187927082,
  94, 0.203099, 0.245687, 0.174713, 0.209325, 0.127952124, 0.167067398, 0.166315, 0.207004637,
  95, 0.223898, 0.267946, 0.189711, 0.224248, 0.141055473, 0.182203121, 0.184831, 0.226995677,
  96, 0.245698, 0.290828, 0.205459, 0.239677, 0.154789490, 0.197762790, 0.204475, 0.247752935,
  97, 0.268334, 0.314124, 0.221906, 0.255555, 0.169050327, 0.213604373, 0.225123, 0.269101888,
  98, 0.291611, 0.337610, 0.238987, 0.271815, 0.183714884, 0.229575027, 0.246617, 0.290846020,
  99, 0.315309, 0.361054, 0.256624, 0.288384, 0.198644587, 0.245516864, 0.268766, 0.312773585,
  100, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000000, 1.000000000, 1.000000, 1.000000000
),nrow = 51, byrow = TRUE)

colnames(mort_tbl) <- c("age", "female_white", "male_white", "female_black", "male_black", "female_asian",  "male_asian", "female_hispanic", "male_hispanic")
mort_tbl <-as.data.frame(mort_tbl)
  


# Hazard ratio mortality due to AD
## Andersen K, Lolk A, Martinussen T, Kragh-Sorensen P. Very mild to severe dementia and mortality: A 14-year follow-up - The Odense study. Dement Geriatr Cogn Disord. 2010;29:61-7.
## Create matrix for coefficients
ADmort_HR <- matrix(c(1, 2.92, 3.85, 9.52,
                      1, 2.43, 2.94, 6.60,
                      1, 3.52, 5.05, 13.74),
                      ncol = 4, nrow = 3, byrow = TRUE)
colnames(ADmort_HR) <- c("MCI", "mild","moderate","severe")
rownames(ADmort_HR) <- c("HR", "95%CI-low", "95%CI-high")

ADmort_coeff <- log(ADmort_HR[1,])


## Create variance-covariance matrix
ADmort_var <- ((log(ADmort_HR[3,])-log(ADmort_HR[2,]))/3.92)^2
ADmort_vcov <- matrix(c(ADmort_var[1], 0, 0, 0, 
                        0, ADmort_var[2], 0, 0,
                        0, 0, ADmort_var[3], 0,  
                        0, 0, 0, ADmort_var[4]
                       ),
                       ncol = 4, nrow = 4, byrow = TRUE)



# Utility
state_names_all <- c(
  "MCI","ADdem111-comm", "ADdem121-comm", "ADdem131-comm", "ADdem112-comm", "ADdem122-comm", "ADdem132-comm", "ADdem113-comm", "ADdem123-comm", "ADdem133-comm", "ADdem211-comm",
  "ADdem221-comm", "ADdem231-comm", "ADdem212-comm", "ADdem222-comm", "ADdem232-comm", "ADdem213-comm", "ADdem223-comm", "ADdem233-comm", "ADdem311-comm", "ADdem321-comm",  
  "ADdem331-comm","ADdem312-comm", "ADdem322-comm", "ADdem332-comm", "ADdem313-comm", "ADdem323-comm", "ADdem333-comm",
  "ADdem111-inst", "ADdem121-inst", "ADdem131-inst", "ADdem112-inst", "ADdem122-inst", "ADdem132-inst", "ADdem113-inst", "ADdem123-inst", "ADdem133-inst", "ADdem211-inst",
  "ADdem221-inst", "ADdem231-inst", "ADdem212-inst", "ADdem222-inst", "ADdem232-inst", "ADdem213-inst", "ADdem223-inst", "ADdem233-inst", "ADdem311-inst", "ADdem321-inst",  
  "ADdem331-inst","ADdem312-inst", "ADdem322-inst", "ADdem332-inst", "ADdem313-inst", "ADdem323-inst", "ADdem333-inst"
  )
utility <- matrix(c(0.73, rep(0.69, 9), rep(0.53, 9), rep(0.38, 9), rep(0.69, 9), rep(0.53, 9), rep(0.38, 9),
                    0.03, rep(0.03, 9), rep(0.03, 9), rep(0.03, 9), rep(0.03, 9), rep(0.03, 9), rep(0.03, 9)   
                    ), ncol = 55, nrow = 2, byrow = TRUE)
colnames(utility) <- state_names_all
rownames(utility) <- c("mean", "se")

u_mean <- utility[1,]

u_se <- utility[2,]




# Community to institution transition
## Data from Green et al 2019 (Supplement 2)
## Exponential time to event distribution
## Create matrix for coefficients
insti_coeff <- matrix(c(-4.26, 0.23, 0.35, 0.25, 0.60, 0.75, 1.43),
                      ncol = 7, nrow = 1, byrow = TRUE)
colnames(insti_coeff) <- c("mild","MMSEmod","MMSEsev","NPImod","NPIsev","FAQmod","FAQsev")


## Create variance-covariance matrix
insti_pvalue <- matrix(c(0.011, 0.008, 0.008, 0.0004, 0.0004, 0.0004),
                       ncol = 6, nrow = 1, byrow = TRUE)
insti_var <- (insti_coeff[,2:7]/qnorm(1-insti_pvalue/2))^2
insti_vcov <- matrix(c(0, 0, 0, 0, 0, 0, 0,
                       0, insti_var[1], 0, 0, 0, 0, 0, 
                       0, 0, insti_var[2], 0, 0, 0, 0,
                       0, 0, 0, insti_var[3], 0, 0, 0,  
                       0, 0, 0, 0, insti_var[4], 0, 0,
                       0, 0, 0, 0, 0, insti_var[5], 0,  
                       0, 0, 0, 0, 0, 0, insti_var[6]
                       ),
                     ncol = 7, nrow = 7, byrow = TRUE)
colnames(insti_vcov) <- rownames(insti_vcov) <- c("mild","MMSEmod","MMSEsev","NPImod","NPIsev","FAQmod","FAQsev")



# Cost of care
cost_care <- matrix(c(13364, rep(26727, 9), rep(31644, 9), rep(40645, 9), rep(111902, 9), rep(111902, 9), rep(113523, 9),
                      1336, rep(2673, 9), rep(3164, 9), rep(4065, 9), rep(11190, 9), rep(11190, 9), rep(11352, 9)  
), ncol = 55, nrow = 2, byrow = TRUE)
colnames(cost_care) <- state_names_all
rownames(cost_care) <- c("mean", "se")

c_care_mean <- cost_care[1,]
#print(c_care_mean)

c_care_se <- cost_care[2,]
#print(c_care_se)



# Treatment cost
c_treat <- matrix(c(tx_cost, rep(0, 54)), ncol = 55, nrow = 1, byrow = TRUE)
colnames(c_treat) <- state_names_all
rownames(c_treat) <- NULL



# Relative treatment effect MCI to AD transition
tx_RR <-   matrix(tx_RR_input,ncol = 4, nrow = 3, byrow = TRUE)
colnames(tx_RR) <- NULL
rownames(tx_RR) <- NULL

tx_coeff <- log(tx_RR[1,])
tx_var <- ((log(tx_RR[3,])-log(tx_RR[2,]))/3.92)^2

tx_vcov <- matrix(c(tx_var[1], 0, 0, 0, 
                        0, tx_var[2], 0, 0,
                        0, 0, tx_var[3], 0,  
                        0, 0, 0, tx_var[4]),
ncol = 4, nrow = 4, byrow = TRUE)




# Combining all parameters
params <- list(
  MCItoADdem_coeff = MCItoADdem_coeff,
  MCItoADdem_vcov = MCItoADdem_vcov,
  alpha_landingstates_ADdem = landingstates_ADdem,
  alpha_ADdem = transitions_ADdem,
  mort = mort_tbl,
  ADmort_coeff = ADmort_coeff,
  ADmort_vcov = ADmort_vcov,
  u_mean = u_mean,
  u_se = u_se,
  insti_coeff = insti_coeff,
  insti_vcov = insti_vcov,
  c_care_mean = c_care_mean,
  c_care_se = c_care_se,
  c_treat = c_treat,
  tx_coeff = tx_coeff,
  tx_vcov = tx_vcov,
  tx_duration = tx_duration
)



# Random number generation
rng_def <- define_rng({
list(
  MCItoADdem_coeff = multi_normal_rng(mu = MCItoADdem_coeff, Sigma = MCItoADdem_vcov),
  p_landingstates_ADdem = dirichlet_rng(alpha_landingstates_ADdem),
  p_ADdem = dirichlet_rng(alpha_ADdem),
  
  mrt_fw = fixed(mort$female_white, names = mort$age),
  mrt_mw = fixed(mort$male_white, names = mort$age),
  mrt_fb = fixed(mort$female_black, names = mort$age),
  mrt_mb = fixed(mort$male_black, names = mort$age),
  mrt_fa = fixed(mort$female_asian, names = mort$age), 
  mrt_ma = fixed(mort$male_asian, names = mort$age),
  mrt_fh = fixed(mort$female_hispanic, names = mort$age), 
  mrt_mh = fixed(mort$male_hispanic, names = mort$age),
  
  ADmort_coeff = multi_normal_rng(mu = ADmort_coeff, Sigma = ADmort_vcov),
  u = beta_rng(mean = u_mean, sd = u_se),
  insti_coeff = multi_normal_rng(mu = insti_coeff, Sigma = insti_vcov),
  c_care = gamma_rng(mean = c_care_mean, sd = c_care_se),
  c_treat = gamma_rng(mean = c_treat, sd = c_treat*0),   #### SOLVE THIS ISSUE; NOT ELEGANT
  log_tx_RR = multi_normal_rng(mu = tx_coeff, Sigma = tx_vcov),
  tx_duration = gamma_rng(mean = tx_duration, sd = tx_duration*0)   #### SOLVE THIS ISSUE; NOT ELEGANT
)  
}, n=PSA)

#evalrng <- eval_rng(x = rng_def, params)



# Transformed parameters
input_data <- expand(hesim_dat, by = c("strategies", "patients"))
#head(input_data)

tparams_def <- define_tparams({
## age and race specific transition from MCI to AD dementia
age6064 <- ifelse(age<65, 1, 0)        
age6569 <- ifelse(age>=65 & age<70, 1, 0) 
age7074 <- ifelse(age>=70 & age<75, 1, 0) 
age7579 <- ifelse(age>=75 & age<80, 1, 0) 
age8084 <- ifelse(age>=80 & age<85, 1, 0) 
age8589 <- ifelse(age>=85 & age<90, 1, 0) 
white <- ifelse(race == "White", 1, 0) 
black <- ifelse(race == "African American", 1, 0)
asian <- ifelse(race == "Asian", 1, 0)
hispanic <- ifelse(race == "Hispanic", 1, 0)
scale <- exp(MCItoADdem_coeff$v2 +    
             MCItoADdem_coeff$v3 * age6569 + 
             MCItoADdem_coeff$v4 * age7074 + 
             MCItoADdem_coeff$v5 * age7579 + 
             MCItoADdem_coeff$v6 * age8084 + 
             MCItoADdem_coeff$v7 * age8589 +
             MCItoADdem_coeff$v8 * black + 
             MCItoADdem_coeff$v9 * asian + 
             MCItoADdem_coeff$v10 * hispanic
             )
shape <- exp(MCItoADdem_coeff$v1)  
p_MCItoADdem <- 1 - exp(scale * (((time-1)^shape) - time^shape))

## mortality rate
age_new <- age + time -1

mrt <- 
  mrt_fw[["50"]] * (sex == "Female" & race == "White" & age_new == 50) + 
  mrt_fw[["51"]] * (sex == "Female" & race == "White" & age_new == 51) +
  mrt_fw[["52"]] * (sex == "Female" & race == "White" & age_new == 52) +
  mrt_fw[["53"]] * (sex == "Female" & race == "White" & age_new == 53) +
  mrt_fw[["54"]] * (sex == "Female" & race == "White" & age_new == 54) +
  mrt_fw[["55"]] * (sex == "Female" & race == "White" & age_new == 55) +
  mrt_fw[["56"]] * (sex == "Female" & race == "White" & age_new == 56) +
  mrt_fw[["57"]] * (sex == "Female" & race == "White" & age_new == 57) +
  mrt_fw[["58"]] * (sex == "Female" & race == "White" & age_new == 58) +
  mrt_fw[["59"]] * (sex == "Female" & race == "White" & age_new == 59) +
  mrt_fw[["60"]] * (sex == "Female" & race == "White" & age_new == 60) +
  mrt_fw[["61"]] * (sex == "Female" & race == "White" & age_new == 61) +
  mrt_fw[["62"]] * (sex == "Female" & race == "White" & age_new == 62) +
  mrt_fw[["63"]] * (sex == "Female" & race == "White" & age_new == 63) +
  mrt_fw[["64"]] * (sex == "Female" & race == "White" & age_new == 64) +
  mrt_fw[["65"]] * (sex == "Female" & race == "White" & age_new == 65) +
  mrt_fw[["66"]] * (sex == "Female" & race == "White" & age_new == 66) +
  mrt_fw[["67"]] * (sex == "Female" & race == "White" & age_new == 67) +
  mrt_fw[["68"]] * (sex == "Female" & race == "White" & age_new == 68) +
  mrt_fw[["69"]] * (sex == "Female" & race == "White" & age_new == 69) +
  mrt_fw[["70"]] * (sex == "Female" & race == "White" & age_new == 70) +
  mrt_fw[["71"]] * (sex == "Female" & race == "White" & age_new == 71) +
  mrt_fw[["72"]] * (sex == "Female" & race == "White" & age_new == 72) +
  mrt_fw[["73"]] * (sex == "Female" & race == "White" & age_new == 73) +
  mrt_fw[["74"]] * (sex == "Female" & race == "White" & age_new == 74) +
  mrt_fw[["75"]] * (sex == "Female" & race == "White" & age_new == 75) +
  mrt_fw[["76"]] * (sex == "Female" & race == "White" & age_new == 76) +
  mrt_fw[["77"]] * (sex == "Female" & race == "White" & age_new == 77) +
  mrt_fw[["78"]] * (sex == "Female" & race == "White" & age_new == 78) +
  mrt_fw[["79"]] * (sex == "Female" & race == "White" & age_new == 79) +
  mrt_fw[["80"]] * (sex == "Female" & race == "White" & age_new == 80) +
  mrt_fw[["81"]] * (sex == "Female" & race == "White" & age_new == 81) +
  mrt_fw[["82"]] * (sex == "Female" & race == "White" & age_new == 82) +
  mrt_fw[["83"]] * (sex == "Female" & race == "White" & age_new == 83) +
  mrt_fw[["84"]] * (sex == "Female" & race == "White" & age_new == 84) +
  mrt_fw[["85"]] * (sex == "Female" & race == "White" & age_new == 85) +
  mrt_fw[["86"]] * (sex == "Female" & race == "White" & age_new == 86) +
  mrt_fw[["87"]] * (sex == "Female" & race == "White" & age_new == 87) +
  mrt_fw[["88"]] * (sex == "Female" & race == "White" & age_new == 88) +
  mrt_fw[["89"]] * (sex == "Female" & race == "White" & age_new == 89) +
  mrt_fw[["90"]] * (sex == "Female" & race == "White" & age_new == 90) +
  mrt_fw[["91"]] * (sex == "Female" & race == "White" & age_new == 91) +
  mrt_fw[["92"]] * (sex == "Female" & race == "White" & age_new == 92) +
  mrt_fw[["93"]] * (sex == "Female" & race == "White" & age_new == 93) +
  mrt_fw[["94"]] * (sex == "Female" & race == "White" & age_new == 94) +
  mrt_fw[["95"]] * (sex == "Female" & race == "White" & age_new == 95) +
  mrt_fw[["96"]] * (sex == "Female" & race == "White" & age_new == 96) +
  mrt_fw[["97"]] * (sex == "Female" & race == "White" & age_new == 97) +
  mrt_fw[["98"]] * (sex == "Female" & race == "White" & age_new == 98) +
  mrt_fw[["99"]] * (sex == "Female" & race == "White" & age_new == 99) +
  mrt_fw[["100"]] * (sex == "Female" & race == "White" & age_new >= 100) +
  mrt_mw[["50"]] * (sex == "Male" & race == "White" & age_new == 50) +
  mrt_mw[["51"]] * (sex == "Male" & race == "White" & age_new == 51) +
  mrt_mw[["52"]] * (sex == "Male" & race == "White" & age_new == 52) +
  mrt_mw[["53"]] * (sex == "Male" & race == "White" & age_new == 53) +
  mrt_mw[["54"]] * (sex == "Male" & race == "White" & age_new == 54) +
  mrt_mw[["55"]] * (sex == "Male" & race == "White" & age_new == 55) +
  mrt_mw[["56"]] * (sex == "Male" & race == "White" & age_new == 56) +
  mrt_mw[["57"]] * (sex == "Male" & race == "White" & age_new == 57) +
  mrt_mw[["58"]] * (sex == "Male" & race == "White" & age_new == 58) +
  mrt_mw[["59"]] * (sex == "Male" & race == "White" & age_new == 59) +
  mrt_mw[["60"]] * (sex == "Male" & race == "White" & age_new == 60) +
  mrt_mw[["61"]] * (sex == "Male" & race == "White" & age_new == 61) +
  mrt_mw[["62"]] * (sex == "Male" & race == "White" & age_new == 62) +
  mrt_mw[["63"]] * (sex == "Male" & race == "White" & age_new == 63) +
  mrt_mw[["64"]] * (sex == "Male" & race == "White" & age_new == 64) +
  mrt_mw[["65"]] * (sex == "Male" & race == "White" & age_new == 65) +
  mrt_mw[["66"]] * (sex == "Male" & race == "White" & age_new == 66) +
  mrt_mw[["67"]] * (sex == "Male" & race == "White" & age_new == 67) +
  mrt_mw[["68"]] * (sex == "Male" & race == "White" & age_new == 68) +
  mrt_mw[["69"]] * (sex == "Male" & race == "White" & age_new == 69) +
  mrt_mw[["70"]] * (sex == "Male" & race == "White" & age_new == 70) +
  mrt_mw[["71"]] * (sex == "Male" & race == "White" & age_new == 71) +
  mrt_mw[["72"]] * (sex == "Male" & race == "White" & age_new == 72) +
  mrt_mw[["73"]] * (sex == "Male" & race == "White" & age_new == 73) +
  mrt_mw[["74"]] * (sex == "Male" & race == "White" & age_new == 74) +
  mrt_mw[["75"]] * (sex == "Male" & race == "White" & age_new == 75) +
  mrt_mw[["76"]] * (sex == "Male" & race == "White" & age_new == 76) +
  mrt_mw[["77"]] * (sex == "Male" & race == "White" & age_new == 77) +
  mrt_mw[["78"]] * (sex == "Male" & race == "White" & age_new == 78) +
  mrt_mw[["79"]] * (sex == "Male" & race == "White" & age_new == 79) +
  mrt_mw[["80"]] * (sex == "Male" & race == "White" & age_new == 80) +
  mrt_mw[["81"]] * (sex == "Male" & race == "White" & age_new == 81) +
  mrt_mw[["82"]] * (sex == "Male" & race == "White" & age_new == 82) +
  mrt_mw[["83"]] * (sex == "Male" & race == "White" & age_new == 83) +
  mrt_mw[["84"]] * (sex == "Male" & race == "White" & age_new == 84) +
  mrt_mw[["85"]] * (sex == "Male" & race == "White" & age_new == 85) +
  mrt_mw[["86"]] * (sex == "Male" & race == "White" & age_new == 86) +
  mrt_mw[["87"]] * (sex == "Male" & race == "White" & age_new == 87) +
  mrt_mw[["88"]] * (sex == "Male" & race == "White" & age_new == 88) +
  mrt_mw[["89"]] * (sex == "Male" & race == "White" & age_new == 89) +
  mrt_mw[["90"]] * (sex == "Male" & race == "White" & age_new == 90) +
  mrt_mw[["91"]] * (sex == "Male" & race == "White" & age_new == 91) +
  mrt_mw[["92"]] * (sex == "Male" & race == "White" & age_new == 92) +
  mrt_mw[["93"]] * (sex == "Male" & race == "White" & age_new == 93) +
  mrt_mw[["94"]] * (sex == "Male" & race == "White" & age_new == 94) +
  mrt_mw[["95"]] * (sex == "Male" & race == "White" & age_new == 95) +
  mrt_mw[["96"]] * (sex == "Male" & race == "White" & age_new == 96) +
  mrt_mw[["97"]] * (sex == "Male" & race == "White" & age_new == 97) +
  mrt_mw[["98"]] * (sex == "Male" & race == "White" & age_new == 98) +
  mrt_mw[["99"]] * (sex == "Male" & race == "White" & age_new == 99) +
  mrt_mw[["100"]] * (sex == "Male" & race == "White" & age_new >= 100) +
  mrt_fb[["50"]] * (sex == "Female" & race == "African American" & age_new == 50) +
  mrt_fb[["51"]] * (sex == "Female" & race == "African American" & age_new == 51) +
  mrt_fb[["52"]] * (sex == "Female" & race == "African American" & age_new == 52) +
  mrt_fb[["53"]] * (sex == "Female" & race == "African American" & age_new == 53) +
  mrt_fb[["54"]] * (sex == "Female" & race == "African American" & age_new == 54) +
  mrt_fb[["55"]] * (sex == "Female" & race == "African American" & age_new == 55) +
  mrt_fb[["56"]] * (sex == "Female" & race == "African American" & age_new == 56) +
  mrt_fb[["57"]] * (sex == "Female" & race == "African American" & age_new == 57) +
  mrt_fb[["58"]] * (sex == "Female" & race == "African American" & age_new == 58) +
  mrt_fb[["59"]] * (sex == "Female" & race == "African American" & age_new == 59) +
  mrt_fb[["60"]] * (sex == "Female" & race == "African American" & age_new == 60) +
  mrt_fb[["61"]] * (sex == "Female" & race == "African American" & age_new == 61) +
  mrt_fb[["62"]] * (sex == "Female" & race == "African American" & age_new == 62) +
  mrt_fb[["63"]] * (sex == "Female" & race == "African American" & age_new == 63) +
  mrt_fb[["64"]] * (sex == "Female" & race == "African American" & age_new == 64) +
  mrt_fb[["65"]] * (sex == "Female" & race == "African American" & age_new == 65) +
  mrt_fb[["66"]] * (sex == "Female" & race == "African American" & age_new == 66) +
  mrt_fb[["67"]] * (sex == "Female" & race == "African American" & age_new == 67) +
  mrt_fb[["68"]] * (sex == "Female" & race == "African American" & age_new == 68) +
  mrt_fb[["69"]] * (sex == "Female" & race == "African American" & age_new == 69) +
  mrt_fb[["70"]] * (sex == "Female" & race == "African American" & age_new == 70) +
  mrt_fb[["71"]] * (sex == "Female" & race == "African American" & age_new == 71) +
  mrt_fb[["72"]] * (sex == "Female" & race == "African American" & age_new == 72) +
  mrt_fb[["73"]] * (sex == "Female" & race == "African American" & age_new == 73) +
  mrt_fb[["74"]] * (sex == "Female" & race == "African American" & age_new == 74) +
  mrt_fb[["75"]] * (sex == "Female" & race == "African American" & age_new == 75) +
  mrt_fb[["76"]] * (sex == "Female" & race == "African American" & age_new == 76) +
  mrt_fb[["77"]] * (sex == "Female" & race == "African American" & age_new == 77) +
  mrt_fb[["78"]] * (sex == "Female" & race == "African American" & age_new == 78) +
  mrt_fb[["79"]] * (sex == "Female" & race == "African American" & age_new == 79) +
  mrt_fb[["80"]] * (sex == "Female" & race == "African American" & age_new == 80) +
  mrt_fb[["81"]] * (sex == "Female" & race == "African American" & age_new == 81) +
  mrt_fb[["82"]] * (sex == "Female" & race == "African American" & age_new == 82) +
  mrt_fb[["83"]] * (sex == "Female" & race == "African American" & age_new == 83) +
  mrt_fb[["84"]] * (sex == "Female" & race == "African American" & age_new == 84) +
  mrt_fb[["85"]] * (sex == "Female" & race == "African American" & age_new == 85) +
  mrt_fb[["86"]] * (sex == "Female" & race == "African American" & age_new == 86) +
  mrt_fb[["87"]] * (sex == "Female" & race == "African American" & age_new == 87) +
  mrt_fb[["88"]] * (sex == "Female" & race == "African American" & age_new == 88) +
  mrt_fb[["89"]] * (sex == "Female" & race == "African American" & age_new == 89) +
  mrt_fb[["90"]] * (sex == "Female" & race == "African American" & age_new == 90) +
  mrt_fb[["91"]] * (sex == "Female" & race == "African American" & age_new == 91) +
  mrt_fb[["92"]] * (sex == "Female" & race == "African American" & age_new == 92) +
  mrt_fb[["93"]] * (sex == "Female" & race == "African American" & age_new == 93) +
  mrt_fb[["94"]] * (sex == "Female" & race == "African American" & age_new == 94) +
  mrt_fb[["95"]] * (sex == "Female" & race == "African American" & age_new == 95) +
  mrt_fb[["96"]] * (sex == "Female" & race == "African American" & age_new == 96) +
  mrt_fb[["97"]] * (sex == "Female" & race == "African American" & age_new == 97) +
  mrt_fb[["98"]] * (sex == "Female" & race == "African American" & age_new == 98) +
  mrt_fb[["99"]] * (sex == "Female" & race == "African American" & age_new == 99) +
  mrt_fb[["100"]] * (sex == "Female" & race == "African American" & age_new >= 100) +
  mrt_mb[["50"]] * (sex == "Male" & race == "African American" & age_new == 50) +
  mrt_mb[["51"]] * (sex == "Male" & race == "African American" & age_new == 51) +
  mrt_mb[["52"]] * (sex == "Male" & race == "African American" & age_new == 52) +
  mrt_mb[["53"]] * (sex == "Male" & race == "African American" & age_new == 53) +
  mrt_mb[["54"]] * (sex == "Male" & race == "African American" & age_new == 54) +
  mrt_mb[["55"]] * (sex == "Male" & race == "African American" & age_new == 55) +
  mrt_mb[["56"]] * (sex == "Male" & race == "African American" & age_new == 56) +
  mrt_mb[["57"]] * (sex == "Male" & race == "African American" & age_new == 57) +
  mrt_mb[["58"]] * (sex == "Male" & race == "African American" & age_new == 58) +
  mrt_mb[["59"]] * (sex == "Male" & race == "African American" & age_new == 59) +
  mrt_mb[["60"]] * (sex == "Male" & race == "African American" & age_new == 60) +
  mrt_mb[["61"]] * (sex == "Male" & race == "African American" & age_new == 61) +
  mrt_mb[["62"]] * (sex == "Male" & race == "African American" & age_new == 62) +
  mrt_mb[["63"]] * (sex == "Male" & race == "African American" & age_new == 63) +
  mrt_mb[["64"]] * (sex == "Male" & race == "African American" & age_new == 64) +
  mrt_mb[["65"]] * (sex == "Male" & race == "African American" & age_new == 65) +
  mrt_mb[["66"]] * (sex == "Male" & race == "African American" & age_new == 66) +
  mrt_mb[["67"]] * (sex == "Male" & race == "African American" & age_new == 67) +
  mrt_mb[["68"]] * (sex == "Male" & race == "African American" & age_new == 68) +
  mrt_mb[["69"]] * (sex == "Male" & race == "African American" & age_new == 69) +
  mrt_mb[["70"]] * (sex == "Male" & race == "African American" & age_new == 70) +
  mrt_mb[["71"]] * (sex == "Male" & race == "African American" & age_new == 71) +
  mrt_mb[["72"]] * (sex == "Male" & race == "African American" & age_new == 72) +
  mrt_mb[["73"]] * (sex == "Male" & race == "African American" & age_new == 73) +
  mrt_mb[["74"]] * (sex == "Male" & race == "African American" & age_new == 74) +
  mrt_mb[["75"]] * (sex == "Male" & race == "African American" & age_new == 75) +
  mrt_mb[["76"]] * (sex == "Male" & race == "African American" & age_new == 76) +
  mrt_mb[["77"]] * (sex == "Male" & race == "African American" & age_new == 77) +
  mrt_mb[["78"]] * (sex == "Male" & race == "African American" & age_new == 78) +
  mrt_mb[["79"]] * (sex == "Male" & race == "African American" & age_new == 79) +
  mrt_mb[["80"]] * (sex == "Male" & race == "African American" & age_new == 80) +
  mrt_mb[["81"]] * (sex == "Male" & race == "African American" & age_new == 81) +
  mrt_mb[["82"]] * (sex == "Male" & race == "African American" & age_new == 82) +
  mrt_mb[["83"]] * (sex == "Male" & race == "African American" & age_new == 83) +
  mrt_mb[["84"]] * (sex == "Male" & race == "African American" & age_new == 84) +
  mrt_mb[["85"]] * (sex == "Male" & race == "African American" & age_new == 85) +
  mrt_mb[["86"]] * (sex == "Male" & race == "African American" & age_new == 86) +
  mrt_mb[["87"]] * (sex == "Male" & race == "African American" & age_new == 87) +
  mrt_mb[["88"]] * (sex == "Male" & race == "African American" & age_new == 88) +
  mrt_mb[["89"]] * (sex == "Male" & race == "African American" & age_new == 89) +
  mrt_mb[["90"]] * (sex == "Male" & race == "African American" & age_new == 90) +
  mrt_mb[["91"]] * (sex == "Male" & race == "African American" & age_new == 91) +
  mrt_mb[["92"]] * (sex == "Male" & race == "African American" & age_new == 92) +
  mrt_mb[["93"]] * (sex == "Male" & race == "African American" & age_new == 93) +
  mrt_mb[["94"]] * (sex == "Male" & race == "African American" & age_new == 94) +
  mrt_mb[["95"]] * (sex == "Male" & race == "African American" & age_new == 95) +
  mrt_mb[["96"]] * (sex == "Male" & race == "African American" & age_new == 96) +
  mrt_mb[["97"]] * (sex == "Male" & race == "African American" & age_new == 97) +
  mrt_mb[["98"]] * (sex == "Male" & race == "African American" & age_new == 98) +
  mrt_mb[["99"]] * (sex == "Male" & race == "African American" & age_new == 99) +
  mrt_mb[["100"]] * (sex == "Male" & race == "African American" & age_new >= 100) +
  mrt_fa[["50"]] * (sex == "Female" & race == "Asian" & age_new == 50) +
  mrt_fa[["51"]] * (sex == "Female" & race == "Asian" & age_new == 51) +
  mrt_fa[["52"]] * (sex == "Female" & race == "Asian" & age_new == 52) +
  mrt_fa[["53"]] * (sex == "Female" & race == "Asian" & age_new == 53) +
  mrt_fa[["54"]] * (sex == "Female" & race == "Asian" & age_new == 54) +
  mrt_fa[["55"]] * (sex == "Female" & race == "Asian" & age_new == 55) +
  mrt_fa[["56"]] * (sex == "Female" & race == "Asian" & age_new == 56) +
  mrt_fa[["57"]] * (sex == "Female" & race == "Asian" & age_new == 57) +
  mrt_fa[["58"]] * (sex == "Female" & race == "Asian" & age_new == 58) +
  mrt_fa[["59"]] * (sex == "Female" & race == "Asian" & age_new == 59) +
  mrt_fa[["60"]] * (sex == "Female" & race == "Asian" & age_new == 60) +
  mrt_fa[["61"]] * (sex == "Female" & race == "Asian" & age_new == 61) +
  mrt_fa[["62"]] * (sex == "Female" & race == "Asian" & age_new == 62) +
  mrt_fa[["63"]] * (sex == "Female" & race == "Asian" & age_new == 63) +
  mrt_fa[["64"]] * (sex == "Female" & race == "Asian" & age_new == 64) +
  mrt_fa[["65"]] * (sex == "Female" & race == "Asian" & age_new == 65) +
  mrt_fa[["66"]] * (sex == "Female" & race == "Asian" & age_new == 66) +
  mrt_fa[["67"]] * (sex == "Female" & race == "Asian" & age_new == 67) +
  mrt_fa[["68"]] * (sex == "Female" & race == "Asian" & age_new == 68) +
  mrt_fa[["69"]] * (sex == "Female" & race == "Asian" & age_new == 69) +
  mrt_fa[["70"]] * (sex == "Female" & race == "Asian" & age_new == 70) +
  mrt_fa[["71"]] * (sex == "Female" & race == "Asian" & age_new == 71) +
  mrt_fa[["72"]] * (sex == "Female" & race == "Asian" & age_new == 72) +
  mrt_fa[["73"]] * (sex == "Female" & race == "Asian" & age_new == 73) +
  mrt_fa[["74"]] * (sex == "Female" & race == "Asian" & age_new == 74) +
  mrt_fa[["75"]] * (sex == "Female" & race == "Asian" & age_new == 75) +
  mrt_fa[["76"]] * (sex == "Female" & race == "Asian" & age_new == 76) +
  mrt_fa[["77"]] * (sex == "Female" & race == "Asian" & age_new == 77) +
  mrt_fa[["78"]] * (sex == "Female" & race == "Asian" & age_new == 78) +
  mrt_fa[["79"]] * (sex == "Female" & race == "Asian" & age_new == 79) +
  mrt_fa[["80"]] * (sex == "Female" & race == "Asian" & age_new == 80) +
  mrt_fa[["81"]] * (sex == "Female" & race == "Asian" & age_new == 81) +
  mrt_fa[["82"]] * (sex == "Female" & race == "Asian" & age_new == 82) +
  mrt_fa[["83"]] * (sex == "Female" & race == "Asian" & age_new == 83) +
  mrt_fa[["84"]] * (sex == "Female" & race == "Asian" & age_new == 84) +
  mrt_fa[["85"]] * (sex == "Female" & race == "Asian" & age_new == 85) +
  mrt_fa[["86"]] * (sex == "Female" & race == "Asian" & age_new == 86) +
  mrt_fa[["87"]] * (sex == "Female" & race == "Asian" & age_new == 87) +
  mrt_fa[["88"]] * (sex == "Female" & race == "Asian" & age_new == 88) +
  mrt_fa[["89"]] * (sex == "Female" & race == "Asian" & age_new == 89) +
  mrt_fa[["90"]] * (sex == "Female" & race == "Asian" & age_new == 90) +
  mrt_fa[["91"]] * (sex == "Female" & race == "Asian" & age_new == 91) +
  mrt_fa[["92"]] * (sex == "Female" & race == "Asian" & age_new == 92) +
  mrt_fa[["93"]] * (sex == "Female" & race == "Asian" & age_new == 93) +
  mrt_fa[["94"]] * (sex == "Female" & race == "Asian" & age_new == 94) +
  mrt_fa[["95"]] * (sex == "Female" & race == "Asian" & age_new == 95) +
  mrt_fa[["96"]] * (sex == "Female" & race == "Asian" & age_new == 96) +
  mrt_fa[["97"]] * (sex == "Female" & race == "Asian" & age_new == 97) +
  mrt_fa[["98"]] * (sex == "Female" & race == "Asian" & age_new == 98) +
  mrt_fa[["99"]] * (sex == "Female" & race == "Asian" & age_new == 99) +
  mrt_fa[["100"]] * (sex == "Female" & race == "Asian" & age_new >= 100) +
  mrt_ma[["50"]] * (sex == "Male" & race == "Asian" & age_new == 50) +
  mrt_ma[["51"]] * (sex == "Male" & race == "Asian" & age_new == 51) +
  mrt_ma[["52"]] * (sex == "Male" & race == "Asian" & age_new == 52) +
  mrt_ma[["53"]] * (sex == "Male" & race == "Asian" & age_new == 53) +
  mrt_ma[["54"]] * (sex == "Male" & race == "Asian" & age_new == 54) +
  mrt_ma[["55"]] * (sex == "Male" & race == "Asian" & age_new == 55) +
  mrt_ma[["56"]] * (sex == "Male" & race == "Asian" & age_new == 56) +
  mrt_ma[["57"]] * (sex == "Male" & race == "Asian" & age_new == 57) +
  mrt_ma[["58"]] * (sex == "Male" & race == "Asian" & age_new == 58) +
  mrt_ma[["59"]] * (sex == "Male" & race == "Asian" & age_new == 59) +
  mrt_ma[["60"]] * (sex == "Male" & race == "Asian" & age_new == 60) +
  mrt_ma[["61"]] * (sex == "Male" & race == "Asian" & age_new == 61) +
  mrt_ma[["62"]] * (sex == "Male" & race == "Asian" & age_new == 62) +
  mrt_ma[["63"]] * (sex == "Male" & race == "Asian" & age_new == 63) +
  mrt_ma[["64"]] * (sex == "Male" & race == "Asian" & age_new == 64) +
  mrt_ma[["65"]] * (sex == "Male" & race == "Asian" & age_new == 65) +
  mrt_ma[["66"]] * (sex == "Male" & race == "Asian" & age_new == 66) +
  mrt_ma[["67"]] * (sex == "Male" & race == "Asian" & age_new == 67) +
  mrt_ma[["68"]] * (sex == "Male" & race == "Asian" & age_new == 68) +
  mrt_ma[["69"]] * (sex == "Male" & race == "Asian" & age_new == 69) +
  mrt_ma[["70"]] * (sex == "Male" & race == "Asian" & age_new == 70) +
  mrt_ma[["71"]] * (sex == "Male" & race == "Asian" & age_new == 71) +
  mrt_ma[["72"]] * (sex == "Male" & race == "Asian" & age_new == 72) +
  mrt_ma[["73"]] * (sex == "Male" & race == "Asian" & age_new == 73) +
  mrt_ma[["74"]] * (sex == "Male" & race == "Asian" & age_new == 74) +
  mrt_ma[["75"]] * (sex == "Male" & race == "Asian" & age_new == 75) +
  mrt_ma[["76"]] * (sex == "Male" & race == "Asian" & age_new == 76) +
  mrt_ma[["77"]] * (sex == "Male" & race == "Asian" & age_new == 77) +
  mrt_ma[["78"]] * (sex == "Male" & race == "Asian" & age_new == 78) +
  mrt_ma[["79"]] * (sex == "Male" & race == "Asian" & age_new == 79) +
  mrt_ma[["80"]] * (sex == "Male" & race == "Asian" & age_new == 80) +
  mrt_ma[["81"]] * (sex == "Male" & race == "Asian" & age_new == 81) +
  mrt_ma[["82"]] * (sex == "Male" & race == "Asian" & age_new == 82) +
  mrt_ma[["83"]] * (sex == "Male" & race == "Asian" & age_new == 83) +
  mrt_ma[["84"]] * (sex == "Male" & race == "Asian" & age_new == 84) +
  mrt_ma[["85"]] * (sex == "Male" & race == "Asian" & age_new == 85) +
  mrt_ma[["86"]] * (sex == "Male" & race == "Asian" & age_new == 86) +
  mrt_ma[["87"]] * (sex == "Male" & race == "Asian" & age_new == 87) +
  mrt_ma[["88"]] * (sex == "Male" & race == "Asian" & age_new == 88) +
  mrt_ma[["89"]] * (sex == "Male" & race == "Asian" & age_new == 89) +
  mrt_ma[["90"]] * (sex == "Male" & race == "Asian" & age_new == 90) +
  mrt_ma[["91"]] * (sex == "Male" & race == "Asian" & age_new == 91) +
  mrt_ma[["92"]] * (sex == "Male" & race == "Asian" & age_new == 92) +
  mrt_ma[["93"]] * (sex == "Male" & race == "Asian" & age_new == 93) +
  mrt_ma[["94"]] * (sex == "Male" & race == "Asian" & age_new == 94) +
  mrt_ma[["95"]] * (sex == "Male" & race == "Asian" & age_new == 95) +
  mrt_ma[["96"]] * (sex == "Male" & race == "Asian" & age_new == 96) +
  mrt_ma[["97"]] * (sex == "Male" & race == "Asian" & age_new == 97) +
  mrt_ma[["98"]] * (sex == "Male" & race == "Asian" & age_new == 98) +
  mrt_ma[["99"]] * (sex == "Male" & race == "Asian" & age_new == 99) +
  mrt_ma[["100"]] * (sex == "Male" & race == "Asian" & age_new >= 100) +
  mrt_fh[["50"]] * (sex == "Female" & race == "Hispanic" & age_new == 50) +
  mrt_fh[["51"]] * (sex == "Female" & race == "Hispanic" & age_new == 51) +
  mrt_fh[["52"]] * (sex == "Female" & race == "Hispanic" & age_new == 52) +
  mrt_fh[["53"]] * (sex == "Female" & race == "Hispanic" & age_new == 53) +
  mrt_fh[["54"]] * (sex == "Female" & race == "Hispanic" & age_new == 54) +
  mrt_fh[["55"]] * (sex == "Female" & race == "Hispanic" & age_new == 55) +
  mrt_fh[["56"]] * (sex == "Female" & race == "Hispanic" & age_new == 56) +
  mrt_fh[["57"]] * (sex == "Female" & race == "Hispanic" & age_new == 57) +
  mrt_fh[["58"]] * (sex == "Female" & race == "Hispanic" & age_new == 58) +
  mrt_fh[["59"]] * (sex == "Female" & race == "Hispanic" & age_new == 59) +
  mrt_fh[["60"]] * (sex == "Female" & race == "Hispanic" & age_new == 60) +
  mrt_fh[["61"]] * (sex == "Female" & race == "Hispanic" & age_new == 61) +
  mrt_fh[["62"]] * (sex == "Female" & race == "Hispanic" & age_new == 62) +
  mrt_fh[["63"]] * (sex == "Female" & race == "Hispanic" & age_new == 63) +
  mrt_fh[["64"]] * (sex == "Female" & race == "Hispanic" & age_new == 64) +
  mrt_fh[["65"]] * (sex == "Female" & race == "Hispanic" & age_new == 65) +
  mrt_fh[["66"]] * (sex == "Female" & race == "Hispanic" & age_new == 66) +
  mrt_fh[["67"]] * (sex == "Female" & race == "Hispanic" & age_new == 67) +
  mrt_fh[["68"]] * (sex == "Female" & race == "Hispanic" & age_new == 68) +
  mrt_fh[["69"]] * (sex == "Female" & race == "Hispanic" & age_new == 69) +
  mrt_fh[["70"]] * (sex == "Female" & race == "Hispanic" & age_new == 70) +
  mrt_fh[["71"]] * (sex == "Female" & race == "Hispanic" & age_new == 71) +
  mrt_fh[["72"]] * (sex == "Female" & race == "Hispanic" & age_new == 72) +
  mrt_fh[["73"]] * (sex == "Female" & race == "Hispanic" & age_new == 73) +
  mrt_fh[["74"]] * (sex == "Female" & race == "Hispanic" & age_new == 74) +
  mrt_fh[["75"]] * (sex == "Female" & race == "Hispanic" & age_new == 75) +
  mrt_fh[["76"]] * (sex == "Female" & race == "Hispanic" & age_new == 76) +
  mrt_fh[["77"]] * (sex == "Female" & race == "Hispanic" & age_new == 77) +
  mrt_fh[["78"]] * (sex == "Female" & race == "Hispanic" & age_new == 78) +
  mrt_fh[["79"]] * (sex == "Female" & race == "Hispanic" & age_new == 79) +
  mrt_fh[["80"]] * (sex == "Female" & race == "Hispanic" & age_new == 80) +
  mrt_fh[["81"]] * (sex == "Female" & race == "Hispanic" & age_new == 81) +
  mrt_fh[["82"]] * (sex == "Female" & race == "Hispanic" & age_new == 82) +
  mrt_fh[["83"]] * (sex == "Female" & race == "Hispanic" & age_new == 83) +
  mrt_fh[["84"]] * (sex == "Female" & race == "Hispanic" & age_new == 84) +
  mrt_fh[["85"]] * (sex == "Female" & race == "Hispanic" & age_new == 85) +
  mrt_fh[["86"]] * (sex == "Female" & race == "Hispanic" & age_new == 86) +
  mrt_fh[["87"]] * (sex == "Female" & race == "Hispanic" & age_new == 87) +
  mrt_fh[["88"]] * (sex == "Female" & race == "Hispanic" & age_new == 88) +
  mrt_fh[["89"]] * (sex == "Female" & race == "Hispanic" & age_new == 89) +
  mrt_fh[["90"]] * (sex == "Female" & race == "Hispanic" & age_new == 90) +
  mrt_fh[["91"]] * (sex == "Female" & race == "Hispanic" & age_new == 91) +
  mrt_fh[["92"]] * (sex == "Female" & race == "Hispanic" & age_new == 92) +
  mrt_fh[["93"]] * (sex == "Female" & race == "Hispanic" & age_new == 93) +
  mrt_fh[["94"]] * (sex == "Female" & race == "Hispanic" & age_new == 94) +
  mrt_fh[["95"]] * (sex == "Female" & race == "Hispanic" & age_new == 95) +
  mrt_fh[["96"]] * (sex == "Female" & race == "Hispanic" & age_new == 96) +
  mrt_fh[["97"]] * (sex == "Female" & race == "Hispanic" & age_new == 97) +
  mrt_fh[["98"]] * (sex == "Female" & race == "Hispanic" & age_new == 98) +
  mrt_fh[["99"]] * (sex == "Female" & race == "Hispanic" & age_new == 99) +
  mrt_fh[["100"]] * (sex == "Female" & race == "Hispanic" & age_new >= 100) +
  mrt_mh[["50"]] * (sex == "Male" & race == "Hispanic" & age_new == 50) +
  mrt_mh[["51"]] * (sex == "Male" & race == "Hispanic" & age_new == 51) +
  mrt_mh[["52"]] * (sex == "Male" & race == "Hispanic" & age_new == 52) +
  mrt_mh[["53"]] * (sex == "Male" & race == "Hispanic" & age_new == 53) +
  mrt_mh[["54"]] * (sex == "Male" & race == "Hispanic" & age_new == 54) +
  mrt_mh[["55"]] * (sex == "Male" & race == "Hispanic" & age_new == 55) +
  mrt_mh[["56"]] * (sex == "Male" & race == "Hispanic" & age_new == 56) +
  mrt_mh[["57"]] * (sex == "Male" & race == "Hispanic" & age_new == 57) +
  mrt_mh[["58"]] * (sex == "Male" & race == "Hispanic" & age_new == 58) +
  mrt_mh[["59"]] * (sex == "Male" & race == "Hispanic" & age_new == 59) +
  mrt_mh[["60"]] * (sex == "Male" & race == "Hispanic" & age_new == 60) +
  mrt_mh[["61"]] * (sex == "Male" & race == "Hispanic" & age_new == 61) +
  mrt_mh[["62"]] * (sex == "Male" & race == "Hispanic" & age_new == 62) +
  mrt_mh[["63"]] * (sex == "Male" & race == "Hispanic" & age_new == 63) +
  mrt_mh[["64"]] * (sex == "Male" & race == "Hispanic" & age_new == 64) +
  mrt_mh[["65"]] * (sex == "Male" & race == "Hispanic" & age_new == 65) +
  mrt_mh[["66"]] * (sex == "Male" & race == "Hispanic" & age_new == 66) +
  mrt_mh[["67"]] * (sex == "Male" & race == "Hispanic" & age_new == 67) +
  mrt_mh[["68"]] * (sex == "Male" & race == "Hispanic" & age_new == 68) +
  mrt_mh[["69"]] * (sex == "Male" & race == "Hispanic" & age_new == 69) +
  mrt_mh[["70"]] * (sex == "Male" & race == "Hispanic" & age_new == 70) +
  mrt_mh[["71"]] * (sex == "Male" & race == "Hispanic" & age_new == 71) +
  mrt_mh[["72"]] * (sex == "Male" & race == "Hispanic" & age_new == 72) +
  mrt_mh[["73"]] * (sex == "Male" & race == "Hispanic" & age_new == 73) +
  mrt_mh[["74"]] * (sex == "Male" & race == "Hispanic" & age_new == 74) +
  mrt_mh[["75"]] * (sex == "Male" & race == "Hispanic" & age_new == 75) +
  mrt_mh[["76"]] * (sex == "Male" & race == "Hispanic" & age_new == 76) +
  mrt_mh[["77"]] * (sex == "Male" & race == "Hispanic" & age_new == 77) +
  mrt_mh[["78"]] * (sex == "Male" & race == "Hispanic" & age_new == 78) +
  mrt_mh[["79"]] * (sex == "Male" & race == "Hispanic" & age_new == 79) +
  mrt_mh[["80"]] * (sex == "Male" & race == "Hispanic" & age_new == 80) +
  mrt_mh[["81"]] * (sex == "Male" & race == "Hispanic" & age_new == 81) +
  mrt_mh[["82"]] * (sex == "Male" & race == "Hispanic" & age_new == 82) +
  mrt_mh[["83"]] * (sex == "Male" & race == "Hispanic" & age_new == 83) +
  mrt_mh[["84"]] * (sex == "Male" & race == "Hispanic" & age_new == 84) +
  mrt_mh[["85"]] * (sex == "Male" & race == "Hispanic" & age_new == 85) +
  mrt_mh[["86"]] * (sex == "Male" & race == "Hispanic" & age_new == 86) +
  mrt_mh[["87"]] * (sex == "Male" & race == "Hispanic" & age_new == 87) +
  mrt_mh[["88"]] * (sex == "Male" & race == "Hispanic" & age_new == 88) +
  mrt_mh[["89"]] * (sex == "Male" & race == "Hispanic" & age_new == 89) +
  mrt_mh[["90"]] * (sex == "Male" & race == "Hispanic" & age_new == 90) +
  mrt_mh[["91"]] * (sex == "Male" & race == "Hispanic" & age_new == 91) +
  mrt_mh[["92"]] * (sex == "Male" & race == "Hispanic" & age_new == 92) +
  mrt_mh[["93"]] * (sex == "Male" & race == "Hispanic" & age_new == 93) +
  mrt_mh[["94"]] * (sex == "Male" & race == "Hispanic" & age_new == 94) +
  mrt_mh[["95"]] * (sex == "Male" & race == "Hispanic" & age_new == 95) +
  mrt_mh[["96"]] * (sex == "Male" & race == "Hispanic" & age_new == 96) +
  mrt_mh[["97"]] * (sex == "Male" & race == "Hispanic" & age_new == 97) +
  mrt_mh[["98"]] * (sex == "Male" & race == "Hispanic" & age_new == 98) +
  mrt_mh[["99"]] * (sex == "Male" & race == "Hispanic" & age_new == 99) +
  mrt_mh[["100"]] * (sex == "Male" & race == "Hispanic" & age_new >= 100) 
  
  
  
## hr mortality AD
hr_MCI <- exp(ADmort_coeff[["MCI"]])
hr_mild <- exp(ADmort_coeff[["mild"]])
hr_moderate <- exp(ADmort_coeff[["moderate"]])
hr_severe <- exp(ADmort_coeff[["severe"]])

## adjusted mortality with AD
mrt_MCI <- 1-exp(-(-log(1-mrt) * hr_MCI))
mrt_mild <- 1-exp(-(-log(1-mrt) * hr_mild))
mrt_moderate <- 1-exp(-(-log(1-mrt) * hr_moderate))
mrt_severe <- 1-exp(-(-log(1-mrt) * hr_severe))

## transitions from community to institution
p_insti_111 <- 1-exp(-exp(insti_coeff$v1))
p_insti_121 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v4))
p_insti_131 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v5))
p_insti_112 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v6))
p_insti_122 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v4+insti_coeff$v6))
p_insti_132 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v5+insti_coeff$v6))
p_insti_113 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v7))
p_insti_123 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v4+insti_coeff$v7))
p_insti_133 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v5+insti_coeff$v7))
p_insti_211 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2))
p_insti_221 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v4))
p_insti_231 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v5))
p_insti_212 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v6))
p_insti_222 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v4+insti_coeff$v6))
p_insti_232 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v5+insti_coeff$v6))
p_insti_213 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v7))
p_insti_223 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v4+insti_coeff$v7))
p_insti_233 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v2+insti_coeff$v5+insti_coeff$v7))
p_insti_311 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3))
p_insti_321 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v4))
p_insti_331 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v5))
p_insti_312 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v6))
p_insti_322 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v4+insti_coeff$v6))
p_insti_332 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v5+insti_coeff$v6))
p_insti_313 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v7))
p_insti_323 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v4+insti_coeff$v7))
p_insti_333 <- 1-exp(-exp(insti_coeff$v1+insti_coeff$v3+insti_coeff$v5+insti_coeff$v7))

## treatment cost
c_tx <- c_treat * (strategy_name == "Tx" & time <= tx_duration)

## treatment effect
rr <- exp(
    log_tx_RR$v1 * white +
    log_tx_RR$v2 * black +
    log_tx_RR$v3 * asian +
    log_tx_RR$v4 * hispanic 
    )
rr <- ifelse((strategy_name == "NoTx" | time > tx_duration), 1, rr)


## define transition matrix    
list(
  tpmatrix = tpmatrix(  
    C,	p_MCItoADdem*p_landingstates_ADdem[,1]*rr,	p_MCItoADdem*p_landingstates_ADdem[,2]*rr,	p_MCItoADdem*p_landingstates_ADdem[,3]*rr,	p_MCItoADdem*p_landingstates_ADdem[,4]*rr,	p_MCItoADdem*p_landingstates_ADdem[,5]*rr,	p_MCItoADdem*p_landingstates_ADdem[,6]*rr,	p_MCItoADdem*p_landingstates_ADdem[,7]*rr,	p_MCItoADdem*p_landingstates_ADdem[,8]*rr,	p_MCItoADdem*p_landingstates_ADdem[,9]*rr,	p_MCItoADdem*p_landingstates_ADdem[,10]*rr,	p_MCItoADdem*p_landingstates_ADdem[,11]*rr,	p_MCItoADdem*p_landingstates_ADdem[,12]*rr,	p_MCItoADdem*p_landingstates_ADdem[,13]*rr,	p_MCItoADdem*p_landingstates_ADdem[,14]*rr,	p_MCItoADdem*p_landingstates_ADdem[,15]*rr,	p_MCItoADdem*p_landingstates_ADdem[,16]*rr,	p_MCItoADdem*p_landingstates_ADdem[,17]*rr,	p_MCItoADdem*p_landingstates_ADdem[,18]*rr,	p_MCItoADdem*p_landingstates_ADdem[,19]*rr,	p_MCItoADdem*p_landingstates_ADdem[,20]*rr,	p_MCItoADdem*p_landingstates_ADdem[,21]*rr,	p_MCItoADdem*p_landingstates_ADdem[,22]*rr,	p_MCItoADdem*p_landingstates_ADdem[,23]*rr,	p_MCItoADdem*p_landingstates_ADdem[,24]*rr,	p_MCItoADdem*p_landingstates_ADdem[,25]*rr,	p_MCItoADdem*p_landingstates_ADdem[,26]*rr,	p_MCItoADdem*p_landingstates_ADdem[,27]*rr,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	mrt_MCI,
    0,	C,	p_ADdem$ADdem111_ADdem121*(1-p_insti_111),	p_ADdem$ADdem111_ADdem131*(1-p_insti_111),	p_ADdem$ADdem111_ADdem112*(1-p_insti_111),	p_ADdem$ADdem111_ADdem122*(1-p_insti_111),	p_ADdem$ADdem111_ADdem132*(1-p_insti_111),	p_ADdem$ADdem111_ADdem113*(1-p_insti_111),	p_ADdem$ADdem111_ADdem123*(1-p_insti_111),	p_ADdem$ADdem111_ADdem133*(1-p_insti_111),	p_ADdem$ADdem111_ADdem211*(1-p_insti_111),	p_ADdem$ADdem111_ADdem221*(1-p_insti_111),	p_ADdem$ADdem111_ADdem231*(1-p_insti_111),	p_ADdem$ADdem111_ADdem212*(1-p_insti_111),	p_ADdem$ADdem111_ADdem222*(1-p_insti_111),	p_ADdem$ADdem111_ADdem232*(1-p_insti_111),	p_ADdem$ADdem111_ADdem213*(1-p_insti_111),	p_ADdem$ADdem111_ADdem223*(1-p_insti_111),	p_ADdem$ADdem111_ADdem233*(1-p_insti_111),	p_ADdem$ADdem111_ADdem311*(1-p_insti_111),	p_ADdem$ADdem111_ADdem321*(1-p_insti_111),	p_ADdem$ADdem111_ADdem331*(1-p_insti_111),	p_ADdem$ADdem111_ADdem312*(1-p_insti_111),	p_ADdem$ADdem111_ADdem322*(1-p_insti_111),	p_ADdem$ADdem111_ADdem332*(1-p_insti_111),	p_ADdem$ADdem111_ADdem313*(1-p_insti_111),	p_ADdem$ADdem111_ADdem323*(1-p_insti_111),	p_ADdem$ADdem111_ADdem333*(1-p_insti_111),	p_ADdem$ADdem111_ADdem111*p_insti_111,	p_ADdem$ADdem111_ADdem121*p_insti_111,	p_ADdem$ADdem111_ADdem131*p_insti_111,	p_ADdem$ADdem111_ADdem112*p_insti_111,	p_ADdem$ADdem111_ADdem122*p_insti_111,	p_ADdem$ADdem111_ADdem132*p_insti_111,	p_ADdem$ADdem111_ADdem113*p_insti_111,	p_ADdem$ADdem111_ADdem123*p_insti_111,	p_ADdem$ADdem111_ADdem133*p_insti_111,	p_ADdem$ADdem111_ADdem211*p_insti_111,	p_ADdem$ADdem111_ADdem221*p_insti_111,	p_ADdem$ADdem111_ADdem231*p_insti_111,	p_ADdem$ADdem111_ADdem212*p_insti_111,	p_ADdem$ADdem111_ADdem222*p_insti_111,	p_ADdem$ADdem111_ADdem232*p_insti_111,	p_ADdem$ADdem111_ADdem213*p_insti_111,	p_ADdem$ADdem111_ADdem223*p_insti_111,	p_ADdem$ADdem111_ADdem233*p_insti_111,	p_ADdem$ADdem111_ADdem311*p_insti_111,	p_ADdem$ADdem111_ADdem321*p_insti_111,	p_ADdem$ADdem111_ADdem331*p_insti_111,	p_ADdem$ADdem111_ADdem312*p_insti_111,	p_ADdem$ADdem111_ADdem322*p_insti_111,	p_ADdem$ADdem111_ADdem332*p_insti_111,	p_ADdem$ADdem111_ADdem313*p_insti_111,	p_ADdem$ADdem111_ADdem323*p_insti_111,	p_ADdem$ADdem111_ADdem333*p_insti_111,	mrt_mild,
    0,	p_ADdem$ADdem121_ADdem111*(1-p_insti_121),	C,	p_ADdem$ADdem121_ADdem131*(1-p_insti_121),	p_ADdem$ADdem121_ADdem112*(1-p_insti_121),	p_ADdem$ADdem121_ADdem122*(1-p_insti_121),	p_ADdem$ADdem121_ADdem132*(1-p_insti_121),	p_ADdem$ADdem121_ADdem113*(1-p_insti_121),	p_ADdem$ADdem121_ADdem123*(1-p_insti_121),	p_ADdem$ADdem121_ADdem133*(1-p_insti_121),	p_ADdem$ADdem121_ADdem211*(1-p_insti_121),	p_ADdem$ADdem121_ADdem221*(1-p_insti_121),	p_ADdem$ADdem121_ADdem231*(1-p_insti_121),	p_ADdem$ADdem121_ADdem212*(1-p_insti_121),	p_ADdem$ADdem121_ADdem222*(1-p_insti_121),	p_ADdem$ADdem121_ADdem232*(1-p_insti_121),	p_ADdem$ADdem121_ADdem213*(1-p_insti_121),	p_ADdem$ADdem121_ADdem223*(1-p_insti_121),	p_ADdem$ADdem121_ADdem233*(1-p_insti_121),	p_ADdem$ADdem121_ADdem311*(1-p_insti_121),	p_ADdem$ADdem121_ADdem321*(1-p_insti_121),	p_ADdem$ADdem121_ADdem331*(1-p_insti_121),	p_ADdem$ADdem121_ADdem312*(1-p_insti_121),	p_ADdem$ADdem121_ADdem322*(1-p_insti_121),	p_ADdem$ADdem121_ADdem332*(1-p_insti_121),	p_ADdem$ADdem121_ADdem313*(1-p_insti_121),	p_ADdem$ADdem121_ADdem323*(1-p_insti_121),	p_ADdem$ADdem121_ADdem333*(1-p_insti_121),	p_ADdem$ADdem121_ADdem111*p_insti_121,	p_ADdem$ADdem121_ADdem121*p_insti_121,	p_ADdem$ADdem121_ADdem131*p_insti_121,	p_ADdem$ADdem121_ADdem112*p_insti_121,	p_ADdem$ADdem121_ADdem122*p_insti_121,	p_ADdem$ADdem121_ADdem132*p_insti_121,	p_ADdem$ADdem121_ADdem113*p_insti_121,	p_ADdem$ADdem121_ADdem123*p_insti_121,	p_ADdem$ADdem121_ADdem133*p_insti_121,	p_ADdem$ADdem121_ADdem211*p_insti_121,	p_ADdem$ADdem121_ADdem221*p_insti_121,	p_ADdem$ADdem121_ADdem231*p_insti_121,	p_ADdem$ADdem121_ADdem212*p_insti_121,	p_ADdem$ADdem121_ADdem222*p_insti_121,	p_ADdem$ADdem121_ADdem232*p_insti_121,	p_ADdem$ADdem121_ADdem213*p_insti_121,	p_ADdem$ADdem121_ADdem223*p_insti_121,	p_ADdem$ADdem121_ADdem233*p_insti_121,	p_ADdem$ADdem121_ADdem311*p_insti_121,	p_ADdem$ADdem121_ADdem321*p_insti_121,	p_ADdem$ADdem121_ADdem331*p_insti_121,	p_ADdem$ADdem121_ADdem312*p_insti_121,	p_ADdem$ADdem121_ADdem322*p_insti_121,	p_ADdem$ADdem121_ADdem332*p_insti_121,	p_ADdem$ADdem121_ADdem313*p_insti_121,	p_ADdem$ADdem121_ADdem323*p_insti_121,	p_ADdem$ADdem121_ADdem333*p_insti_121,	mrt_mild,
    0,	p_ADdem$ADdem131_ADdem111*(1-p_insti_131),	p_ADdem$ADdem131_ADdem121*(1-p_insti_131),	C,	p_ADdem$ADdem131_ADdem112*(1-p_insti_131),	p_ADdem$ADdem131_ADdem122*(1-p_insti_131),	p_ADdem$ADdem131_ADdem132*(1-p_insti_131),	p_ADdem$ADdem131_ADdem113*(1-p_insti_131),	p_ADdem$ADdem131_ADdem123*(1-p_insti_131),	p_ADdem$ADdem131_ADdem133*(1-p_insti_131),	p_ADdem$ADdem131_ADdem211*(1-p_insti_131),	p_ADdem$ADdem131_ADdem221*(1-p_insti_131),	p_ADdem$ADdem131_ADdem231*(1-p_insti_131),	p_ADdem$ADdem131_ADdem212*(1-p_insti_131),	p_ADdem$ADdem131_ADdem222*(1-p_insti_131),	p_ADdem$ADdem131_ADdem232*(1-p_insti_131),	p_ADdem$ADdem131_ADdem213*(1-p_insti_131),	p_ADdem$ADdem131_ADdem223*(1-p_insti_131),	p_ADdem$ADdem131_ADdem233*(1-p_insti_131),	p_ADdem$ADdem131_ADdem311*(1-p_insti_131),	p_ADdem$ADdem131_ADdem321*(1-p_insti_131),	p_ADdem$ADdem131_ADdem331*(1-p_insti_131),	p_ADdem$ADdem131_ADdem312*(1-p_insti_131),	p_ADdem$ADdem131_ADdem322*(1-p_insti_131),	p_ADdem$ADdem131_ADdem332*(1-p_insti_131),	p_ADdem$ADdem131_ADdem313*(1-p_insti_131),	p_ADdem$ADdem131_ADdem323*(1-p_insti_131),	p_ADdem$ADdem131_ADdem333*(1-p_insti_131),	p_ADdem$ADdem131_ADdem111*p_insti_131,	p_ADdem$ADdem131_ADdem121*p_insti_131,	p_ADdem$ADdem131_ADdem131*p_insti_131,	p_ADdem$ADdem131_ADdem112*p_insti_131,	p_ADdem$ADdem131_ADdem122*p_insti_131,	p_ADdem$ADdem131_ADdem132*p_insti_131,	p_ADdem$ADdem131_ADdem113*p_insti_131,	p_ADdem$ADdem131_ADdem123*p_insti_131,	p_ADdem$ADdem131_ADdem133*p_insti_131,	p_ADdem$ADdem131_ADdem211*p_insti_131,	p_ADdem$ADdem131_ADdem221*p_insti_131,	p_ADdem$ADdem131_ADdem231*p_insti_131,	p_ADdem$ADdem131_ADdem212*p_insti_131,	p_ADdem$ADdem131_ADdem222*p_insti_131,	p_ADdem$ADdem131_ADdem232*p_insti_131,	p_ADdem$ADdem131_ADdem213*p_insti_131,	p_ADdem$ADdem131_ADdem223*p_insti_131,	p_ADdem$ADdem131_ADdem233*p_insti_131,	p_ADdem$ADdem131_ADdem311*p_insti_131,	p_ADdem$ADdem131_ADdem321*p_insti_131,	p_ADdem$ADdem131_ADdem331*p_insti_131,	p_ADdem$ADdem131_ADdem312*p_insti_131,	p_ADdem$ADdem131_ADdem322*p_insti_131,	p_ADdem$ADdem131_ADdem332*p_insti_131,	p_ADdem$ADdem131_ADdem313*p_insti_131,	p_ADdem$ADdem131_ADdem323*p_insti_131,	p_ADdem$ADdem131_ADdem333*p_insti_131,	mrt_moderate,
    0,	p_ADdem$ADdem112_ADdem111*(1-p_insti_112),	p_ADdem$ADdem112_ADdem121*(1-p_insti_112),	p_ADdem$ADdem112_ADdem131*(1-p_insti_112),	C,	p_ADdem$ADdem112_ADdem122*(1-p_insti_112),	p_ADdem$ADdem112_ADdem132*(1-p_insti_112),	p_ADdem$ADdem112_ADdem113*(1-p_insti_112),	p_ADdem$ADdem112_ADdem123*(1-p_insti_112),	p_ADdem$ADdem112_ADdem133*(1-p_insti_112),	p_ADdem$ADdem112_ADdem211*(1-p_insti_112),	p_ADdem$ADdem112_ADdem221*(1-p_insti_112),	p_ADdem$ADdem112_ADdem231*(1-p_insti_112),	p_ADdem$ADdem112_ADdem212*(1-p_insti_112),	p_ADdem$ADdem112_ADdem222*(1-p_insti_112),	p_ADdem$ADdem112_ADdem232*(1-p_insti_112),	p_ADdem$ADdem112_ADdem213*(1-p_insti_112),	p_ADdem$ADdem112_ADdem223*(1-p_insti_112),	p_ADdem$ADdem112_ADdem233*(1-p_insti_112),	p_ADdem$ADdem112_ADdem311*(1-p_insti_112),	p_ADdem$ADdem112_ADdem321*(1-p_insti_112),	p_ADdem$ADdem112_ADdem331*(1-p_insti_112),	p_ADdem$ADdem112_ADdem312*(1-p_insti_112),	p_ADdem$ADdem112_ADdem322*(1-p_insti_112),	p_ADdem$ADdem112_ADdem332*(1-p_insti_112),	p_ADdem$ADdem112_ADdem313*(1-p_insti_112),	p_ADdem$ADdem112_ADdem323*(1-p_insti_112),	p_ADdem$ADdem112_ADdem333*(1-p_insti_112),	p_ADdem$ADdem112_ADdem111*p_insti_112,	p_ADdem$ADdem112_ADdem121*p_insti_112,	p_ADdem$ADdem112_ADdem131*p_insti_112,	p_ADdem$ADdem112_ADdem112*p_insti_112,	p_ADdem$ADdem112_ADdem122*p_insti_112,	p_ADdem$ADdem112_ADdem132*p_insti_112,	p_ADdem$ADdem112_ADdem113*p_insti_112,	p_ADdem$ADdem112_ADdem123*p_insti_112,	p_ADdem$ADdem112_ADdem133*p_insti_112,	p_ADdem$ADdem112_ADdem211*p_insti_112,	p_ADdem$ADdem112_ADdem221*p_insti_112,	p_ADdem$ADdem112_ADdem231*p_insti_112,	p_ADdem$ADdem112_ADdem212*p_insti_112,	p_ADdem$ADdem112_ADdem222*p_insti_112,	p_ADdem$ADdem112_ADdem232*p_insti_112,	p_ADdem$ADdem112_ADdem213*p_insti_112,	p_ADdem$ADdem112_ADdem223*p_insti_112,	p_ADdem$ADdem112_ADdem233*p_insti_112,	p_ADdem$ADdem112_ADdem311*p_insti_112,	p_ADdem$ADdem112_ADdem321*p_insti_112,	p_ADdem$ADdem112_ADdem331*p_insti_112,	p_ADdem$ADdem112_ADdem312*p_insti_112,	p_ADdem$ADdem112_ADdem322*p_insti_112,	p_ADdem$ADdem112_ADdem332*p_insti_112,	p_ADdem$ADdem112_ADdem313*p_insti_112,	p_ADdem$ADdem112_ADdem323*p_insti_112,	p_ADdem$ADdem112_ADdem333*p_insti_112,	mrt_mild,
    0,	p_ADdem$ADdem122_ADdem111*(1-p_insti_122),	p_ADdem$ADdem122_ADdem121*(1-p_insti_122),	p_ADdem$ADdem122_ADdem131*(1-p_insti_122),	p_ADdem$ADdem122_ADdem112*(1-p_insti_122),	C,	p_ADdem$ADdem122_ADdem132*(1-p_insti_122),	p_ADdem$ADdem122_ADdem113*(1-p_insti_122),	p_ADdem$ADdem122_ADdem123*(1-p_insti_122),	p_ADdem$ADdem122_ADdem133*(1-p_insti_122),	p_ADdem$ADdem122_ADdem211*(1-p_insti_122),	p_ADdem$ADdem122_ADdem221*(1-p_insti_122),	p_ADdem$ADdem122_ADdem231*(1-p_insti_122),	p_ADdem$ADdem122_ADdem212*(1-p_insti_122),	p_ADdem$ADdem122_ADdem222*(1-p_insti_122),	p_ADdem$ADdem122_ADdem232*(1-p_insti_122),	p_ADdem$ADdem122_ADdem213*(1-p_insti_122),	p_ADdem$ADdem122_ADdem223*(1-p_insti_122),	p_ADdem$ADdem122_ADdem233*(1-p_insti_122),	p_ADdem$ADdem122_ADdem311*(1-p_insti_122),	p_ADdem$ADdem122_ADdem321*(1-p_insti_122),	p_ADdem$ADdem122_ADdem331*(1-p_insti_122),	p_ADdem$ADdem122_ADdem312*(1-p_insti_122),	p_ADdem$ADdem122_ADdem322*(1-p_insti_122),	p_ADdem$ADdem122_ADdem332*(1-p_insti_122),	p_ADdem$ADdem122_ADdem313*(1-p_insti_122),	p_ADdem$ADdem122_ADdem323*(1-p_insti_122),	p_ADdem$ADdem122_ADdem333*(1-p_insti_122),	p_ADdem$ADdem122_ADdem111*p_insti_122,	p_ADdem$ADdem122_ADdem121*p_insti_122,	p_ADdem$ADdem122_ADdem131*p_insti_122,	p_ADdem$ADdem122_ADdem112*p_insti_122,	p_ADdem$ADdem122_ADdem122*p_insti_122,	p_ADdem$ADdem122_ADdem132*p_insti_122,	p_ADdem$ADdem122_ADdem113*p_insti_122,	p_ADdem$ADdem122_ADdem123*p_insti_122,	p_ADdem$ADdem122_ADdem133*p_insti_122,	p_ADdem$ADdem122_ADdem211*p_insti_122,	p_ADdem$ADdem122_ADdem221*p_insti_122,	p_ADdem$ADdem122_ADdem231*p_insti_122,	p_ADdem$ADdem122_ADdem212*p_insti_122,	p_ADdem$ADdem122_ADdem222*p_insti_122,	p_ADdem$ADdem122_ADdem232*p_insti_122,	p_ADdem$ADdem122_ADdem213*p_insti_122,	p_ADdem$ADdem122_ADdem223*p_insti_122,	p_ADdem$ADdem122_ADdem233*p_insti_122,	p_ADdem$ADdem122_ADdem311*p_insti_122,	p_ADdem$ADdem122_ADdem321*p_insti_122,	p_ADdem$ADdem122_ADdem331*p_insti_122,	p_ADdem$ADdem122_ADdem312*p_insti_122,	p_ADdem$ADdem122_ADdem322*p_insti_122,	p_ADdem$ADdem122_ADdem332*p_insti_122,	p_ADdem$ADdem122_ADdem313*p_insti_122,	p_ADdem$ADdem122_ADdem323*p_insti_122,	p_ADdem$ADdem122_ADdem333*p_insti_122,	mrt_mild,
    0,	p_ADdem$ADdem132_ADdem111*(1-p_insti_132),	p_ADdem$ADdem132_ADdem121*(1-p_insti_132),	p_ADdem$ADdem132_ADdem131*(1-p_insti_132),	p_ADdem$ADdem132_ADdem112*(1-p_insti_132),	p_ADdem$ADdem132_ADdem122*(1-p_insti_132),	C,	p_ADdem$ADdem132_ADdem113*(1-p_insti_132),	p_ADdem$ADdem132_ADdem123*(1-p_insti_132),	p_ADdem$ADdem132_ADdem133*(1-p_insti_132),	p_ADdem$ADdem132_ADdem211*(1-p_insti_132),	p_ADdem$ADdem132_ADdem221*(1-p_insti_132),	p_ADdem$ADdem132_ADdem231*(1-p_insti_132),	p_ADdem$ADdem132_ADdem212*(1-p_insti_132),	p_ADdem$ADdem132_ADdem222*(1-p_insti_132),	p_ADdem$ADdem132_ADdem232*(1-p_insti_132),	p_ADdem$ADdem132_ADdem213*(1-p_insti_132),	p_ADdem$ADdem132_ADdem223*(1-p_insti_132),	p_ADdem$ADdem132_ADdem233*(1-p_insti_132),	p_ADdem$ADdem132_ADdem311*(1-p_insti_132),	p_ADdem$ADdem132_ADdem321*(1-p_insti_132),	p_ADdem$ADdem132_ADdem331*(1-p_insti_132),	p_ADdem$ADdem132_ADdem312*(1-p_insti_132),	p_ADdem$ADdem132_ADdem322*(1-p_insti_132),	p_ADdem$ADdem132_ADdem332*(1-p_insti_132),	p_ADdem$ADdem132_ADdem313*(1-p_insti_132),	p_ADdem$ADdem132_ADdem323*(1-p_insti_132),	p_ADdem$ADdem132_ADdem333*(1-p_insti_132),	p_ADdem$ADdem132_ADdem111*p_insti_132,	p_ADdem$ADdem132_ADdem121*p_insti_132,	p_ADdem$ADdem132_ADdem131*p_insti_132,	p_ADdem$ADdem132_ADdem112*p_insti_132,	p_ADdem$ADdem132_ADdem122*p_insti_132,	p_ADdem$ADdem132_ADdem132*p_insti_132,	p_ADdem$ADdem132_ADdem113*p_insti_132,	p_ADdem$ADdem132_ADdem123*p_insti_132,	p_ADdem$ADdem132_ADdem133*p_insti_132,	p_ADdem$ADdem132_ADdem211*p_insti_132,	p_ADdem$ADdem132_ADdem221*p_insti_132,	p_ADdem$ADdem132_ADdem231*p_insti_132,	p_ADdem$ADdem132_ADdem212*p_insti_132,	p_ADdem$ADdem132_ADdem222*p_insti_132,	p_ADdem$ADdem132_ADdem232*p_insti_132,	p_ADdem$ADdem132_ADdem213*p_insti_132,	p_ADdem$ADdem132_ADdem223*p_insti_132,	p_ADdem$ADdem132_ADdem233*p_insti_132,	p_ADdem$ADdem132_ADdem311*p_insti_132,	p_ADdem$ADdem132_ADdem321*p_insti_132,	p_ADdem$ADdem132_ADdem331*p_insti_132,	p_ADdem$ADdem132_ADdem312*p_insti_132,	p_ADdem$ADdem132_ADdem322*p_insti_132,	p_ADdem$ADdem132_ADdem332*p_insti_132,	p_ADdem$ADdem132_ADdem313*p_insti_132,	p_ADdem$ADdem132_ADdem323*p_insti_132,	p_ADdem$ADdem132_ADdem333*p_insti_132,	mrt_moderate,
    0,	p_ADdem$ADdem113_ADdem111*(1-p_insti_113),	p_ADdem$ADdem113_ADdem121*(1-p_insti_113),	p_ADdem$ADdem113_ADdem131*(1-p_insti_113),	p_ADdem$ADdem113_ADdem112*(1-p_insti_113),	p_ADdem$ADdem113_ADdem122*(1-p_insti_113),	p_ADdem$ADdem113_ADdem132*(1-p_insti_113),	C,	p_ADdem$ADdem113_ADdem123*(1-p_insti_113),	p_ADdem$ADdem113_ADdem133*(1-p_insti_113),	p_ADdem$ADdem113_ADdem211*(1-p_insti_113),	p_ADdem$ADdem113_ADdem221*(1-p_insti_113),	p_ADdem$ADdem113_ADdem231*(1-p_insti_113),	p_ADdem$ADdem113_ADdem212*(1-p_insti_113),	p_ADdem$ADdem113_ADdem222*(1-p_insti_113),	p_ADdem$ADdem113_ADdem232*(1-p_insti_113),	p_ADdem$ADdem113_ADdem213*(1-p_insti_113),	p_ADdem$ADdem113_ADdem223*(1-p_insti_113),	p_ADdem$ADdem113_ADdem233*(1-p_insti_113),	p_ADdem$ADdem113_ADdem311*(1-p_insti_113),	p_ADdem$ADdem113_ADdem321*(1-p_insti_113),	p_ADdem$ADdem113_ADdem331*(1-p_insti_113),	p_ADdem$ADdem113_ADdem312*(1-p_insti_113),	p_ADdem$ADdem113_ADdem322*(1-p_insti_113),	p_ADdem$ADdem113_ADdem332*(1-p_insti_113),	p_ADdem$ADdem113_ADdem313*(1-p_insti_113),	p_ADdem$ADdem113_ADdem323*(1-p_insti_113),	p_ADdem$ADdem113_ADdem333*(1-p_insti_113),	p_ADdem$ADdem113_ADdem111*p_insti_113,	p_ADdem$ADdem113_ADdem121*p_insti_113,	p_ADdem$ADdem113_ADdem131*p_insti_113,	p_ADdem$ADdem113_ADdem112*p_insti_113,	p_ADdem$ADdem113_ADdem122*p_insti_113,	p_ADdem$ADdem113_ADdem132*p_insti_113,	p_ADdem$ADdem113_ADdem113*p_insti_113,	p_ADdem$ADdem113_ADdem123*p_insti_113,	p_ADdem$ADdem113_ADdem133*p_insti_113,	p_ADdem$ADdem113_ADdem211*p_insti_113,	p_ADdem$ADdem113_ADdem221*p_insti_113,	p_ADdem$ADdem113_ADdem231*p_insti_113,	p_ADdem$ADdem113_ADdem212*p_insti_113,	p_ADdem$ADdem113_ADdem222*p_insti_113,	p_ADdem$ADdem113_ADdem232*p_insti_113,	p_ADdem$ADdem113_ADdem213*p_insti_113,	p_ADdem$ADdem113_ADdem223*p_insti_113,	p_ADdem$ADdem113_ADdem233*p_insti_113,	p_ADdem$ADdem113_ADdem311*p_insti_113,	p_ADdem$ADdem113_ADdem321*p_insti_113,	p_ADdem$ADdem113_ADdem331*p_insti_113,	p_ADdem$ADdem113_ADdem312*p_insti_113,	p_ADdem$ADdem113_ADdem322*p_insti_113,	p_ADdem$ADdem113_ADdem332*p_insti_113,	p_ADdem$ADdem113_ADdem313*p_insti_113,	p_ADdem$ADdem113_ADdem323*p_insti_113,	p_ADdem$ADdem113_ADdem333*p_insti_113,	mrt_moderate,
    0,	p_ADdem$ADdem123_ADdem111*(1-p_insti_123),	p_ADdem$ADdem123_ADdem121*(1-p_insti_123),	p_ADdem$ADdem123_ADdem131*(1-p_insti_123),	p_ADdem$ADdem123_ADdem112*(1-p_insti_123),	p_ADdem$ADdem123_ADdem122*(1-p_insti_123),	p_ADdem$ADdem123_ADdem132*(1-p_insti_123),	p_ADdem$ADdem123_ADdem113*(1-p_insti_123),	C,	p_ADdem$ADdem123_ADdem133*(1-p_insti_123),	p_ADdem$ADdem123_ADdem211*(1-p_insti_123),	p_ADdem$ADdem123_ADdem221*(1-p_insti_123),	p_ADdem$ADdem123_ADdem231*(1-p_insti_123),	p_ADdem$ADdem123_ADdem212*(1-p_insti_123),	p_ADdem$ADdem123_ADdem222*(1-p_insti_123),	p_ADdem$ADdem123_ADdem232*(1-p_insti_123),	p_ADdem$ADdem123_ADdem213*(1-p_insti_123),	p_ADdem$ADdem123_ADdem223*(1-p_insti_123),	p_ADdem$ADdem123_ADdem233*(1-p_insti_123),	p_ADdem$ADdem123_ADdem311*(1-p_insti_123),	p_ADdem$ADdem123_ADdem321*(1-p_insti_123),	p_ADdem$ADdem123_ADdem331*(1-p_insti_123),	p_ADdem$ADdem123_ADdem312*(1-p_insti_123),	p_ADdem$ADdem123_ADdem322*(1-p_insti_123),	p_ADdem$ADdem123_ADdem332*(1-p_insti_123),	p_ADdem$ADdem123_ADdem313*(1-p_insti_123),	p_ADdem$ADdem123_ADdem323*(1-p_insti_123),	p_ADdem$ADdem123_ADdem333*(1-p_insti_123),	p_ADdem$ADdem123_ADdem111*p_insti_123,	p_ADdem$ADdem123_ADdem121*p_insti_123,	p_ADdem$ADdem123_ADdem131*p_insti_123,	p_ADdem$ADdem123_ADdem112*p_insti_123,	p_ADdem$ADdem123_ADdem122*p_insti_123,	p_ADdem$ADdem123_ADdem132*p_insti_123,	p_ADdem$ADdem123_ADdem113*p_insti_123,	p_ADdem$ADdem123_ADdem123*p_insti_123,	p_ADdem$ADdem123_ADdem133*p_insti_123,	p_ADdem$ADdem123_ADdem211*p_insti_123,	p_ADdem$ADdem123_ADdem221*p_insti_123,	p_ADdem$ADdem123_ADdem231*p_insti_123,	p_ADdem$ADdem123_ADdem212*p_insti_123,	p_ADdem$ADdem123_ADdem222*p_insti_123,	p_ADdem$ADdem123_ADdem232*p_insti_123,	p_ADdem$ADdem123_ADdem213*p_insti_123,	p_ADdem$ADdem123_ADdem223*p_insti_123,	p_ADdem$ADdem123_ADdem233*p_insti_123,	p_ADdem$ADdem123_ADdem311*p_insti_123,	p_ADdem$ADdem123_ADdem321*p_insti_123,	p_ADdem$ADdem123_ADdem331*p_insti_123,	p_ADdem$ADdem123_ADdem312*p_insti_123,	p_ADdem$ADdem123_ADdem322*p_insti_123,	p_ADdem$ADdem123_ADdem332*p_insti_123,	p_ADdem$ADdem123_ADdem313*p_insti_123,	p_ADdem$ADdem123_ADdem323*p_insti_123,	p_ADdem$ADdem123_ADdem333*p_insti_123,	mrt_moderate,
    0,	p_ADdem$ADdem133_ADdem111*(1-p_insti_133),	p_ADdem$ADdem133_ADdem121*(1-p_insti_133),	p_ADdem$ADdem133_ADdem131*(1-p_insti_133),	p_ADdem$ADdem133_ADdem112*(1-p_insti_133),	p_ADdem$ADdem133_ADdem122*(1-p_insti_133),	p_ADdem$ADdem133_ADdem132*(1-p_insti_133),	p_ADdem$ADdem133_ADdem113*(1-p_insti_133),	p_ADdem$ADdem133_ADdem123*(1-p_insti_133),	C,	p_ADdem$ADdem133_ADdem211*(1-p_insti_133),	p_ADdem$ADdem133_ADdem221*(1-p_insti_133),	p_ADdem$ADdem133_ADdem231*(1-p_insti_133),	p_ADdem$ADdem133_ADdem212*(1-p_insti_133),	p_ADdem$ADdem133_ADdem222*(1-p_insti_133),	p_ADdem$ADdem133_ADdem232*(1-p_insti_133),	p_ADdem$ADdem133_ADdem213*(1-p_insti_133),	p_ADdem$ADdem133_ADdem223*(1-p_insti_133),	p_ADdem$ADdem133_ADdem233*(1-p_insti_133),	p_ADdem$ADdem133_ADdem311*(1-p_insti_133),	p_ADdem$ADdem133_ADdem321*(1-p_insti_133),	p_ADdem$ADdem133_ADdem331*(1-p_insti_133),	p_ADdem$ADdem133_ADdem312*(1-p_insti_133),	p_ADdem$ADdem133_ADdem322*(1-p_insti_133),	p_ADdem$ADdem133_ADdem332*(1-p_insti_133),	p_ADdem$ADdem133_ADdem313*(1-p_insti_133),	p_ADdem$ADdem133_ADdem323*(1-p_insti_133),	p_ADdem$ADdem133_ADdem333*(1-p_insti_133),	p_ADdem$ADdem133_ADdem111*p_insti_133,	p_ADdem$ADdem133_ADdem121*p_insti_133,	p_ADdem$ADdem133_ADdem131*p_insti_133,	p_ADdem$ADdem133_ADdem112*p_insti_133,	p_ADdem$ADdem133_ADdem122*p_insti_133,	p_ADdem$ADdem133_ADdem132*p_insti_133,	p_ADdem$ADdem133_ADdem113*p_insti_133,	p_ADdem$ADdem133_ADdem123*p_insti_133,	p_ADdem$ADdem133_ADdem133*p_insti_133,	p_ADdem$ADdem133_ADdem211*p_insti_133,	p_ADdem$ADdem133_ADdem221*p_insti_133,	p_ADdem$ADdem133_ADdem231*p_insti_133,	p_ADdem$ADdem133_ADdem212*p_insti_133,	p_ADdem$ADdem133_ADdem222*p_insti_133,	p_ADdem$ADdem133_ADdem232*p_insti_133,	p_ADdem$ADdem133_ADdem213*p_insti_133,	p_ADdem$ADdem133_ADdem223*p_insti_133,	p_ADdem$ADdem133_ADdem233*p_insti_133,	p_ADdem$ADdem133_ADdem311*p_insti_133,	p_ADdem$ADdem133_ADdem321*p_insti_133,	p_ADdem$ADdem133_ADdem331*p_insti_133,	p_ADdem$ADdem133_ADdem312*p_insti_133,	p_ADdem$ADdem133_ADdem322*p_insti_133,	p_ADdem$ADdem133_ADdem332*p_insti_133,	p_ADdem$ADdem133_ADdem313*p_insti_133,	p_ADdem$ADdem133_ADdem323*p_insti_133,	p_ADdem$ADdem133_ADdem333*p_insti_133,	mrt_moderate,
    0,	p_ADdem$ADdem211_ADdem111*(1-p_insti_211),	p_ADdem$ADdem211_ADdem121*(1-p_insti_211),	p_ADdem$ADdem211_ADdem131*(1-p_insti_211),	p_ADdem$ADdem211_ADdem112*(1-p_insti_211),	p_ADdem$ADdem211_ADdem122*(1-p_insti_211),	p_ADdem$ADdem211_ADdem132*(1-p_insti_211),	p_ADdem$ADdem211_ADdem113*(1-p_insti_211),	p_ADdem$ADdem211_ADdem123*(1-p_insti_211),	p_ADdem$ADdem211_ADdem133*(1-p_insti_211),	C,	p_ADdem$ADdem211_ADdem221*(1-p_insti_211),	p_ADdem$ADdem211_ADdem231*(1-p_insti_211),	p_ADdem$ADdem211_ADdem212*(1-p_insti_211),	p_ADdem$ADdem211_ADdem222*(1-p_insti_211),	p_ADdem$ADdem211_ADdem232*(1-p_insti_211),	p_ADdem$ADdem211_ADdem213*(1-p_insti_211),	p_ADdem$ADdem211_ADdem223*(1-p_insti_211),	p_ADdem$ADdem211_ADdem233*(1-p_insti_211),	p_ADdem$ADdem211_ADdem311*(1-p_insti_211),	p_ADdem$ADdem211_ADdem321*(1-p_insti_211),	p_ADdem$ADdem211_ADdem331*(1-p_insti_211),	p_ADdem$ADdem211_ADdem312*(1-p_insti_211),	p_ADdem$ADdem211_ADdem322*(1-p_insti_211),	p_ADdem$ADdem211_ADdem332*(1-p_insti_211),	p_ADdem$ADdem211_ADdem313*(1-p_insti_211),	p_ADdem$ADdem211_ADdem323*(1-p_insti_211),	p_ADdem$ADdem211_ADdem333*(1-p_insti_211),	p_ADdem$ADdem211_ADdem111*p_insti_211,	p_ADdem$ADdem211_ADdem121*p_insti_211,	p_ADdem$ADdem211_ADdem131*p_insti_211,	p_ADdem$ADdem211_ADdem112*p_insti_211,	p_ADdem$ADdem211_ADdem122*p_insti_211,	p_ADdem$ADdem211_ADdem132*p_insti_211,	p_ADdem$ADdem211_ADdem113*p_insti_211,	p_ADdem$ADdem211_ADdem123*p_insti_211,	p_ADdem$ADdem211_ADdem133*p_insti_211,	p_ADdem$ADdem211_ADdem211*p_insti_211,	p_ADdem$ADdem211_ADdem221*p_insti_211,	p_ADdem$ADdem211_ADdem231*p_insti_211,	p_ADdem$ADdem211_ADdem212*p_insti_211,	p_ADdem$ADdem211_ADdem222*p_insti_211,	p_ADdem$ADdem211_ADdem232*p_insti_211,	p_ADdem$ADdem211_ADdem213*p_insti_211,	p_ADdem$ADdem211_ADdem223*p_insti_211,	p_ADdem$ADdem211_ADdem233*p_insti_211,	p_ADdem$ADdem211_ADdem311*p_insti_211,	p_ADdem$ADdem211_ADdem321*p_insti_211,	p_ADdem$ADdem211_ADdem331*p_insti_211,	p_ADdem$ADdem211_ADdem312*p_insti_211,	p_ADdem$ADdem211_ADdem322*p_insti_211,	p_ADdem$ADdem211_ADdem332*p_insti_211,	p_ADdem$ADdem211_ADdem313*p_insti_211,	p_ADdem$ADdem211_ADdem323*p_insti_211,	p_ADdem$ADdem211_ADdem333*p_insti_211,	mrt_moderate,
    0,	p_ADdem$ADdem221_ADdem111*(1-p_insti_221),	p_ADdem$ADdem221_ADdem121*(1-p_insti_221),	p_ADdem$ADdem221_ADdem131*(1-p_insti_221),	p_ADdem$ADdem221_ADdem112*(1-p_insti_221),	p_ADdem$ADdem221_ADdem122*(1-p_insti_221),	p_ADdem$ADdem221_ADdem132*(1-p_insti_221),	p_ADdem$ADdem221_ADdem113*(1-p_insti_221),	p_ADdem$ADdem221_ADdem123*(1-p_insti_221),	p_ADdem$ADdem221_ADdem133*(1-p_insti_221),	p_ADdem$ADdem221_ADdem211*(1-p_insti_221),	C,	p_ADdem$ADdem221_ADdem231*(1-p_insti_221),	p_ADdem$ADdem221_ADdem212*(1-p_insti_221),	p_ADdem$ADdem221_ADdem222*(1-p_insti_221),	p_ADdem$ADdem221_ADdem232*(1-p_insti_221),	p_ADdem$ADdem221_ADdem213*(1-p_insti_221),	p_ADdem$ADdem221_ADdem223*(1-p_insti_221),	p_ADdem$ADdem221_ADdem233*(1-p_insti_221),	p_ADdem$ADdem221_ADdem311*(1-p_insti_221),	p_ADdem$ADdem221_ADdem321*(1-p_insti_221),	p_ADdem$ADdem221_ADdem331*(1-p_insti_221),	p_ADdem$ADdem221_ADdem312*(1-p_insti_221),	p_ADdem$ADdem221_ADdem322*(1-p_insti_221),	p_ADdem$ADdem221_ADdem332*(1-p_insti_221),	p_ADdem$ADdem221_ADdem313*(1-p_insti_221),	p_ADdem$ADdem221_ADdem323*(1-p_insti_221),	p_ADdem$ADdem221_ADdem333*(1-p_insti_221),	p_ADdem$ADdem221_ADdem111*p_insti_221,	p_ADdem$ADdem221_ADdem121*p_insti_221,	p_ADdem$ADdem221_ADdem131*p_insti_221,	p_ADdem$ADdem221_ADdem112*p_insti_221,	p_ADdem$ADdem221_ADdem122*p_insti_221,	p_ADdem$ADdem221_ADdem132*p_insti_221,	p_ADdem$ADdem221_ADdem113*p_insti_221,	p_ADdem$ADdem221_ADdem123*p_insti_221,	p_ADdem$ADdem221_ADdem133*p_insti_221,	p_ADdem$ADdem221_ADdem211*p_insti_221,	p_ADdem$ADdem221_ADdem221*p_insti_221,	p_ADdem$ADdem221_ADdem231*p_insti_221,	p_ADdem$ADdem221_ADdem212*p_insti_221,	p_ADdem$ADdem221_ADdem222*p_insti_221,	p_ADdem$ADdem221_ADdem232*p_insti_221,	p_ADdem$ADdem221_ADdem213*p_insti_221,	p_ADdem$ADdem221_ADdem223*p_insti_221,	p_ADdem$ADdem221_ADdem233*p_insti_221,	p_ADdem$ADdem221_ADdem311*p_insti_221,	p_ADdem$ADdem221_ADdem321*p_insti_221,	p_ADdem$ADdem221_ADdem331*p_insti_221,	p_ADdem$ADdem221_ADdem312*p_insti_221,	p_ADdem$ADdem221_ADdem322*p_insti_221,	p_ADdem$ADdem221_ADdem332*p_insti_221,	p_ADdem$ADdem221_ADdem313*p_insti_221,	p_ADdem$ADdem221_ADdem323*p_insti_221,	p_ADdem$ADdem221_ADdem333*p_insti_221,	mrt_moderate,
    0,	p_ADdem$ADdem231_ADdem111*(1-p_insti_231),	p_ADdem$ADdem231_ADdem121*(1-p_insti_231),	p_ADdem$ADdem231_ADdem131*(1-p_insti_231),	p_ADdem$ADdem231_ADdem112*(1-p_insti_231),	p_ADdem$ADdem231_ADdem122*(1-p_insti_231),	p_ADdem$ADdem231_ADdem132*(1-p_insti_231),	p_ADdem$ADdem231_ADdem113*(1-p_insti_231),	p_ADdem$ADdem231_ADdem123*(1-p_insti_231),	p_ADdem$ADdem231_ADdem133*(1-p_insti_231),	p_ADdem$ADdem231_ADdem211*(1-p_insti_231),	p_ADdem$ADdem231_ADdem221*(1-p_insti_231),	C,	p_ADdem$ADdem231_ADdem212*(1-p_insti_231),	p_ADdem$ADdem231_ADdem222*(1-p_insti_231),	p_ADdem$ADdem231_ADdem232*(1-p_insti_231),	p_ADdem$ADdem231_ADdem213*(1-p_insti_231),	p_ADdem$ADdem231_ADdem223*(1-p_insti_231),	p_ADdem$ADdem231_ADdem233*(1-p_insti_231),	p_ADdem$ADdem231_ADdem311*(1-p_insti_231),	p_ADdem$ADdem231_ADdem321*(1-p_insti_231),	p_ADdem$ADdem231_ADdem331*(1-p_insti_231),	p_ADdem$ADdem231_ADdem312*(1-p_insti_231),	p_ADdem$ADdem231_ADdem322*(1-p_insti_231),	p_ADdem$ADdem231_ADdem332*(1-p_insti_231),	p_ADdem$ADdem231_ADdem313*(1-p_insti_231),	p_ADdem$ADdem231_ADdem323*(1-p_insti_231),	p_ADdem$ADdem231_ADdem333*(1-p_insti_231),	p_ADdem$ADdem231_ADdem111*p_insti_231,	p_ADdem$ADdem231_ADdem121*p_insti_231,	p_ADdem$ADdem231_ADdem131*p_insti_231,	p_ADdem$ADdem231_ADdem112*p_insti_231,	p_ADdem$ADdem231_ADdem122*p_insti_231,	p_ADdem$ADdem231_ADdem132*p_insti_231,	p_ADdem$ADdem231_ADdem113*p_insti_231,	p_ADdem$ADdem231_ADdem123*p_insti_231,	p_ADdem$ADdem231_ADdem133*p_insti_231,	p_ADdem$ADdem231_ADdem211*p_insti_231,	p_ADdem$ADdem231_ADdem221*p_insti_231,	p_ADdem$ADdem231_ADdem231*p_insti_231,	p_ADdem$ADdem231_ADdem212*p_insti_231,	p_ADdem$ADdem231_ADdem222*p_insti_231,	p_ADdem$ADdem231_ADdem232*p_insti_231,	p_ADdem$ADdem231_ADdem213*p_insti_231,	p_ADdem$ADdem231_ADdem223*p_insti_231,	p_ADdem$ADdem231_ADdem233*p_insti_231,	p_ADdem$ADdem231_ADdem311*p_insti_231,	p_ADdem$ADdem231_ADdem321*p_insti_231,	p_ADdem$ADdem231_ADdem331*p_insti_231,	p_ADdem$ADdem231_ADdem312*p_insti_231,	p_ADdem$ADdem231_ADdem322*p_insti_231,	p_ADdem$ADdem231_ADdem332*p_insti_231,	p_ADdem$ADdem231_ADdem313*p_insti_231,	p_ADdem$ADdem231_ADdem323*p_insti_231,	p_ADdem$ADdem231_ADdem333*p_insti_231,	mrt_moderate,
    0,	p_ADdem$ADdem212_ADdem111*(1-p_insti_212),	p_ADdem$ADdem212_ADdem121*(1-p_insti_212),	p_ADdem$ADdem212_ADdem131*(1-p_insti_212),	p_ADdem$ADdem212_ADdem112*(1-p_insti_212),	p_ADdem$ADdem212_ADdem122*(1-p_insti_212),	p_ADdem$ADdem212_ADdem132*(1-p_insti_212),	p_ADdem$ADdem212_ADdem113*(1-p_insti_212),	p_ADdem$ADdem212_ADdem123*(1-p_insti_212),	p_ADdem$ADdem212_ADdem133*(1-p_insti_212),	p_ADdem$ADdem212_ADdem211*(1-p_insti_212),	p_ADdem$ADdem212_ADdem221*(1-p_insti_212),	p_ADdem$ADdem212_ADdem231*(1-p_insti_212),	C,	p_ADdem$ADdem212_ADdem222*(1-p_insti_212),	p_ADdem$ADdem212_ADdem232*(1-p_insti_212),	p_ADdem$ADdem212_ADdem213*(1-p_insti_212),	p_ADdem$ADdem212_ADdem223*(1-p_insti_212),	p_ADdem$ADdem212_ADdem233*(1-p_insti_212),	p_ADdem$ADdem212_ADdem311*(1-p_insti_212),	p_ADdem$ADdem212_ADdem321*(1-p_insti_212),	p_ADdem$ADdem212_ADdem331*(1-p_insti_212),	p_ADdem$ADdem212_ADdem312*(1-p_insti_212),	p_ADdem$ADdem212_ADdem322*(1-p_insti_212),	p_ADdem$ADdem212_ADdem332*(1-p_insti_212),	p_ADdem$ADdem212_ADdem313*(1-p_insti_212),	p_ADdem$ADdem212_ADdem323*(1-p_insti_212),	p_ADdem$ADdem212_ADdem333*(1-p_insti_212),	p_ADdem$ADdem212_ADdem111*p_insti_212,	p_ADdem$ADdem212_ADdem121*p_insti_212,	p_ADdem$ADdem212_ADdem131*p_insti_212,	p_ADdem$ADdem212_ADdem112*p_insti_212,	p_ADdem$ADdem212_ADdem122*p_insti_212,	p_ADdem$ADdem212_ADdem132*p_insti_212,	p_ADdem$ADdem212_ADdem113*p_insti_212,	p_ADdem$ADdem212_ADdem123*p_insti_212,	p_ADdem$ADdem212_ADdem133*p_insti_212,	p_ADdem$ADdem212_ADdem211*p_insti_212,	p_ADdem$ADdem212_ADdem221*p_insti_212,	p_ADdem$ADdem212_ADdem231*p_insti_212,	p_ADdem$ADdem212_ADdem212*p_insti_212,	p_ADdem$ADdem212_ADdem222*p_insti_212,	p_ADdem$ADdem212_ADdem232*p_insti_212,	p_ADdem$ADdem212_ADdem213*p_insti_212,	p_ADdem$ADdem212_ADdem223*p_insti_212,	p_ADdem$ADdem212_ADdem233*p_insti_212,	p_ADdem$ADdem212_ADdem311*p_insti_212,	p_ADdem$ADdem212_ADdem321*p_insti_212,	p_ADdem$ADdem212_ADdem331*p_insti_212,	p_ADdem$ADdem212_ADdem312*p_insti_212,	p_ADdem$ADdem212_ADdem322*p_insti_212,	p_ADdem$ADdem212_ADdem332*p_insti_212,	p_ADdem$ADdem212_ADdem313*p_insti_212,	p_ADdem$ADdem212_ADdem323*p_insti_212,	p_ADdem$ADdem212_ADdem333*p_insti_212,	mrt_moderate,
    0,	p_ADdem$ADdem222_ADdem111*(1-p_insti_222),	p_ADdem$ADdem222_ADdem121*(1-p_insti_222),	p_ADdem$ADdem222_ADdem131*(1-p_insti_222),	p_ADdem$ADdem222_ADdem112*(1-p_insti_222),	p_ADdem$ADdem222_ADdem122*(1-p_insti_222),	p_ADdem$ADdem222_ADdem132*(1-p_insti_222),	p_ADdem$ADdem222_ADdem113*(1-p_insti_222),	p_ADdem$ADdem222_ADdem123*(1-p_insti_222),	p_ADdem$ADdem222_ADdem133*(1-p_insti_222),	p_ADdem$ADdem222_ADdem211*(1-p_insti_222),	p_ADdem$ADdem222_ADdem221*(1-p_insti_222),	p_ADdem$ADdem222_ADdem231*(1-p_insti_222),	p_ADdem$ADdem222_ADdem212*(1-p_insti_222),	C,	p_ADdem$ADdem222_ADdem232*(1-p_insti_222),	p_ADdem$ADdem222_ADdem213*(1-p_insti_222),	p_ADdem$ADdem222_ADdem223*(1-p_insti_222),	p_ADdem$ADdem222_ADdem233*(1-p_insti_222),	p_ADdem$ADdem222_ADdem311*(1-p_insti_222),	p_ADdem$ADdem222_ADdem321*(1-p_insti_222),	p_ADdem$ADdem222_ADdem331*(1-p_insti_222),	p_ADdem$ADdem222_ADdem312*(1-p_insti_222),	p_ADdem$ADdem222_ADdem322*(1-p_insti_222),	p_ADdem$ADdem222_ADdem332*(1-p_insti_222),	p_ADdem$ADdem222_ADdem313*(1-p_insti_222),	p_ADdem$ADdem222_ADdem323*(1-p_insti_222),	p_ADdem$ADdem222_ADdem333*(1-p_insti_222),	p_ADdem$ADdem222_ADdem111*p_insti_222,	p_ADdem$ADdem222_ADdem121*p_insti_222,	p_ADdem$ADdem222_ADdem131*p_insti_222,	p_ADdem$ADdem222_ADdem112*p_insti_222,	p_ADdem$ADdem222_ADdem122*p_insti_222,	p_ADdem$ADdem222_ADdem132*p_insti_222,	p_ADdem$ADdem222_ADdem113*p_insti_222,	p_ADdem$ADdem222_ADdem123*p_insti_222,	p_ADdem$ADdem222_ADdem133*p_insti_222,	p_ADdem$ADdem222_ADdem211*p_insti_222,	p_ADdem$ADdem222_ADdem221*p_insti_222,	p_ADdem$ADdem222_ADdem231*p_insti_222,	p_ADdem$ADdem222_ADdem212*p_insti_222,	p_ADdem$ADdem222_ADdem222*p_insti_222,	p_ADdem$ADdem222_ADdem232*p_insti_222,	p_ADdem$ADdem222_ADdem213*p_insti_222,	p_ADdem$ADdem222_ADdem223*p_insti_222,	p_ADdem$ADdem222_ADdem233*p_insti_222,	p_ADdem$ADdem222_ADdem311*p_insti_222,	p_ADdem$ADdem222_ADdem321*p_insti_222,	p_ADdem$ADdem222_ADdem331*p_insti_222,	p_ADdem$ADdem222_ADdem312*p_insti_222,	p_ADdem$ADdem222_ADdem322*p_insti_222,	p_ADdem$ADdem222_ADdem332*p_insti_222,	p_ADdem$ADdem222_ADdem313*p_insti_222,	p_ADdem$ADdem222_ADdem323*p_insti_222,	p_ADdem$ADdem222_ADdem333*p_insti_222,	mrt_moderate,
    0,	p_ADdem$ADdem232_ADdem111*(1-p_insti_232),	p_ADdem$ADdem232_ADdem121*(1-p_insti_232),	p_ADdem$ADdem232_ADdem131*(1-p_insti_232),	p_ADdem$ADdem232_ADdem112*(1-p_insti_232),	p_ADdem$ADdem232_ADdem122*(1-p_insti_232),	p_ADdem$ADdem232_ADdem132*(1-p_insti_232),	p_ADdem$ADdem232_ADdem113*(1-p_insti_232),	p_ADdem$ADdem232_ADdem123*(1-p_insti_232),	p_ADdem$ADdem232_ADdem133*(1-p_insti_232),	p_ADdem$ADdem232_ADdem211*(1-p_insti_232),	p_ADdem$ADdem232_ADdem221*(1-p_insti_232),	p_ADdem$ADdem232_ADdem231*(1-p_insti_232),	p_ADdem$ADdem232_ADdem212*(1-p_insti_232),	p_ADdem$ADdem232_ADdem222*(1-p_insti_232),	C,	p_ADdem$ADdem232_ADdem213*(1-p_insti_232),	p_ADdem$ADdem232_ADdem223*(1-p_insti_232),	p_ADdem$ADdem232_ADdem233*(1-p_insti_232),	p_ADdem$ADdem232_ADdem311*(1-p_insti_232),	p_ADdem$ADdem232_ADdem321*(1-p_insti_232),	p_ADdem$ADdem232_ADdem331*(1-p_insti_232),	p_ADdem$ADdem232_ADdem312*(1-p_insti_232),	p_ADdem$ADdem232_ADdem322*(1-p_insti_232),	p_ADdem$ADdem232_ADdem332*(1-p_insti_232),	p_ADdem$ADdem232_ADdem313*(1-p_insti_232),	p_ADdem$ADdem232_ADdem323*(1-p_insti_232),	p_ADdem$ADdem232_ADdem333*(1-p_insti_232),	p_ADdem$ADdem232_ADdem111*p_insti_232,	p_ADdem$ADdem232_ADdem121*p_insti_232,	p_ADdem$ADdem232_ADdem131*p_insti_232,	p_ADdem$ADdem232_ADdem112*p_insti_232,	p_ADdem$ADdem232_ADdem122*p_insti_232,	p_ADdem$ADdem232_ADdem132*p_insti_232,	p_ADdem$ADdem232_ADdem113*p_insti_232,	p_ADdem$ADdem232_ADdem123*p_insti_232,	p_ADdem$ADdem232_ADdem133*p_insti_232,	p_ADdem$ADdem232_ADdem211*p_insti_232,	p_ADdem$ADdem232_ADdem221*p_insti_232,	p_ADdem$ADdem232_ADdem231*p_insti_232,	p_ADdem$ADdem232_ADdem212*p_insti_232,	p_ADdem$ADdem232_ADdem222*p_insti_232,	p_ADdem$ADdem232_ADdem232*p_insti_232,	p_ADdem$ADdem232_ADdem213*p_insti_232,	p_ADdem$ADdem232_ADdem223*p_insti_232,	p_ADdem$ADdem232_ADdem233*p_insti_232,	p_ADdem$ADdem232_ADdem311*p_insti_232,	p_ADdem$ADdem232_ADdem321*p_insti_232,	p_ADdem$ADdem232_ADdem331*p_insti_232,	p_ADdem$ADdem232_ADdem312*p_insti_232,	p_ADdem$ADdem232_ADdem322*p_insti_232,	p_ADdem$ADdem232_ADdem332*p_insti_232,	p_ADdem$ADdem232_ADdem313*p_insti_232,	p_ADdem$ADdem232_ADdem323*p_insti_232,	p_ADdem$ADdem232_ADdem333*p_insti_232,	mrt_moderate,
    0,	p_ADdem$ADdem213_ADdem111*(1-p_insti_213),	p_ADdem$ADdem213_ADdem121*(1-p_insti_213),	p_ADdem$ADdem213_ADdem131*(1-p_insti_213),	p_ADdem$ADdem213_ADdem112*(1-p_insti_213),	p_ADdem$ADdem213_ADdem122*(1-p_insti_213),	p_ADdem$ADdem213_ADdem132*(1-p_insti_213),	p_ADdem$ADdem213_ADdem113*(1-p_insti_213),	p_ADdem$ADdem213_ADdem123*(1-p_insti_213),	p_ADdem$ADdem213_ADdem133*(1-p_insti_213),	p_ADdem$ADdem213_ADdem211*(1-p_insti_213),	p_ADdem$ADdem213_ADdem221*(1-p_insti_213),	p_ADdem$ADdem213_ADdem231*(1-p_insti_213),	p_ADdem$ADdem213_ADdem212*(1-p_insti_213),	p_ADdem$ADdem213_ADdem222*(1-p_insti_213),	p_ADdem$ADdem213_ADdem232*(1-p_insti_213),	C,	p_ADdem$ADdem213_ADdem223*(1-p_insti_213),	p_ADdem$ADdem213_ADdem233*(1-p_insti_213),	p_ADdem$ADdem213_ADdem311*(1-p_insti_213),	p_ADdem$ADdem213_ADdem321*(1-p_insti_213),	p_ADdem$ADdem213_ADdem331*(1-p_insti_213),	p_ADdem$ADdem213_ADdem312*(1-p_insti_213),	p_ADdem$ADdem213_ADdem322*(1-p_insti_213),	p_ADdem$ADdem213_ADdem332*(1-p_insti_213),	p_ADdem$ADdem213_ADdem313*(1-p_insti_213),	p_ADdem$ADdem213_ADdem323*(1-p_insti_213),	p_ADdem$ADdem213_ADdem333*(1-p_insti_213),	p_ADdem$ADdem213_ADdem111*p_insti_213,	p_ADdem$ADdem213_ADdem121*p_insti_213,	p_ADdem$ADdem213_ADdem131*p_insti_213,	p_ADdem$ADdem213_ADdem112*p_insti_213,	p_ADdem$ADdem213_ADdem122*p_insti_213,	p_ADdem$ADdem213_ADdem132*p_insti_213,	p_ADdem$ADdem213_ADdem113*p_insti_213,	p_ADdem$ADdem213_ADdem123*p_insti_213,	p_ADdem$ADdem213_ADdem133*p_insti_213,	p_ADdem$ADdem213_ADdem211*p_insti_213,	p_ADdem$ADdem213_ADdem221*p_insti_213,	p_ADdem$ADdem213_ADdem231*p_insti_213,	p_ADdem$ADdem213_ADdem212*p_insti_213,	p_ADdem$ADdem213_ADdem222*p_insti_213,	p_ADdem$ADdem213_ADdem232*p_insti_213,	p_ADdem$ADdem213_ADdem213*p_insti_213,	p_ADdem$ADdem213_ADdem223*p_insti_213,	p_ADdem$ADdem213_ADdem233*p_insti_213,	p_ADdem$ADdem213_ADdem311*p_insti_213,	p_ADdem$ADdem213_ADdem321*p_insti_213,	p_ADdem$ADdem213_ADdem331*p_insti_213,	p_ADdem$ADdem213_ADdem312*p_insti_213,	p_ADdem$ADdem213_ADdem322*p_insti_213,	p_ADdem$ADdem213_ADdem332*p_insti_213,	p_ADdem$ADdem213_ADdem313*p_insti_213,	p_ADdem$ADdem213_ADdem323*p_insti_213,	p_ADdem$ADdem213_ADdem333*p_insti_213,	mrt_moderate,
    0,	p_ADdem$ADdem223_ADdem111*(1-p_insti_223),	p_ADdem$ADdem223_ADdem121*(1-p_insti_223),	p_ADdem$ADdem223_ADdem131*(1-p_insti_223),	p_ADdem$ADdem223_ADdem112*(1-p_insti_223),	p_ADdem$ADdem223_ADdem122*(1-p_insti_223),	p_ADdem$ADdem223_ADdem132*(1-p_insti_223),	p_ADdem$ADdem223_ADdem113*(1-p_insti_223),	p_ADdem$ADdem223_ADdem123*(1-p_insti_223),	p_ADdem$ADdem223_ADdem133*(1-p_insti_223),	p_ADdem$ADdem223_ADdem211*(1-p_insti_223),	p_ADdem$ADdem223_ADdem221*(1-p_insti_223),	p_ADdem$ADdem223_ADdem231*(1-p_insti_223),	p_ADdem$ADdem223_ADdem212*(1-p_insti_223),	p_ADdem$ADdem223_ADdem222*(1-p_insti_223),	p_ADdem$ADdem223_ADdem232*(1-p_insti_223),	p_ADdem$ADdem223_ADdem213*(1-p_insti_223),	C,	p_ADdem$ADdem223_ADdem233*(1-p_insti_223),	p_ADdem$ADdem223_ADdem311*(1-p_insti_223),	p_ADdem$ADdem223_ADdem321*(1-p_insti_223),	p_ADdem$ADdem223_ADdem331*(1-p_insti_223),	p_ADdem$ADdem223_ADdem312*(1-p_insti_223),	p_ADdem$ADdem223_ADdem322*(1-p_insti_223),	p_ADdem$ADdem223_ADdem332*(1-p_insti_223),	p_ADdem$ADdem223_ADdem313*(1-p_insti_223),	p_ADdem$ADdem223_ADdem323*(1-p_insti_223),	p_ADdem$ADdem223_ADdem333*(1-p_insti_223),	p_ADdem$ADdem223_ADdem111*p_insti_223,	p_ADdem$ADdem223_ADdem121*p_insti_223,	p_ADdem$ADdem223_ADdem131*p_insti_223,	p_ADdem$ADdem223_ADdem112*p_insti_223,	p_ADdem$ADdem223_ADdem122*p_insti_223,	p_ADdem$ADdem223_ADdem132*p_insti_223,	p_ADdem$ADdem223_ADdem113*p_insti_223,	p_ADdem$ADdem223_ADdem123*p_insti_223,	p_ADdem$ADdem223_ADdem133*p_insti_223,	p_ADdem$ADdem223_ADdem211*p_insti_223,	p_ADdem$ADdem223_ADdem221*p_insti_223,	p_ADdem$ADdem223_ADdem231*p_insti_223,	p_ADdem$ADdem223_ADdem212*p_insti_223,	p_ADdem$ADdem223_ADdem222*p_insti_223,	p_ADdem$ADdem223_ADdem232*p_insti_223,	p_ADdem$ADdem223_ADdem213*p_insti_223,	p_ADdem$ADdem223_ADdem223*p_insti_223,	p_ADdem$ADdem223_ADdem233*p_insti_223,	p_ADdem$ADdem223_ADdem311*p_insti_223,	p_ADdem$ADdem223_ADdem321*p_insti_223,	p_ADdem$ADdem223_ADdem331*p_insti_223,	p_ADdem$ADdem223_ADdem312*p_insti_223,	p_ADdem$ADdem223_ADdem322*p_insti_223,	p_ADdem$ADdem223_ADdem332*p_insti_223,	p_ADdem$ADdem223_ADdem313*p_insti_223,	p_ADdem$ADdem223_ADdem323*p_insti_223,	p_ADdem$ADdem223_ADdem333*p_insti_223,	mrt_moderate,
    0,	p_ADdem$ADdem233_ADdem111*(1-p_insti_233),	p_ADdem$ADdem233_ADdem121*(1-p_insti_233),	p_ADdem$ADdem233_ADdem131*(1-p_insti_233),	p_ADdem$ADdem233_ADdem112*(1-p_insti_233),	p_ADdem$ADdem233_ADdem122*(1-p_insti_233),	p_ADdem$ADdem233_ADdem132*(1-p_insti_233),	p_ADdem$ADdem233_ADdem113*(1-p_insti_233),	p_ADdem$ADdem233_ADdem123*(1-p_insti_233),	p_ADdem$ADdem233_ADdem133*(1-p_insti_233),	p_ADdem$ADdem233_ADdem211*(1-p_insti_233),	p_ADdem$ADdem233_ADdem221*(1-p_insti_233),	p_ADdem$ADdem233_ADdem231*(1-p_insti_233),	p_ADdem$ADdem233_ADdem212*(1-p_insti_233),	p_ADdem$ADdem233_ADdem222*(1-p_insti_233),	p_ADdem$ADdem233_ADdem232*(1-p_insti_233),	p_ADdem$ADdem233_ADdem213*(1-p_insti_233),	p_ADdem$ADdem233_ADdem223*(1-p_insti_233),	C,	p_ADdem$ADdem233_ADdem311*(1-p_insti_233),	p_ADdem$ADdem233_ADdem321*(1-p_insti_233),	p_ADdem$ADdem233_ADdem331*(1-p_insti_233),	p_ADdem$ADdem233_ADdem312*(1-p_insti_233),	p_ADdem$ADdem233_ADdem322*(1-p_insti_233),	p_ADdem$ADdem233_ADdem332*(1-p_insti_233),	p_ADdem$ADdem233_ADdem313*(1-p_insti_233),	p_ADdem$ADdem233_ADdem323*(1-p_insti_233),	p_ADdem$ADdem233_ADdem333*(1-p_insti_233),	p_ADdem$ADdem233_ADdem111*p_insti_233,	p_ADdem$ADdem233_ADdem121*p_insti_233,	p_ADdem$ADdem233_ADdem131*p_insti_233,	p_ADdem$ADdem233_ADdem112*p_insti_233,	p_ADdem$ADdem233_ADdem122*p_insti_233,	p_ADdem$ADdem233_ADdem132*p_insti_233,	p_ADdem$ADdem233_ADdem113*p_insti_233,	p_ADdem$ADdem233_ADdem123*p_insti_233,	p_ADdem$ADdem233_ADdem133*p_insti_233,	p_ADdem$ADdem233_ADdem211*p_insti_233,	p_ADdem$ADdem233_ADdem221*p_insti_233,	p_ADdem$ADdem233_ADdem231*p_insti_233,	p_ADdem$ADdem233_ADdem212*p_insti_233,	p_ADdem$ADdem233_ADdem222*p_insti_233,	p_ADdem$ADdem233_ADdem232*p_insti_233,	p_ADdem$ADdem233_ADdem213*p_insti_233,	p_ADdem$ADdem233_ADdem223*p_insti_233,	p_ADdem$ADdem233_ADdem233*p_insti_233,	p_ADdem$ADdem233_ADdem311*p_insti_233,	p_ADdem$ADdem233_ADdem321*p_insti_233,	p_ADdem$ADdem233_ADdem331*p_insti_233,	p_ADdem$ADdem233_ADdem312*p_insti_233,	p_ADdem$ADdem233_ADdem322*p_insti_233,	p_ADdem$ADdem233_ADdem332*p_insti_233,	p_ADdem$ADdem233_ADdem313*p_insti_233,	p_ADdem$ADdem233_ADdem323*p_insti_233,	p_ADdem$ADdem233_ADdem333*p_insti_233,	mrt_severe,
    0,	p_ADdem$ADdem311_ADdem111*(1-p_insti_311),	p_ADdem$ADdem311_ADdem121*(1-p_insti_311),	p_ADdem$ADdem311_ADdem131*(1-p_insti_311),	p_ADdem$ADdem311_ADdem112*(1-p_insti_311),	p_ADdem$ADdem311_ADdem122*(1-p_insti_311),	p_ADdem$ADdem311_ADdem132*(1-p_insti_311),	p_ADdem$ADdem311_ADdem113*(1-p_insti_311),	p_ADdem$ADdem311_ADdem123*(1-p_insti_311),	p_ADdem$ADdem311_ADdem133*(1-p_insti_311),	p_ADdem$ADdem311_ADdem211*(1-p_insti_311),	p_ADdem$ADdem311_ADdem221*(1-p_insti_311),	p_ADdem$ADdem311_ADdem231*(1-p_insti_311),	p_ADdem$ADdem311_ADdem212*(1-p_insti_311),	p_ADdem$ADdem311_ADdem222*(1-p_insti_311),	p_ADdem$ADdem311_ADdem232*(1-p_insti_311),	p_ADdem$ADdem311_ADdem213*(1-p_insti_311),	p_ADdem$ADdem311_ADdem223*(1-p_insti_311),	p_ADdem$ADdem311_ADdem233*(1-p_insti_311),	C,	p_ADdem$ADdem311_ADdem321*(1-p_insti_311),	p_ADdem$ADdem311_ADdem331*(1-p_insti_311),	p_ADdem$ADdem311_ADdem312*(1-p_insti_311),	p_ADdem$ADdem311_ADdem322*(1-p_insti_311),	p_ADdem$ADdem311_ADdem332*(1-p_insti_311),	p_ADdem$ADdem311_ADdem313*(1-p_insti_311),	p_ADdem$ADdem311_ADdem323*(1-p_insti_311),	p_ADdem$ADdem311_ADdem333*(1-p_insti_311),	p_ADdem$ADdem311_ADdem111*p_insti_311,	p_ADdem$ADdem311_ADdem121*p_insti_311,	p_ADdem$ADdem311_ADdem131*p_insti_311,	p_ADdem$ADdem311_ADdem112*p_insti_311,	p_ADdem$ADdem311_ADdem122*p_insti_311,	p_ADdem$ADdem311_ADdem132*p_insti_311,	p_ADdem$ADdem311_ADdem113*p_insti_311,	p_ADdem$ADdem311_ADdem123*p_insti_311,	p_ADdem$ADdem311_ADdem133*p_insti_311,	p_ADdem$ADdem311_ADdem211*p_insti_311,	p_ADdem$ADdem311_ADdem221*p_insti_311,	p_ADdem$ADdem311_ADdem231*p_insti_311,	p_ADdem$ADdem311_ADdem212*p_insti_311,	p_ADdem$ADdem311_ADdem222*p_insti_311,	p_ADdem$ADdem311_ADdem232*p_insti_311,	p_ADdem$ADdem311_ADdem213*p_insti_311,	p_ADdem$ADdem311_ADdem223*p_insti_311,	p_ADdem$ADdem311_ADdem233*p_insti_311,	p_ADdem$ADdem311_ADdem311*p_insti_311,	p_ADdem$ADdem311_ADdem321*p_insti_311,	p_ADdem$ADdem311_ADdem331*p_insti_311,	p_ADdem$ADdem311_ADdem312*p_insti_311,	p_ADdem$ADdem311_ADdem322*p_insti_311,	p_ADdem$ADdem311_ADdem332*p_insti_311,	p_ADdem$ADdem311_ADdem313*p_insti_311,	p_ADdem$ADdem311_ADdem323*p_insti_311,	p_ADdem$ADdem311_ADdem333*p_insti_311,	mrt_severe,
    0,	p_ADdem$ADdem321_ADdem111*(1-p_insti_321),	p_ADdem$ADdem321_ADdem121*(1-p_insti_321),	p_ADdem$ADdem321_ADdem131*(1-p_insti_321),	p_ADdem$ADdem321_ADdem112*(1-p_insti_321),	p_ADdem$ADdem321_ADdem122*(1-p_insti_321),	p_ADdem$ADdem321_ADdem132*(1-p_insti_321),	p_ADdem$ADdem321_ADdem113*(1-p_insti_321),	p_ADdem$ADdem321_ADdem123*(1-p_insti_321),	p_ADdem$ADdem321_ADdem133*(1-p_insti_321),	p_ADdem$ADdem321_ADdem211*(1-p_insti_321),	p_ADdem$ADdem321_ADdem221*(1-p_insti_321),	p_ADdem$ADdem321_ADdem231*(1-p_insti_321),	p_ADdem$ADdem321_ADdem212*(1-p_insti_321),	p_ADdem$ADdem321_ADdem222*(1-p_insti_321),	p_ADdem$ADdem321_ADdem232*(1-p_insti_321),	p_ADdem$ADdem321_ADdem213*(1-p_insti_321),	p_ADdem$ADdem321_ADdem223*(1-p_insti_321),	p_ADdem$ADdem321_ADdem233*(1-p_insti_321),	p_ADdem$ADdem321_ADdem311*(1-p_insti_321),	C,	p_ADdem$ADdem321_ADdem331*(1-p_insti_321),	p_ADdem$ADdem321_ADdem312*(1-p_insti_321),	p_ADdem$ADdem321_ADdem322*(1-p_insti_321),	p_ADdem$ADdem321_ADdem332*(1-p_insti_321),	p_ADdem$ADdem321_ADdem313*(1-p_insti_321),	p_ADdem$ADdem321_ADdem323*(1-p_insti_321),	p_ADdem$ADdem321_ADdem333*(1-p_insti_321),	p_ADdem$ADdem321_ADdem111*p_insti_321,	p_ADdem$ADdem321_ADdem121*p_insti_321,	p_ADdem$ADdem321_ADdem131*p_insti_321,	p_ADdem$ADdem321_ADdem112*p_insti_321,	p_ADdem$ADdem321_ADdem122*p_insti_321,	p_ADdem$ADdem321_ADdem132*p_insti_321,	p_ADdem$ADdem321_ADdem113*p_insti_321,	p_ADdem$ADdem321_ADdem123*p_insti_321,	p_ADdem$ADdem321_ADdem133*p_insti_321,	p_ADdem$ADdem321_ADdem211*p_insti_321,	p_ADdem$ADdem321_ADdem221*p_insti_321,	p_ADdem$ADdem321_ADdem231*p_insti_321,	p_ADdem$ADdem321_ADdem212*p_insti_321,	p_ADdem$ADdem321_ADdem222*p_insti_321,	p_ADdem$ADdem321_ADdem232*p_insti_321,	p_ADdem$ADdem321_ADdem213*p_insti_321,	p_ADdem$ADdem321_ADdem223*p_insti_321,	p_ADdem$ADdem321_ADdem233*p_insti_321,	p_ADdem$ADdem321_ADdem311*p_insti_321,	p_ADdem$ADdem321_ADdem321*p_insti_321,	p_ADdem$ADdem321_ADdem331*p_insti_321,	p_ADdem$ADdem321_ADdem312*p_insti_321,	p_ADdem$ADdem321_ADdem322*p_insti_321,	p_ADdem$ADdem321_ADdem332*p_insti_321,	p_ADdem$ADdem321_ADdem313*p_insti_321,	p_ADdem$ADdem321_ADdem323*p_insti_321,	p_ADdem$ADdem321_ADdem333*p_insti_321,	mrt_severe,
    0,	p_ADdem$ADdem331_ADdem111*(1-p_insti_331),	p_ADdem$ADdem331_ADdem121*(1-p_insti_331),	p_ADdem$ADdem331_ADdem131*(1-p_insti_331),	p_ADdem$ADdem331_ADdem112*(1-p_insti_331),	p_ADdem$ADdem331_ADdem122*(1-p_insti_331),	p_ADdem$ADdem331_ADdem132*(1-p_insti_331),	p_ADdem$ADdem331_ADdem113*(1-p_insti_331),	p_ADdem$ADdem331_ADdem123*(1-p_insti_331),	p_ADdem$ADdem331_ADdem133*(1-p_insti_331),	p_ADdem$ADdem331_ADdem211*(1-p_insti_331),	p_ADdem$ADdem331_ADdem221*(1-p_insti_331),	p_ADdem$ADdem331_ADdem231*(1-p_insti_331),	p_ADdem$ADdem331_ADdem212*(1-p_insti_331),	p_ADdem$ADdem331_ADdem222*(1-p_insti_331),	p_ADdem$ADdem331_ADdem232*(1-p_insti_331),	p_ADdem$ADdem331_ADdem213*(1-p_insti_331),	p_ADdem$ADdem331_ADdem223*(1-p_insti_331),	p_ADdem$ADdem331_ADdem233*(1-p_insti_331),	p_ADdem$ADdem331_ADdem311*(1-p_insti_331),	p_ADdem$ADdem331_ADdem321*(1-p_insti_331),	C,	p_ADdem$ADdem331_ADdem312*(1-p_insti_331),	p_ADdem$ADdem331_ADdem322*(1-p_insti_331),	p_ADdem$ADdem331_ADdem332*(1-p_insti_331),	p_ADdem$ADdem331_ADdem313*(1-p_insti_331),	p_ADdem$ADdem331_ADdem323*(1-p_insti_331),	p_ADdem$ADdem331_ADdem333*(1-p_insti_331),	p_ADdem$ADdem331_ADdem111*p_insti_331,	p_ADdem$ADdem331_ADdem121*p_insti_331,	p_ADdem$ADdem331_ADdem131*p_insti_331,	p_ADdem$ADdem331_ADdem112*p_insti_331,	p_ADdem$ADdem331_ADdem122*p_insti_331,	p_ADdem$ADdem331_ADdem132*p_insti_331,	p_ADdem$ADdem331_ADdem113*p_insti_331,	p_ADdem$ADdem331_ADdem123*p_insti_331,	p_ADdem$ADdem331_ADdem133*p_insti_331,	p_ADdem$ADdem331_ADdem211*p_insti_331,	p_ADdem$ADdem331_ADdem221*p_insti_331,	p_ADdem$ADdem331_ADdem231*p_insti_331,	p_ADdem$ADdem331_ADdem212*p_insti_331,	p_ADdem$ADdem331_ADdem222*p_insti_331,	p_ADdem$ADdem331_ADdem232*p_insti_331,	p_ADdem$ADdem331_ADdem213*p_insti_331,	p_ADdem$ADdem331_ADdem223*p_insti_331,	p_ADdem$ADdem331_ADdem233*p_insti_331,	p_ADdem$ADdem331_ADdem311*p_insti_331,	p_ADdem$ADdem331_ADdem321*p_insti_331,	p_ADdem$ADdem331_ADdem331*p_insti_331,	p_ADdem$ADdem331_ADdem312*p_insti_331,	p_ADdem$ADdem331_ADdem322*p_insti_331,	p_ADdem$ADdem331_ADdem332*p_insti_331,	p_ADdem$ADdem331_ADdem313*p_insti_331,	p_ADdem$ADdem331_ADdem323*p_insti_331,	p_ADdem$ADdem331_ADdem333*p_insti_331,	mrt_severe,
    0,	p_ADdem$ADdem312_ADdem111*(1-p_insti_312),	p_ADdem$ADdem312_ADdem121*(1-p_insti_312),	p_ADdem$ADdem312_ADdem131*(1-p_insti_312),	p_ADdem$ADdem312_ADdem112*(1-p_insti_312),	p_ADdem$ADdem312_ADdem122*(1-p_insti_312),	p_ADdem$ADdem312_ADdem132*(1-p_insti_312),	p_ADdem$ADdem312_ADdem113*(1-p_insti_312),	p_ADdem$ADdem312_ADdem123*(1-p_insti_312),	p_ADdem$ADdem312_ADdem133*(1-p_insti_312),	p_ADdem$ADdem312_ADdem211*(1-p_insti_312),	p_ADdem$ADdem312_ADdem221*(1-p_insti_312),	p_ADdem$ADdem312_ADdem231*(1-p_insti_312),	p_ADdem$ADdem312_ADdem212*(1-p_insti_312),	p_ADdem$ADdem312_ADdem222*(1-p_insti_312),	p_ADdem$ADdem312_ADdem232*(1-p_insti_312),	p_ADdem$ADdem312_ADdem213*(1-p_insti_312),	p_ADdem$ADdem312_ADdem223*(1-p_insti_312),	p_ADdem$ADdem312_ADdem233*(1-p_insti_312),	p_ADdem$ADdem312_ADdem311*(1-p_insti_312),	p_ADdem$ADdem312_ADdem321*(1-p_insti_312),	p_ADdem$ADdem312_ADdem331*(1-p_insti_312),	C,	p_ADdem$ADdem312_ADdem322*(1-p_insti_312),	p_ADdem$ADdem312_ADdem332*(1-p_insti_312),	p_ADdem$ADdem312_ADdem313*(1-p_insti_312),	p_ADdem$ADdem312_ADdem323*(1-p_insti_312),	p_ADdem$ADdem312_ADdem333*(1-p_insti_312),	p_ADdem$ADdem312_ADdem111*p_insti_312,	p_ADdem$ADdem312_ADdem121*p_insti_312,	p_ADdem$ADdem312_ADdem131*p_insti_312,	p_ADdem$ADdem312_ADdem112*p_insti_312,	p_ADdem$ADdem312_ADdem122*p_insti_312,	p_ADdem$ADdem312_ADdem132*p_insti_312,	p_ADdem$ADdem312_ADdem113*p_insti_312,	p_ADdem$ADdem312_ADdem123*p_insti_312,	p_ADdem$ADdem312_ADdem133*p_insti_312,	p_ADdem$ADdem312_ADdem211*p_insti_312,	p_ADdem$ADdem312_ADdem221*p_insti_312,	p_ADdem$ADdem312_ADdem231*p_insti_312,	p_ADdem$ADdem312_ADdem212*p_insti_312,	p_ADdem$ADdem312_ADdem222*p_insti_312,	p_ADdem$ADdem312_ADdem232*p_insti_312,	p_ADdem$ADdem312_ADdem213*p_insti_312,	p_ADdem$ADdem312_ADdem223*p_insti_312,	p_ADdem$ADdem312_ADdem233*p_insti_312,	p_ADdem$ADdem312_ADdem311*p_insti_312,	p_ADdem$ADdem312_ADdem321*p_insti_312,	p_ADdem$ADdem312_ADdem331*p_insti_312,	p_ADdem$ADdem312_ADdem312*p_insti_312,	p_ADdem$ADdem312_ADdem322*p_insti_312,	p_ADdem$ADdem312_ADdem332*p_insti_312,	p_ADdem$ADdem312_ADdem313*p_insti_312,	p_ADdem$ADdem312_ADdem323*p_insti_312,	p_ADdem$ADdem312_ADdem333*p_insti_312,	mrt_severe,
    0,	p_ADdem$ADdem322_ADdem111*(1-p_insti_322),	p_ADdem$ADdem322_ADdem121*(1-p_insti_322),	p_ADdem$ADdem322_ADdem131*(1-p_insti_322),	p_ADdem$ADdem322_ADdem112*(1-p_insti_322),	p_ADdem$ADdem322_ADdem122*(1-p_insti_322),	p_ADdem$ADdem322_ADdem132*(1-p_insti_322),	p_ADdem$ADdem322_ADdem113*(1-p_insti_322),	p_ADdem$ADdem322_ADdem123*(1-p_insti_322),	p_ADdem$ADdem322_ADdem133*(1-p_insti_322),	p_ADdem$ADdem322_ADdem211*(1-p_insti_322),	p_ADdem$ADdem322_ADdem221*(1-p_insti_322),	p_ADdem$ADdem322_ADdem231*(1-p_insti_322),	p_ADdem$ADdem322_ADdem212*(1-p_insti_322),	p_ADdem$ADdem322_ADdem222*(1-p_insti_322),	p_ADdem$ADdem322_ADdem232*(1-p_insti_322),	p_ADdem$ADdem322_ADdem213*(1-p_insti_322),	p_ADdem$ADdem322_ADdem223*(1-p_insti_322),	p_ADdem$ADdem322_ADdem233*(1-p_insti_322),	p_ADdem$ADdem322_ADdem311*(1-p_insti_322),	p_ADdem$ADdem322_ADdem321*(1-p_insti_322),	p_ADdem$ADdem322_ADdem331*(1-p_insti_322),	p_ADdem$ADdem322_ADdem312*(1-p_insti_322),	C,	p_ADdem$ADdem322_ADdem332*(1-p_insti_322),	p_ADdem$ADdem322_ADdem313*(1-p_insti_322),	p_ADdem$ADdem322_ADdem323*(1-p_insti_322),	p_ADdem$ADdem322_ADdem333*(1-p_insti_322),	p_ADdem$ADdem322_ADdem111*p_insti_322,	p_ADdem$ADdem322_ADdem121*p_insti_322,	p_ADdem$ADdem322_ADdem131*p_insti_322,	p_ADdem$ADdem322_ADdem112*p_insti_322,	p_ADdem$ADdem322_ADdem122*p_insti_322,	p_ADdem$ADdem322_ADdem132*p_insti_322,	p_ADdem$ADdem322_ADdem113*p_insti_322,	p_ADdem$ADdem322_ADdem123*p_insti_322,	p_ADdem$ADdem322_ADdem133*p_insti_322,	p_ADdem$ADdem322_ADdem211*p_insti_322,	p_ADdem$ADdem322_ADdem221*p_insti_322,	p_ADdem$ADdem322_ADdem231*p_insti_322,	p_ADdem$ADdem322_ADdem212*p_insti_322,	p_ADdem$ADdem322_ADdem222*p_insti_322,	p_ADdem$ADdem322_ADdem232*p_insti_322,	p_ADdem$ADdem322_ADdem213*p_insti_322,	p_ADdem$ADdem322_ADdem223*p_insti_322,	p_ADdem$ADdem322_ADdem233*p_insti_322,	p_ADdem$ADdem322_ADdem311*p_insti_322,	p_ADdem$ADdem322_ADdem321*p_insti_322,	p_ADdem$ADdem322_ADdem331*p_insti_322,	p_ADdem$ADdem322_ADdem312*p_insti_322,	p_ADdem$ADdem322_ADdem322*p_insti_322,	p_ADdem$ADdem322_ADdem332*p_insti_322,	p_ADdem$ADdem322_ADdem313*p_insti_322,	p_ADdem$ADdem322_ADdem323*p_insti_322,	p_ADdem$ADdem322_ADdem333*p_insti_322,	mrt_severe,
    0,	p_ADdem$ADdem332_ADdem111*(1-p_insti_332),	p_ADdem$ADdem332_ADdem121*(1-p_insti_332),	p_ADdem$ADdem332_ADdem131*(1-p_insti_332),	p_ADdem$ADdem332_ADdem112*(1-p_insti_332),	p_ADdem$ADdem332_ADdem122*(1-p_insti_332),	p_ADdem$ADdem332_ADdem132*(1-p_insti_332),	p_ADdem$ADdem332_ADdem113*(1-p_insti_332),	p_ADdem$ADdem332_ADdem123*(1-p_insti_332),	p_ADdem$ADdem332_ADdem133*(1-p_insti_332),	p_ADdem$ADdem332_ADdem211*(1-p_insti_332),	p_ADdem$ADdem332_ADdem221*(1-p_insti_332),	p_ADdem$ADdem332_ADdem231*(1-p_insti_332),	p_ADdem$ADdem332_ADdem212*(1-p_insti_332),	p_ADdem$ADdem332_ADdem222*(1-p_insti_332),	p_ADdem$ADdem332_ADdem232*(1-p_insti_332),	p_ADdem$ADdem332_ADdem213*(1-p_insti_332),	p_ADdem$ADdem332_ADdem223*(1-p_insti_332),	p_ADdem$ADdem332_ADdem233*(1-p_insti_332),	p_ADdem$ADdem332_ADdem311*(1-p_insti_332),	p_ADdem$ADdem332_ADdem321*(1-p_insti_332),	p_ADdem$ADdem332_ADdem331*(1-p_insti_332),	p_ADdem$ADdem332_ADdem312*(1-p_insti_332),	p_ADdem$ADdem332_ADdem322*(1-p_insti_332),	C,	p_ADdem$ADdem332_ADdem313*(1-p_insti_332),	p_ADdem$ADdem332_ADdem323*(1-p_insti_332),	p_ADdem$ADdem332_ADdem333*(1-p_insti_332),	p_ADdem$ADdem332_ADdem111*p_insti_332,	p_ADdem$ADdem332_ADdem121*p_insti_332,	p_ADdem$ADdem332_ADdem131*p_insti_332,	p_ADdem$ADdem332_ADdem112*p_insti_332,	p_ADdem$ADdem332_ADdem122*p_insti_332,	p_ADdem$ADdem332_ADdem132*p_insti_332,	p_ADdem$ADdem332_ADdem113*p_insti_332,	p_ADdem$ADdem332_ADdem123*p_insti_332,	p_ADdem$ADdem332_ADdem133*p_insti_332,	p_ADdem$ADdem332_ADdem211*p_insti_332,	p_ADdem$ADdem332_ADdem221*p_insti_332,	p_ADdem$ADdem332_ADdem231*p_insti_332,	p_ADdem$ADdem332_ADdem212*p_insti_332,	p_ADdem$ADdem332_ADdem222*p_insti_332,	p_ADdem$ADdem332_ADdem232*p_insti_332,	p_ADdem$ADdem332_ADdem213*p_insti_332,	p_ADdem$ADdem332_ADdem223*p_insti_332,	p_ADdem$ADdem332_ADdem233*p_insti_332,	p_ADdem$ADdem332_ADdem311*p_insti_332,	p_ADdem$ADdem332_ADdem321*p_insti_332,	p_ADdem$ADdem332_ADdem331*p_insti_332,	p_ADdem$ADdem332_ADdem312*p_insti_332,	p_ADdem$ADdem332_ADdem322*p_insti_332,	p_ADdem$ADdem332_ADdem332*p_insti_332,	p_ADdem$ADdem332_ADdem313*p_insti_332,	p_ADdem$ADdem332_ADdem323*p_insti_332,	p_ADdem$ADdem332_ADdem333*p_insti_332,	mrt_severe,
    0,	p_ADdem$ADdem313_ADdem111*(1-p_insti_313),	p_ADdem$ADdem313_ADdem121*(1-p_insti_313),	p_ADdem$ADdem313_ADdem131*(1-p_insti_313),	p_ADdem$ADdem313_ADdem112*(1-p_insti_313),	p_ADdem$ADdem313_ADdem122*(1-p_insti_313),	p_ADdem$ADdem313_ADdem132*(1-p_insti_313),	p_ADdem$ADdem313_ADdem113*(1-p_insti_313),	p_ADdem$ADdem313_ADdem123*(1-p_insti_313),	p_ADdem$ADdem313_ADdem133*(1-p_insti_313),	p_ADdem$ADdem313_ADdem211*(1-p_insti_313),	p_ADdem$ADdem313_ADdem221*(1-p_insti_313),	p_ADdem$ADdem313_ADdem231*(1-p_insti_313),	p_ADdem$ADdem313_ADdem212*(1-p_insti_313),	p_ADdem$ADdem313_ADdem222*(1-p_insti_313),	p_ADdem$ADdem313_ADdem232*(1-p_insti_313),	p_ADdem$ADdem313_ADdem213*(1-p_insti_313),	p_ADdem$ADdem313_ADdem223*(1-p_insti_313),	p_ADdem$ADdem313_ADdem233*(1-p_insti_313),	p_ADdem$ADdem313_ADdem311*(1-p_insti_313),	p_ADdem$ADdem313_ADdem321*(1-p_insti_313),	p_ADdem$ADdem313_ADdem331*(1-p_insti_313),	p_ADdem$ADdem313_ADdem312*(1-p_insti_313),	p_ADdem$ADdem313_ADdem322*(1-p_insti_313),	p_ADdem$ADdem313_ADdem332*(1-p_insti_313),	C,	p_ADdem$ADdem313_ADdem323*(1-p_insti_313),	p_ADdem$ADdem313_ADdem333*(1-p_insti_313),	p_ADdem$ADdem313_ADdem111*p_insti_313,	p_ADdem$ADdem313_ADdem121*p_insti_313,	p_ADdem$ADdem313_ADdem131*p_insti_313,	p_ADdem$ADdem313_ADdem112*p_insti_313,	p_ADdem$ADdem313_ADdem122*p_insti_313,	p_ADdem$ADdem313_ADdem132*p_insti_313,	p_ADdem$ADdem313_ADdem113*p_insti_313,	p_ADdem$ADdem313_ADdem123*p_insti_313,	p_ADdem$ADdem313_ADdem133*p_insti_313,	p_ADdem$ADdem313_ADdem211*p_insti_313,	p_ADdem$ADdem313_ADdem221*p_insti_313,	p_ADdem$ADdem313_ADdem231*p_insti_313,	p_ADdem$ADdem313_ADdem212*p_insti_313,	p_ADdem$ADdem313_ADdem222*p_insti_313,	p_ADdem$ADdem313_ADdem232*p_insti_313,	p_ADdem$ADdem313_ADdem213*p_insti_313,	p_ADdem$ADdem313_ADdem223*p_insti_313,	p_ADdem$ADdem313_ADdem233*p_insti_313,	p_ADdem$ADdem313_ADdem311*p_insti_313,	p_ADdem$ADdem313_ADdem321*p_insti_313,	p_ADdem$ADdem313_ADdem331*p_insti_313,	p_ADdem$ADdem313_ADdem312*p_insti_313,	p_ADdem$ADdem313_ADdem322*p_insti_313,	p_ADdem$ADdem313_ADdem332*p_insti_313,	p_ADdem$ADdem313_ADdem313*p_insti_313,	p_ADdem$ADdem313_ADdem323*p_insti_313,	p_ADdem$ADdem313_ADdem333*p_insti_313,	mrt_severe,
    0,	p_ADdem$ADdem323_ADdem111*(1-p_insti_323),	p_ADdem$ADdem323_ADdem121*(1-p_insti_323),	p_ADdem$ADdem323_ADdem131*(1-p_insti_323),	p_ADdem$ADdem323_ADdem112*(1-p_insti_323),	p_ADdem$ADdem323_ADdem122*(1-p_insti_323),	p_ADdem$ADdem323_ADdem132*(1-p_insti_323),	p_ADdem$ADdem323_ADdem113*(1-p_insti_323),	p_ADdem$ADdem323_ADdem123*(1-p_insti_323),	p_ADdem$ADdem323_ADdem133*(1-p_insti_323),	p_ADdem$ADdem323_ADdem211*(1-p_insti_323),	p_ADdem$ADdem323_ADdem221*(1-p_insti_323),	p_ADdem$ADdem323_ADdem231*(1-p_insti_323),	p_ADdem$ADdem323_ADdem212*(1-p_insti_323),	p_ADdem$ADdem323_ADdem222*(1-p_insti_323),	p_ADdem$ADdem323_ADdem232*(1-p_insti_323),	p_ADdem$ADdem323_ADdem213*(1-p_insti_323),	p_ADdem$ADdem323_ADdem223*(1-p_insti_323),	p_ADdem$ADdem323_ADdem233*(1-p_insti_323),	p_ADdem$ADdem323_ADdem311*(1-p_insti_323),	p_ADdem$ADdem323_ADdem321*(1-p_insti_323),	p_ADdem$ADdem323_ADdem331*(1-p_insti_323),	p_ADdem$ADdem323_ADdem312*(1-p_insti_323),	p_ADdem$ADdem323_ADdem322*(1-p_insti_323),	p_ADdem$ADdem323_ADdem332*(1-p_insti_323),	p_ADdem$ADdem323_ADdem313*(1-p_insti_323),	C,	p_ADdem$ADdem323_ADdem333*(1-p_insti_323),	p_ADdem$ADdem323_ADdem111*p_insti_323,	p_ADdem$ADdem323_ADdem121*p_insti_323,	p_ADdem$ADdem323_ADdem131*p_insti_323,	p_ADdem$ADdem323_ADdem112*p_insti_323,	p_ADdem$ADdem323_ADdem122*p_insti_323,	p_ADdem$ADdem323_ADdem132*p_insti_323,	p_ADdem$ADdem323_ADdem113*p_insti_323,	p_ADdem$ADdem323_ADdem123*p_insti_323,	p_ADdem$ADdem323_ADdem133*p_insti_323,	p_ADdem$ADdem323_ADdem211*p_insti_323,	p_ADdem$ADdem323_ADdem221*p_insti_323,	p_ADdem$ADdem323_ADdem231*p_insti_323,	p_ADdem$ADdem323_ADdem212*p_insti_323,	p_ADdem$ADdem323_ADdem222*p_insti_323,	p_ADdem$ADdem323_ADdem232*p_insti_323,	p_ADdem$ADdem323_ADdem213*p_insti_323,	p_ADdem$ADdem323_ADdem223*p_insti_323,	p_ADdem$ADdem323_ADdem233*p_insti_323,	p_ADdem$ADdem323_ADdem311*p_insti_323,	p_ADdem$ADdem323_ADdem321*p_insti_323,	p_ADdem$ADdem323_ADdem331*p_insti_323,	p_ADdem$ADdem323_ADdem312*p_insti_323,	p_ADdem$ADdem323_ADdem322*p_insti_323,	p_ADdem$ADdem323_ADdem332*p_insti_323,	p_ADdem$ADdem323_ADdem313*p_insti_323,	p_ADdem$ADdem323_ADdem323*p_insti_323,	p_ADdem$ADdem323_ADdem333*p_insti_323,	mrt_severe,
    0,	p_ADdem$ADdem333_ADdem111*(1-p_insti_333),	p_ADdem$ADdem333_ADdem121*(1-p_insti_333),	p_ADdem$ADdem333_ADdem131*(1-p_insti_333),	p_ADdem$ADdem333_ADdem112*(1-p_insti_333),	p_ADdem$ADdem333_ADdem122*(1-p_insti_333),	p_ADdem$ADdem333_ADdem132*(1-p_insti_333),	p_ADdem$ADdem333_ADdem113*(1-p_insti_333),	p_ADdem$ADdem333_ADdem123*(1-p_insti_333),	p_ADdem$ADdem333_ADdem133*(1-p_insti_333),	p_ADdem$ADdem333_ADdem211*(1-p_insti_333),	p_ADdem$ADdem333_ADdem221*(1-p_insti_333),	p_ADdem$ADdem333_ADdem231*(1-p_insti_333),	p_ADdem$ADdem333_ADdem212*(1-p_insti_333),	p_ADdem$ADdem333_ADdem222*(1-p_insti_333),	p_ADdem$ADdem333_ADdem232*(1-p_insti_333),	p_ADdem$ADdem333_ADdem213*(1-p_insti_333),	p_ADdem$ADdem333_ADdem223*(1-p_insti_333),	p_ADdem$ADdem333_ADdem233*(1-p_insti_333),	p_ADdem$ADdem333_ADdem311*(1-p_insti_333),	p_ADdem$ADdem333_ADdem321*(1-p_insti_333),	p_ADdem$ADdem333_ADdem331*(1-p_insti_333),	p_ADdem$ADdem333_ADdem312*(1-p_insti_333),	p_ADdem$ADdem333_ADdem322*(1-p_insti_333),	p_ADdem$ADdem333_ADdem332*(1-p_insti_333),	p_ADdem$ADdem333_ADdem313*(1-p_insti_333),	p_ADdem$ADdem333_ADdem323*(1-p_insti_333),	C,	p_ADdem$ADdem333_ADdem111*p_insti_333,	p_ADdem$ADdem333_ADdem121*p_insti_333,	p_ADdem$ADdem333_ADdem131*p_insti_333,	p_ADdem$ADdem333_ADdem112*p_insti_333,	p_ADdem$ADdem333_ADdem122*p_insti_333,	p_ADdem$ADdem333_ADdem132*p_insti_333,	p_ADdem$ADdem333_ADdem113*p_insti_333,	p_ADdem$ADdem333_ADdem123*p_insti_333,	p_ADdem$ADdem333_ADdem133*p_insti_333,	p_ADdem$ADdem333_ADdem211*p_insti_333,	p_ADdem$ADdem333_ADdem221*p_insti_333,	p_ADdem$ADdem333_ADdem231*p_insti_333,	p_ADdem$ADdem333_ADdem212*p_insti_333,	p_ADdem$ADdem333_ADdem222*p_insti_333,	p_ADdem$ADdem333_ADdem232*p_insti_333,	p_ADdem$ADdem333_ADdem213*p_insti_333,	p_ADdem$ADdem333_ADdem223*p_insti_333,	p_ADdem$ADdem333_ADdem233*p_insti_333,	p_ADdem$ADdem333_ADdem311*p_insti_333,	p_ADdem$ADdem333_ADdem321*p_insti_333,	p_ADdem$ADdem333_ADdem331*p_insti_333,	p_ADdem$ADdem333_ADdem312*p_insti_333,	p_ADdem$ADdem333_ADdem322*p_insti_333,	p_ADdem$ADdem333_ADdem332*p_insti_333,	p_ADdem$ADdem333_ADdem313*p_insti_333,	p_ADdem$ADdem333_ADdem323*p_insti_333,	p_ADdem$ADdem333_ADdem333*p_insti_333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	C,	p_ADdem$ADdem111_ADdem121,	p_ADdem$ADdem111_ADdem131,	p_ADdem$ADdem111_ADdem112,	p_ADdem$ADdem111_ADdem122,	p_ADdem$ADdem111_ADdem132,	p_ADdem$ADdem111_ADdem113,	p_ADdem$ADdem111_ADdem123,	p_ADdem$ADdem111_ADdem133,	p_ADdem$ADdem111_ADdem211,	p_ADdem$ADdem111_ADdem221,	p_ADdem$ADdem111_ADdem231,	p_ADdem$ADdem111_ADdem212,	p_ADdem$ADdem111_ADdem222,	p_ADdem$ADdem111_ADdem232,	p_ADdem$ADdem111_ADdem213,	p_ADdem$ADdem111_ADdem223,	p_ADdem$ADdem111_ADdem233,	p_ADdem$ADdem111_ADdem311,	p_ADdem$ADdem111_ADdem321,	p_ADdem$ADdem111_ADdem331,	p_ADdem$ADdem111_ADdem312,	p_ADdem$ADdem111_ADdem322,	p_ADdem$ADdem111_ADdem332,	p_ADdem$ADdem111_ADdem313,	p_ADdem$ADdem111_ADdem323,	p_ADdem$ADdem111_ADdem333,	mrt_mild,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem121_ADdem111,	C,	p_ADdem$ADdem121_ADdem131,	p_ADdem$ADdem121_ADdem112,	p_ADdem$ADdem121_ADdem122,	p_ADdem$ADdem121_ADdem132,	p_ADdem$ADdem121_ADdem113,	p_ADdem$ADdem121_ADdem123,	p_ADdem$ADdem121_ADdem133,	p_ADdem$ADdem121_ADdem211,	p_ADdem$ADdem121_ADdem221,	p_ADdem$ADdem121_ADdem231,	p_ADdem$ADdem121_ADdem212,	p_ADdem$ADdem121_ADdem222,	p_ADdem$ADdem121_ADdem232,	p_ADdem$ADdem121_ADdem213,	p_ADdem$ADdem121_ADdem223,	p_ADdem$ADdem121_ADdem233,	p_ADdem$ADdem121_ADdem311,	p_ADdem$ADdem121_ADdem321,	p_ADdem$ADdem121_ADdem331,	p_ADdem$ADdem121_ADdem312,	p_ADdem$ADdem121_ADdem322,	p_ADdem$ADdem121_ADdem332,	p_ADdem$ADdem121_ADdem313,	p_ADdem$ADdem121_ADdem323,	p_ADdem$ADdem121_ADdem333,	mrt_mild,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem131_ADdem111,	p_ADdem$ADdem131_ADdem121,	C,	p_ADdem$ADdem131_ADdem112,	p_ADdem$ADdem131_ADdem122,	p_ADdem$ADdem131_ADdem132,	p_ADdem$ADdem131_ADdem113,	p_ADdem$ADdem131_ADdem123,	p_ADdem$ADdem131_ADdem133,	p_ADdem$ADdem131_ADdem211,	p_ADdem$ADdem131_ADdem221,	p_ADdem$ADdem131_ADdem231,	p_ADdem$ADdem131_ADdem212,	p_ADdem$ADdem131_ADdem222,	p_ADdem$ADdem131_ADdem232,	p_ADdem$ADdem131_ADdem213,	p_ADdem$ADdem131_ADdem223,	p_ADdem$ADdem131_ADdem233,	p_ADdem$ADdem131_ADdem311,	p_ADdem$ADdem131_ADdem321,	p_ADdem$ADdem131_ADdem331,	p_ADdem$ADdem131_ADdem312,	p_ADdem$ADdem131_ADdem322,	p_ADdem$ADdem131_ADdem332,	p_ADdem$ADdem131_ADdem313,	p_ADdem$ADdem131_ADdem323,	p_ADdem$ADdem131_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem112_ADdem111,	p_ADdem$ADdem112_ADdem121,	p_ADdem$ADdem112_ADdem131,	C,	p_ADdem$ADdem112_ADdem122,	p_ADdem$ADdem112_ADdem132,	p_ADdem$ADdem112_ADdem113,	p_ADdem$ADdem112_ADdem123,	p_ADdem$ADdem112_ADdem133,	p_ADdem$ADdem112_ADdem211,	p_ADdem$ADdem112_ADdem221,	p_ADdem$ADdem112_ADdem231,	p_ADdem$ADdem112_ADdem212,	p_ADdem$ADdem112_ADdem222,	p_ADdem$ADdem112_ADdem232,	p_ADdem$ADdem112_ADdem213,	p_ADdem$ADdem112_ADdem223,	p_ADdem$ADdem112_ADdem233,	p_ADdem$ADdem112_ADdem311,	p_ADdem$ADdem112_ADdem321,	p_ADdem$ADdem112_ADdem331,	p_ADdem$ADdem112_ADdem312,	p_ADdem$ADdem112_ADdem322,	p_ADdem$ADdem112_ADdem332,	p_ADdem$ADdem112_ADdem313,	p_ADdem$ADdem112_ADdem323,	p_ADdem$ADdem112_ADdem333,	mrt_mild,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem122_ADdem111,	p_ADdem$ADdem122_ADdem121,	p_ADdem$ADdem122_ADdem131,	p_ADdem$ADdem122_ADdem112,	C,	p_ADdem$ADdem122_ADdem132,	p_ADdem$ADdem122_ADdem113,	p_ADdem$ADdem122_ADdem123,	p_ADdem$ADdem122_ADdem133,	p_ADdem$ADdem122_ADdem211,	p_ADdem$ADdem122_ADdem221,	p_ADdem$ADdem122_ADdem231,	p_ADdem$ADdem122_ADdem212,	p_ADdem$ADdem122_ADdem222,	p_ADdem$ADdem122_ADdem232,	p_ADdem$ADdem122_ADdem213,	p_ADdem$ADdem122_ADdem223,	p_ADdem$ADdem122_ADdem233,	p_ADdem$ADdem122_ADdem311,	p_ADdem$ADdem122_ADdem321,	p_ADdem$ADdem122_ADdem331,	p_ADdem$ADdem122_ADdem312,	p_ADdem$ADdem122_ADdem322,	p_ADdem$ADdem122_ADdem332,	p_ADdem$ADdem122_ADdem313,	p_ADdem$ADdem122_ADdem323,	p_ADdem$ADdem122_ADdem333,	mrt_mild,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem132_ADdem111,	p_ADdem$ADdem132_ADdem121,	p_ADdem$ADdem132_ADdem131,	p_ADdem$ADdem132_ADdem112,	p_ADdem$ADdem132_ADdem122,	C,	p_ADdem$ADdem132_ADdem113,	p_ADdem$ADdem132_ADdem123,	p_ADdem$ADdem132_ADdem133,	p_ADdem$ADdem132_ADdem211,	p_ADdem$ADdem132_ADdem221,	p_ADdem$ADdem132_ADdem231,	p_ADdem$ADdem132_ADdem212,	p_ADdem$ADdem132_ADdem222,	p_ADdem$ADdem132_ADdem232,	p_ADdem$ADdem132_ADdem213,	p_ADdem$ADdem132_ADdem223,	p_ADdem$ADdem132_ADdem233,	p_ADdem$ADdem132_ADdem311,	p_ADdem$ADdem132_ADdem321,	p_ADdem$ADdem132_ADdem331,	p_ADdem$ADdem132_ADdem312,	p_ADdem$ADdem132_ADdem322,	p_ADdem$ADdem132_ADdem332,	p_ADdem$ADdem132_ADdem313,	p_ADdem$ADdem132_ADdem323,	p_ADdem$ADdem132_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem113_ADdem111,	p_ADdem$ADdem113_ADdem121,	p_ADdem$ADdem113_ADdem131,	p_ADdem$ADdem113_ADdem112,	p_ADdem$ADdem113_ADdem122,	p_ADdem$ADdem113_ADdem132,	C,	p_ADdem$ADdem113_ADdem123,	p_ADdem$ADdem113_ADdem133,	p_ADdem$ADdem113_ADdem211,	p_ADdem$ADdem113_ADdem221,	p_ADdem$ADdem113_ADdem231,	p_ADdem$ADdem113_ADdem212,	p_ADdem$ADdem113_ADdem222,	p_ADdem$ADdem113_ADdem232,	p_ADdem$ADdem113_ADdem213,	p_ADdem$ADdem113_ADdem223,	p_ADdem$ADdem113_ADdem233,	p_ADdem$ADdem113_ADdem311,	p_ADdem$ADdem113_ADdem321,	p_ADdem$ADdem113_ADdem331,	p_ADdem$ADdem113_ADdem312,	p_ADdem$ADdem113_ADdem322,	p_ADdem$ADdem113_ADdem332,	p_ADdem$ADdem113_ADdem313,	p_ADdem$ADdem113_ADdem323,	p_ADdem$ADdem113_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem123_ADdem111,	p_ADdem$ADdem123_ADdem121,	p_ADdem$ADdem123_ADdem131,	p_ADdem$ADdem123_ADdem112,	p_ADdem$ADdem123_ADdem122,	p_ADdem$ADdem123_ADdem132,	p_ADdem$ADdem123_ADdem113,	C,	p_ADdem$ADdem123_ADdem133,	p_ADdem$ADdem123_ADdem211,	p_ADdem$ADdem123_ADdem221,	p_ADdem$ADdem123_ADdem231,	p_ADdem$ADdem123_ADdem212,	p_ADdem$ADdem123_ADdem222,	p_ADdem$ADdem123_ADdem232,	p_ADdem$ADdem123_ADdem213,	p_ADdem$ADdem123_ADdem223,	p_ADdem$ADdem123_ADdem233,	p_ADdem$ADdem123_ADdem311,	p_ADdem$ADdem123_ADdem321,	p_ADdem$ADdem123_ADdem331,	p_ADdem$ADdem123_ADdem312,	p_ADdem$ADdem123_ADdem322,	p_ADdem$ADdem123_ADdem332,	p_ADdem$ADdem123_ADdem313,	p_ADdem$ADdem123_ADdem323,	p_ADdem$ADdem123_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem133_ADdem111,	p_ADdem$ADdem133_ADdem121,	p_ADdem$ADdem133_ADdem131,	p_ADdem$ADdem133_ADdem112,	p_ADdem$ADdem133_ADdem122,	p_ADdem$ADdem133_ADdem132,	p_ADdem$ADdem133_ADdem113,	p_ADdem$ADdem133_ADdem123,	C,	p_ADdem$ADdem133_ADdem211,	p_ADdem$ADdem133_ADdem221,	p_ADdem$ADdem133_ADdem231,	p_ADdem$ADdem133_ADdem212,	p_ADdem$ADdem133_ADdem222,	p_ADdem$ADdem133_ADdem232,	p_ADdem$ADdem133_ADdem213,	p_ADdem$ADdem133_ADdem223,	p_ADdem$ADdem133_ADdem233,	p_ADdem$ADdem133_ADdem311,	p_ADdem$ADdem133_ADdem321,	p_ADdem$ADdem133_ADdem331,	p_ADdem$ADdem133_ADdem312,	p_ADdem$ADdem133_ADdem322,	p_ADdem$ADdem133_ADdem332,	p_ADdem$ADdem133_ADdem313,	p_ADdem$ADdem133_ADdem323,	p_ADdem$ADdem133_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem211_ADdem111,	p_ADdem$ADdem211_ADdem121,	p_ADdem$ADdem211_ADdem131,	p_ADdem$ADdem211_ADdem112,	p_ADdem$ADdem211_ADdem122,	p_ADdem$ADdem211_ADdem132,	p_ADdem$ADdem211_ADdem113,	p_ADdem$ADdem211_ADdem123,	p_ADdem$ADdem211_ADdem133,	C,	p_ADdem$ADdem211_ADdem221,	p_ADdem$ADdem211_ADdem231,	p_ADdem$ADdem211_ADdem212,	p_ADdem$ADdem211_ADdem222,	p_ADdem$ADdem211_ADdem232,	p_ADdem$ADdem211_ADdem213,	p_ADdem$ADdem211_ADdem223,	p_ADdem$ADdem211_ADdem233,	p_ADdem$ADdem211_ADdem311,	p_ADdem$ADdem211_ADdem321,	p_ADdem$ADdem211_ADdem331,	p_ADdem$ADdem211_ADdem312,	p_ADdem$ADdem211_ADdem322,	p_ADdem$ADdem211_ADdem332,	p_ADdem$ADdem211_ADdem313,	p_ADdem$ADdem211_ADdem323,	p_ADdem$ADdem211_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem221_ADdem111,	p_ADdem$ADdem221_ADdem121,	p_ADdem$ADdem221_ADdem131,	p_ADdem$ADdem221_ADdem112,	p_ADdem$ADdem221_ADdem122,	p_ADdem$ADdem221_ADdem132,	p_ADdem$ADdem221_ADdem113,	p_ADdem$ADdem221_ADdem123,	p_ADdem$ADdem221_ADdem133,	p_ADdem$ADdem221_ADdem211,	C,	p_ADdem$ADdem221_ADdem231,	p_ADdem$ADdem221_ADdem212,	p_ADdem$ADdem221_ADdem222,	p_ADdem$ADdem221_ADdem232,	p_ADdem$ADdem221_ADdem213,	p_ADdem$ADdem221_ADdem223,	p_ADdem$ADdem221_ADdem233,	p_ADdem$ADdem221_ADdem311,	p_ADdem$ADdem221_ADdem321,	p_ADdem$ADdem221_ADdem331,	p_ADdem$ADdem221_ADdem312,	p_ADdem$ADdem221_ADdem322,	p_ADdem$ADdem221_ADdem332,	p_ADdem$ADdem221_ADdem313,	p_ADdem$ADdem221_ADdem323,	p_ADdem$ADdem221_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem231_ADdem111,	p_ADdem$ADdem231_ADdem121,	p_ADdem$ADdem231_ADdem131,	p_ADdem$ADdem231_ADdem112,	p_ADdem$ADdem231_ADdem122,	p_ADdem$ADdem231_ADdem132,	p_ADdem$ADdem231_ADdem113,	p_ADdem$ADdem231_ADdem123,	p_ADdem$ADdem231_ADdem133,	p_ADdem$ADdem231_ADdem211,	p_ADdem$ADdem231_ADdem221,	C,	p_ADdem$ADdem231_ADdem212,	p_ADdem$ADdem231_ADdem222,	p_ADdem$ADdem231_ADdem232,	p_ADdem$ADdem231_ADdem213,	p_ADdem$ADdem231_ADdem223,	p_ADdem$ADdem231_ADdem233,	p_ADdem$ADdem231_ADdem311,	p_ADdem$ADdem231_ADdem321,	p_ADdem$ADdem231_ADdem331,	p_ADdem$ADdem231_ADdem312,	p_ADdem$ADdem231_ADdem322,	p_ADdem$ADdem231_ADdem332,	p_ADdem$ADdem231_ADdem313,	p_ADdem$ADdem231_ADdem323,	p_ADdem$ADdem231_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem212_ADdem111,	p_ADdem$ADdem212_ADdem121,	p_ADdem$ADdem212_ADdem131,	p_ADdem$ADdem212_ADdem112,	p_ADdem$ADdem212_ADdem122,	p_ADdem$ADdem212_ADdem132,	p_ADdem$ADdem212_ADdem113,	p_ADdem$ADdem212_ADdem123,	p_ADdem$ADdem212_ADdem133,	p_ADdem$ADdem212_ADdem211,	p_ADdem$ADdem212_ADdem221,	p_ADdem$ADdem212_ADdem231,	C,	p_ADdem$ADdem212_ADdem222,	p_ADdem$ADdem212_ADdem232,	p_ADdem$ADdem212_ADdem213,	p_ADdem$ADdem212_ADdem223,	p_ADdem$ADdem212_ADdem233,	p_ADdem$ADdem212_ADdem311,	p_ADdem$ADdem212_ADdem321,	p_ADdem$ADdem212_ADdem331,	p_ADdem$ADdem212_ADdem312,	p_ADdem$ADdem212_ADdem322,	p_ADdem$ADdem212_ADdem332,	p_ADdem$ADdem212_ADdem313,	p_ADdem$ADdem212_ADdem323,	p_ADdem$ADdem212_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem222_ADdem111,	p_ADdem$ADdem222_ADdem121,	p_ADdem$ADdem222_ADdem131,	p_ADdem$ADdem222_ADdem112,	p_ADdem$ADdem222_ADdem122,	p_ADdem$ADdem222_ADdem132,	p_ADdem$ADdem222_ADdem113,	p_ADdem$ADdem222_ADdem123,	p_ADdem$ADdem222_ADdem133,	p_ADdem$ADdem222_ADdem211,	p_ADdem$ADdem222_ADdem221,	p_ADdem$ADdem222_ADdem231,	p_ADdem$ADdem222_ADdem212,	C,	p_ADdem$ADdem222_ADdem232,	p_ADdem$ADdem222_ADdem213,	p_ADdem$ADdem222_ADdem223,	p_ADdem$ADdem222_ADdem233,	p_ADdem$ADdem222_ADdem311,	p_ADdem$ADdem222_ADdem321,	p_ADdem$ADdem222_ADdem331,	p_ADdem$ADdem222_ADdem312,	p_ADdem$ADdem222_ADdem322,	p_ADdem$ADdem222_ADdem332,	p_ADdem$ADdem222_ADdem313,	p_ADdem$ADdem222_ADdem323,	p_ADdem$ADdem222_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem232_ADdem111,	p_ADdem$ADdem232_ADdem121,	p_ADdem$ADdem232_ADdem131,	p_ADdem$ADdem232_ADdem112,	p_ADdem$ADdem232_ADdem122,	p_ADdem$ADdem232_ADdem132,	p_ADdem$ADdem232_ADdem113,	p_ADdem$ADdem232_ADdem123,	p_ADdem$ADdem232_ADdem133,	p_ADdem$ADdem232_ADdem211,	p_ADdem$ADdem232_ADdem221,	p_ADdem$ADdem232_ADdem231,	p_ADdem$ADdem232_ADdem212,	p_ADdem$ADdem232_ADdem222,	C,	p_ADdem$ADdem232_ADdem213,	p_ADdem$ADdem232_ADdem223,	p_ADdem$ADdem232_ADdem233,	p_ADdem$ADdem232_ADdem311,	p_ADdem$ADdem232_ADdem321,	p_ADdem$ADdem232_ADdem331,	p_ADdem$ADdem232_ADdem312,	p_ADdem$ADdem232_ADdem322,	p_ADdem$ADdem232_ADdem332,	p_ADdem$ADdem232_ADdem313,	p_ADdem$ADdem232_ADdem323,	p_ADdem$ADdem232_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem213_ADdem111,	p_ADdem$ADdem213_ADdem121,	p_ADdem$ADdem213_ADdem131,	p_ADdem$ADdem213_ADdem112,	p_ADdem$ADdem213_ADdem122,	p_ADdem$ADdem213_ADdem132,	p_ADdem$ADdem213_ADdem113,	p_ADdem$ADdem213_ADdem123,	p_ADdem$ADdem213_ADdem133,	p_ADdem$ADdem213_ADdem211,	p_ADdem$ADdem213_ADdem221,	p_ADdem$ADdem213_ADdem231,	p_ADdem$ADdem213_ADdem212,	p_ADdem$ADdem213_ADdem222,	p_ADdem$ADdem213_ADdem232,	C,	p_ADdem$ADdem213_ADdem223,	p_ADdem$ADdem213_ADdem233,	p_ADdem$ADdem213_ADdem311,	p_ADdem$ADdem213_ADdem321,	p_ADdem$ADdem213_ADdem331,	p_ADdem$ADdem213_ADdem312,	p_ADdem$ADdem213_ADdem322,	p_ADdem$ADdem213_ADdem332,	p_ADdem$ADdem213_ADdem313,	p_ADdem$ADdem213_ADdem323,	p_ADdem$ADdem213_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem223_ADdem111,	p_ADdem$ADdem223_ADdem121,	p_ADdem$ADdem223_ADdem131,	p_ADdem$ADdem223_ADdem112,	p_ADdem$ADdem223_ADdem122,	p_ADdem$ADdem223_ADdem132,	p_ADdem$ADdem223_ADdem113,	p_ADdem$ADdem223_ADdem123,	p_ADdem$ADdem223_ADdem133,	p_ADdem$ADdem223_ADdem211,	p_ADdem$ADdem223_ADdem221,	p_ADdem$ADdem223_ADdem231,	p_ADdem$ADdem223_ADdem212,	p_ADdem$ADdem223_ADdem222,	p_ADdem$ADdem223_ADdem232,	p_ADdem$ADdem223_ADdem213,	C,	p_ADdem$ADdem223_ADdem233,	p_ADdem$ADdem223_ADdem311,	p_ADdem$ADdem223_ADdem321,	p_ADdem$ADdem223_ADdem331,	p_ADdem$ADdem223_ADdem312,	p_ADdem$ADdem223_ADdem322,	p_ADdem$ADdem223_ADdem332,	p_ADdem$ADdem223_ADdem313,	p_ADdem$ADdem223_ADdem323,	p_ADdem$ADdem223_ADdem333,	mrt_moderate,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem233_ADdem111,	p_ADdem$ADdem233_ADdem121,	p_ADdem$ADdem233_ADdem131,	p_ADdem$ADdem233_ADdem112,	p_ADdem$ADdem233_ADdem122,	p_ADdem$ADdem233_ADdem132,	p_ADdem$ADdem233_ADdem113,	p_ADdem$ADdem233_ADdem123,	p_ADdem$ADdem233_ADdem133,	p_ADdem$ADdem233_ADdem211,	p_ADdem$ADdem233_ADdem221,	p_ADdem$ADdem233_ADdem231,	p_ADdem$ADdem233_ADdem212,	p_ADdem$ADdem233_ADdem222,	p_ADdem$ADdem233_ADdem232,	p_ADdem$ADdem233_ADdem213,	p_ADdem$ADdem233_ADdem223,	C,	p_ADdem$ADdem233_ADdem311,	p_ADdem$ADdem233_ADdem321,	p_ADdem$ADdem233_ADdem331,	p_ADdem$ADdem233_ADdem312,	p_ADdem$ADdem233_ADdem322,	p_ADdem$ADdem233_ADdem332,	p_ADdem$ADdem233_ADdem313,	p_ADdem$ADdem233_ADdem323,	p_ADdem$ADdem233_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem311_ADdem111,	p_ADdem$ADdem311_ADdem121,	p_ADdem$ADdem311_ADdem131,	p_ADdem$ADdem311_ADdem112,	p_ADdem$ADdem311_ADdem122,	p_ADdem$ADdem311_ADdem132,	p_ADdem$ADdem311_ADdem113,	p_ADdem$ADdem311_ADdem123,	p_ADdem$ADdem311_ADdem133,	p_ADdem$ADdem311_ADdem211,	p_ADdem$ADdem311_ADdem221,	p_ADdem$ADdem311_ADdem231,	p_ADdem$ADdem311_ADdem212,	p_ADdem$ADdem311_ADdem222,	p_ADdem$ADdem311_ADdem232,	p_ADdem$ADdem311_ADdem213,	p_ADdem$ADdem311_ADdem223,	p_ADdem$ADdem311_ADdem233,	C,	p_ADdem$ADdem311_ADdem321,	p_ADdem$ADdem311_ADdem331,	p_ADdem$ADdem311_ADdem312,	p_ADdem$ADdem311_ADdem322,	p_ADdem$ADdem311_ADdem332,	p_ADdem$ADdem311_ADdem313,	p_ADdem$ADdem311_ADdem323,	p_ADdem$ADdem311_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem321_ADdem111,	p_ADdem$ADdem321_ADdem121,	p_ADdem$ADdem321_ADdem131,	p_ADdem$ADdem321_ADdem112,	p_ADdem$ADdem321_ADdem122,	p_ADdem$ADdem321_ADdem132,	p_ADdem$ADdem321_ADdem113,	p_ADdem$ADdem321_ADdem123,	p_ADdem$ADdem321_ADdem133,	p_ADdem$ADdem321_ADdem211,	p_ADdem$ADdem321_ADdem221,	p_ADdem$ADdem321_ADdem231,	p_ADdem$ADdem321_ADdem212,	p_ADdem$ADdem321_ADdem222,	p_ADdem$ADdem321_ADdem232,	p_ADdem$ADdem321_ADdem213,	p_ADdem$ADdem321_ADdem223,	p_ADdem$ADdem321_ADdem233,	p_ADdem$ADdem321_ADdem311,	C,	p_ADdem$ADdem321_ADdem331,	p_ADdem$ADdem321_ADdem312,	p_ADdem$ADdem321_ADdem322,	p_ADdem$ADdem321_ADdem332,	p_ADdem$ADdem321_ADdem313,	p_ADdem$ADdem321_ADdem323,	p_ADdem$ADdem321_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem331_ADdem111,	p_ADdem$ADdem331_ADdem121,	p_ADdem$ADdem331_ADdem131,	p_ADdem$ADdem331_ADdem112,	p_ADdem$ADdem331_ADdem122,	p_ADdem$ADdem331_ADdem132,	p_ADdem$ADdem331_ADdem113,	p_ADdem$ADdem331_ADdem123,	p_ADdem$ADdem331_ADdem133,	p_ADdem$ADdem331_ADdem211,	p_ADdem$ADdem331_ADdem221,	p_ADdem$ADdem331_ADdem231,	p_ADdem$ADdem331_ADdem212,	p_ADdem$ADdem331_ADdem222,	p_ADdem$ADdem331_ADdem232,	p_ADdem$ADdem331_ADdem213,	p_ADdem$ADdem331_ADdem223,	p_ADdem$ADdem331_ADdem233,	p_ADdem$ADdem331_ADdem311,	p_ADdem$ADdem331_ADdem321,	C,	p_ADdem$ADdem331_ADdem312,	p_ADdem$ADdem331_ADdem322,	p_ADdem$ADdem331_ADdem332,	p_ADdem$ADdem331_ADdem313,	p_ADdem$ADdem331_ADdem323,	p_ADdem$ADdem331_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem312_ADdem111,	p_ADdem$ADdem312_ADdem121,	p_ADdem$ADdem312_ADdem131,	p_ADdem$ADdem312_ADdem112,	p_ADdem$ADdem312_ADdem122,	p_ADdem$ADdem312_ADdem132,	p_ADdem$ADdem312_ADdem113,	p_ADdem$ADdem312_ADdem123,	p_ADdem$ADdem312_ADdem133,	p_ADdem$ADdem312_ADdem211,	p_ADdem$ADdem312_ADdem221,	p_ADdem$ADdem312_ADdem231,	p_ADdem$ADdem312_ADdem212,	p_ADdem$ADdem312_ADdem222,	p_ADdem$ADdem312_ADdem232,	p_ADdem$ADdem312_ADdem213,	p_ADdem$ADdem312_ADdem223,	p_ADdem$ADdem312_ADdem233,	p_ADdem$ADdem312_ADdem311,	p_ADdem$ADdem312_ADdem321,	p_ADdem$ADdem312_ADdem331,	C,	p_ADdem$ADdem312_ADdem322,	p_ADdem$ADdem312_ADdem332,	p_ADdem$ADdem312_ADdem313,	p_ADdem$ADdem312_ADdem323,	p_ADdem$ADdem312_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem322_ADdem111,	p_ADdem$ADdem322_ADdem121,	p_ADdem$ADdem322_ADdem131,	p_ADdem$ADdem322_ADdem112,	p_ADdem$ADdem322_ADdem122,	p_ADdem$ADdem322_ADdem132,	p_ADdem$ADdem322_ADdem113,	p_ADdem$ADdem322_ADdem123,	p_ADdem$ADdem322_ADdem133,	p_ADdem$ADdem322_ADdem211,	p_ADdem$ADdem322_ADdem221,	p_ADdem$ADdem322_ADdem231,	p_ADdem$ADdem322_ADdem212,	p_ADdem$ADdem322_ADdem222,	p_ADdem$ADdem322_ADdem232,	p_ADdem$ADdem322_ADdem213,	p_ADdem$ADdem322_ADdem223,	p_ADdem$ADdem322_ADdem233,	p_ADdem$ADdem322_ADdem311,	p_ADdem$ADdem322_ADdem321,	p_ADdem$ADdem322_ADdem331,	p_ADdem$ADdem322_ADdem312,	C,	p_ADdem$ADdem322_ADdem332,	p_ADdem$ADdem322_ADdem313,	p_ADdem$ADdem322_ADdem323,	p_ADdem$ADdem322_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem332_ADdem111,	p_ADdem$ADdem332_ADdem121,	p_ADdem$ADdem332_ADdem131,	p_ADdem$ADdem332_ADdem112,	p_ADdem$ADdem332_ADdem122,	p_ADdem$ADdem332_ADdem132,	p_ADdem$ADdem332_ADdem113,	p_ADdem$ADdem332_ADdem123,	p_ADdem$ADdem332_ADdem133,	p_ADdem$ADdem332_ADdem211,	p_ADdem$ADdem332_ADdem221,	p_ADdem$ADdem332_ADdem231,	p_ADdem$ADdem332_ADdem212,	p_ADdem$ADdem332_ADdem222,	p_ADdem$ADdem332_ADdem232,	p_ADdem$ADdem332_ADdem213,	p_ADdem$ADdem332_ADdem223,	p_ADdem$ADdem332_ADdem233,	p_ADdem$ADdem332_ADdem311,	p_ADdem$ADdem332_ADdem321,	p_ADdem$ADdem332_ADdem331,	p_ADdem$ADdem332_ADdem312,	p_ADdem$ADdem332_ADdem322,	C,	p_ADdem$ADdem332_ADdem313,	p_ADdem$ADdem332_ADdem323,	p_ADdem$ADdem332_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem313_ADdem111,	p_ADdem$ADdem313_ADdem121,	p_ADdem$ADdem313_ADdem131,	p_ADdem$ADdem313_ADdem112,	p_ADdem$ADdem313_ADdem122,	p_ADdem$ADdem313_ADdem132,	p_ADdem$ADdem313_ADdem113,	p_ADdem$ADdem313_ADdem123,	p_ADdem$ADdem313_ADdem133,	p_ADdem$ADdem313_ADdem211,	p_ADdem$ADdem313_ADdem221,	p_ADdem$ADdem313_ADdem231,	p_ADdem$ADdem313_ADdem212,	p_ADdem$ADdem313_ADdem222,	p_ADdem$ADdem313_ADdem232,	p_ADdem$ADdem313_ADdem213,	p_ADdem$ADdem313_ADdem223,	p_ADdem$ADdem313_ADdem233,	p_ADdem$ADdem313_ADdem311,	p_ADdem$ADdem313_ADdem321,	p_ADdem$ADdem313_ADdem331,	p_ADdem$ADdem313_ADdem312,	p_ADdem$ADdem313_ADdem322,	p_ADdem$ADdem313_ADdem332,	C,	p_ADdem$ADdem313_ADdem323,	p_ADdem$ADdem313_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem323_ADdem111,	p_ADdem$ADdem323_ADdem121,	p_ADdem$ADdem323_ADdem131,	p_ADdem$ADdem323_ADdem112,	p_ADdem$ADdem323_ADdem122,	p_ADdem$ADdem323_ADdem132,	p_ADdem$ADdem323_ADdem113,	p_ADdem$ADdem323_ADdem123,	p_ADdem$ADdem323_ADdem133,	p_ADdem$ADdem323_ADdem211,	p_ADdem$ADdem323_ADdem221,	p_ADdem$ADdem323_ADdem231,	p_ADdem$ADdem323_ADdem212,	p_ADdem$ADdem323_ADdem222,	p_ADdem$ADdem323_ADdem232,	p_ADdem$ADdem323_ADdem213,	p_ADdem$ADdem323_ADdem223,	p_ADdem$ADdem323_ADdem233,	p_ADdem$ADdem323_ADdem311,	p_ADdem$ADdem323_ADdem321,	p_ADdem$ADdem323_ADdem331,	p_ADdem$ADdem323_ADdem312,	p_ADdem$ADdem323_ADdem322,	p_ADdem$ADdem323_ADdem332,	p_ADdem$ADdem323_ADdem313,	C,	p_ADdem$ADdem323_ADdem333,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	p_ADdem$ADdem333_ADdem111,	p_ADdem$ADdem333_ADdem121,	p_ADdem$ADdem333_ADdem131,	p_ADdem$ADdem333_ADdem112,	p_ADdem$ADdem333_ADdem122,	p_ADdem$ADdem333_ADdem132,	p_ADdem$ADdem333_ADdem113,	p_ADdem$ADdem333_ADdem123,	p_ADdem$ADdem333_ADdem133,	p_ADdem$ADdem333_ADdem211,	p_ADdem$ADdem333_ADdem221,	p_ADdem$ADdem333_ADdem231,	p_ADdem$ADdem333_ADdem212,	p_ADdem$ADdem333_ADdem222,	p_ADdem$ADdem333_ADdem232,	p_ADdem$ADdem333_ADdem213,	p_ADdem$ADdem333_ADdem223,	p_ADdem$ADdem333_ADdem233,	p_ADdem$ADdem333_ADdem311,	p_ADdem$ADdem333_ADdem321,	p_ADdem$ADdem333_ADdem331,	p_ADdem$ADdem333_ADdem312,	p_ADdem$ADdem333_ADdem322,	p_ADdem$ADdem333_ADdem332,	p_ADdem$ADdem333_ADdem313,	p_ADdem$ADdem333_ADdem323,	C,	mrt_severe,
    0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1
    ),
  
  utility=u,
  
  costs=list(
    care = c_care,
    treatment = c_tx
  )
)  
}, times = 1:follow_up)





##############
# SIMULATION #
##############

# Constructing the economic model
mod_def <- define_model(tparams_def = tparams_def,
                        rng_def = rng_def,
                        params = params)

#evalmod <- eval_model(mod_def, input_data)

econmod <- create_CohortDtstm(mod_def, input_data)



# Simulating outcomes
## Health state probabilities
econmod$sim_stateprobs(n_cycles = follow_up)


## Costs and QALYs
econmod$sim_qalys(dr = .03, integrate_method = "trapz",lys = TRUE)
econmod$sim_costs(dr = .03, integrate_method = "trapz")





##############################################
# DISTRIBUTIONAL COST-EFFECTIVENESS ANALYSIS #
##############################################

ce_sim_grp <- econmod$summarize(by_grp = TRUE)
ce_sim_dcea<-left_join(ce_sim_grp$qalys,filter(ce_sim_grp$costs,category == "total"),by=c("sample", "strategy_id","grp_id"))
patient_wt_dcea <- rep(patient_wt, PSA*max(strategies$strategy_id)) 
ce_sim_dcea <-bind_cols(ce_sim_dcea,patient_wt_dcea)
ce_sim_dcea<-rename(ce_sim_dcea,"weight"="...9")
ce_sim_dcea[, grp_id_new := 
                            factor(
                              1 * (grp_id == 1  | grp_id == 2 | grp_id == 3 | grp_id == 13 | grp_id == 14 | grp_id == 15 ) +
                              2 * (grp_id == 4  | grp_id == 5 | grp_id == 6 | grp_id == 16 | grp_id == 17 | grp_id == 18 ) +
                              3 * (grp_id == 7  | grp_id == 8 | grp_id == 9 | grp_id == 19 | grp_id == 20 | grp_id == 21 ) +
                              4 * (grp_id == 10  | grp_id == 11 | grp_id == 12 | grp_id == 22 | grp_id == 23 | grp_id == 24 )    
                            )
                          ]

ce_sim_dcea_overall <- ce_sim_dcea[, lapply(.SD, stats::weighted.mean, w = weight),
                           by = c("sample","strategy_id"), .SDcols = c("qalys","costs")]
ce_sim_dcea_overall[, grp_id_new :=5]
ce_sim_dcea_overall[, weight_new := 1]
ce_sim_dcea_subgroups <- ce_sim_dcea[, lapply(.SD, stats::weighted.mean, w = weight),
                       by = c("grp_id_new","sample","strategy_id"), .SDcols = c("qalys","costs")]
ce_sim_dcea_subgroups_weights <- ce_sim_dcea[, .(weight_new = sum(weight)),
                                      by = c("grp_id_new","sample","strategy_id")]
ce_sim_dcea_subgroups <-left_join(ce_sim_dcea_subgroups,ce_sim_dcea_subgroups_weights,by = c("grp_id_new","sample","strategy_id"))

ce_sim_dcea_long <- rbind(ce_sim_dcea_subgroups,ce_sim_dcea_overall)

ce_sim_dcea_wide = ce_sim_dcea_long %>% 
  tidyr::pivot_wider(names_from = grp_id_new, values_from = c(qalys, costs, weight_new)) %>%
  dplyr::rename("qalys_overall"="qalys_5","costs_overall"="costs_5") %>%
  as.data.table()

ce_sim_dcea_wide[, NHB_1 := qalys_1 - costs_overall/wtp_dcea] #opportunity costs equally divided
ce_sim_dcea_wide[, NHB_2 := qalys_2 - costs_overall/wtp_dcea]
ce_sim_dcea_wide[, NHB_3 := qalys_3 - costs_overall/wtp_dcea]
ce_sim_dcea_wide[, NHB_4 := qalys_4 - costs_overall/wtp_dcea]
ce_sim_dcea_wide[, NHB_overall := qalys_overall - costs_overall/wtp_dcea]

ce_sim_dcea_wide[, NHB_ede := (((((NHB_1/NHB_overall)^(1-atkinson_dcea))*weight_new_1) + 
                                  (((NHB_2/NHB_overall)^(1-atkinson_dcea))*weight_new_2) +
                                  (((NHB_3/NHB_overall)^(1-atkinson_dcea))*weight_new_3) +
                                  (((NHB_4/NHB_overall)^(1-atkinson_dcea))*weight_new_4))^(1/(1-atkinson_dcea)))*NHB_overall
                 ]
ce_sim_dcea_wide[, inequality :=  1-(NHB_ede/NHB_overall) ] 

ce_sim_dcea_wide[, NHB_ede_kolm := NHB_overall-((1/kolm_dcea)*log(
                                                                   weight_new_1*exp(kolm_dcea*(NHB_overall-NHB_1)) +
                                                                   weight_new_2*exp(kolm_dcea*(NHB_overall-NHB_2)) +
                                                                   weight_new_3*exp(kolm_dcea*(NHB_overall-NHB_3)) + 
                                                                   weight_new_4*exp(kolm_dcea*(NHB_overall-NHB_4)) ))
                         ]                               
                                              
ce_sim_dcea_wide[, inequality_kolm :=  NHB_overall-NHB_ede_kolm ] 
                 

QALY_table <- ce_sim_dcea_wide[,.(
  qalys_overall = mean(qalys_overall),
  "95%CI-low" = quantile(qalys_overall, .025),
  "95%CI-high" = quantile(qalys_overall, .975),
  qalys_1 = mean(qalys_1),
  "95%CI-low" = quantile(qalys_1, .025),
  "95%CI-high" = quantile(qalys_1, .975),
  qalys_2 = mean(qalys_2),
  "95%CI-low" = quantile(qalys_2, .025),
  "95%CI-high" = quantile(qalys_2, .975),
  qalys_3 = mean(qalys_3),
  "95%CI-low" = quantile(qalys_3, .025),
  "95%CI-high" = quantile(qalys_3, .975),
  qalys_4 = mean(qalys_4),
  "95%CI-low" = quantile(qalys_4, .025),
  "95%CI-high" = quantile(qalys_4, .975)
),
by = c("strategy_id")]



DCEA_table <- ce_sim_dcea_wide[,.(
  NHB_overall = mean(NHB_overall),
  "95%CI-low" = quantile(NHB_overall, .025),
  "95%CI-high" = quantile(NHB_overall, .975),
  NHB_1 = mean(NHB_1),
  "95%CI-low" = quantile(NHB_1, .025),
  "95%CI-high" = quantile(NHB_1, .975),
  NHB_2 = mean(NHB_2),
  "95%CI-low" = quantile(NHB_2, .025),
  "95%CI-high" = quantile(NHB_2, .975),
  NHB_3 = mean(NHB_3),
  "95%CI-low" = quantile(NHB_3, .025),
  "95%CI-high" = quantile(NHB_3, .975),
  NHB_4 = mean(NHB_4),
  "95%CI-low" = quantile(NHB_4, .025),
  "95%CI-high" = quantile(NHB_4, .975)
),
by = c("strategy_id")]




ce_sim_dcea_pairwise = ce_sim_dcea_wide[, c("sample","strategy_id","NHB_overall", "NHB_ede", "inequality", "NHB_ede_kolm", "inequality_kolm")] %>% 
  tidyr::pivot_wider(names_from = strategy_id, values_from = c(NHB_overall, NHB_ede, inequality, NHB_ede_kolm, inequality_kolm)) %>%
  as.data.table()

ce_sim_dcea_pairwise[, dNHB_overall := NHB_overall_2 - NHB_overall_1 ]
ce_sim_dcea_pairwise[, hei := 0-(inequality_2 - inequality_1) ]
ce_sim_dcea_pairwise[, dNHB_ede := NHB_ede_2 - NHB_ede_1 ]
ce_sim_dcea_pairwise[, hei_kolm := 0-(inequality_kolm_2 - inequality_kolm_1) ]
ce_sim_dcea_pairwise[, dNHB_ede_kolm := NHB_ede_kolm_2 - NHB_ede_kolm_1 ]


ce_sim_dcea_pairwise <- ce_sim_dcea_pairwise[,c("sample","dNHB_overall","hei","dNHB_ede","hei_kolm","dNHB_ede_kolm")]


DCEA_pairwise_table <- ce_sim_dcea_pairwise[,.(
  dNHB_overall = mean(dNHB_overall),
  "95%CI-low" = quantile(dNHB_overall, .025),
  "95%CI-high" = quantile(dNHB_overall, .975),
  hei = mean(hei, na.rm=TRUE),
  "95%CI-low" = quantile(hei, .025, na.rm=TRUE),
  "95%CI-high" = quantile(hei, .975, na.rm=TRUE),
  dNHB_ede = mean(dNHB_ede, na.rm=TRUE),
  "95%CI-low" = quantile(dNHB_ede, .025, na.rm=TRUE),
  "95%CI-high" = quantile(dNHB_ede, .975, na.rm=TRUE),
  hei_kolm = mean(hei_kolm, na.rm=TRUE),
  "95%CI-low" = quantile(hei_kolm, .025, na.rm=TRUE),
  "95%CI-high" = quantile(hei_kolm, .975, na.rm=TRUE),
  dNHB_ede_kolm = mean(dNHB_ede_kolm, na.rm=TRUE),
  "95%CI-low" = quantile(dNHB_ede_kolm, .025, na.rm=TRUE),
  "95%CI-high" = quantile(dNHB_ede_kolm, .975, na.rm=TRUE)
)]


## HDE_ede by inequality parameter
### Atkinson
atkinson_range <- c(0,0.9,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
n_atkinson <- length(atkinson_range)
out <- array(NA, dim = c(PSA, 3, n_atkinson))
dimnames(out) <- list(NULL, 
                      c("sample", "dNHB_ede","atkinson_value"), 
                      NULL)
i <- 1
for (i in 1:n_atkinson){
  NHB_ede_atkinson <- ce_sim_dcea_wide[,c("sample","strategy_id","NHB_1","NHB_2","NHB_3","NHB_4","NHB_overall","weight_new_1","weight_new_2","weight_new_3","weight_new_4")]
  atkinson_value <- atkinson_range[i]
 
  NHB_ede_atkinson <- NHB_ede_atkinson[, NHB_ede := (((((NHB_1/NHB_overall)^(1-atkinson_value))*weight_new_1) + 
                                                        (((NHB_2/NHB_overall)^(1-atkinson_value))*weight_new_2) +
                                                        (((NHB_3/NHB_overall)^(1-atkinson_value))*weight_new_3) +
                                                        (((NHB_4/NHB_overall)^(1-atkinson_value))*weight_new_4))^(1/(1-atkinson_value)))*NHB_overall
                                       ]
  NHB_ede_atkinson <- NHB_ede_atkinson[,c("sample","strategy_id","NHB_ede")]
  NHB_ede_atkinson[, atkinson_value := atkinson_value]
  
  NHB_ede_atkinson_pairwise <- NHB_ede_atkinson %>% 
    tidyr::pivot_wider(names_from = strategy_id, values_from = NHB_ede) %>%
    dplyr::rename("NHB_ede_1"="1","NHB_ede_2"="2") %>%
    as.data.table()
  NHB_ede_atkinson_pairwise[, dNHB_ede := NHB_ede_2 - NHB_ede_1 ]
  NHB_ede_atkinson_pairwise <- NHB_ede_atkinson_pairwise[,c("sample","dNHB_ede","atkinson_value")]
  out[,,i] <- cbind(NHB_ede_atkinson_pairwise$sample, NHB_ede_atkinson_pairwise$dNHB_ede,NHB_ede_atkinson_pairwise$atkinson_value)
  
  i <- i + 1
}

#### rbind an array of matrices into a single matrix
rbind_array <- function(x){
  n_rows <- dim(x)[3] * dim(x)[1]
  x_mat <- matrix(c(aperm(x, perm = c(2, 1, 3))),
                  nrow = n_rows, byrow = TRUE)
  colnames(x_mat) <- dimnames(x)[[2]]
  return(x_mat)
}
dNHB_ede_atkinson <- as.data.table(rbind_array(out))

dNHB_ede_table <- dNHB_ede_atkinson[,
                                    .(dNHB_ede = mean(dNHB_ede,na.rm=TRUE),
                                      "95%CI-low" = quantile(dNHB_ede, .025, na.rm = TRUE),
                                      "95%CI-high" = quantile(dNHB_ede, .975,na.rm = TRUE)),
                                    by = atkinson_value]

### Kolm
kolm_range <- c(0.0001,0.005,0.01,seq(0.025, 0.30, 0.025))
n_kolm <- length(kolm_range)
out <- array(NA, dim = c(PSA, 3, n_kolm))
dimnames(out) <- list(NULL, 
                      c("sample", "dNHB_ede_kolm","kolm_value"), 
                      NULL)
i <- 1
for (i in 1:n_kolm){
  NHB_ede_kolm <- ce_sim_dcea_wide[,c("sample","strategy_id","NHB_1","NHB_2","NHB_3","NHB_4","NHB_overall","weight_new_1","weight_new_2","weight_new_3","weight_new_4")]
  kolm_value <- kolm_range[i]
  
  NHB_ede_kolm <- NHB_ede_kolm[, NHB_ede_kolm := 
                                 NHB_overall-((1/kolm_value)*log(
                                     weight_new_1*exp(kolm_value*(NHB_overall-NHB_1)) +
                                     weight_new_2*exp(kolm_value*(NHB_overall-NHB_2)) +
                                     weight_new_3*exp(kolm_value*(NHB_overall-NHB_3)) + 
                                     weight_new_4*exp(kolm_value*(NHB_overall-NHB_4)) ))
                               ]
  NHB_ede_kolm <- NHB_ede_kolm[,c("sample","strategy_id","NHB_ede_kolm")]
  NHB_ede_kolm[, kolm_value := kolm_value]
  
  NHB_ede_kolm_pairwise <- NHB_ede_kolm %>% 
    tidyr::pivot_wider(names_from = strategy_id, values_from = NHB_ede_kolm) %>%
    dplyr::rename("NHB_ede_kolm_1"="1","NHB_ede_kolm_2"="2") %>%
    as.data.table()
  NHB_ede_kolm_pairwise[, dNHB_ede_kolm := NHB_ede_kolm_2 - NHB_ede_kolm_1 ]
  NHB_ede_kolm_pairwise <- NHB_ede_kolm_pairwise[,c("sample","dNHB_ede_kolm","kolm_value")]
  out[,,i] <- cbind(NHB_ede_kolm_pairwise$sample, NHB_ede_kolm_pairwise$dNHB_ede, NHB_ede_kolm_pairwise$dNHB_ede_alt, NHB_ede_kolm_pairwise$kolm_value)
  
  i <- i + 1
}

#### rbind an array of matrices into a single matrix
rbind_array <- function(x){
  n_rows <- dim(x)[3] * dim(x)[1]
  x_mat <- matrix(c(aperm(x, perm = c(2, 1, 3))),
                  nrow = n_rows, byrow = TRUE)
  colnames(x_mat) <- dimnames(x)[[2]]
  return(x_mat)
}
dNHB_ede_kolm <- as.data.table(rbind_array(out))

dNHB_ede_kolm_table <- dNHB_ede_kolm[,
                                     .(dNHB_ede_kolm = mean(dNHB_ede_kolm,na.rm=TRUE),
                                       "95%CI-low" = quantile(dNHB_ede_kolm, .025, na.rm = TRUE),
                                       "95%CI-high" = quantile(dNHB_ede_kolm, .975,na.rm = TRUE)),
                                     by = kolm_value]



###########
# OUTPUTS #
###########

output <- list(
         QALY_table, 
         DCEA_table,
         DCEA_pairwise_table,
         dNHB_ede_table,dNHB_ede_kolm_table)

return(output)
}




################
# RUN ANALYSIS #
################

# $10,000
ad_model(  
  tx_cost = 10000,
  PSA = 125)  # depending on computer memory multiple sets of PSA simulations may need to be performed to obtain 1000 simulations


# $56,000
ad_model(  
  tx_cost = 56000,
  PSA = 125)

