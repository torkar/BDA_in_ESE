library(dagitty)
library(dplyr)
library(HDInterval)
library(bayesplot)
color_scheme_set("darkgray")
library(latex2exp)
library(ggplot2)
library(ggthemes)
library(brms)
library(mice)
library(gridExtra)
library(sjPlot)

options(mc.cores = parallel::detectCores())

################################################################################
# Build DAGs! AFP is a mediator for several response variables!
# It is routine to worry about mistaken inferences that arise from omitting 
# predictor variables. Such mistakes are often called omitted variable bias. It 
# is much less routine to worry about mistaken inferences arising from including 
# variables that are consequences of other variables, i.e., post-treatment bias.

# AFP is built from IFPUG in our case, but IFPUG uses input, output, 
# enquiry, etc. It looks like the graph would then be a DAG graph representing 
# The Pipe confound. Conditioning on AFP induces d-separation.
# In short, AFP mediates association between the other predictors and Effort.
# Conditioning on Z removes dependency between our predictors and Effort, e.g.,
# Input _||_ Effort | AFP
dag <- dagitty( "dag {
                      Input -> AFP
                      Output -> AFP
                      Enquiry -> AFP
                      File -> AFP
                      Interface -> AFP
                      Added -> AFP
                      Changed -> AFP
                      AFP -> Effort
                      }")

coordinates(dag) <- list(x=c(Effort=2, AFP = 1, Input = 0.5, Output = 0.25, 
                             Enquiry = 0.15, File = 0.1, Interface = 0.15, 
                             Added = 0.25,Changed = 0.5),
                         y=c(Effort=0, AFP = 0, Input = -1.5, Output = -1, 
                             Enquiry = -0.5, File = 0, Interface = 0.5, 
                             Added = 1, Changed = 1.5))

impliedConditionalIndependencies(dag)
dseparated(dag, "Input", "Effort") # should eval to FALSE
plot(dag)
# In short, AFP mediates association between the other predictors and Effort.
# Conditioning on AFP removes dependency between our predictors and Effort, 
# e.g., Input _||_ Effort | AFP
################################################################################
# IFPUG v10 data
d <- read.csv2("~/Documents/cth/Research Projects, Studies & Data/Bayesian/missing data/data.csv", sep=";", na.strings=c("","NA"))

d$Effort <- as.numeric(d$Summary.Work.Effort)
d$AFP <- as.numeric(d$Adjusted.Function.Points)
d$Input <- as.numeric(d$Input.count)
d$Output <- as.numeric(d$Output.count)
d$Enquiry <- as.numeric(d$Enquiry.count)
d$File <- as.numeric(d$File.count)
d$Interface <- as.numeric(d$Interface.count)
d$Added <- as.numeric(d$Added.count)
d$Changed <- as.numeric(d$Changed.count)
d$Deleted <- as.numeric(d$Deleted.count)
d$DQR <- as.factor(d$Data.Quality.Rating)


# Columns w/ DV/IV to keep according to original study. We've added a few here,
# DQR which we will use as a varying intercept and DQR_idx which are index 
# variables for DQR. In addition we have Count.Approach and FP.Standard to 
# distinguish between count approaches and IFPUG versions, respectively, when
# setting up our data
keepColumns <- as.vector(c("Effort", "Input","Output", "Enquiry", "File", 
                           "Interface", "Added", "Changed", "Deleted", "DQR", 
                           "Count.Approach", "FP.Standard"))

d <- d[ , keepColumns]

# And a set where we don't have Count.Approach and FP.Standard so we don't
# condition on irrelevant columns when removing rows w/ NAs.
keepColumns_o <- as.vector(c("Effort", "Input","Output", "Enquiry", "File", 
                             "Interface", "Added", "Changed", "Deleted", "DQR"))

# Now we create four datasets (DQR == A || A-D) but only w/ IFPUG >= 4, 
# i.e., according to ISBSG, IFPUG >= 4 can not be compared with IFPUG <4. 
# And we either keep NAs or remove them (_clean).
d_ifpug_AD <- d[d$FP.Standard == "IFPUG 4" | d$FP.Standard == "IFPUG 4.0" | 
               d$FP.Standard == "IFPUG 4.1" | d$FP.Standard == "IFPUG 4.2", ]
d_ifpug_AD <- d_ifpug_AD[ , keepColumns_o]

# _clean w/ removed NAs but contains zeros
d_ifpug_AD_clean <- d_ifpug_AD[complete.cases(d_ifpug_AD) , ]
d_ifpug_AD <- droplevels(d_ifpug_AD)
d_ifpug_AD_clean <- droplevels(d_ifpug_AD_clean)

# Only keep DQR=='A'
d_ifpug_A <- d_ifpug_AD[d_ifpug_AD$DQR == 'A' , ]

# _clean w/ removed NAs but contains zeros
d_ifpug_A_clean <- d_ifpug_A[complete.cases(d_ifpug_A) , ]
d_ifpug_A <- droplevels(d_ifpug_A)
d_ifpug_A_clean <- droplevels(d_ifpug_A_clean)

################################################################################
# Check to see so we haven't messed up the data sets...
stopifnot(dim(d)[1] == 4106)
stopifnot(dim(d_ifpug_A)[1] == 501)
stopifnot(dim(d_ifpug_A_clean)[1] == 211)
stopifnot(dim(d_ifpug_AD)[1] == 1689)
stopifnot(dim(d_ifpug_AD_clean)[1] == 490)
################################################################################

# Strong correlations between variables is generally speaking a bad thing when 
# building statistical models. The model will make good predictions, but 
# understandability will decrease since the impression will be that variables 
# that have strong correlation, do not seem to have much predictive power, while 
# in fact they have strong associations with the 
# outcome (McElreath, 2015).

# Multiple ways to handle this issue exist: manually check all 
# variables for correlations (e.g., using cor() in R), simply examine a pairs 
# plot where all combinations of variables and their correlations are visualized 
# (using alias() in R is also an option of course), or check if the matrix is a 
# full rank matrix.

# Check if we have a full rank matrix, i.e., no multi-collinearity
# https://www.cds.caltech.edu/~murray/wiki/What_is_matrix_rank_and_how_do_i_calculate_it%3F

# Thanks to Max Mantei for input on this
# https://discourse.mc-stan.org/t/blog-post-identifying-non-identifiability/4201/3  
# -1 for no intercept
X_full <- model.matrix(Effort ~ Input + Output + Enquiry + File + Interface +
                         Added + Changed + Deleted -1, data=d_ifpug_A)

# Compute the QR decomposition of the matrix
qrX_full <- qr(X_full)

# Returns the original matrix from which the object was constructed or the 
# components of the decomposition, abs and sum.
# This should evaluate to TRUE (otherwise check the "0.1" value), when compared
# to ncol minus the rank

ncol(X_full) - qrX_full$rank == sum(abs(diag(qr.R(qrX_full))) < 0.1)

# Check which columns are problematic
colnames(qr.R(qrX_full))[abs(diag(qr.R(qrX_full))) < 0.1]
# so... 'Deleted' should not be used as a predictor
# What we've done is calculate the QR decomposition \mathsf{Q} of the matrix 
# \mathsf{X}. Completethe matrix \mathsf{R} by binding zero-value rows beneath 
# the square uppertriangle. Extract the diagonal of the matrix.

################################################################################
# Take II
#
# Let's create new datasets where we don't use 'Deleted' and 'AFP'. 
# Perhaps more non-NAs will give us even more data to play with...
# Load the original data set
d <- read.csv2("~/bar/data.csv", sep=";", na.strings=c("","NA"))

# Same as before but we don't use AFP and Deleteds
d$Effort <- as.integer(d$Summary.Work.Effort)
d$Input <- as.integer(d$Input.count)
d$Output <- as.integer(d$Output.count)
d$Enquiry <- as.integer(d$Enquiry.count)
d$File <- as.integer(d$File.count)
d$Interface <- as.integer(d$Interface.count)
d$Added <- as.integer(d$Added.count)
d$Changed <- as.integer(d$Changed.count)
d$DQR <- as.factor(d$Data.Quality.Rating)

# Columns w/ DV/IV to keep according to original study (w/o AFP and Deleted).
# We have Count.Approach and FP.Standard to distinguish between count 
# approaches and IFPUG versions, respectively, when setting up our data
keepColumns <- as.vector(c("Effort", "Input","Output", "Enquiry", "File", 
                           "Interface", "Added", "Changed", "DQR", 
                           "Count.Approach", "FP.Standard"))

d <- d[ , keepColumns]

# And a set where we don't have Count.Approach and FP.Standard so we don't
# condition on irrelevant columns when removing rows w/ NAs.
keepColumns_o <- as.vector(c("Effort", "Input","Output", "Enquiry", "File", 
                             "Interface", "Added", "Changed", "DQR"))

# Now we create datasets (DQR == A | A-D) but only w/ IFPUG >= 4, 
# i.e., according to ISBSG, IFPUG >= 4 can't be compared with IFPUG <4. 
# And we either keep NAs or remove them (_clean).
d_ifpug_AD <- d[d$FP.Standard == "IFPUG 4" | d$FP.Standard == "IFPUG 4.0" | 
                  d$FP.Standard == "IFPUG 4.1" | d$FP.Standard == "IFPUG 4.2",]
d_ifpug_AD <- d_ifpug_AD[ , keepColumns_o]

# _clean, i.e., removed NAs but still contains zeros
d_ifpug_AD_clean <- d_ifpug_AD[complete.cases(d_ifpug_AD) , ]
d_ifpug_AD <- droplevels(d_ifpug_AD)
d_ifpug_AD_clean <- droplevels(d_ifpug_AD_clean)

# Only keep DQR=='A'
d_ifpug_A <- d_ifpug_AD[d_ifpug_AD$DQR == 'A' , ]
d_ifpug_A$DQR <- NULL # not needed now
d_ifpug_A <- droplevels(d_ifpug_A)

# _clean w/ removed NAs but contains zeros
d_ifpug_A_clean <- d_ifpug_A[complete.cases(d_ifpug_A) , ]
d_ifpug_A_clean <- droplevels(d_ifpug_A_clean)

d_ifpug_A_clean <- d_ifpug_A_clean %>% 
  mutate_at(c(2:8), funs(c(scale(.))))

d_ifpug_AD_clean <- d_ifpug_AD_clean %>% 
  mutate_at(c(2:8), funs(c(scale(.))))

# d_ifpug_A and d_ifpug_AD we scale later after imputing values! 
# Never scale before imputation!

################################################################################
# Check to see so we haven't messed up the data sets...
stopifnot(dim(d)[1] == 4106) # original dataset

stopifnot(dim(d_ifpug_AD)[1] == 1689) # w/ 
stopifnot(dim(d_ifpug_A)[1] == 501)

stopifnot(dim(d_ifpug_AD_clean)[1] == 494) # +4 rows compared to before
stopifnot(dim(d_ifpug_A_clean)[1] == 214) # +3 rows 

################################################################################
#
# Sensitivity analysis of priors
# W/o looking at the data, concerning Effort there are 9.7% projects 
# with >20 people in the team. If we calculate, roughly 1700 h/year for an 
# individual, having 60 people in a dev team sums up to ~100,000 hours. Have a 
# hard time seeing that the data has Effort>100,000 hours so let's assume that 
# is the approximate max for now.
# For good examples see Statistical Rethinking (McElreath, 2015)

# Set sample size
N <- 100

# Having N(5, 4) provides us a mean equal to
a <- rnorm(1e4, 5, 4)
lambda <- exp(a)
mean(lambda) # ~260,000 h, i.e., approx. 150 FTE/year hours in the project
max(lambda)
b1 <- rnorm(N, 0, 10) # very broad priors
b2 <- rnorm(N, 0, 10)
b3 <- rnorm(N, 0, 10)
b4 <- rnorm(N, 0, 10)
b5 <- rnorm(N, 0, 10)
b6 <- rnorm(N, 0, 10)
b7 <- rnorm(N, 0, 10)

plot(NULL, xlim=c(-2,2), ylim=c(0, 5000000), xlab="", ylab="", axes=F, pch=16, 
     type="b")

axis(1, at=seq(-2, 2, 1), label=c(-2, -1, 0, 1, 2), tick=F)
axis(2, at=c(0, 1e6, 3e6), label=c(0, expression(1 %*% 10^6), 
                                   expression(3 %*% 10^6)), 
     tick=F, las=2)
abline(h=1e6,lty=2)

for ( i in 1:N ) curve( exp( a[i] + b1[i]*x + b2[i]*x + b3[i]*x + b4[i]*x + 
                               b5[i]*x + b6[i]*x + b7[i]*x) , add=TRUE)

b1 <- rnorm(N, 0, 0.25)
b2 <- rnorm(N, 0, 0.25)
b3 <- rnorm(N, 0, 0.25)
b4 <- rnorm(N, 0, 0.25)
b5 <- rnorm(N, 0, 0.25)
b6 <- rnorm(N, 0, 0.25)
b7 <- rnorm(N, 0, 0.25)

plot( NULL , xlim=c(-2,2), ylim=c(0, 5000000), xlab="", ylab="", axes=F, pch=16, 
      type="b" )

axis(1, at=seq(-2, 2, 1), label=c(-2, -1, 0, 1, 2), tick=F)
axis(2, at=c(0, 1e6, 3e6), label=c(0, expression(1 %*% 10^6), 
                                   expression(3 %*% 10^6)), 
     tick=F, las=2)
abline(h=1e6,lty=2)

for ( i in 1:N ) curve( exp( a[i] + b1[i]*x + b2[i]*x + b3[i]*x + b4[i]*x + 
                               b5[i]*x + b6[i]*x + b7[i]*x) , add=TRUE)

################################################################################
# Missing data analysis
#pdf("md_pattern_A.pdf", width=14, height=21)
md.pattern(d_ifpug_A)
#dev.off()

foo <- as.data.frame(flux(d_ifpug_A))
set.seed(1)
ggplot(foo, aes(x = influx, y = outflux)) + 
  scale_x_continuous(limits = c(-0.01, 1.0)) +
  scale_y_continuous(limits = c(-0.01, 1.01)) +
  geom_abline(aes(intercept=1, slope=-1), size=0.1) +
  geom_jitter(alpha=0.5, size=5, width=0.01, height=0.01) +
  theme_tufte(ticks = FALSE) +
  xlab(TeX('$I_j$')) + ylab(TeX('$O_j$')) +
  annotate("text", x = 0.05, y = 0.37, label = "Effort", family="serif") +
  annotate("text", x = 0.20, y = 0.03, label = "Changed", family="serif")

################################################################################

# Need to use negative binomial (Gamma Poisson) since the variance is *much* 
# larger than the mean, however we are talking about "counting", i.e., Effort in 
# hours 0 -> theoretical \inf. Predictor variables adjust the shape of the 
# Gamma-Poisson, not the expected value of each observation, so Gamma-Poisson is 
# a nice alternative, a.k.a. a negative binomial.

# Set for reproducibility
SEED = 061215

# AD_clean
AD_clean <- brm(bf(Effort ~ Input + Output + Enquiry + File + Interface + 
                     Added + Changed + (1 | DQR)), # shape ~ DQR can also be used,
                data=d_ifpug_AD_clean, seed=SEED,
                family = negbinomial,
                prior = c(prior(normal(5, 4), class = Intercept),
                          prior(normal(0, 0.25), class = b),
                          prior(gamma(0.5, 0.5), class = shape),
                          prior(cauchy(0, 1), class=sd)),
                cores=4, chains=4, control=list(adapt_delta=0.99))

posterior <- as.matrix(AD_clean)
color_scheme_set("gray")
p4 <- mcmc_areas_ridges(posterior, pars = c("b_Input", "b_Output", "b_Enquiry", 
                                         "b_File", "b_Interface", "b_Added", 
                                         "b_Changed"),
                     prob=0.95, prob_outer=1) + 
  scale_x_continuous(limits=c(-0.9, 0.9)) +
  theme_tufte(base_size = 14) +
  theme(axis.text.y=element_blank(),
         axis.ticks.y=element_blank())

color_scheme_set("gray")
p_alphas_ad_clean <- mcmc_intervals(posterior, regex_pars = c("^r"),
                                    prob_outer=0.95) + 
  scale_y_discrete(limits=c("r_DQR[D,Intercept]","r_DQR[C,Intercept]",
                            "r_DQR[B,Intercept]","r_DQR[A,Intercept]"),
                   labels=c("r_DQR[A,Intercept]" = bquote(alpha[A]),
                            "r_DQR[B,Intercept]" = bquote(alpha[B]),
                            "r_DQR[C,Intercept]" = bquote(alpha[C]),
                            "r_DQR[D,Intercept]" = bquote(alpha[D]))) +  
  scale_x_continuous(limits=c(-1.5, 1.5)) +
  theme_tufte(base_size = 16) 

A_clean <- brm(bf(Effort ~ Input + Output + Enquiry + File + 
                         Interface + Added + Changed),
                    data=d_ifpug_A_clean, seed=SEED,
                    family = negbinomial,
                    prior = c(prior(normal(5, 4), class = Intercept),
                              prior(normal(0, 0.25), class = b),
                              prior(normal(0.5, 0.5), class = shape)),
                    cores=4, chains=4)

posterior <- as.matrix(A_clean)
color_scheme_set("gray")
p3 <- mcmc_areas_ridges(posterior, pars = c("b_Input", "b_Output", "b_Enquiry", 
                                         "b_File", "b_Interface", "b_Added", 
                                         "b_Changed"), 
                     prob=0.95, 
                     prob_outer=1) +
  scale_x_continuous(limits=c(-0.9, 0.9)) +
  scale_y_discrete(limits=c("b_Changed","b_Added","b_Interface","b_File",
                            "b_Enquiry","b_Output","b_Input"),
                   labels=c("b_Changed" = bquote(beta[Changed]),
                            "b_Added" = bquote(beta[Added]),
                            "b_Interface" = bquote(beta[Interface]),
                            "b_File" = bquote(beta[File]),
                            "b_Enquiry" = bquote(beta[Enquiry]),
                            "b_Output" = bquote(beta[Output]),
                            "b_Input" = bquote(beta[Input]))) +
  theme_tufte(base_size = 14) #+
  #theme(axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank())


# So, easy peasy to fit with brms. The above datasets had NAs removed. 
# The below datasets have NAs...
################################################################################

# Number of imputations
m = 25

# Imputation
imp_AD <- mice(d_ifpug_AD, m, print = FALSE, seed=SEED)

# Extract our imputed datasets
df_list <- list()
for (i in 1:m) {
  df_list[[i]] <- complete(imp_AD, i)
}

# For each dataset, standardize our predictors if they're not categorical
for (i in 1:m) {
  df_list[[i]] <- df_list[[i]] %>% mutate_at(.vars = vars(Input, Output,
                                                          Enquiry, File,
                                                          Interface, Added,
                                                          Changed),
                                             .funs = funs(scale(.)))
}

AD <- brm_multiple(bf(Effort ~ Input + Output + Enquiry + File + Interface + 
                        Added + Changed + (1 | DQR)),
                   data = df_list, 
                   family = negbinomial,
                   prior = c(prior(normal(5, 4), class = Intercept),
                             prior(normal(0, 0.25), class = b),
                             prior(gamma(0.5, 0.5), class=shape),
                             prior(cauchy(0, 1), class=sd)),
                   cores = 4, chains = 4, control = list(adapt_delta=0.999999, 
                                                         max_treedepth=13))

# If we do summary(AD) we'll see large \hat{R} values and low n_eff.
# \hat{R} is the estimated between-chains and within-chain variances for each 
# model parameter. Large differences between these variances indicate 
# nonconvergence. BUT, if one sampling is ok then we can assume the others are 
# also ok. \hat{R} are useful for one sampling, (i.e., the 8 chains), but if 
# we sample with the same model 25 times (m=25) then we will have 200 chains to 
# calculate \hat{R}, i.e., leading to bias.
# See Brooks, S. P., and A. Gelman. 1997. General Methods for Monitoring 
# Convergence of Iterative Simulations. Journal of Computational and Graphical 
# Statistics 7: 434–455.
# Check for max value among all and each sampling:
# > max(AD$rhats)
# [1] 1.008099

posterior <- as.matrix(AD)

color_scheme_set("darkgray")
p2 <- mcmc_areas_ridges(posterior, pars = c("b_Input", "b_Output", "b_Enquiry", 
                                         "b_File", "b_Interface", "b_Added", 
                                         "b_Changed"),
                     prob=0.95, prob_outer=1) + 
  scale_x_continuous(limits=c(-0.9, 0.9)) +
  theme_tufte(base_size = 14) +
  ggtitle("AD") +
  theme(axis.text.y=element_blank(),
       axis.ticks.y=element_blank(),
       axis.text.x=element_blank(), 
       axis.ticks.x=element_blank(),
       plot.title = element_text(hjust = 0.5))

color_scheme_set("darkgray")
p_alphas_ad <- mcmc_intervals(posterior, regex_pars = c("^r"),
                                    prob_outer=0.95) + 
  scale_y_discrete(limits=c("r_DQR[D,Intercept]","r_DQR[C,Intercept]",
                            "r_DQR[B,Intercept]","r_DQR[A,Intercept]"),
                   labels=c("r_DQR[A,Intercept]" = bquote(alpha[A]),
                            "r_DQR[B,Intercept]" = bquote(alpha[B]),
                            "r_DQR[C,Intercept]" = bquote(alpha[C]),
                            "r_DQR[D,Intercept]" = bquote(alpha[D]))) +
  theme_tufte(base_size = 16) +
  scale_x_continuous(limits=c(-1, 1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# theme(axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())
################################################################################

# Imputation
imp_A <- mice(d_ifpug_A, m, print = FALSE, seed=SEED) #impute

# Extract our imputed datasets
df_list <- list()
for (i in 1:m) {
  df_list[[i]] <- complete(imp_A, i)
}

# For each dataset, scale our predictors if they're not categorical
for (i in 1:m) {
  df_list[[i]] <- df_list[[i]] %>% mutate_at(.vars = vars(Input, Output,
                                                          Enquiry, File,
                                                          Interface, Added,
                                                          Changed),
                                             .funs = funs(scale(.)))
}

A <- brm_multiple(bf(Effort ~ Input + Output + Enquiry + File + Interface + 
                       Added + Changed),
                  data = df_list,
                  family = negbinomial,
                  prior = c(prior(normal(5, 4), class = Intercept),
                            prior(normal(0, 0.25), class = b), 
                            prior(gamma(0.5, 0.5), class=shape)),
                  cores=4, chains=4, control=list(adapt_delta=0.95))

# Once again, when we have m=25 we disregard \hat{R} and n_eff, and check each 
# sampling separately.
# > max(A$rhats)
# [1] 1.006118

posterior <- as.matrix(A)
color_scheme_set("darkgray")
p1 <- mcmc_areas_ridges(posterior, pars = c("b_Input", "b_Output", "b_Enquiry", 
                                         "b_File", "b_Interface", "b_Added", 
                                         "b_Changed"), 
                     prob=0.95, 
                     prob_outer=1) +
  scale_x_continuous(limits=c(-0.9, 0.9)) +
  scale_y_discrete(limits=c("b_Changed","b_Added","b_Interface","b_File",
                            "b_Enquiry","b_Output","b_Input"),
                   labels=c("b_Changed" = bquote(beta[Changed]),
                            "b_Added" = bquote(beta[Added]),
                            "b_Interface" = bquote(beta[Interface]),
                            "b_File" = bquote(beta[File]),
                            "b_Enquiry" = bquote(beta[Enquiry]),
                            "b_Output" = bquote(beta[Output]),
                            "b_Input" = bquote(beta[Input]))) +
  ggtitle("A") +
  theme_tufte(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
  

################################################################################

# Interval plots of each of the seven estimated parameters
pdf("~/betas.pdf", height=7, width=7)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

# Interval plots of the intercepts for DQR
grid.arrange(p_alphas_ad, p_alphas_ad_clean, ncol = 1)

# # Check out level A data
# data_A <- d_ifpug_AD_clean[d_ifpug_AD_clean$DQR == 
#                              levels(d_ifpug_AD_clean$DQR)[1],]
# data_B <- d_ifpug_AD_clean[d_ifpug_AD_clean$DQR == 
#                              levels(d_ifpug_AD_clean$DQR)[2],]
# data_C <- d_ifpug_AD_clean[d_ifpug_AD_clean$DQR == 
#                              levels(d_ifpug_AD_clean$DQR)[4],]
# data_D <- d_ifpug_AD_clean[d_ifpug_AD_clean$DQR == 
#                              levels(d_ifpug_AD_clean$DQR)[4],]
# # condition on `A'
# PPD_A <- posterior_predict(AD_clean, newdata = data_A)
# PPD_B <- posterior_predict(AD_clean, newdata = data_B)
# PPD_C <- posterior_predict(AD_clean, newdata = data_C)
# PPD_D <- posterior_predict(AD_clean, newdata = data_D)
# 
# AD_clean %>%
#   spread_draws(r_DQR[condition,]) %>%
#   compare_levels(r_DQR, by = condition) %>%
#   ggplot(aes(y = condition, x = r_DQR)) +
#   geom_halfeyeh() +
#   theme_tufte()

# Include random effects for the posterior distributions of interest
ran_a_clean <- predict(A_clean, nsamples = 4000)
ran_ad <- predict(AD, nsamples = 4000)

med_a_clean <- log10(median(ran_a_clean[,1])) # 2937
# hdi(ran_a_clean[,1]) # 2045;15225

med_ad <- log10(median(ran_ad[,1])) # 3276
# hdi(ran_ad[,1]) # 1925;17937

df1 <- data.frame("dens" = log10(ran_a_clean[,1]))
df1$set <- rep("A_clean", times=214)

df2 <- data.frame("dens" = log10(ran_ad[,1]))
df2$set <- rep("AD", times=1689)

df <- rbind(df1, df2)
df$set <- as.factor(df$set)

ggplot(df,aes(x=dens, fill=set)) + 
  geom_density(alpha=0.5) + 
  scale_fill_manual(values=c("gray", "gray20")) +
  coord_cartesian(xlim = c(3,5), ylim=c(0,4)) +
  theme(legend.position = "none", axis.title.y=element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey90"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
        #text=element_text(size=10,  family="serif")
        ) +
  geom_vline(xintercept=med_a_clean, linetype="dashed") +
  geom_vline(xintercept=med_ad, linetype="dashed") +
  annotate("text", x=3.58, y=2.8, label="AD") +
  annotate("text", x=3.25, y=2.8, label="A_clean") +
  xlab(expression(paste(log[10], " (Effort)"))) 

# posterior predictive check
color_scheme_set("darkgray")

ppc_dens_overlay(y = d_ifpug_AD_clean$Effort,
                 yrep = posterior_predict(AD_clean, nsamples = 50)) +
  scale_x_log10() +
  theme_tufte() +
  theme(legend.position="none") +
  xlab(expression(paste(log[10], " (Effort)")))

np <- rhat(AD_clean)
mcmc_rhat(np) + theme_tufte() + theme(legend.position="none") + 
  geom_vline(xintercept = 1.0035, linetype="dashed") +
  coord_cartesian(xlim = c(1,1.005)) +
  scale_x_continuous(breaks=c(1,1.0035), labels=c(1,1.0035))

nr <- neff_ratio(AD_clean)
mcmc_neff(nr) + theme_tufte() + theme(legend.position="none") 

mcmc_trace(as.array(AD_clean), pars = c("b_Intercept", "b_Input", "b_Changed", "b_Enquiry")) + 
  theme_tufte() + theme(legend.position = "none")

mon <- monitor(as.array(AD_clean))
samp <- as.array(AD_clean)

# Pick the param with the lowest neff in the tail
which_min_ess <- which.min(mon[1:14, 'Tail_ESS'])
xmin <- as.integer(paste0(which_min_ess))

mcmc_hist_r_scale(samp[, , xmin])

# nuts energy levels next (Bayesian fraction of missing information)
np <- nuts_params(AD_clean)
mcmc_nuts_energy(np) + theme_tufte() + theme(legend.position = "none")

# PPC
yrep <- posterior_predict(AD_clean)
ppc_intervals(y=AD_clean$data$Effort, yrep=yrep) + 
  coord_cartesian(xlim = c(1,10), ylim = c(0,50000)) +
  scale_x_continuous(breaks = seq(1:10)) +
  xlab("Project #") +
  ylab("Effort") +
  theme_tufte() +
  theme(legend.position = "none")

mean(as.data.frame(AD_clean)$b_Input > 0) # prob that Input is >0
################################################################################
# TODO
# * MCMC diag. http://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
# * Remove cruft?
# exclude <- c("lp__")
# par_names <- setdiff(dimnames(x)[[3]], exclude)
# mcmc_parcoord(x, pars = par_names)