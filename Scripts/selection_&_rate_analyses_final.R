
## ----source_functions_&_libraries, results='hide', echo=FALSE, message=FALSE----
rm(list=ls())

require(MASS)
require(lme4)
require(MCMCglmm)
require( plyr )
require( lme4 )
require( VGAM )

# function to mean-standardize covariance matrices
StanCovMat <- function( X, X_means ) {
  X / ( X_means %*% t(X_means) )
}

# function for including correlations on the pairs plots
panel.cor <- function(x, y, digits=2) {
	r <- cor(x, y, use="complete")
	par( usr =c(0, 1, 0, 1))
	Cor <- format(c(r, 0.123456789), digits=digits)[1]
 	text(0.5, 0.5, paste("r=", Cor), cex=1.5)
}

# function to apply the folded normal distribution
FoldedNormal <- function( X ) {
	X.means <- apply( X, 2, mean )
	X.sd 	<- apply( X, 2, sd )
	X.dim	<- apply( X, 2, length )
	fn.X 	<- matrix( NA, nrow=dim(X)[1], ncol=dim(X)[2] ) 
	for( i in 1:dim(X)[2] ) {
		fn.X[,i] <- rfnorm( X.dim[i], X.means[i], X.sd[i]) }
	fn.X
}

# 'apply'-able HPDinterval function
hpd <- function( X ) {
  HPDinterval( as.mcmc( X ), probs=c(0.025, 0.975) )[1:2]
}

# mode (central tendency, not data-storage type) function
Mode <- function( x ) {
  ux <- unique( x )
  ux[ which.max( tabulate( match( x, ux))) ]
}

# univariate Rsquared calculator
Rsq <- function( model ){
  fitted.variance <- var(model$fitted)
	total.variance	<- var(model$fitted) + var(model$resid)
	fitted.variance / total.variance
}

# geometric mean function
g_mean <- function( X ) {
    log_data <- log( X )
    gm <- exp( mean( log_data[ is.finite( log_data )] ))
    return( gm )
}



## ----read_in_data, echo=TRUE, results='asis', message=TRUE---------------
response <- read.csv("../Data/All_rates_Haldanes.csv")

response <- response[response$Trait.type!="P",]

response$Trait.type <- factor(response$Trait.type)

response$absRate <- abs(response$Rate)

# str(response)
# head(response)


## ----data_hygiene, results='hide', fig.keep='none', echo=FALSE-----------
hist( response$Rate )		# no troublesome outliers by eye
pairs( response[, c(2, 5, 9)] ) # hmmm, some oddly high rates
response[ response$absRate >1.5 ,]  # both from study s105 (exp. evo.)

plot( response$absRate ~ response$NatExp )
Rsq( lm( response$absRate ~ response$NatExp ) )	# 0.016 - not huge, but needs to be in the model


## ----summary_information_for_table, echo=FALSE---------------------------
ddply( response, .(Taxon2, Trait.type), summarise, count=length( absRate ), Mean=round( mean( absRate ),3 ), 
       Gmean=round( g_mean( absRate ), 3 ), Median=round(median( absRate ), 3), Mode=Mode( absRate ) )

# plant estimates are too few and are unbalancing our models
responseP <- response[ response$Taxon2 == 'P', ]
responseA <- response[ response$Taxon2 == 'A', ]


## ----model_selection_for_rates_Animals, echo=FALSE, message=FALSE, results='hide', warning=FALSE, fig.keep='none'----
# let's fit some models
DZ.reml0 <- lmer( absRate ~ 1 + (1|Study), data=responseA )

DZ.reml1 <- lmer( absRate ~ 1 + Trait.type + NatExp + (1|Study), data=responseA )

DZ.reml2 <- lmer( absRate ~ 1 + Trait.type + NatExp + (1|Study) + (1|Species), data=responseA )

DZ.reml3 <- lmer( absRate ~ 1 + Trait.type + NatExp + (Trait.type|Species) + (1|Study), data=responseA )

bbmle::BICtab( DZ.reml0, DZ.reml1, DZ.reml2, DZ.reml3, nobs=nrow(responseA), weights=T )
    #          dBIC df weight
    # DZ.reml0  0.0 3  1     
    # DZ.reml1 36.4 6  <0.001
    # DZ.reml2 44.2 7  <0.001
    # DZ.reml3 54.0 12 <0.001

bbmle::AICctab( DZ.reml0, DZ.reml1, DZ.reml2, DZ.reml3, nobs=nrow(responseA), weights=T )
    #          dAICc df weight
    # DZ.reml0  0.0  3  0.691 
    # DZ.reml3  1.6  12 0.309 
    # DZ.reml1 18.9  6  <0.001
    # DZ.reml2 20.9  7  <0.001

# BIC & AIC both prefer the simplest model...

# parametric bootstrap

# model 0 vs model 1
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
	model.lower <- lmer(y.sim$sim_1 ~ 1 + (1|Study), data=responseA )
	model.full  <- lmer(y.sim$sim_1 ~ 1 + Trait.type + NatExp + (1|Study), data=responseA )
	LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( DZ.reml0 ))
LR.model <- -as.numeric(deviance( DZ.reml1 ) - deviance( DZ.reml0 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.2207792 this suggests sticking with model0


# model 1 vs model 2
LikRatioSim <- function( mod ) {
	y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
	model.lower <- lmer(y.sim$sim_1 ~ 1 + Trait.type + NatExp + (1|Study), data=responseA )
	model.full  <- lmer(y.sim$sim_1 ~ 1 + Trait.type + NatExp + (1|Study) + (1|Species), data=responseA )
	LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( DZ.reml1 ))
LR.model <- -as.numeric(deviance( DZ.reml2 ) - deviance( DZ.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.8251748 this suggests sticking with model1


# model 1 vs model 3
LikRatioSim <- function( mod ) {
	y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
	model.lower <- lmer(y.sim$sim_1 ~ 1 + Trait.type + NatExp + (1|Study), data=responseA )
	model.full  <- lmer(y.sim$sim_1 ~ 1 + Trait.type + NatExp + (Trait.type|Species) + (1|Study), data=responseA )
	LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( DZ.reml1 ))
LR.model <- -as.numeric(deviance( DZ.reml3 ) - deviance( DZ.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.000999001 support for model3 over simpler models in disagreement with BIC...




# ...and now, the bayesian version:
myprior1 <- list( list(V=1, nu=0.002), list( G1=list(V=1, nu=0.002)) )

myprior2 <- list( list(V=1, nu=0.002), list( G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002) ) )

myprior3 <- list( list(V=1, nu=0.002), list( G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002) ) )

# NB: remember to put in supp. that we tried priors with V=1 & nu=0 as well without changing the results

DZ.Bayes0 <- MCMCglmm( absRate ~ 1, random=~Study+Species, nitt=10000, burnin=500, thin=10 , data=responseA, verbose=F, family="gaussian")
summary( DZ.Bayes0 )$solutions
            # post.mean   l-95% CI  u-95% CI eff.samp        pMCMC
# (Intercept) 0.1299152 0.08540443 0.1739915     9500 0.0001052632

resp <- FoldedNormal( DZ.Bayes0$Sol )
mean( resp )   ;   hpd( resp )       ;   Mode( resp )
# 0.1324803    0.08826479 - 0.17283631   0.1456772


DZ.Bayes1 <- MCMCglmm( absRate ~ 1 + Trait.type + NatExp, random=~ Study, nitt=50000, burnin=1000, thin=15 , data=responseA, verbose=F, family="gaussian", prior=myprior1 )

DZ.Bayes2 <- MCMCglmm( absRate ~ 1 + Trait.type + NatExp, random=~ Species+Study, nitt=50000, burnin=1000, thin=15 , data=responseA, verbose=F, family="gaussian", prior=myprior2 )

DZ.Bayes3 <- MCMCglmm( absRate ~ 1 + Trait.type + NatExp, random=~ idh(Trait.type):Species + Study, nitt=50000, burnin=1000, thin=25 , data=responseA, verbose=F, family="gaussian", prior=myprior2 )

# DIC prefers the more complex model like LR but unlike BIC... we shall go with model 3
DZ.Bayes0$DIC # -5740.594
DZ.Bayes1$DIC	# -5739.073
DZ.Bayes2$DIC	# -5739.149
DZ.Bayes3$DIC	# -5785.765


plot(DZ.Bayes3$VCV)
acf(DZ.Bayes3$VCV)
plot(DZ.Bayes3$Sol)
acf(DZ.Bayes3$Sol)
summary(DZ.Bayes3)$solutions
              # post.mean    l-95% CI  u-95% CI eff.samp     pMCMC
# (Intercept) 0.018889154 -0.17412625 0.2107472 19600.00 0.8424490
# Trait.typeM 0.004775399 -0.10479113 0.1013593 19600.00 0.8676531
# Trait.typeS 0.063132204 -0.05431215 0.1796105 19221.73 0.2743878
# NatExpN     0.106405569 -0.06432589 0.2786911 19600.00 0.2144898

# weight estimates by numbers of natural/experimental estimates (2483/88)
weightedintercept <- ((DZ.Bayes3$Sol[,1] * length(responseA$Rate[ responseA$NatExp=="E" ]) ) / dim(responseA)[1]) +  
              ((DZ.Bayes3$Sol[,1] + DZ.Bayes3$Sol[,4]) * length(responseA$Rate[ responseA$NatExp=="N" ]) ) / dim(responseA)[1]

estimates <- FoldedNormal( cbind( weightedintercept,
  	                      ( weightedintercept + DZ.Bayes3$Sol[,2] ),
		                      ( weightedintercept + DZ.Bayes3$Sol[,3] ) ))

colnames( estimates ) <- c( 'A.LH', 'A.M', 'A.S' )

# tabulate estimates
apply( estimates, 2, mean )  ;	apply( estimates, 2, hpd )
		      # A.LH        A.M        A.S 
	# mean	0.1209507	0.1260378 0.1845940
	# -95% 	0.01470875 0.08102147 0.1135921
	# +95%  0.22481614 0.17068093 0.2523750


## ----plot_results_Animals, echo=FALSE, results='hide', fig.keep='high'----
Rfixed <- FoldedNormal( DZ.Bayes3$Sol )
plot(Rfixed[1:3,1], type="n", xaxt="n", ylab="Rate of evolution (haldanes)", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(-0.05,0.45), main="Rate of Evolution")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(estimates[,3])), lwd=2 )

A.est <- c( mean( estimates[,1] ), mean( estimates[,2] ), mean( estimates[,3] ) )
points( 1:3, A.est, pch=16, cex=2)  #animals - filled points


## ----model_selection_for_Plants, echo=FALSE, results='hide', message=FALSE, fig.keep='none', warning=FALSE----
DZ.reml.P.0 <- lmer( absRate ~ 1 + (1|Species), data=responseP )

DZ.reml.P.1 <- lmer( absRate ~ 1 + Trait.type + (1|Species), data=responseP )

DZ.reml.P.2 <- lmer( absRate ~ 1 + Trait.type + Study + (1|Species), data=responseP )

bbmle::BICtab( DZ.reml.P.0, DZ.reml.P.1, DZ.reml.P.2, nobs=nrow(responseP), weights=T )
      #             dBIC df weight
      # DZ.reml.P.1 0.0  4  0.575 
      # DZ.reml.P.0 1.3  3  0.293 
      # DZ.reml.P.2 3.0  5  0.131

bbmle::AICctab( DZ.reml.P.0, DZ.reml.P.1, DZ.reml.P.2, nobs=nrow(responseP), weights=T )
      #             dAICc df weight
      # DZ.reml.P.1 0.0   4  0.606 
      # DZ.reml.P.0 2.2   3  0.197 
      # DZ.reml.P.2 2.3   5  0.197


# model 0 vs model 1
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
	model.lower <- lmer(y.sim$sim_1 ~ 1 + (1|Species), data=responseP )
	model.full  <- lmer(y.sim$sim_1 ~ 1 + Trait.type + (1|Species), data=responseP )
	LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( DZ.reml.P.0 ))
LR.model <- -as.numeric(deviance( DZ.reml.P.1 ) - deviance( DZ.reml.P.0 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.003996004 support here for moving up to model 1


# model 1 vs model 2
LikRatioSim <- function( mod ) {
	y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
	model.lower <- lmer(y.sim$sim_1 ~ 1 + Trait.type + (1|Species), data=responseP )
	model.full  <- lmer(y.sim$sim_1 ~ 1 + Trait.type + Study + (1|Species), data=responseP )
	LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
	return(LRSim)
	rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( DZ.reml.P.1 ))
LR.model <- -as.numeric(deviance( DZ.reml.P.2 ) - deviance( DZ.reml.P.1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.1948052 no support for model 2 over model 1 though...


DZ.Bayes.P.0 <- MCMCglmm( absRate ~ 1, random=~Species, nitt=500000, burnin=10000, thin=20 , data=responseP, verbose=F, family="gaussian")

respP <- FoldedNormal( DZ.Bayes.P.0$Sol )
mean( respP )     ;   hpd( respP )    ;   median( respP )
# 0.1891405     0.08628083 - 0.28893978     0.1894191

DZ.Bayes.P.1 <- MCMCglmm( absRate ~ 1 + Trait.type, random=~Species, nitt=500000, burnin=10000, thin=20 , data=responseP, verbose=F, family="gaussian")

DZ.Bayes.P.2 <- MCMCglmm( absRate ~ 1 + Trait.type + Species, random=~Species, nitt=500000, burnin=10000, thin=20 , data=responseP, verbose=F, family="gaussian")

# DIC gives a very close call between models 1 & 2, but overall everything points to model 1
DZ.Bayes.P.0$DIC  ;   DZ.Bayes.P.1$DIC	;	DZ.Bayes.P.2$DIC
# -67.59389          -79.80158	    				-79.46739


plot(DZ.Bayes.P.1$VCV)
plot(DZ.Bayes.P.1$Sol)
summary( DZ.Bayes.P.1 )$solutions

# NB: this is based on only 33 estimates and so we suggest interpteting with caution. 
             # post.mean   l-95% CI    u-95% CI eff.samp        pMCMC
# (Intercept)  0.2983510  0.1785164  0.41695026 99500.00 6.030151e-05
# Trait.typeM -0.1434201 -0.2370725 -0.04935445 65338.46 3.035176e-03

plant.estimates <- FoldedNormal( cbind( DZ.Bayes.P.1$Sol[,1], DZ.Bayes.P.1$Sol[,1] + DZ.Bayes.P.1$Sol[,2] ))

apply( plant.estimates, 2, mean )	;	apply( plant.estimates, 2, hpd )
	# mean 	0.2983510 	0.1549309
	# -95% 	0.1785164 	0.05338595
	# +95%  0.4169503 	0.25250416


## ----results_plants_&_animals, results='hide', echo=FALSE, fig.keep='high'----
# remake rates plot with plants added
plot(Rfixed[1:3,1], type="n", xaxt="n", ylab="Rate of evolution (haldanes)", xlab="Trait type", xlim=c(0.7,3.3), 
      ylim=c(-0.05,0.45), main="Rate of Evolution")
      axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(estimates[,3])), lwd=2 )

lines( c(1.1,1.1), c(hpd(plant.estimates[,1])), lwd=2 )
lines( c(2.1,2.1), c(hpd(plant.estimates[,2])), lwd=2 )

points( 1:3, A.est, pch=16, cex=2)  #animals - filled points

P.est <- c( mean( plant.estimates[,1] ), mean( plant.estimates[,2] ) )
points( 1.1:2.1, P.est, pch=21, bg="white", cex=2)  # plants - empty points


## ----read_in_selection_data, echo=FALSE, results='hide'------------------
select <- read.csv("../Data/combined_selection_data.csv", na.strings="NA")

str(select)
names(select)
summary(select)
levels( select$Species )

select$Taxon2 <- factor( sub( "V", "A", sub( "I", "A", select$Taxon )))


## ----selection_data_hygiene, echo=FALSE, message=FALSE-------------------
names( select )

pairs( select[, c(5,7,8,10,11)] )
plot( select$GamPos )	;	# identify( select$GamPos )	# check 459 466
select <- select[ -c(459, 466),]

# Summary information for table?

print( "counts of gradients" )
ddply( select, .(Taxon2, Trait.type), summarise, n_abs_Beta= sum(!is.na(abs_Beta)), n_GamPos= sum(!is.na(GamPos)), n_GamNeg= sum(!is.na(GamNeg)) )

print( "means of gradients" )
ddply( select, .(Taxon2, Trait.type), summarise, m_abs_Beta= mean(!is.na(abs_Beta)), m_GamPos= mean(!is.na(GamPos)), m_GamNeg= mean(!is.na(GamNeg)) )

print( "geometric means of gradients" )
ddply( select, .(Taxon2, Trait.type), summarise, gm_abs_Beta= g_mean((abs_Beta)), gm_GamPos= g_mean((GamPos)), gm_GamNeg= g_mean(abs(GamNeg)) )

print( "medians of gradients" )
ddply( select, .(Taxon2, Trait.type), summarise, me_abs_Beta= median(na.omit(abs_Beta)), me_GamPos= median(na.omit(GamPos)), me_GamNeg= median(na.omit(GamNeg)) )


## ----modelling_beta, results='hide', echo=FALSE, fig.keep='none'---------
# A more formal meta-analysis, after Kingsolver et al. 2002
myprior <- list(R = list(V=1, nu=0.002), G= list( G1= list (V=1, nu=0.002), G2= list (V=1, nu=0.002)))

# for this analysis we use only those gradients for which standard errors are available.
select_B <- subset( select, select$Beta_SE > 0 )
measurement_error_variance =  select_B$Beta_SE^2  

with( select_B, by( abs_Beta, list(Taxon2, Trait.type), mean, na.rm=T ) )

ddply( select_B, .(Taxon2, Trait.type), summarise, n_abs_Beta= sum(!is.na(abs_Beta)), n_GamPos= sum(!is.na(GamPos)), n_GamNeg= sum(!is.na(GamNeg)) )
dim( select_B )

Beta_Bayes <- MCMCglmm( Beta ~ 1 + Trait.type * Taxon2, mev= measurement_error_variance, data= select_B, random= ~ StudyID + Species, prior= myprior, nitt=500000, burnin=10000, thin=25 )

plot(Beta_Bayes$VCV)
acf(Beta_Bayes$VCV)
plot(Beta_Bayes$Sol)
acf(Beta_Bayes$Sol)

summary(Beta_Bayes)$solutions
                       # post.mean    l-95% CI  u-95% CI eff.samp pMCMC
# (Intercept)          0.003630710 -0.13526199 0.1388960 20006.55 0.96
# Trait.typeM          0.058247419 -0.09131960 0.2071727 19600.00 0.44
# Trait.typeS          0.091705159 -0.06637597 0.2534982 20054.42 0.25
# Taxon2P              0.002583934 -0.17788472 0.1807363 19141.49 0.96
# Trait.typeM:Taxon2P  0.060256293 -0.18036532 0.3147182 19600.00 0.63
# Trait.typeS:Taxon2P -0.100812658 -0.44168627 0.2751529 19114.98 0.57

estimatesB <- cbind( Beta_Bayes$Sol[,1],
  			( Beta_Bayes$Sol[,1] + Beta_Bayes$Sol[,2] ),
				( Beta_Bayes$Sol[,1] + Beta_Bayes$Sol[,3] ),
				( Beta_Bayes$Sol[,1] + Beta_Bayes$Sol[,4] ),
				( Beta_Bayes$Sol[,1] + Beta_Bayes$Sol[,2] + Beta_Bayes$Sol[,4] + Beta_Bayes$Sol[,5] ),
				 ( Beta_Bayes$Sol[,1] + Beta_Bayes$Sol[,3] + Beta_Bayes$Sol[,4] + Beta_Bayes$Sol[,6] ) )

colnames( estimatesB ) <- c( 'A.LH', 'A.M', 'A.S', 'P.LH', 'P.M', 'P.S' )

# tabulate estimates
apply( estimatesB, 2, mean )  ;	apply( estimatesB, 2, hpd )
     # A.LH          A.M          A.S         P.LH          P.M       P.S
# mean 0.00293137  0.059443773 0.09534618  0.06382877  0.13169555 -0.00588
# -95% -0.1269036 -0.004199983 0.01557039 -0.04471114 -0.01632193 -0.27560
# +95%  0.1288168  0.125664502 0.17454144  0.16498689  0.28107967  0.25879



## ----plotting_results_for_beta, echo=FALSE, message=FALSE, fig.keep='high'----
Bfixed <- summary(Beta_Bayes)$solutions
plot(Bfixed[c(1,3,4),1], type="n", xaxt="n", ylab="Beta values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(-0.3,0.5), main="Beta")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(estimatesB[,1])), lwd=2 )
lines( c(2,2), c(hpd(estimatesB[,2])), lwd=2 )
lines( c(3,3), c(hpd(estimatesB[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(estimatesB[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(estimatesB[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(estimatesB[,6])), lwd=2 )

A.est <- c( mean(estimatesB[,1]), mean(estimatesB[,2]), mean(estimatesB[,3]) )
points( 1:3, A.est, pch=16, cex=2)  #animals - filled points
P.est <- c( mean(estimatesB[,4]), mean(estimatesB[,5]), mean(estimatesB[,6]) )
points( 1.1:3.1, P.est, pch=21, bg="white", cex=2)


## ----modelling_absolute_beta, results='hide', echo=FALSE, fig.keep='none'----
absBeta_Bayes <- MCMCglmm( abs_Beta ~ 1 + Trait.type * Taxon2, mev= measurement_error_variance, data= select_B, random= ~ StudyID + Species,  nitt=500000, burnin=10000, thin=25 )

plot(absBeta_Bayes$VCV)
acf(absBeta_Bayes$VCV)
plot(absBeta_Bayes$Sol)
acf(absBeta_Bayes$Sol)

summary(absBeta_Bayes)$solutions
                      # post.mean    l-95% CI   u-95% CI eff.samp  pMCMC
# (Intercept)          0.09107684 -0.02856087 0.20983240 19600.00 0.137
# Trait.typeM          0.12764927 -0.00347855 0.26598995 19600.00 0.061
# Trait.typeS          0.09434900 -0.04357560 0.23587526 19600.00 0.186
# Taxon2P              0.18650175  0.03538563 0.33824004 18348.24 0.017
# Trait.typeM:Taxon2P -0.18440815 -0.39395161 0.02579394 19600.00 0.082
# Trait.typeS:Taxon2P -0.23177637 -0.53649314 0.07403859 19600.00 0.132

estimatesAB <- FoldedNormal( cbind( absBeta_Bayes$Sol[,1],
  			                  ( absBeta_Bayes$Sol[,1] + absBeta_Bayes$Sol[,2] ),
				                  ( absBeta_Bayes$Sol[,1] + absBeta_Bayes$Sol[,3] ),
				                  ( absBeta_Bayes$Sol[,1] + absBeta_Bayes$Sol[,4] ),
				                  ( absBeta_Bayes$Sol[,1] + absBeta_Bayes$Sol[,2] + absBeta_Bayes$Sol[,4] + absBeta_Bayes$Sol[,5] ),
				                  ( absBeta_Bayes$Sol[,1] + absBeta_Bayes$Sol[,3] + absBeta_Bayes$Sol[,4] + absBeta_Bayes$Sol[,6] ) ))

colnames( estimatesAB ) <- c( 'A.LH', 'A.M', 'A.S', 'P.LH', 'P.M', 'P.S' )


# tabulate estimates
apply( estimatesAB, 2, mean )  ;	apply( estimatesAB, 2, hpd )
      # A.LH        A.M        A.S       P.LH        P.M        P.S 
# mean  0.09195647 0.2193015 0.1868024 0.3084369 0.21890802  0.1388694 
# -95% -0.03255756 0.1564934 0.1113526 0.2091150 0.07847622 -0.1369955
# +95%  0.21634752 0.2846237 0.2642353 0.4030195 0.36491964  0.4040568


## ----plotting_results_for_absolute_beta, echo=FALSE, message=FALSE, fig.keep='high'----

aBfixed <- summary(absBeta_Bayes)$solutions
plot(aBfixed[c(1,3,5),1], type="n", xaxt="n", ylab="(absolute) Beta values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(-0.1,0.5), main="Beta (absolute)")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(estimatesAB[,1])), lwd=2 )
lines( c(2,2), c(hpd(estimatesAB[,2])), lwd=2 )
lines( c(3,3), c(hpd(estimatesAB[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(estimatesAB[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(estimatesAB[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(estimatesAB[,6])), lwd=2 )

A.est <- c( mean(estimatesAB[,1]), mean(estimatesAB[,2]), mean(estimatesAB[,3]) )
points( 1:3, A.est, pch=16, cex=2)  #animals - filled points
P.est <- c( mean(estimatesAB[,4]), mean(estimatesAB[,5]), mean(estimatesAB[,6]) )
points( 1.1:3.1, P.est, pch=21, bg="white", cex=2)


## ----modelling_gamma+, results='hide', echo=FALSE, fig.keep='last'-------

select_G <- subset( select, select$Gamma_SE > 0 )

ddply( select_G, .(Taxon2, Trait.type), summarise, n_abs_Beta= sum(!is.na(abs_Beta)), n_GamPos= sum(!is.na(GamPos)), n_GamNeg= sum(!is.na(GamNeg)) )
dim( select_G )

ddply( select_G, .(Taxon2, Trait.type), summarise, m_GamPos= mean(na.omit(GamPos)), med_GamPos= median(na.omit(GamPos)), m_GamNeg= mean(na.omit(GamNeg)), med_GamNeg= median(na.omit(GamNeg)) )

measurement_error_variance <- select_G$Gamma_SE^2

GamPos_1_Bayes <- MCMCglmm( GamPos ~ 1, mev= measurement_error_variance, data= select_G, random= ~ StudyID + Species, prior= myprior, nitt=500000 )
summary( GamPos_1_Bayes )$solutions
            # post.mean   l-95% CI  u-95% CI eff.samp       pMCMC
# (Intercept) 0.1035137 0.04507204 0.1630329 1932.665 0.001276596

GP1B<- FoldedNormal( GamPos_1_Bayes$Sol )
mean( GP1B )  	;	  hpd( GP1B )   ;     Mode( GP1B )
# 0.104019   0.03829851 -- 0.17745633         0.1048835


GamPos_Bayes <- MCMCglmm( GamPos ~ 1 + Trait.type * Taxon2, mev= measurement_error_variance, data= select_G, random= ~ StudyID + Species, prior= myprior, nitt=500000, burnin=10000, thin=25 )

plot(GamPos_Bayes)
acf(GamPos_Bayes$VCV)
summary(GamPos_Bayes)$solutions
                      # post.mean    l-95% CI  u-95% CI eff.samp  pMCMC
# (Intercept)          0.01164284 -0.15621767 0.1984974 1492.240 0.9223
# Trait.typeM          0.08137383 -0.12917569 0.2672039 1146.888 0.4208
# Trait.typeS          0.16081478 -0.06034077 0.3639002 1246.503 0.1451
# Taxon2P              0.16123892 -0.06709793 0.3892172 2666.024 0.1560
# Trait.typeM:Taxon2P -0.20538125 -0.54001117 0.1085685 2482.749 0.2036
# Trait.typeS:Taxon2P -0.31377590 -0.77041684 0.1119918 5689.849 0.1450

estimatesGP <- FoldedNormal( cbind( GamPos_Bayes$Sol[,1],
				                  ( GamPos_Bayes$Sol[,1] + GamPos_Bayes$Sol[,2] ),
                  				( GamPos_Bayes$Sol[,1] + GamPos_Bayes$Sol[,3] ),
                  				( GamPos_Bayes$Sol[,1] + GamPos_Bayes$Sol[,4] ),
                  				( GamPos_Bayes$Sol[,1] + GamPos_Bayes$Sol[,2] + GamPos_Bayes$Sol[,4] + GamPos_Bayes$Sol[,5] ),
                  				( GamPos_Bayes$Sol[,1] + GamPos_Bayes$Sol[,3] + GamPos_Bayes$Sol[,4] + GamPos_Bayes$Sol[,6] ) ))

colnames( estimatesGP ) <- c( 'A.LH', 'A.M', 'A.S', 'P.LH', 'P.M', 'P.S' )

# tabulate estimates
apply( estimatesGP, 2, mean )  ;	apply( estimatesGP, 2, hpd )
      # A.LH        A.M        A.S       P.LH        P.M        P.S 
# mean  0.0116428  0.09301666 0.1724576 0.17288175  0.0488743  0.0199206 
# -95% -0.1562177 -0.01971555 0.0500691 0.02210829 -0.1666227 -0.3495664
# +95%  0.1984974  0.20876162 0.2964155 0.32136885  0.2645039  0.3624161


Gpfixed <- apply( estimatesGP, 2, mean )
plot(Gpfixed[c(1,3,5)], type="n", xaxt="n", ylab="gamma values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,0.5), main="Gamma (+ve)")
axis(side=1, c("LH","M","S"), at=c(1,2,3))


lines( c(1,1), c(hpd(estimatesGP[,1])), lwd=2 )
lines( c(2,2), c(hpd(estimatesGP[,2])), lwd=2 )
lines( c(3,3), c(hpd(estimatesGP[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(estimatesGP[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(estimatesGP[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(estimatesGP[,6])), lwd=2 )

A.est <- c( mean(estimatesGP[,1]), mean(estimatesGP[,2]), mean(estimatesGP[,3]) )
points( 1:3, A.est, pch=16, cex=2)  #animals - filled points
P.est <- c( mean(estimatesGP[,4]), mean(estimatesGP[,5]), mean(estimatesGP[,6]) )
points( 1.1:3.1, P.est, pch=21, bg="white", cex=2)


## ----modelling_gamma-, results='hide', echo=FALSE, fig.keep='last', warning=FALSE----
select_G <- subset( select, select$Gamma_SE > 0 )

measurement_error_variance <- select_G$Gamma_SE^2

GamNeg_1_Bayes <- MCMCglmm( GamNeg ~ 1, mev= measurement_error_variance, data= select_G, random= ~ StudyID + Species, prior= myprior, nitt=500000 )
summary( GamNeg_1_Bayes )$solutions
#               post.mean    l-95% CI   u-95% CI eff.samp       pMCMC
# (Intercept) -0.05176697 -0.07517503 -0.0304819 780.3846 0.000212766

GN1B<- FoldedNormal( GamNeg_1_Bayes$Sol )
mean( GN1B )	  ; 	hpd( GN1B )     ;   Mode( GN1B )
# 0.06047757    0.03698929 -- 0.08779185     0.06761388


GamNeg_Bayes <- MCMCglmm( abs(GamNeg) ~ 1 + Trait.type * Taxon2, mev= measurement_error_variance, data= select_G, random= ~ StudyID + Species, prior= myprior, nitt=500000, burnin=10000, thin=25 )

plot(GamNeg_Bayes)
acf(GamNeg_Bayes$VCV)
summary(GamNeg_Bayes)$solutions
                       # post.mean     l-95% CI   u-95% CI  eff.samp       pMCMC
# (Intercept)          0.020858281 -0.004494951 0.04822164 14622.373 0.090408163
# Trait.typeM          0.029973906 -0.008515954 0.07572058  1694.011 0.086020408
# Trait.typeS          0.062272146  0.024058278 0.10019579  7955.334 0.003571429
# Taxon2P              0.002935027 -0.030167158 0.04098947  5975.385 0.912244898
# Trait.typeM:Taxon2P -0.023121864 -0.077721519 0.02412212  1579.405 0.326122449


# here I'm fitting a -ve folded normal distribution for the sake of visualizaiton
estimatesGN <- -FoldedNormal( cbind( GamNeg_Bayes$Sol[,1],
                				( GamNeg_Bayes$Sol[,1] + GamNeg_Bayes$Sol[,2] ),
                				( GamNeg_Bayes$Sol[,1] + GamNeg_Bayes$Sol[,3] ),
                				( GamNeg_Bayes$Sol[,1] + GamNeg_Bayes$Sol[,4] ),
                				( GamNeg_Bayes$Sol[,1] + GamNeg_Bayes$Sol[,2] + GamNeg_Bayes$Sol[,4] + GamNeg_Bayes$Sol[,5] ) ))

colnames( estimatesGN ) <- c( 'A.LH', 'A.M', 'A.S', 'P.LH', 'P.M' )

# tabulate estimates
apply( estimatesGN, 2, mean )  ;	apply( estimatesGN, 2, hpd )
       # A.LH         A.M         A.S        P.LH         P.M       
# mean -0.03020001 -0.07364298 -0.09690765 -0.03567132 -0.04385565 
# -95% -0.0690654007 -0.11646180 -0.13642832 -0.073781694 -0.090643135
# +95% -0.0006772858 -0.03229998 -0.05374732 -0.001617268 -0.003421427



Gnfixed <- apply( estimatesGN, 2, mean )
plot(Gnfixed[c(1,3,5)], type="n", xaxt="n", ylab="gamma values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(-0.2,0), main="(abs) Gamma (-ve)")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(estimatesGN[,1])), lwd=2 )
lines( c(2,2), c(hpd(estimatesGN[,2])), lwd=2 )
lines( c(3,3), c(hpd(estimatesGN[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(estimatesGN[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(estimatesGN[,5])), lwd=2 )

A.est <- c( mean(estimatesGN[,1]), mean(estimatesGN[,2]), mean(estimatesGN[,3]) )
points( 1:3, A.est, pch=16, cex=2)  #animals - filled points
P.est <- c( mean(estimatesGN[,4]), mean(estimatesGN[,5]) )
points( 1.1:2.1, P.est, pch=21, bg="white", cex=2)


