
## ----source_functions_&_libraries, echo=FALSE, results='hide', message=FALSE----

rm( list= ls())

# load required packages
require( plyr )   # version 1.8
require( lme4 )   # version 1.16
require( bbmle )  # version 1.0.5.2
require( MCMCglmm )   # version 2.17
require( VGAM )   # version 0.9-1
require( effects )  # version 2.2-4

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


## ----read_in_cov_matrices_to_list, echo=TRUE, eval=FALSE-----------------
## # read_matrices&means
## setwd("..//Data/Gmats_&_means_as_CSVs")
##   matrices <- dir()
##   no.mats <- length(matrices)
##   matrix_list <- list()
## 
## # this loop reads in each matrix from the folder of .csv files and writes them into a list
## for ( i in 1:no.mats ) {
##     matrix_list[[i]] <- read.csv(matrices[i])
##   }
##   names( matrix_list ) <- matrices
## 
## # I shall write this list out as an object in order to make it as simple as possible for readers
## save( matrix_list, file="covariance_matrix_list.R" )


## ----read_in_matrix_list, echo=TRUE, results='asis'----------------------
setwd( "..//Data" )

load( file="covariance_matrix_list.R" )

no.mats <- length( matrix_list )    # building empty variables for each metric
END <- numeric(no.mats)
TGV <- numeric(no.mats)
Emax <- numeric(no.mats)
gmax <- numeric(no.mats)
AEvo <- numeric(no.mats)
Even <- numeric(no.mats)

# this loop reads each matrix, calculates each metric and writes them to the empty vectors above
for (i in 1:no.mats) {
  matrix.and.means <- as.matrix( matrix_list[[i]] )
  means <- matrix.and.means[,1]
  matrix <- matrix.and.means[,2:ncol(matrix.and.means)]
  if ( sum(means)==0 ) { st.matrix <- matrix }
  else { st.matrix <- StanCovMat( as.matrix(matrix), means ) }
  diag.matrix <- svd( st.matrix )
  eigen.values <- diag.matrix$d
  TGV[i] <- sum(eigen.values)   # 'total genetic varaince' after Kirkpatrick 2008
  END[i] <- sum(eigen.values / eigen.values[1])   # 'effective no. dimensions' after Kirkpatrick 2008
  Emax[i] <- sqrt(eigen.values[1])    # 'maximum evolvability' after Kirkpatrick 2008
  gmax[i] <- eigen.values[1]      # storing gmax; the first eigenvalue
  AEvo[i] <- mean( eigen.values )     # 'average evolvability' after Hansen & Houle 2008
  lamda_tilde <- abs(eigen.values) / sum( abs(eigen.values) )
  Even[i] <- - sum( lamda_tilde*log(lamda_tilde) / log( length(eigen.values) ) )  # 'eigenvalue eveness' after Agrawal & Stinchcombe 2008 (though note that they intended this metric to be calculated from correlation matrices. see below)
}

Metrics <- data.frame( TGV, END, Emax, gmax, AEvo, Even ) 
matrices_names <- sub( ".csv", "", names( matrix_list ))
Metrics$filename <- factor( matrices_names )
# str(Metrics)

# now to read in the csv that holds the ID table
Index <- read.csv("MatrixIndexFinal.csv")
# str(Index)

# merge the metrics with the index
gmats <- merge(Metrics, Index, "filename")
# str(gmats)

# dump columns we don't need for analysis
gmats <- gmats[, c(-10, -17)]


## ----data_hygiene, echo=FALSE, results='hide', fig.keep='none'-----------
# remove 2 matrices from Sherrard.et.al.2009 & 1 from House&Simmons because values seem unbelievable
gmats <- gmats[ -c(59,83,84),]

plot( gmats$AEvo ) ; # identify( gmats$AEvo )
gmats[c(27,28,29,30,68,69,70),]
# these unusually large AEvo/Emax values are from 2 studies: Blo2003 & Mer1996
# Mer1996 values are correct as published
# Blo2003 excluded because I can't work out exactly what they did (Blows!)
gmats <- gmats[ gmats$study.code !="Blo2003",]

# here we're going to exclude those matrices whose component traits comprise awkward units e.g. mm^2 or cm^3
# NB: after all the previous data fisking this removes precisely nothing
gmats[ gmats$problem.units =="Y" ,]

dim( gmats )  # 81 matrices after PDF


## ----tabulate_data_for_tables, echo=FALSE--------------------------------
ddply( gmats, .(taxon2, trait.type), summarise, m_TGV= round(mean(TGV), 3), m_END= round(mean(END), 3), m_Emax= round(mean(Emax), 3), m_Even= round(mean(Even), 3), m_AEvo= round(mean(AEvo), 3), m_Trait.no= round( mean(trait.no), 2), no.rows= length(TGV), m_no.families= round( mean(na.omit( no.families)), 1) )

# information for metric correlation table
round( cor( gmats[,c(2:7, 14)] ), 2)

round( cor( na.omit( gmats[,c(2:7, 14, 18)] )), 2)


## ----cov_pairs_plot, echo=FALSE, results='hide', message=FALSE, warning=FALSE, fig.height=7, fig.height=7----
# pairs figure for covariance matrices
names( gmats[,c(2:7, 14)] )

metric.names <- expression( italic("tgv"),
                            italic("n"["D"]),
                            italic("e"["max"]),
                            italic(bolditalic("g")["max"]),
                            italic("\u0113"),
                            italic("E"["\u03bb"]),
                            italic("n") )

pairs( gmats[,c(2:7, 14)], cex.labels=2, labels=metric.names, lower.panel=panel.cor )


## ----model_selection_for_END, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

# first we're looking at the Effective Number of Dimensions

END.reml1 <- lmer(END ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)

END.reml2 <- lmer(END ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)

END.reml3 <- lmer(END ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

END.reml4 <- lmer(END ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

BICtab( END.reml1, END.reml2, END.reml3, END.reml4, weights=T, nobs=nrow(gmats) )

# this suggests that we don't need the interaction!?
# model 1 is the best fitting one to include trait.no
# dBIC df weight 
# END.reml1  0.0 7  0.83534
# END.reml3  4.4 8  0.09282
# END.reml2  5.1 9  0.06466
# END.reml4  9.5 10 0.00718


ENDt <- with( gmats, END/trait.no )
END.reml5 <- lmer( ENDt ~ 1 + trait.type + taxon2 + (1|study.code), data=gmats)
END.reml5.1 <- lmer(ENDt ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
END.reml5.2 <- lmer(ENDt ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
END.reml5.3 <- lmer(ENDt ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
	BICtab( END.reml5, END.reml5.1, END.reml5.2, END.reml5.3, weights=T, nobs=nrow(gmats) )
#		            dBIC df weight 
#		END.reml5.2  0.0 8  0.79911
#		END.reml5.1  3.1 9  0.17286
#		END.reml5.3  7.5 10 0.01921
#		END.reml5    9.0 6  0.00882


ENDt2 <- with( gmats, END/(trait.no^2) )
END.reml6 <- lmer( ENDt2 ~ 1 + trait.type + taxon2 + (1|study.code), data=gmats)
END.reml6.1 <- lmer(ENDt2 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
END.reml6.2 <- lmer(ENDt2 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
END.reml6.3 <- lmer(ENDt2 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
	BICtab( END.reml6, END.reml6.1, END.reml6.2, END.reml6.3, weights=T, nobs=nrow(gmats) )
#	            dBIC df weight
#	END.reml6.2  0.0 8  0.918 
#	END.reml6.1  4.9 9  0.081 
#	END.reml6   14.4 6  <0.001
#	END.reml6.3 14.9 10 <0.001


# to what extent do we need the trait.no standardization as opposed to modelling trait.no as a covaraite?
cor( summary( END.reml5 )$coefficients[,1], summary( END.reml6 )$coefficients[,1] )
cor( summary( END.reml1 )$coefficients[1:4,1], summary( END.reml6 )$coefficients[,1] )
cor( summary( END.reml1 )$coefficients[1:4,1], summary( END.reml5 )$coefficients[,1] )
# estimates from all these 3 models are correlated at >0.97 so I think we are good with any of them.

# now I'll use the parametric bootstrap to confirm...

# model 1 vs model 2
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( END.reml1 ))
LR.model <- -as.numeric(deviance( END.reml2 ) - deviance( END.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.1558 this does not suggest stepping up to model 2...

# model 1 vs model 3
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( END.reml1 ))
LR.model <- -as.numeric(deviance( END.reml3 ) - deviance( END.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.7332 supports sticking with model 1 over model 3





myprior1 <- list(R = list(V =1, nu = 0.002), G= list( G1= list (V=1, nu=0.002)))

myprior2 <- list(R = list(V =1, nu = 0.002), G= list( G1= list (V=1, nu=0.002), G2= list (V=1, nu=0.002)))

myprior3 <- list( list(V=1, nu=0.002), list( G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002) ) )



END.Bayes0 <- MCMCglmm(END ~ 1, nitt=1000000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )
summary(END.Bayes0)$solutions
# post.mean l-95% CI u-95% CI eff.samp        pMCMC
# (Intercept)  1.532492 1.391107 1.673519 34391.33 3.015045e-05

ENDB0 <- rfnorm( dim(END.Bayes0$Sol)[1], mean(END.Bayes0$Sol[,1]), sd(END.Bayes0$Sol[,1]) )
mean( ENDB0 )  ;  hpd( ENDB0 ) ; Mode( ENDB0 )
# 1.531906    1.388636 -- 1.673171		1.503468


END.Bayes1 <- MCMCglmm(END ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

END.Bayes2 <- MCMCglmm(END ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

END.Bayes3 <- MCMCglmm(END ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

END.Bayes4 <- MCMCglmm(END ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# DIC agrees with BIC agrees with LRsim on model 1
END.Bayes1$DIC	;	END.Bayes2$DIC ; END.Bayes3$DIC ; END.Bayes4$DIC
# 	50.77			51.80				61.75			64.36



END.Bayes5 <- MCMCglmm(ENDt ~ 1 + trait.type + taxon2, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )
END.Bayes5.1 <- MCMCglmm(ENDt ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )
END.Bayes5.2 <- MCMCglmm(ENDt ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
END.Bayes5.3 <- MCMCglmm(ENDt ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
	END.Bayes5$DIC	;	END.Bayes5.1$DIC	;	END.Bayes5.2$DIC	;	END.Bayes5.3$DIC


END.Bayes6 <- MCMCglmm(ENDt2 ~ 1 + trait.type + taxon2, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )
END.Bayes6.1 <- MCMCglmm(ENDt2 ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )
END.Bayes6.2 <- MCMCglmm(ENDt2 ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
END.Bayes6.3 <- MCMCglmm(ENDt2 ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
	END.Bayes6$DIC	;	END.Bayes6.1$DIC	;	END.Bayes6.2$DIC	;	END.Bayes6.3$DIC


# as with the REML estimates, these fit better but...
END.Bayes5$DIC  ;	END.Bayes6$DIC
# -185.51       ;   -357.03
# ...once again the estimates are all correlated at >0.97 so using model1 makes sense for interpretability
cor( colMeans( END.Bayes5$Sol ), colMeans( END.Bayes6$Sol ) )
cor( colMeans( END.Bayes5$Sol ), colMeans( END.Bayes1$Sol )[1:4] )
cor( colMeans( END.Bayes6$Sol ), colMeans( END.Bayes1$Sol )[1:4] )


plot(END.Bayes1$VCV)
acf(END.Bayes1$VCV)
plot(END.Bayes1$Sol)
acf(END.Bayes1$Sol)
summary(END.Bayes1)$solutions
# post.mean     l-95% CI  u-95% CI eff.samp        pMCMC
# (Intercept) 1.14205496  0.765377638 1.5201990     3167 0.0003157562
# trait.typeM 0.06254670 -0.287247203 0.4062500     3167 0.7224502684
# trait.typeS 0.16213894 -0.207232412 0.5740124     3167 0.4180612567
# taxon2P     0.08734098 -0.273186713 0.4561255     3167 0.6226712978
# trait.no    0.05129427 -0.005488227 0.1102218     3167 0.0991474582


## ----END_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----
END.estimates <- cbind( END.Bayes1$Sol[,1],
                        ( END.Bayes1$Sol[,1] + END.Bayes1$Sol[,2] ),
                        ( END.Bayes1$Sol[,1] + END.Bayes1$Sol[,3] ),
                        ( END.Bayes1$Sol[,1] + END.Bayes1$Sol[,4] ),
                        ( END.Bayes1$Sol[,1] + END.Bayes1$Sol[,2] + END.Bayes1$Sol[,4] ),
                        ( END.Bayes1$Sol[,1] + END.Bayes1$Sol[,3] + END.Bayes1$Sol[,4] ) )

# END cannot be -ve, so we need to apply a folded normal distribution to our posteriors
END.estimates <- FoldedNormal( END.estimates )

# Tabulate estimates
apply( END.estimates, 2, mean )
apply( END.estimates, 2, Mode )
apply( END.estimates, 2, hpd )
# A.LH       A.M       A.S      P.LH       P.M       P.S 
# mean 1.133488 1.200328 1.308200 1.227267 1.293960 1.388670
# mode 1.3978042 1.1902635 2.0576041 1.2795317 0.9810789 1.8061362
# -95% 0.7560626 0.8177292 0.8472309 0.7720113 0.8598855 0.9942053
# +95% 1.5045455 1.5379631 1.8155311 1.6481765 1.7099160 1.8107346


# END plot
ENDfixed <- apply( END.estimates, 2, mean )
plot(ENDfixed, type="n", xaxt="n", ylab="END values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,2), main="Effective Dimensionality")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(END.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(END.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(END.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(END.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(END.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(END.estimates[,6])), lwd=2 )

points( 1:3, ENDfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, ENDfixed[4:6], pch=21, bg="white", cex=2)



## ----comparing_nD_scaling_options, echo=FALSE, fig.keep='last', message=FALSE----

END.est.5 <- FoldedNormal( cbind( END.Bayes5.1$Sol[,1],
	                	( END.Bayes5.1$Sol[,1] + END.Bayes5.1$Sol[,2] ),
    	                ( END.Bayes5.1$Sol[,1] + END.Bayes5.1$Sol[,3] ),
        	            ( END.Bayes5.1$Sol[,1] + END.Bayes5.1$Sol[,4] ),
            	        ( END.Bayes5.1$Sol[,1] + END.Bayes5.1$Sol[,2] + END.Bayes5.1$Sol[,4] + END.Bayes5.1$Sol[,6] ),
                	    ( END.Bayes5.1$Sol[,1] + END.Bayes5.1$Sol[,3] + END.Bayes5.1$Sol[,4] + END.Bayes5.1$Sol[,7] ) ) )

END.est.6 <- FoldedNormal( cbind( END.Bayes6.1$Sol[,1],
	                	( END.Bayes6.1$Sol[,1] + END.Bayes6.1$Sol[,2] ),
    	                ( END.Bayes6.1$Sol[,1] + END.Bayes6.1$Sol[,3] ),
        	            ( END.Bayes6.1$Sol[,1] + END.Bayes6.1$Sol[,4] ),
            	        ( END.Bayes6.1$Sol[,1] + END.Bayes6.1$Sol[,2] + END.Bayes6.1$Sol[,4] + END.Bayes6.1$Sol[,6] ),
                	    ( END.Bayes6.1$Sol[,1] + END.Bayes6.1$Sol[,3] + END.Bayes6.1$Sol[,4] + END.Bayes6.1$Sol[,7] ) ) )

par( mfrow=c(1, 3) )

plot( apply( END.estimates, 2, mean ), type="n", xaxt="n", ylab=expression(italic("n"["D"])), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,2), main="END with k as covariate")
		axis(side=1, c("LH","M","S"), at=c(1,2,3))
		lines( c(1,1), c(hpd(END.estimates[,1])), lwd=2 )
		lines( c(2,2), c(hpd(END.estimates[,2])), lwd=2 )
		lines( c(3,3), c(hpd(END.estimates[,3])), lwd=2 )
		lines( c(1.1,1.1), c(hpd(END.estimates[,4])), lwd=2 )
		lines( c(2.1,2.1), c(hpd(END.estimates[,5])), lwd=2 )
		lines( c(3.1,3.1), c(hpd(END.estimates[,6])), lwd=2 )
		points( 1:3, ENDfixed[1:3], pch=16, cex=2)
		#animals - filled points ^ #plants - open points v
		points( 1.1:3.1, ENDfixed[4:6], pch=21, bg="white", cex=2)

plot( apply( END.est.5, 2, mean ), type="n", xaxt="n", ylab=expression(italic("scaled n"["D"])), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,0.7), main="END / k")
		axis(side=1, c("LH","M","S"), at=c(1,2,3))
		lines( c(1,1), c(hpd(END.est.5[,1])), lwd=2 )
		lines( c(2,2), c(hpd(END.est.5[,2])), lwd=2 )
		lines( c(3,3), c(hpd(END.est.5[,3])), lwd=2 )
		lines( c(1.1,1.1), c(hpd(END.est.5[,4])), lwd=2 )
		lines( c(2.1,2.1), c(hpd(END.est.5[,5])), lwd=2 )
		lines( c(3.1,3.1), c(hpd(END.est.5[,6])), lwd=2 )
		points( 1:3, apply(END.est.5[,1:3], 2, mean), pch=16, cex=2)
		#animals - filled points ^ #plants - open points v
		points( 1.1:3.1, apply(END.est.5[,4:6], 2, mean), pch=21, bg="white", cex=2)

plot( apply( END.est.6, 2, mean ), type="n", xaxt="n", ylab=expression(italic("scaled n"["D"])), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,0.4), main="END / k^2")
		axis(side=1, c("LH","M","S"), at=c(1,2,3))
		lines( c(1,1), c(hpd(END.est.6[,1])), lwd=2 )
		lines( c(2,2), c(hpd(END.est.6[,2])), lwd=2 )
		lines( c(3,3), c(hpd(END.est.6[,3])), lwd=2 )
		lines( c(1.1,1.1), c(hpd(END.est.6[,4])), lwd=2 )
		lines( c(2.1,2.1), c(hpd(END.est.6[,5])), lwd=2 )
		lines( c(3.1,3.1), c(hpd(END.est.6[,6])), lwd=2 )
		points( 1:3, apply(END.est.6[,1:3], 2, mean), pch=16, cex=2)
		#animals - filled points ^ #plants - open points v
		points( 1.1:3.1, apply(END.est.6[,4:6], 2, mean), pch=21, bg="white", cex=2)

par( mfrow=c(1,1))


## ----model_selection_for_Emax, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----
# now to analyse Max. Evolvability

EMAX.reml1 <- lmer(Emax ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)

EMAX.reml2 <- lmer(Emax ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)

EMAX.reml3 <- lmer(Emax ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

EMAX.reml4 <- lmer(Emax ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

BICtab( EMAX.reml1, EMAX.reml2, EMAX.reml3, EMAX.reml4, weights=T, nobs=nrow(gmats) )

# this suggests that we don't need the interaction!?
# dBIC df weight
# EMAX.reml1  0.0 7  0.6194
# EMAX.reml2  1.6 9  0.2806
# EMAX.reml3  4.4 8  0.0688
# EMAX.reml4  6.0 10 0.0312

# now I'll use the parametric bootstrap to confirm...

# model 1 vs model 2
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( EMAX.reml1 ))
LR.model <- -as.numeric(deviance( EMAX.reml2 ) - deviance( EMAX.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.705 supports sticking with model 1 over model 2


# model 1 vs model 3
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( EMAX.reml1 ))
LR.model <- -as.numeric(deviance( EMAX.reml3 ) - deviance( EMAX.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.423 supports sticking with model 1 over model 3






EMAX.Bayes0 <- MCMCglmm(Emax ~ 1, nitt=100000, burnin=10000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
summary( EMAX.Bayes0 )$solutions
# post.mean  l-95% CI u-95% CI eff.samp        pMCMC
# (Intercept) 0.6067565 0.2490954 0.962691     3000 0.0003333333

EMAXB0 <- rfnorm(3000, mean(EMAX.Bayes0$Sol[,1]), sd(EMAX.Bayes0$Sol[,1]))
mean( EMAXB0 )	;	hpd( EMAXB0 ) 	 ;	  Mode( EMAXB0 )
# 0.6112283		0.2587635 -- 0.9674873		0.4653126


EMAX.Bayes1 <- MCMCglmm(Emax ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

EMAX.Bayes2 <- MCMCglmm(Emax ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

EMAX.Bayes3 <- MCMCglmm(Emax ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

EMAX.Bayes4 <- MCMCglmm(Emax ~ 1 + trait.type * taxon2 + trait.no, nitt=500000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )


# DIC is in disagreement with BIC and LRsim... I going to err on the side of caution and go with model3
EMAX.Bayes1$DIC ; EMAX.Bayes2$DIC ; EMAX.Bayes3$DIC	; EMAX.Bayes4$DIC
# 287.4724			287.9171			268.0664		268.3305


plot(EMAX.Bayes3$VCV)
acf(EMAX.Bayes3$VCV)
plot(EMAX.Bayes3$Sol)
acf(EMAX.Bayes3$Sol)
summary(EMAX.Bayes3)$solutions
# post.mean   l-95% CI  u-95% CI eff.samp     pMCMC
# (Intercept)  0.09065319 -0.9488144 1.0829023 3167.000 0.8601200
# trait.typeM  1.13577784 -0.5354211 3.0375432 3167.000 0.1957689
# trait.typeS  0.57123654 -0.6798971 1.7049119 3167.000 0.3328071
# taxon2P     -0.64535373 -1.6896664 0.3090051 2981.949 0.1951374
# trait.no     0.06964617 -0.1448010 0.2947385 3167.000 0.5424692


## ----Emax_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----
EMAX.estimates <- cbind( EMAX.Bayes3$Sol[,1],
                         ( EMAX.Bayes3$Sol[,1] + EMAX.Bayes3$Sol[,2] ),
                         ( EMAX.Bayes3$Sol[,1] + EMAX.Bayes3$Sol[,3] ),
                         ( EMAX.Bayes3$Sol[,1] + EMAX.Bayes3$Sol[,4] ),
                         ( EMAX.Bayes3$Sol[,1] + EMAX.Bayes3$Sol[,2] + EMAX.Bayes3$Sol[,4] ),
                         ( EMAX.Bayes3$Sol[,1] + EMAX.Bayes3$Sol[,3] + EMAX.Bayes3$Sol[,4] ) )

colnames(EMAX.estimates) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

EMAX.estimates <- FoldedNormal( EMAX.estimates )

apply( EMAX.estimates, 2, mean )
apply( EMAX.estimates, 2, Mode )
apply( EMAX.estimates, 2, hpd )
# A.LH       A.M       A.S      P.LH       P.M        P.S 
# mean 0.4334804 1.2491385 0.7823528 0.5888987 0.9102743 0.5250986
# mode 0.2579401 0.0065128 0.8597651 0.8885794 0.5685091 0.1098872
# -95% 0.0003335 0.0007561 0.0011306 0.0001440 0.0028559 0.0002239
# +95% 1.0429934 2.8338943 1.8755535 1.3894593 2.2157249 1.2813908


# EMAX plot
EMAXfixed <- apply( EMAX.estimates, 2, mean )
plot(EMAXfixed, type="n", xaxt="n", ylab="EMAX values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 5), main="Max. Evolvability")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(EMAX.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(EMAX.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(EMAX.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(EMAX.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(EMAX.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(EMAX.estimates[,6])), lwd=2 )

points( 1:3, EMAXfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, EMAXfixed[4:6], pch=21, bg="white", cex=2)



## ----model_selection_for_TGV, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

TGV.reml1 <- lmer(TGV ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)

TGV.reml2 <- lmer(TGV ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)

TGV.reml3 <- lmer(TGV ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

TGV.reml4 <- lmer(TGV ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

# model2 is the best-fitting with REML and BIC
BICtab( TGV.reml1, TGV.reml2, TGV.reml3, TGV.reml4, weights=T, nobs=nrow(gmats) )
# dBIC df weight 
# TGV.reml2  0.0 9  0.89568
# TGV.reml4  4.4 10 0.09952
# TGV.reml1 10.7 7  0.00432
# TGV.reml3 15.1 8  < 0.001

# now the parametric bootstrap

# model 1 vs model 2
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( TGV.reml1 ))
LR.model <- -as.numeric(deviance( TGV.reml2 ) - deviance( TGV.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.776 bootstrap doesn't support model2 over model1 in contrast with BIC...


# model 2 vs model 4
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( TGV.reml2 ))
LR.model <- -as.numeric(deviance( TGV.reml4 ) - deviance( TGV.reml2 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.642 stay with model 2





TGV.Bayes0 <- MCMCglmm(TGV ~ 1, nitt=500000, burnin=1000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
summary( TGV.Bayes0 )$solutions

TGVB0 <- rfnorm( length( TGV.Bayes0$Sol[,1] ), mean( TGV.Bayes0$Sol[,1] ), sd( TGV.Bayes0$Sol[,1] ) )
mean( TGVB0 )	;	hpd( TGVB0 ) 		 ;	 Mode( TGVB0 )
# 3.135937		0.001847855 -to- 7.566797628		2.940357


TGV.Bayes1 <- MCMCglmm(TGV ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

TGV.Bayes2 <- MCMCglmm(TGV ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

TGV.Bayes3 <- MCMCglmm(TGV ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

TGV.Bayes4 <- MCMCglmm(TGV ~ 1 + trait.type * taxon2 + trait.no, nitt=500000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# DIC suggests that model3 is a better fit than model2
TGV.Bayes1$DIC ; TGV.Bayes2$DIC ; TGV.Bayes3$DIC	;	TGV.Bayes4$DIC
# 776.5418			777.4574			740.8967			741.07

# the shared coef.s are highly correlated (0.95), so I'm going to err on the side of complexity and go with model3
cor( c(colMeans( TGV.Bayes2$Sol )[1:5]), c(colMeans( TGV.Bayes3$Sol )) )


plot(TGV.Bayes3$VCV)
acf(TGV.Bayes3$VCV)
plot(TGV.Bayes3$Sol)
acf(TGV.Bayes3$Sol)
summary(TGV.Bayes3)$solutions
# post.mean   l-95% CI  u-95% CI eff.samp     pMCMC
# (Intercept) -0.4570179 -19.554800 20.442918 2791.130 0.9573729
# trait.typeM 20.2472241 -17.633004 63.580113 2981.970 0.3151247
# trait.typeS  2.7124723 -18.610255 27.785068 3598.117 0.8317019
# taxon2P     -2.9255618 -22.614285 16.022980 3004.807 0.7786549
# trait.no     0.2375523  -4.358424  4.294138 3167.000 0.9131670

TGV.estimates <- cbind( TGV.Bayes3$Sol[,1],
                        ( TGV.Bayes3$Sol[,1] + TGV.Bayes3$Sol[,2] ),
                        ( TGV.Bayes3$Sol[,1] + TGV.Bayes3$Sol[,3] ),
                        ( TGV.Bayes3$Sol[,1] + TGV.Bayes3$Sol[,4] ),
                        ( TGV.Bayes3$Sol[,1] +TGV.Bayes3$Sol[,2] +TGV.Bayes3$Sol[,4] ),
                        ( TGV.Bayes3$Sol[,1] +TGV.Bayes3$Sol[,3] +TGV.Bayes3$Sol[,4] ) )

colnames(TGV.estimates) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

TGV.estimates <- FoldedNormal( TGV.estimates )

# Tabulate estimates
apply( TGV.estimates, 2, mean )
apply( TGV.estimates, 2, Mode )
apply( TGV.estimates, 2, hpd )
# A.LH       A.M       A.S      P.LH       P.M       P.S 
# mean 8.622775 25.408110 14.317302  9.383316 24.088053 11.043411
# mode 17.15973 23.57304 21.19176   17.36984  39.67697  10.85037
# -95%  0.006315  0.2961  0.003593  0.008477  0.008882  0.0068360
# +95% 18.113191 52.4658 30.889906 19.636975 50.480431 22.8630851



## ----TGV_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----

TGVfixed <- apply( TGV.estimates, 2, mean )
plot(TGVfixed, type="n", xaxt="n", ylab="TGV values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,75), main="Total Genetic Variance")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(TGV.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(TGV.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(TGV.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(TGV.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(TGV.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(TGV.estimates[,6])), lwd=2 )

points( 1:3, TGVfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, TGVfixed[4:6], pch=21, bg="white", cex=2)


## ----model_selection_Gmax, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

gmax.reml1 <- lmer(gmax ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)

gmax.reml2 <- lmer(gmax ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)

gmax.reml3 <- lmer(gmax ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

gmax.reml4 <- lmer(gmax ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

# BIC reccommends model2
BICtab( gmax.reml1, gmax.reml2, gmax.reml3, gmax.reml4, weights=T, nobs=nrow(gmats) )
# dBIC df weight 
# gmats.reml2  0.0 9  0.89211
# gmats.reml4  4.4 10 0.09912
# gmats.reml1  9.5 7  0.00789
# gmats.reml3 13.8 8  < 0.001


# model 2 vs model 4
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( gmax.reml2 ))
LR.model <- -as.numeric(deviance( gmax.reml4 ) - deviance( gmax.reml2 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.735 stay with model 2



gmax.Bayes1 <- MCMCglmm(gmax ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

gmax.Bayes2 <- MCMCglmm(gmax ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

gmax.Bayes3 <- MCMCglmm(gmax ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

gmax.Bayes4 <- MCMCglmm(gmax ~ 1 + trait.type * taxon2 + trait.no, nitt=500000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# Bayesian approach supports model 3 
gmax.Bayes1$DIC ; gmax.Bayes2$DIC ; gmax.Bayes3$DIC  ;	gmax.Bayes4$DIC
# 801.9224			804.683				774.3294			777.2449

# using model3 gives estimates that are correlated at >0.94 with model2, so I'm going with model3
cor( c(colMeans( gmax.Bayes2$Sol )[1:5]), c(colMeans( gmax.Bayes3$Sol )) )


plot(gmax.Bayes3$VCV)
acf(gmax.Bayes3$VCV)
plot(gmax.Bayes3$Sol)
acf(gmax.Bayes3$Sol)
summary(gmax.Bayes3)$solutions
#              post.mean   l-95% CI  u-95% CI eff.samp     pMCMC
# (Intercept) -2.1772503 -24.356400 19.024930 3167.000 0.8487528
# trait.typeM 12.2396237 -15.523725 40.112034 3562.468 0.3808020
# trait.typeS  0.9493931 -23.160816 24.021480 3262.508 0.9340069
# taxon2P     -1.7168857 -20.590525 18.684758 2818.922 0.8658036
# trait.no     0.6474424  -3.782605  5.165388 2925.330 0.7843385


gmax.estimates <- cbind( gmax.Bayes3$Sol[,1],
                         ( gmax.Bayes3$Sol[,1] + gmax.Bayes3$Sol[,2] ),
                         ( gmax.Bayes3$Sol[,1] + gmax.Bayes3$Sol[,3] ),
                         ( gmax.Bayes3$Sol[,1] + gmax.Bayes3$Sol[,4] ),
                         ( gmax.Bayes3$Sol[,1] +gmax.Bayes3$Sol[,2] +gmax.Bayes3$Sol[,4] ),
                         ( gmax.Bayes3$Sol[,1] +gmax.Bayes3$Sol[,3] +gmax.Bayes3$Sol[,4] ) )

gmax.estimates <- FoldedNormal( gmax.estimates )

colnames(gmax.estimates) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

# Tabulate estimates
apply( gmax.estimates, 2, mean )
apply( gmax.estimates, 2, Mode )
apply( gmax.estimates, 2, hpd )
# A.LH       A.M       A.S      P.LH       P.M       P.S 
# 9.071529 17.573916 14.774674 22.417751 29.892646 11.165447 
# 23.761995 11.810164 19.365173  6.093903 57.117965  4.838458 
# 3.2708e-04  0.003470 2.0487e-04 4.5605e-04  0.00434 4.8361e-04
# 2.2386e+01 42.566804 3.6349e+01 5.4547e+01 73.60713 2.7470e+01


## ----gmax_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----

gmaxfixed <- apply( gmax.estimates, 2, mean )
plot(gmaxfixed, type="n", xaxt="n", ylab="gmax values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,50), main="Gmax (1st eigenvalue)")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(gmax.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(gmax.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(gmax.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(gmax.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(gmax.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(gmax.estimates[,6])), lwd=2 )

points( 1:3, gmaxfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, gmaxfixed[4:6], pch=21, bg="white", cex=2)


## ----average_evolvability_model_selection, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

AEvo.reml1 <- lmer(AEvo ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)

AEvo.reml2 <- lmer(AEvo ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)

AEvo.reml3 <- lmer(AEvo ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

AEvo.reml4 <- lmer(AEvo ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

# BIC comes out in favour of model2 again...
BICtab( AEvo.reml1, AEvo.reml2, AEvo.reml3, AEvo.reml4, weights=T, nobs=81 )
# dBIC df weight
# AEvo.reml2  0.0 9  0.7300
# AEvo.reml1  2.9 7  0.1700
# AEvo.reml4  4.4 10 0.0811
# AEvo.reml3  7.3 8  0.0189


# model 1 vs model 2
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( AEvo.reml1 ))
LR.model <- -as.numeric(deviance( AEvo.reml2 ) - deviance( AEvo.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.769 nothing in it between model 2 & model 1  <-  this is not in agreement with BIC...


# model 2 vs model 4
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( AEvo.reml2 ))
LR.model <- -as.numeric(deviance( AEvo.reml4 ) - deviance( AEvo.reml2 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.505 nothing in it between model 2 and 4 either...


# Bayesian version
AEvo.Bayes0 <- MCMCglmm(AEvo ~ 1, nitt=500000, burnin=1000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )
summary( AEvo.Bayes0 )$solutions
# post.mean   l-95% CI u-95% CI eff.samp    pMCMC
# (Intercept) 0.3355777 -0.6093231 1.283184    16634 0.488758

AEvoB0 <- rfnorm( length( AEvo.Bayes0$Sol[,1] ), mean( AEvo.Bayes0$Sol[,1] ), sd( AEvo.Bayes0$Sol[,1] ) )
mean( AEvoB0 )	;	hpd( AEvoB0 ) 		 ;		 Mode( AEvoB0 )
# 0.4775428		5.870392e-05 -to- 1.148672e+00		0.65104


AEvo.Bayes1 <- MCMCglmm(AEvo ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

AEvo.Bayes2 <- MCMCglmm(AEvo ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

AEvo.Bayes3 <- MCMCglmm(AEvo ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

AEvo.Bayes4 <- MCMCglmm(AEvo ~ 1 + trait.type * taxon2 + trait.no, nitt=500000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# Bayesian approach indicates model3 is the best fitting
AEvo.Bayes1$DIC ; AEvo.Bayes2$DIC ; AEvo.Bayes3$DIC	;	AEvo.Bayes4$DIC
# 471.9372			473.4084			434.9815			437.3735

# the choice between model2 & model3 again makes little quantitative difference <- estimates correlated at >0.94
cor( c(colMeans( AEvo.Bayes2$Sol )[1:5]), c(colMeans( AEvo.Bayes3$Sol )) )

# going with model3, erring on the side of the more complex models with additional covariates. 
summary( AEvo.Bayes3 )$solutions
# post.mean   l-95% CI  u-95% CI eff.samp     pMCMC
# (Intercept)  0.02891164 -2.8527564 2.9017820     3167 0.9838964
# trait.typeM  3.08145461 -2.5856689 8.9037925     3167 0.2854436
# trait.typeS  0.48952770 -2.6076508 3.4216634     3167 0.7641301
# taxon2P     -0.48971460 -3.0229550 1.9963982     3167 0.7199242
# trait.no     0.01663702 -0.6234972 0.6099593     3167 0.9718977

AEvo.estimates <- cbind( AEvo.Bayes3$Sol[,1],
                         ( AEvo.Bayes3$Sol[,1] + AEvo.Bayes3$Sol[,2] ),
                         ( AEvo.Bayes3$Sol[,1] + AEvo.Bayes3$Sol[,3] ),
                         ( AEvo.Bayes3$Sol[,1] + AEvo.Bayes3$Sol[,4] ),
                         ( AEvo.Bayes3$Sol[,1] +AEvo.Bayes3$Sol[,2] +AEvo.Bayes3$Sol[,4] ),
                         ( AEvo.Bayes3$Sol[,1] +AEvo.Bayes3$Sol[,3] +AEvo.Bayes3$Sol[,4] ) )

AEvo.estimates <- FoldedNormal( AEvo.estimates )

colnames(AEvo.estimates) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

# Tabulate estimates
apply( AEvo.estimates, 2, mean )
apply( AEvo.estimates, 2, Mode )
apply( AEvo.estimates, 2, hpd )
# A.LH       A.M       A.S      P.LH       P.M       P.S 
# 1.171787 3.692067 1.948575 1.292362 3.402012 1.468000 
# 0.5468495 1.7414221 1.2911786 0.7985577 0.6695847 1.6241952 
# 0.0001642 0.0062499 0.0007516 0.002241 0.001609879 0.000541
# 2.9434666 8.5480446 4.7517676 3.198912 7.825105986 3.618291


## ----AEvo_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----
# AEvo plot
AEvofixed <- apply( AEvo.estimates, 2, mean )
plot(AEvofixed, type="n", xaxt="n", ylab="AEvo values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 12), main="Average Evolvability")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(AEvo.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(AEvo.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(AEvo.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(AEvo.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(AEvo.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(AEvo.estimates[,6])), lwd=2 )

points( 1:3, AEvofixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, AEvofixed[4:6], pch=21, bg="white", cex=2)


## ----AEvo_vs_TGV_plot, results='hide', fig.height=4, fig.width=9, echo=FALSE----
# AEvo versus TGV
par( mfrow= c(1,2) )

mycor <- cor( log10( gmats$TGV), log10( gmats$AEvo) )	# 0.987
plot( log10( gmats$TGV), log10( gmats$AEvo), xlab="(log) total gen. var.", ylab="(log) average evo." )
abline( lm( log10( gmats$AEvo) ~ log10( gmats$TGV) ), col="red" )
text( -1, 1, paste( "cor=", round(mycor, 3), sep='', coll='' ) )

mycor2 <- cor( gmats$TGV, gmats$AEvo)	# 0.999
plot( gmats$TGV, gmats$AEvo, xlab="total gen. var.", ylab="average evo." )
abline( lm( gmats$AEvo ~ gmats$TGV ), col="red" )
text( 100, 55, paste( "cor=", round(mycor2, 3), sep='', coll='' ) )

mtext( text="TGV & AEvo are highly correlated", side=3, line=2, at=-200 )

par( mfrow=c(1,1))


## ----model_selection_for_eveness, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

Even.reml1 <- lmer(Even ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)

Even.reml2 <- lmer(Even ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=gmats)

Even.reml3 <- lmer(Even ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

Even.reml4 <- lmer(Even ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)

# clear support from BIC for model1
BICtab( Even.reml1, Even.reml2, Even.reml3, Even.reml4, weights=T, nobs=nrow(gmats) )
# dBIC df weight 
# Even.reml1  0.0 7  0.85936
# Even.reml3  4.4 8  0.09557
# Even.reml2  6.1 9  0.04057
# Even.reml4 10.5 10 0.00451


# model 1 vs model 3
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=gmats)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=gmats)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( Even.reml1 ))
LR.model <- -as.numeric(deviance( Even.reml3 ) - deviance( Even.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.425 no support for stepping up to model 3


Even.Bayes1 <- MCMCglmm(Even ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

Even.Bayes2 <- MCMCglmm(Even ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ study.code, prior= myprior1 )

Even.Bayes3 <- MCMCglmm(Even ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

Even.Bayes4 <- MCMCglmm(Even ~ 1 + trait.type * taxon2 + trait.no, nitt=500000, burnin=5000, thin=30, data=gmats, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# model 1 is supported across the board
Even.Bayes1$DIC	; Even.Bayes2$DIC	; Even.Bayes3$DIC ;	Even.Bayes4$DIC
# -43.71037			-43.40248		-35.102		-32.53772



summary( Even.Bayes1 )$solutions
# post.mean   l-95% CI   u-95% CI eff.samp        pMCMC
# (Intercept)  0.517137942  0.3248217 0.72888772 2976.144 0.0003157562
# trait.typeM -0.032967906 -0.2072599 0.15582019 2672.422 0.7098200189
# trait.typeS  0.121768106 -0.0847587 0.34980164 3167.000 0.2797600253
# taxon2P     -0.012270879 -0.2096049 0.18496107 3167.000 0.9125355226
# trait.no    -0.005475639 -0.0374283 0.02381878 3003.118 0.7218187559

Even.estimates <- cbind( Even.Bayes1$Sol[,1],
                         ( Even.Bayes1$Sol[,1] + Even.Bayes1$Sol[,2] ),
                         ( Even.Bayes1$Sol[,1] + Even.Bayes1$Sol[,3] ),
                         ( Even.Bayes1$Sol[,1] + Even.Bayes1$Sol[,4] ),
                         ( Even.Bayes1$Sol[,1] +Even.Bayes1$Sol[,2] +Even.Bayes1$Sol[,4]),
                         ( Even.Bayes1$Sol[,1] +Even.Bayes1$Sol[,3] +Even.Bayes1$Sol[,4] ) )

Even.estimates <- FoldedNormal( Even.estimates )

colnames(Even.estimates) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

# Tabulate estimates
apply( Even.estimates, 2, mean )
apply( Even.estimates, 2, Mode )
apply( Even.estimates, 2, hpd )

# A.LH       A.M       A.S      P.LH       P.M       P.S 
# 0.5117348 0.4863453 0.6360569 0.4987666 0.4731764 0.6242820 
# 0.4289423 0.3762670 0.8149248 0.5339408 0.5757301 0.5496561 
# 0.3203705 0.2968269 0.3823017 0.2498971 0.2466708 0.4255370
# 0.7088060 0.6679766 0.8793206 0.7469892 0.7048321 0.8508209


## ----eveness_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----
Evenfixed <- apply( Even.estimates, 2, mean )
plot(Evenfixed, type="n", xaxt="n", ylab="Evenness values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,1), main="Evenness")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(Even.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(Even.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(Even.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(Even.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(Even.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(Even.estimates[,6])), lwd=2 )

points( 1:3, Evenfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, Evenfixed[4:6], pch=21, bg="white", cex=2)


## ----read_in_correlation_matrices, echo=TRUE, eval=FALSE-----------------
## setwd("..//Data//Gmats_Cor_as_CSVs")
## matrices <- dir()
## no.mats <- length(matrices)
## cor_matrix_list <- list()
## 
## # this loop reads in each matrix from the folder of .csv files and writes them into a list
##   for ( i in 1:no.mats ) {
##     cor_matrix_list[[i]] <- read.csv( matrices[i] )
##   }
## names( cor_matrix_list ) <- matrices
## 
## # I shall write this list out as an object in order to make it as simple as possible for readers
## save( cor_matrix_list, file="correlation_matrix_list.R" )


## ----read_in_cor_matrix_list, echo=TRUE, results='asis'------------------

setwd( "..//Data" )
load( file="correlation_matrix_list.R" )

no.mats <- length( cor_matrix_list )
REV <- numeric(no.mats)   # allocating empty variables for the 
Even <- numeric(no.mats)
Evar <- numeric(no.mats)
gmax <- numeric(no.mats)

for (i in 1:no.mats) {
cor.matrix <- as.matrix( cor_matrix_list[[i]])
diag.matrix <- svd( cor.matrix )
eigen.values <- diag.matrix$d
gmax[i] <- eigen.values[1]    # gmax - the principal eigenvalue
Evar[i] <- sum( (eigen.values - 1 )^2 ) / length(eigen.values)    # 'eigenvalue variance' following Pavlicev et al 2009
REV[i] <- Evar[i] / ( length(eigen.values) -1 )   # 'relative eigenvalue variance' following Pavlicev et al 2009
lamda_tilde <- abs(eigen.values) / sum( abs(eigen.values) )
Even[i] <- - sum( lamda_tilde*log(lamda_tilde) / log( length(eigen.values) ) )    # 'eigenvalue eveness' after Agrawal & Stinchcombe 2008
}

cor_matrices_names <- sub( ".csv", "", names( cor_matrix_list ))
REV.mats <- data.frame( REV, Even, Evar, gmax, cor_matrices_names )
REV.mats$filename <- factor( cor_matrices_names )

# reading in the csv that holds the ID table (this may still be in memory, but just in case you didn't start from line 1 )
Index <- read.csv("MatrixIndexFinal.csv")
# str(Index)

# merge the metrics with the index
g.cor <- merge( Index, REV.mats, "filename" )
# str(g.cor)

# remove uneeded columns
g.cor <- g.cor[, -c(11, 12, 19)]


## ----data_hygiene_cor_matrices, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----
# check for weird values  #  122, 145 are visual outliers
plot( REV, Evar ) ;# identify( REV, Evar )
plot( g.cor$trait.no, Evar ) ;# identify( g.cor$trait.no, Evar )
plot( g.cor$trait.no, g.cor$gmax )  ;#   identify( g.cor$trait.no, g.cor$gmax )

# exclude McGuigan.et.al.2005, McGuigan.2006.MolEco, House.Simmons.2005.JEB and Petfield.et.al.2005.PNAS for >1 correlations
g.cor <- g.cor[ g.cor$study.code != "McG2005" ,]
g.cor <- g.cor[ g.cor$study.code != "McG2006" ,]
g.cor <- g.cor[ g.cor$study.code != "Pet2005" ,]
g.cor <- g.cor[ g.cor$study.code != "Hou2005" ,]

# exclude matrices with problematic units (e.g. mm^2)
g.cor <- g.cor[ g.cor$problem.units == "N", ]

dim( g.cor )  # 221 rows after data hygiene


## ----tabulate_descriptive_stats, echo=FALSE------------------------------
ddply( g.cor, .(taxon2, trait.type), summarise, m_REV= round(mean(REV), 3), m_Even= round(mean(Even), 3), m_Evar= round(mean(Evar), 3), m_gmax= round(mean(gmax), 3), m_Trait.no= round( mean(trait.no), 3), no.rows= length(gmax), m_no.families= round( mean(na.omit( no.families)), 1) )

# how do the metrics relate to each other?
round( cor( g.cor[,c( 13:16, 9)]), 3)

round( cor( na.omit(g.cor[,c( 13:16, 9, 12)]) ), 3)


## ----cor_pairs_plot, echo=FALSE, results='hide', message=FALSE, warning=FALSE, fig.height=7, fig.height=7----
metric.names <- expression( italic(paste("Var"["rel"], "(\u03BB)")),
                            italic("E"["\u03BB"]),
                            italic("Var(\u03BB)"),
                            italic(bolditalic("g")["max"]),
                            italic("n") )

pairs( g.cor[,c(13:16, 9)], labels=metric.names, lower.panel=panel.cor )


## ----more_supplementary_plots, message=FALSE, echo=FALSE, fig.keep='high'----
par( mfrow=c(1,2) )
plot( density( gmats$trait.no ), main="n traits from cov. mats." )
plot( density( g.cor$trait.no ), main="n traits from cor. mats" )
par( mfrow=c(1,1) )

table( gmats$trait.no )
table( g.cor$trait.no )


par( mfrow=c(1,2) )
plot( density( na.omit( gmats$no.families )), main="no. families from cov. mats." )
plot( density( na.omit( g.cor$no.families )), main="no. families from cor. mats" )
par( mfrow=c(1,1) )


names( gmats )[c(3, 6, 7, 18)]

pairs( gmats[,c(3, 6, 7, 18)], lower.panel=panel.cor )


names( g.cor )[c(13, 14, 15, 12)]

pairs( g.cor[c(13, 14, 15, 12)], lower.panel=panel.cor )



## ----model_selection_using_REV, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

REV.reml1 <- lmer(REV ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=g.cor)

REV.reml2 <- lmer(REV ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=g.cor)

REV.reml3 <- lmer(REV ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)

REV.reml4 <- lmer(REV ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)

BICtab( REV.reml1, REV.reml2, REV.reml3, REV.reml4, weights=T, nobs=nrow(g.cor) )
# model1 is looking good
                # dBIC df weight 
      # REV.reml1  0.0 7  0.93566
      # REV.reml3  5.4 8  0.06308
      # REV.reml2 13.3 9  0.00118
      # REV.reml4 18.7 10 < 0.001


# model 1 vs model 3
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=g.cor)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( REV.reml1 ))
LR.model <- -as.numeric(deviance( REV.reml3 ) - deviance( REV.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.675 no support for stepping up to model 3




myprior1 <- list(R = list(V =1, nu = 0.002), G= list( G1= list (V=1, nu=0.002)))

myprior2 <- list(R = list(V =1, nu = 0.002), G= list( G1= list (V=1, nu=0.002), G2= list (V=1, nu=0.002)))

myprior3 <- list( list(V=1, nu=0.002), list( G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002) ) )



# the bayesian approach
REV.Bayes0 <- MCMCglmm(REV ~ 1, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )
summary( REV.Bayes0 )$solutions

mean( REV.Bayes0$Sol ) ; Mode( REV.Bayes0$Sol )	; hpd( REV.Bayes0$Sol )
# 0.3242037				0.3507863			0.2342975 -- 0.4132790



REV.Bayes1 <- MCMCglmm(REV ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

REV.Bayes2 <- MCMCglmm(REV ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

REV.Bayes3 <- MCMCglmm(REV ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

REV.Bayes4 <- MCMCglmm(REV ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# support for model 1 all around
REV.Bayes1$DIC ; REV.Bayes2$DIC ; REV.Bayes3$DIC ; REV.Bayes4$DIC
# 16.74648			18.96781		35.18126		34.33423


plot(REV.Bayes1$VCV)
acf(REV.Bayes1$VCV)
plot(REV.Bayes1$Sol)
acf(REV.Bayes1$Sol)
summary(REV.Bayes1)$solutions
# post.mean    l-95% CI    u-95% CI eff.samp        pMCMC
# (Intercept)  0.357175556  0.13123310  0.55197790 3467.055 0.0003157562
# trait.typeM  0.079971791 -0.08056580  0.23336105 3888.525 0.3195453110
# trait.typeS  0.136484340 -0.06267555  0.31821710 3167.000 0.1654562678
# taxon2P     -0.206129518 -0.39344862 -0.01365675 2832.467 0.0359962109
# matdim      -0.007201403 -0.03740919  0.02333252 3167.000 0.6308809599

REV.estimates <- cbind( REV.Bayes1$Sol[,1],
                        ( REV.Bayes1$Sol[,1] + REV.Bayes1$Sol[,2] ),
                        ( REV.Bayes1$Sol[,1] + REV.Bayes1$Sol[,3] ),
                        ( REV.Bayes1$Sol[,1] + REV.Bayes1$Sol[,4] ),
                        ( REV.Bayes1$Sol[,1] +REV.Bayes1$Sol[,2] +REV.Bayes1$Sol[,4] ),
                        ( REV.Bayes1$Sol[,1] +REV.Bayes1$Sol[,3] +REV.Bayes1$Sol[,4] ))

REV.estimates <- FoldedNormal( REV.estimates )

colnames( REV.estimates ) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

apply( REV.estimates, 2, mean )
apply( REV.estimates, 2, Mode )
apply( REV.estimates, 2, hpd )
# A.LH       A.M    	 A.S     P.LH  	   P.M  	    P.S 
# mean	0.3582150 0.4336948 0.4962601 0.1572021 0.2341973 0.2897198 
# mode	0.45188673 0.45670191 0.49335045 0.06834723 0.24025053 0.12569836 
# -95%	0.1500290 0.2248246 0.2465886 0.000581226 0.002818709 0.06351661
# +95%	0.5595551 0.6341115 0.7452494 0.333986606 0.419158291 0.53022921


## ----REV_results, results='hide', fig.height=4, fig.width=4, echo=FALSE----
# REV plot
REVfixed <- apply( REV.estimates, 2, mean )
plot(REVfixed, type="n", xaxt="n", ylab="REV values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 0.9), main="Relative Eigenvalue Variance")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(REV.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(REV.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(REV.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(REV.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(REV.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(REV.estimates[,6])), lwd=2 )

points( 1:3, REVfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, REVfixed[4:6], pch=21, bg="white", cex=2)


## ----model_selection_for_Eigenvalue_Variance, echo=FALSE, results='hide', message=FALSE, fig.keep='none'----

Evar.reml1 <- lmer(Evar ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=g.cor)

Evar.reml2 <- lmer(Evar ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=g.cor)

Evar.reml3 <- lmer(Evar ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)

Evar.reml4 <- lmer(Evar ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)

BICtab( Evar.reml1, Evar.reml2, Evar.reml3, Evar.reml4, weights=T, nobs=nrow(g.cor) )
# model1 supported here
                # dBIC df weight
      # Evar.reml1  0.0 7  0.9206
      # Evar.reml3  5.4 8  0.0619
      # Evar.reml2  8.1 9  0.0164
      # Evar.reml4 13.5 10 0.0011


# model 1 vs model 3
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=g.cor)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( Evar.reml1 ))
LR.model <- -as.numeric(deviance( Evar.reml3 ) - deviance( Evar.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.44955 no support for stepping up to model 3



# the bayesian approach
Evar.Bayes0 <- MCMCglmm(Evar ~ 1, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )
summary( Evar.Bayes0 )$solutions

mean( Evar.Bayes0$Sol) ; Mode( Evar.Bayes0$Sol ) ; hpd( Evar.Bayes0$Sol)
# 1.395017					1.111799			0.9189918 -- 1.8634107


Evar.Bayes1 <- MCMCglmm(Evar ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

Evar.Bayes2 <- MCMCglmm(Evar ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

Evar.Bayes3 <- MCMCglmm(Evar ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

Evar.Bayes4 <- MCMCglmm(Evar ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# support for model3 here...
Evar.Bayes1$DIC ; Evar.Bayes2$DIC ; Evar.Bayes3$DIC ; Evar.Bayes4$DIC
# 493.9349			494.7357		    481.8635		      481.9546

# correlation is strong, but less so than in previous cases: r=0.67, however ∆BIC from model1->model3 was ~5 whereas ∆DIC from model3->model1 was ~10... I therefore am going to err on the side of the more complicated model
cor( colMeans(Evar.Bayes1$Sol), colMeans(Evar.Bayes3$Sol) )


plot(Evar.Bayes3$VCV)
acf(Evar.Bayes3$VCV)
plot(Evar.Bayes3$Sol)
acf(Evar.Bayes3$Sol)
summary(Evar.Bayes3)$solutions
#               post.mean   l-95% CI   u-95% CI eff.samp        pMCMC
# (Intercept) -0.34041740 -1.1542601 0.40912513 3167.000 0.3763814335
# trait.typeM  0.17250118 -0.8994697 1.19582255 3167.000 0.7420271550
# trait.typeS  0.09206386 -0.5237765 0.76050889 3261.582 0.7742342911
# taxon2P     -0.48873006 -1.0704802 0.07634632 3167.000 0.0972529207
# trait.no     0.33841058  0.2270651 0.45882085 3167.000 0.0003157562

Evar.estimates <- cbind( Evar.Bayes3$Sol[,1],
                         ( Evar.Bayes3$Sol[,1] + Evar.Bayes3$Sol[,2] ),
                         ( Evar.Bayes3$Sol[,1] + Evar.Bayes3$Sol[,3] ),
                         ( Evar.Bayes3$Sol[,1] + Evar.Bayes3$Sol[,4] ),
                         ( Evar.Bayes3$Sol[,1] + Evar.Bayes3$Sol[,2] + Evar.Bayes3$Sol[,4] ),
                         ( Evar.Bayes3$Sol[,1] + Evar.Bayes3$Sol[,3] + Evar.Bayes3$Sol[,4] ))

Evar.estimates <- FoldedNormal( Evar.estimates )

colnames( Evar.estimates ) <- c('A.LH','A.M','A.S','P.LH','P.M','P.S')

apply( Evar.estimates, 2, mean )
apply( Evar.estimates, 2, Mode )
apply( Evar.estimates, 2, hpd )
# A.LH       A.M    	 A.S     P.LH  	   P.M  	    P.S 
# mean	0.4260990 0.4840164 0.4126727 0.8348012 0.7256184 0.7326746
# mode	0.0398283 0.8957608 0.2789196 1.0728697 1.2359638 0.4357133 
# -95%	0.0001375 3.095e-05 0.0001118 0.0284474 7.427e-05 0.0073272
# +95%	0.9860031 1.183e+00 1.0339848 1.5183364 1.604e+00 1.3428949


## ----Eigenvalue_variance_results, results='hide', echo=FALSE, fig.height=4, fig.width=4----
Evarfixed <- apply( Evar.estimates, 2, mean )
plot(Evarfixed, type="n", xaxt="n", ylab="Evar values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 2), main="Eigenvalue Variance")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(Evar.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(Evar.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(Evar.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(Evar.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(Evar.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(Evar.estimates[,6])), lwd=2 )

points( 1:3, Evarfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, Evarfixed[4:6], pch=21, bg="white", cex=2)


## ----model_selection_for_Eveness, echo=FALSE, results='hide', fig.keep='none', message=FALSE----

Even.reml1 <- lmer(Even ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=g.cor)

Even.reml2 <- lmer(Even ~ 1 + trait.type * taxon2 + trait.no + (1|study.code), data=g.cor)

Even.reml3 <- lmer(Even ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)

Even.reml4 <- lmer(Even ~ 1 + trait.type * taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)

BICtab( Even.reml1, Even.reml2, Even.reml3, Even.reml4, weights=T, nobs=nrow(g.cor) )
# support for model1 again here
                  # dBIC df weight
        # Even.reml1  0.0 7  0.9313
        # Even.reml3  5.2 8  0.0685
        # Even.reml2 17.1 9  <0.001
        # Even.reml4 22.2 10 <0.001

# model 1 vs model 3
LikRatioSim <- function( mod ) {
  y.sim <- simulate(mod, nsim = 1)  # one set of simulated data
  model.lower <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code), data=g.cor)
  model.full  <- lmer(y.sim$sim_1 ~ 1 + trait.type + taxon2 + trait.no + (1|study.code) + (1|species), data=g.cor)
  LRSim <-  -as.numeric(deviance(model.full) - deviance(model.lower))
  return(LRSim)
  rm(y.sim, model.lower, model.full, LRSim)
}

LikRatParBoot <- replicate(n = 1000, LikRatioSim( Even.reml1 ))
LR.model <- -as.numeric(deviance( Even.reml3 ) - deviance( Even.reml1 ))
mean(c((LikRatParBoot > LR.model), 1))
# 0.2487512 no support for stepping up to model3


Even.Bayes0 <- MCMCglmm(Even ~ 1, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

mean(Even.Bayes0$Sol) ; Mode(Even.Bayes0$Sol) ; hpd(Even.Bayes0$Sol)
# 0.7347764				0.7570196			0.7045816 -- 0.7661564


Even.Bayes1 <- MCMCglmm(Even ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

Even.Bayes2 <- MCMCglmm(Even ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ study.code, prior= myprior1 )

Even.Bayes3 <- MCMCglmm(Even ~ 1 + trait.type + taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

Even.Bayes4 <- MCMCglmm(Even ~ 1 + trait.type * taxon2 + trait.no, nitt=100000, burnin=5000, thin=30, data=g.cor, verbose=F, random= ~ idh(trait.type):study.code + species, prior= myprior3 )

# DIC is in agreement in this case <- model1 is supported
Even.Bayes1$DIC ; Even.Bayes2$DIC ; Even.Bayes3$DIC ; Even.Bayes4$DIC
# -174.24			-171.81				-156.10			-152.75

plot(Even.Bayes1$VCV)
acf(Even.Bayes1$VCV)
plot(Even.Bayes1$Sol)
acf(Even.Bayes1$Sol)
summary(Even.Bayes1)$solutions
# post.mean    l-95% CI     u-95% CI eff.samp     pMCMC
# (Intercept)  0.7582261090  0.68046150  0.843512765     3167 0.0003157
# trait.typeM -0.0624660677 -0.13150820  0.004840475     3167 0.0776760
# trait.typeS -0.0908189412 -0.17028441 -0.012640371     3167 0.0296810
# taxon2P      0.0965238982  0.03635891  0.165503796     3167 0.0063151
# matdim      -0.0002662439 -0.01247044  0.011770649     3167 0.9674771


Even.estimates <- cbind( Even.Bayes1$Sol[,1],
                         ( Even.Bayes1$Sol[,1] + Even.Bayes1$Sol[,2] ),
                         ( Even.Bayes1$Sol[,1] + Even.Bayes1$Sol[,3] ),
                         ( Even.Bayes1$Sol[,1] + Even.Bayes1$Sol[,4] ),
                         ( Even.Bayes1$Sol[,1] +Even.Bayes1$Sol[,2] +Even.Bayes1$Sol[,4] ),
                         ( Even.Bayes1$Sol[,1] +Even.Bayes1$Sol[,3] +Even.Bayes1$Sol[,4] ) )

Even.estimates <- FoldedNormal( Even.estimates )

colnames( Even.estimates )<-c('A.LH','A.M','A.S','P.LH','P.M','P.S')

apply( Even.estimates, 2, mean )
apply( Even.estimates, 2, Mode )
apply( Even.estimates, 2, hpd )
# A.LH       A.M       A.S      P.LH       P.M       P.S 
# 0.7573074 0.6951077 0.6677275 0.8559004 0.7919052 0.7657368 
# 0.7817484 0.7708938 0.6880481 0.7980492 0.8138261 0.8059144 
# 0.6742629 0.6136220 0.5717078 0.7691409 0.7080655 0.6790868
# 0.8353725 0.7712633 0.7651096 0.9442544 0.8722907 0.8568488


## ----Eveness_results, results='hide', echo=FALSE, fig.height=4, fig.width=4----
Evenfixed <- apply( Even.estimates, 2, mean )
plot(Evenfixed, type="n", xaxt="n", ylab="Evenness values", xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0.4, 1), main="Evenness")
axis(side=1, c("LH","M","S"), at=c(1,2,3))

lines( c(1,1), c(hpd(Even.estimates[,1])), lwd=2 )
lines( c(2,2), c(hpd(Even.estimates[,2])), lwd=2 )
lines( c(3,3), c(hpd(Even.estimates[,3])), lwd=2 )
lines( c(1.1,1.1), c(hpd(Even.estimates[,4])), lwd=2 )
lines( c(2.1,2.1), c(hpd(Even.estimates[,5])), lwd=2 )
lines( c(3.1,3.1), c(hpd(Even.estimates[,6])), lwd=2 )

points( 1:3, Evenfixed[1:3], pch=16, cex=2)
#animals - filled points ^ #plants - open points v
points( 1.1:3.1, Evenfixed[4:6], pch=21, bg="white", cex=2)


## ----ouputting_plots_for_the_MS, eval=FALSE, echo=TRUE-------------------
## 
## ##################################################
## ######## This block makes MS-quality plots #######
## ##################################################
## 
## ### this plot is about covariance matrices  <- Figure 6 in the paper.
## par( mfrow=c( 2, 2 ), mar=c( 2, 4, 2, 2 ))
## #animals - filled points ^ #plants - open points v
## 
## # top left
## EMAXfixed <- apply( EMAX.estimates, 2, mean )
## plot(EMAXfixed, type="n", xaxt="n", ylab=expression(italic("e"["max"])), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,5) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(EMAX.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(EMAX.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(EMAX.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(EMAX.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(EMAX.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(EMAX.estimates[,6])), lwd=2 )
## points( 1:3, EMAXfixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, EMAXfixed[4:6], pch=21, bg="white", cex=2)
## 
## # top right
## TGVfixed <- apply( TGV.estimates, 2, mean )
## plot(TGVfixed, type="n", xaxt="n", ylab=expression(italic("tgv")), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,75) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(TGV.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(TGV.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(TGV.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(TGV.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(TGV.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(TGV.estimates[,6])), lwd=2 )
## points( 1:3, TGVfixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, TGVfixed[4:6], pch=21, bg="white", cex=2)
## 
## #bottom left
## AEvofixed <- apply( AEvo.estimates, 2, mean )
## plot(AEvofixed, type="n", xaxt="n", ylab=expression(italic("\u0113")), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 12) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(AEvo.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(AEvo.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(AEvo.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(AEvo.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(AEvo.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(AEvo.estimates[,6])), lwd=2 )
## points( 1:3, AEvofixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, AEvofixed[4:6], pch=21, bg="white", cex=2)
## 
## # bottom right
## ENDfixed <- apply( END.estimates, 2, mean )
## plot(ENDfixed, type="n", xaxt="n", ylab=expression(italic("n"["D"])), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0,2) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(END.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(END.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(END.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(END.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(END.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(END.estimates[,6])), lwd=2 )
## points( 1:3, ENDfixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, ENDfixed[4:6], pch=21, bg="white", cex=2)
## 
## ###############################
## 
## ### this plot is about correlation matrices  <- Figure 7 in the paper.
## par( mfrow=c( 1, 3 ), mar=c( 2, 4, 2, 2 ))
## #animals - filled points ^ #plants - open points v
## 
## # left panel
## REVfixed <- apply( REV.estimates, 2, mean )
## plot(REVfixed, type="n", xaxt="n", ylab=expression(italic(paste("Var"["rel"], "(\u03BB)"))), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 0.9) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(REV.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(REV.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(REV.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(REV.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(REV.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(REV.estimates[,6])), lwd=2 )
## points( 1:3, REVfixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, REVfixed[4:6], pch=21, bg="white", cex=2)
## 
## # middle panel
## Evarfixed <- apply( Evar.estimates, 2, mean )
## plot(Evarfixed, type="n", xaxt="n", ylab=expression(italic(paste("Var(\u03BB)"))), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0, 2) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(Evar.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(Evar.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(Evar.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(Evar.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(Evar.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(Evar.estimates[,6])), lwd=2 )
## points( 1:3, Evarfixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, Evarfixed[4:6], pch=21, bg="white", cex=2)
## 
## # right panel
## Evenfixed <- apply( Even.estimates, 2, mean )
## plot(Evenfixed, type="n", xaxt="n", ylab=expression(italic("E"["\u03BB"])), xlab="Trait type", xlim=c(0.7,3.3), ylim=c(0.4, 1) )
## axis(side=1, c("LH","M","S"), at=c(1,2,3))
## 
## lines( c(1,1), c(hpd(Even.estimates[,1])), lwd=2 )
## lines( c(2,2), c(hpd(Even.estimates[,2])), lwd=2 )
## lines( c(3,3), c(hpd(Even.estimates[,3])), lwd=2 )
## lines( c(1.1,1.1), c(hpd(Even.estimates[,4])), lwd=2 )
## lines( c(2.1,2.1), c(hpd(Even.estimates[,5])), lwd=2 )
## lines( c(3.1,3.1), c(hpd(Even.estimates[,6])), lwd=2 )
## points( 1:3, Evenfixed[1:3], pch=16, cex=2)
## points( 1.1:3.1, Evenfixed[4:6], pch=21, bg="white", cex=2)
## 
## #####


