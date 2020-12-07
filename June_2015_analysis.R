##################set the working directory and import the respective packages required######

setwd("C:/Users/newscripts")
library(rjags)
library(geosphere)
library(mcmcplots)
library(ggplot2)
library(runjags)
library(ggmcmc)
library(mcmcplots)
library(RColorBrewer)
library(gridExtra)

#####################importing the file containing the isotope data##########

isotope = read.csv("C:/Users/June_2015_distribution_sites.csv")

isite = subset(isotope, type=="Distribution_site")   ######This  contains all the distribution site data
isource = subset(isotope, type=="source_water")  #######This contains all the source data

##########################################importing source data 

Sources = read.csv("C:/Users/Sources_June_2015.csv")
SRC1 = subset(Sources, source== "SRC1")
SRC2 = subset(Sources, source== "SRC2")
SRC3 = subset(Sources, source== "SRC3")
SRC4 = subset(Sources, source== "SRC4")
SRC5 = subset(Sources, source== "SRC5")
SRC6 = subset(Sources, source== "SRC6")
SRC7 = subset(Sources, source== "SRC7")
SRC8 = subset(Sources, source== "SRC8")
SRC9 = subset(Sources, source== "SRC9")
SRC10 = subset(Sources, source== "SRC10")
SRC11 = subset(Sources, source== "SRC11")


SRC1 <- cbind(SRC1$mean_d18O,SRC1$mean_d2H)
SRC2 <- cbind(SRC2$mean_d18O,SRC2$mean_d2H)
SRC3 <- cbind(SRC3$mean_d18O,SRC3$mean_d2H)
SRC4 <- cbind(SRC4$mean_d18O,SRC4$mean_d2H)
SRC5 <- cbind(SRC5$mean_d18O,SRC5$mean_d2H)
SRC6 <- cbind(SRC6$mean_d18O,SRC6$mean_d2H)
SRC7 <- cbind(SRC7$mean_d18O,SRC7$mean_d2H)
SRC8 <- cbind(SRC8$mean_d18O,SRC8$mean_d2H)
SRC9 <- cbind(SRC9$mean_d18O,SRC9$mean_d2H)
SRC10 <- cbind(SRC10$mean_d18O,SRC10$mean_d2H)
SRC11 <- cbind(SRC11$mean_d18O,SRC11$mean_d2H)


L1 = nrow(SRC1)
L2 = nrow(SRC2)
L3 = nrow(SRC3)
L4 = nrow(SRC4)
L5 = nrow(SRC5)
L6 = nrow(SRC6)
L7 = nrow(SRC7)
L8 = nrow(SRC8)
L9 = nrow(SRC9)
L10 = nrow(SRC10)
L11 = nrow(SRC11)

#####plotting mean source and distribution site isotope values ##########

s_mean = cbind(isource$mean_d18O,isource$mean_d2H)
S_SD = cbind(isource$stdev_d18O, isource$stdev_d2H)
y = cbind(isite$mean_d18O,isite$mean_d2H)
y_SD = cbind(isite$stdev_d18O, isite$stdev_d2H)
plot.df = rbind(isite,isource)
ggplot(plot.df, aes(x=mean_d18O, y=mean_d2H, shape=type, col=type)) + geom_point(size=3)

############################################

################# assigning projection to lat long data ( this is done to calculate the distance between the supply sites and respective 
## sources)

isite.sp = isite
coordinates(isite.sp) <- ~lon+lat
proj4string(isite.sp) <- CRS("+init=epsg:4326")

isource.sp = isource
coordinates(isource.sp) <- ~lon+lat
proj4string(isource.sp) <- CRS("+init=epsg:4326")


distM = cbind(
  distGeo(p1=isite.sp, p2=isource.sp[1,]),   ########calculating distance bewen source 1 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[2,]),   ########calculating distance bewen source 2 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[3,]),   ########calculating distance bewen source 3 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[4,]),   ########calculating distance bewen source 4 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[5,]),   ########calculating distance bewen source 5 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[6,]),   ########calculating distance bewen source 6 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[7,]),   ########calculating distance bewen source 7 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[8,]),   ########calculating distance bewen source 8 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[9,]),   ########calculating distance bewen source 9 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[10,]),  ########calculating distance bewen source 10 and distribution sites
  distGeo(p1=isite.sp, p2=isource.sp[11,])   ########calculating distance bewen source 11 and distribution sites
)/1000


#### no. of supply sites, sources and isotope used
Ntotal = dim(y)[1]
Jtotal = dim(y)[2]
Ktotal = dim(s_mean)[1]


###### the next part is used for calculating the priors for each source supplying a given site based upon its distance from the source
### and the volume of water supplied by the source

##### the matrix distM1 gives normalized distance weight. i.e. the sum of the distance of a given supply sites from all the sources will be
### equal to 1. The source that will be further away from the supply site will have higher weight. 

distM1= rep(1, Ktotal)
distM1 = matrix(rep(distM1,each=Ntotal),nrow=Ntotal)
q <- matrix((rep(1, Ntotal)), nrow  = Ntotal)
for (i in 1: Ntotal){
  for(j in 1:Ktotal){
    q[i,] <- sum(distM[i,])
    distM1[i,j] <- 1/(distM[i,j]/q[i,])   #####matrix containing normalize distance
  }
}

## alpha 1 is  a 1X11 matrix (11 sources) with each value representing the fraction of water provided by the respective source
#for example source 4 with a value of 0.61 supplies 61% of water
alpha1 = matrix(1, nrow=Ntotal, ncol=Ktotal)
alpha1 = isource$volumes/sum(isource$volumes)


alpha_volume =matrix(rep(alpha1,each=Ntotal),nrow=Ntotal)###alpha 1 is repeated Ntotal no. of times. where Ntotal is the total no. of distribution sites
alpha_volume
alpha2 = alpha_volume*distM1 
alpha2

##calculate the final  priors

alpha= rep(1, Ktotal)
alpha = matrix(1, nrow=Ntotal, ncol=Ktotal)
q1 <- matrix((rep(1, Ntotal)), nrow  = Ntotal)
for (i in 1: Ntotal){
  for(j in 1:Ktotal){
    q1[i,] <- sum(alpha2[i,])
    alpha[i,j] <- alpha2[i,j]/q1[i,]
  }
}

#alpha gives the final prior  probability of a given source supplying a site.
# This depends upon the distanc of the source from the supply site and volume of water supplied by the respective sources.

alpha

######################### calculation of prior finished
#########################


##############Getting ready to run JAGS now

## S3 is the identity matrix for the Wishart
S3 = diag(Jtotal)

#######this datalist is again requried for JAGS. All the variables needed to run the model should be here.
dataList = list(
  S3 = S3,
  y = y,
  SRC1 = SRC1,
  SRC2 = SRC2,
  SRC3 = SRC3,
  SRC4 = SRC4,
  SRC5 = SRC5,
  SRC6 = SRC6,
  SRC7 = SRC7,
  SRC8 = SRC8,
  SRC9 = SRC9,
  SRC10 = SRC10,
  SRC11 = SRC11,
  L1 = L1,
  L2 = L2,
  L3 = L3,
  L4 = L4,
  L5 = L5,
  L6 = L6,
  L7 = L7,
  L8 = L8,
  L9 = L9,
  L10 = L10,
  L11 = L11,
  alpha = alpha,
  K = Ktotal,
  N = Ntotal,
  J = dim(y)[2]
)

parList = c("p", "mu") ##list of parameters to be observed in the posterior distribution


####This model string contains the JAGS model to calculate the posterior 

modelString = "
  model {
    # Likelihood
      for (k in 1:11){
      for (j in 1:2){
         s_mean1[k,j] ~ dnorm(0, 0.000000001)
      }
    }
     for ( i in 1:L1){
       SRC1[i,1:2] ~ dmnorm(s_mean1[1,1:2], tau1 )
     }
    for ( i in 1:L2){
      SRC2[i,1:2] ~ dmnorm(s_mean1[2,1:2], tau2 )
    }
    for ( i in 1:L3){
      SRC3[i,1:2] ~ dmnorm(s_mean1[3,1:2], tau3 )
    }
    for ( i in 1:L4){
       SRC4[i,1:2] ~ dmnorm(s_mean1[4,1:2], tau4 )
     }
    for ( i in 1:L5){
      SRC5[i,1:2] ~ dmnorm(s_mean1[5,1:2], tau5 )
    }
    for ( i in 1:L6){
      SRC6[i,1:2] ~ dmnorm(s_mean1[6,1:2], tau6 )
    } 
    for ( i in 1:L7){
      SRC7[i,1:2] ~ dmnorm(s_mean1[7,1:2], tau7 )
    }
    for ( i in 1:L8){
       SRC8[i,1:2] ~ dmnorm(s_mean1[8,1:2], tau8 )
    }
    for ( i in 1:L9){
       SRC9[i,1:2] ~ dmnorm(s_mean1[9,1:2], tau9 )
    }
    for ( i in 1:L10){
       SRC10[i,1:2] ~ dmnorm(s_mean1[10,1:2], tau10 )
    }
    for ( i in 1:L11){
       SRC11[i,1:2] ~ dmnorm(s_mean1[11,1:2], tau11 )
    }
    # Now p must become p[i,1:k] row per obs, col for each source
    for (i in 1:N) {
      # Generate mu
      for (j in 1:J) {
        mu[i,j] = inprod(f[i,1:K],s_mean1[1:K,j]) ##### mu is the mean of the supply site as calculated by multyplying p's with respective
                                                #####  source isotope ratio and summing them 
      }
     
      y[i,1:J] ~ dmnorm( mu[i,1:J], tau )
    }
    # 
    # Dirichlet prior on p
    # Wrap this in N loop
    for (i in 1:N) {
        f[i,1:K] ~ ddirch(alpha[i,])      #####  dirchlet distribution to calculate the f's which will be used to calculate the p's 
      # Now estimate p from sum of alpha
     for(k in 1:K) {
        p[i,k] <- f[i,k] / sum(f[i,1:K])    ## proportion of water supplied by different sources 
      }
    }
    # Tau - precision on sample variance
    # Wishart prior: J d.f.
    tau[1:J,1:J] ~  dwish(S3[,],J)
    tau1[1:J,1:J] ~  dwish(S3[,],J)
    tau2[1:J,1:J] ~  dwish(S3[,],J)
    tau3[1:J,1:J] ~  dwish(S3[,],J)
    tau4[1:J,1:J] ~  dwish(S3[,],J)
    tau5[1:J,1:J] ~  dwish(S3[,],J)
    tau6[1:J,1:J] ~  dwish(S3[,],J)
    tau7[1:J,1:J] ~  dwish(S3[,],J)
    tau8[1:J,1:J] ~  dwish(S3[,],J)
    tau9[1:J,1:J] ~  dwish(S3[,],J)
    tau10[1:J,1:J] ~  dwish(S3[,],J)
    tau11[1:J,1:J] ~  dwish(S3[,],J)
  }
"


writeLines( modelString, con ="jags_isotopes_v12.txt") ### this writes the model in txt file to be called by JAGS


####running the model in jags, defining burnin  and no. of samples to be run

runJagsOut30 <-run.jags(method = "parallel", model = "jags_isotopes_v12.txt", monitor = c("p", "mu"), 
                        inits = NA, data = dataList,
                        n.chains = 3, adapt = 50000, burnin = 40000, sample = 6000, thin = 50, 
                        summarise = FALSE, plots = FALSE )

codaSamples =as.mcmc.list(runJagsOut30)
#################
########################################analysis of JAGS output in R
#################################
options(max.print=1000000)
#### analyzing the robustness of the posterior distribution. 
gelman.diag(codaSamples, multivariate = FALSE)

x <- summary(codaSamples)   ### summarzing the credible intervals, mean SD, for the posteriors

#### median values of the different sources for the different sites#################

isite.sp$pSRC1_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",1]"),][,3]
isite.sp$pSRC2_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",2]"),][,3]
isite.sp$pSRC3_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",3]"),][,3]
isite.sp$pSRC4_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",4]"),][,3]
isite.sp$pSRC5_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",5]"),][,3]
isite.sp$pSRC6_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",6]"),][,3]
isite.sp$pSRC7_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",7]"),][,3]
isite.sp$pSRC8_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",8]"),][,3]
isite.sp$pSRC9_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",9]"),][,3]
isite.sp$pSRC10_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",10]"),][,3]
isite.sp$pSRC11_median <- x$quantiles[paste0("p[",seq(1,Ntotal),",11]"),][,3]

##################credible intervals of the different sources##########################

isite.sp$pSRC1_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",1]"),][,1]
isite.sp$pSRC2_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",2]"),][,1]
isite.sp$pSRC3_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",3]"),][,1]
isite.sp$pSRC4_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",4]"),][,1]
isite.sp$pSRC5_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",5]"),][,1]
isite.sp$pSRC6_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",6]"),][,1]
isite.sp$pSRC7_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",7]"),][,1]
isite.sp$pSRC8_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",8]"),][,1]
isite.sp$pSRC9_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",9]"),][,1]
isite.sp$pSRC10_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",10]"),][,1]
isite.sp$pSRC11_lower <- x$quantiles[paste0("p[",seq(1,Ntotal),",11]"),][,1]

isite.sp$pSRC1_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",1]"),][,5]
isite.sp$pSRC2_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",2]"),][,5]
isite.sp$pSRC3_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",3]"),][,5]
isite.sp$pSRC4_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",4]"),][,5]
isite.sp$pSRC5_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",5]"),][,5]
isite.sp$pSRC6_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",6]"),][,5]
isite.sp$pSRC7_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",7]"),][,5]
isite.sp$pSRC8_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",8]"),][,5]
isite.sp$pSRC9_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",9]"),][,5]
isite.sp$pSRC10_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",10]"),][,5]
isite.sp$pSRC11_upper <- x$quantiles[paste0("p[",seq(1,Ntotal),",11]"),][,5]

####################### mean and standard deviation of the contribtuion from different sources#######################

isite.sp$pSRC1_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",1]"),][,1]
isite.sp$pSRC2_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",2]"),][,1]
isite.sp$pSRC3_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",3]"),][,1]
isite.sp$pSRC4_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",4]"),][,1]
isite.sp$pSRC5_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",5]"),][,1]
isite.sp$pSRC6_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",6]"),][,1]
isite.sp$pSRC7_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",7]"),][,1]
isite.sp$pSRC8_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",8]"),][,1]
isite.sp$pSRC9_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",9]"),][,1]
isite.sp$pSRC10_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",10]"),][,1]
isite.sp$pSRC11_mean <- x$statistics[paste0("p[",seq(1,Ntotal),",11]"),][,1]

isite.sp$pSRC1_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",1]"),][,2]
isite.sp$pSRC2_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",2]"),][,2]
isite.sp$pSRC3_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",3]"),][,2]
isite.sp$pSRC4_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",4]"),][,2]
isite.sp$pSRC5_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",5]"),][,2]
isite.sp$pSRC6_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",6]"),][,2]
isite.sp$pSRC7_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",7]"),][,2]
isite.sp$pSRC8_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",8]"),][,2]
isite.sp$pSRC9_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",9]"),][,2]
isite.sp$pSRC10_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",10]"),][,2]
isite.sp$pSRC11_SD <- x$statistics[paste0("p[",seq(1,Ntotal),",11]"),][,2]

########################################### We export the output as a csv file for post processing.
###This file contains mean, median, standard deviation and 95% CI of the contribution of different sources at the distribution sites.
source.prop <- data.frame(isite.sp)
write.csv(source.prop, file = "June_2015.csv")
