
# Note: Results produced by this script may differ slightly from those presented in the main text.
#  Results in the main text were from the website, for which we had used rounded values as inputs.
#  Here we use exact inputs. Also, confidence intervals produced by bootstrapping can differ slightly 
#  across multiple runs. 

rm( list = ls() )

####################### SET UP ####################### 

library(MetaUtility)
library(ggplot2)
library(dplyr)
library(metafor)
library(testthat)
library(boot)
#library(EValue)
library(here)

# code.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Code (git)"
private.data.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Data (private)"
# results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Results from R"
# overleaf.dir = "~/Dropbox/Apps/Overleaf/Metasens tutorial"

# setwd(code.dir)
# source("helper.R")

setwd(private.data.dir)
dfs = list( read.csv("kodama_prepped.csv"),
            read.csv("flegal_prepped.csv"),
            read.csv("gbc_prepped.csv") )

meta.names = c( "KODAMA", "FLEGAL", "GBMC" )

# rounding digits
digits = 2

# use local versions of EValue code
detach(package:EValue, unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/R packages/EValue package (git)/evalue_package/EValue/R")
source("meta-analysis.R")
source("utils.R")
source("effect_measures.R")
source("EValue.R")

####################### NAIVE META-ANALYSIS FOR ALL THREE METAS ####################### 

# got through each meta-analysis and fit naive analysis
# also check normality, relevant only for parametric sensitivity analyses that use
#   heterogeneous bias

for ( i in 1:length(dfs) ){
  
  dat = dfs[[i]]
  cat("\n\n***************** ", meta.names[i], " *****************\n")
  
  ##### simple meta-analysis #####
  ( meta = rma.uni( yi = dat$yi, 
                    vi = dat$vi, 
                    method = "REML", 
                    knha = TRUE ) )
  
  mu = meta$b
  t2 = meta$tau2
  mu.lo = meta$ci.lb
  mu.hi = meta$ci.ub
  mu.se = meta$se
  mu.pval = meta$pval
  
  print(meta)
  
  # to report
  cat("\n***Naive estimate (RR or HR) and CI:", paste( round( exp(mu), digits ),
                                               " [",
                                               round( exp(mu.lo), digits ),
                                               ", ",
                                               round( exp(mu.hi), digits ),
                                               "]",
                                               sep = "") )  

  ##### normality #####
  # all seem nicely normal :)
  dat$calib = calib_ests( yi = dat$yi,
                          sei = sqrt(dat$vi) )
  
  cat("\n\nShapiro test: ")
  print( shapiro.test(dat$calib) )
  
  qqnorm(dat$calib)
  qqline(dat$calib, col = "red")
}


####################### ANALYZE KODAMA ####################### 

dat = dfs[[1]]

# refit naive meta-analysis
( meta = rma.uni( yi = dat$yi, 
                  vi = dat$vi, 
                  method = "REML", 
                  knha = TRUE ) )

# E-value for point estimate
evalue( est = RR( exp(meta$b) ),
        lo = RR( exp(meta$ci.lb) ),
        hi = RR( exp(meta$ci.ub) ) )

# uncorrected estimate of effects with RR > 1.1
# Tmin and Gmin here also give amount of bias and of confounding required
#   to reduce to less than 15% the percentage of meaningfully large effects
confounded_meta( method = "calibrated",
                 q = log(1.1),
                 r = 0.15,
                 tail = "above",
                 muB = 0,
                 sigB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )

# heterogeneous bias (Supplement)
confounded_meta( method = "parametric",
                q = log(1.1),
                r = 0.15,
                tail = "above",
                muB = log(1.5),
                sigB = sqrt( 0.8 * meta$tau2 ),
                yr = meta$b,
                vyr = meta$se^2,
                t2 = meta$tau2,
                vt2 = meta$se.tau2^2 )


####################### ANALYZE FLEGAL ####################### 


dat = dfs[[2]]

# refit naive meta-analysis
( meta = rma.uni( yi = dat$yi, 
                  vi = dat$vi, 
                  method = "REML", 
                  knha = TRUE ) )

# E-value for point estimate
# HR for a rare outcome approximates an RR
evalue( est = RR( exp(meta$b) ),
        lo = RR( exp(meta$ci.lb) ),
        hi = RR( exp(meta$ci.ub) ) )

# uncorrected estimate of effects with RR < 0.90
# Tmin and Gmin here also give amount of bias and of confounding required
#   to reduce to less than 15% the percentage of meaningfully large effects
confounded_meta( method = "calibrated",
                 q = log(0.90),
                 r = 0.15,
                 tail = "below",
                 muB = 0,
                 sigB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )

# uncorrected estimate of effects with RR > 1.1
confounded_meta( method = "calibrated",
                 q = log(1.1),
                 r = 0.15,
                 tail = "above",
                 muB = 0,
                 sigB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )

# heterogeneous bias (Supplement)
# the 1.02 is chosen to avoid going below 15% and hence compromising parametric estimation
confounded_meta( method = "parametric",
                 q = log(0.90),
                 tail = "below",
                 muB = log(1.02),
                 sigB = sqrt( 0.8 * meta$tau2 ),
                 yr = meta$b,
                 vyr = meta$se^2,
                 t2 = meta$tau2,
                 vt2 = meta$se.tau2^2 )


####################### ANALYZE GBMC ####################### 

dat = dfs[[3]]

# refit naive meta-analysis
( meta = rma.uni( yi = dat$yi, 
                  vi = dat$vi, 
                  method = "REML", 
                  knha = TRUE ) )

# E-value for point estimate
# HR for a rare outcome approximates an RR
evalue( est = RR( exp(meta$b) ),
        lo = RR( exp(meta$ci.lb) ),
        hi = RR( exp(meta$ci.ub) ) )

# uncorrected estimate of effects with RR < 0.90
# Tmin and Gmin here also give amount of bias and of confounding required
#   to reduce to less than 15% the percentage of meaningfully large effects
confounded_meta( method = "calibrated",
                 q = log(0.90),
                 r = 0.15,
                 tail = "below",
                 muB = 0,
                 sigB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )

# uncorrected estimate of effects with RR > 1.1
confounded_meta( method = "calibrated",
                 q = log(1.1),
                 r = 0.15,
                 tail = "above",
                 muB = 0,
                 sigB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )

# heterogeneous bias (Supplement)
# the 1.01 is chosen to avoid going below 15% and hence compromising parametric estimation
confounded_meta( method = "parametric",
                 q = log(1.1),
                 tail = "above",
                 muB = log(1.01),
                 sigB = sqrt( 0.8 * meta$tau2 ),
                 yr = meta$b,
                 vyr = meta$se^2,
                 t2 = meta$tau2,
                 vt2 = meta$se.tau2^2 )
