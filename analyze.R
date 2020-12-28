
rm( list = ls() )

####################### SET UP ####################### 

code.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Code (git)"
private.data.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Data (private)"
results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Results from R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/Metasens tutorial"

setwd(code.dir)
source("helper.R")

setwd(private.data.dir)
dfs = list( read.csv("kodama_prepped.csv"),
            read.csv("flegal_prepped.csv"),
            read.csv("gbc_prepped.csv") )

meta.names = c("KODAMA", "FLEGAL", "GBMC")


library(MetaUtility)
library(ggplot2)
library(dplyr)
library(metafor)
library(testthat)
library(boot)

digits = 2

# # use local versions of EValue code
# detach(package:EValue, unload = TRUE)
# setwd("~/Dropbox/Personal computer/Independent studies/EValue package/evalue_package_git/EValue/R")
# source("EValue.R")

####################### NORMALITY ASSUMPTION ####################### 


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
  cat("\nEstimate (RR or HR) and CI:", paste( round( exp(mu), digits ),
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


####################### SIMPLE META-ANALYSIS - KODAMA ####################### 

dat = dfs[[1]]

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

# to report
exp(mu)
exp(mu.lo)
exp(mu.hi)
sqrt(t2)

# example of using R package for the Supplement
library(EValue)

# E-value
evalue( est = RR( exp(meta$b) ),
        lo = RR( exp(meta$ci.lb) ),
        hi = RR( exp(meta$ci.ub) ) )

# homogeneous bias
confounded_meta( method = "calibrated",
                 q = log(1.1),
                 r = 0.15,
                 tail = "above",
                 muB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )

# heterogeneous bias
confounded_meta( method = "parametric",
                q = log(1.1),
                r = 0.15,
                muB = log(1.5),
                sigB = sqrt( 0.8 * t2 ),
                yr = meta$b,
                vyr = meta$se^2,
                t2 = meta$tau2,
                vt2 = meta$se.tau2^2 )


####################### SIMPLE ANALYSES - FLEGAL ####################### 


dat = dfs[[1]]
meta.name = "Flegal"
q = log(.9)
tail = c("below")
Bmax = 1.3
muB.phat = log(1.03)

##### Naive Meta-Analysis #####
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

sigB = varB.ratio * sqrt(t2)

# to report
exp(mu)
exp(mu.lo)
exp(mu.hi)
sqrt(t2)

####################### SIMPLE ANALYSES - GBC ####################### 

dat = dfs[[2]]
meta.name = "GBC"
tail = c( "above" )
muB.phat = log(1.03)
q = log(1.1)



##### Naive Meta-Analysis #####
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

sigB = varB.ratio * sqrt(t2)

# to report
exp(mu)
exp(mu.lo)
exp(mu.hi)
sqrt(t2)
