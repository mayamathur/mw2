
# Note: Results produced by this script may differ slightly from those presented in the main text.
#  Results in the main text were from the website, for which we had used rounded values as inputs.
#  Here we use exact inputs. Also, confidence intervals produced by bootstrapping can differ slightly 
#  across multiple runs. 

rm( list = ls() )



# PRELIMINARIES ---------------------------------------------------

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
renv::restore()

# data-wrangling packages
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
library(MetaUtility)
library(EValue)
# other
library(here)
library(xtable)
library(testthat)

# run this only if you want to update the R environment specs
# renv::snapshot()


# set working directories
code.dir = here()

# go up a level in parent directory
results.dir = str_replace_all( string = here(),
                               pattern = "Code \\(git\\)",
                               replacement = "Results from R" ) 

private.data.dir = str_replace_all( string = here(),
                                    pattern = "Linked to OSF \\(MW2\\)/Code \\(git\\)",
                                    replacement = "Data (private)" )

#@not yet synced with DropBox?
overleaf.dir = results.dir


# get helper fns
setwd(code.dir)
source("helper.R")

# no sci notation
options(scipen=999)
# rounding digits
digits = 2


# get prepped data
setwd(private.data.dir)
dfs = list( read.csv("flegal_prepped.csv"),
            read.csv("gbc_prepped.csv") )

meta.names = c( "FLEGAL", "GBMC" )



# NAIVE META-ANALYSIS FOR BOTH METAS ---------------------------------------------------

# got through each meta-analysis and fit naive analysis
# also check normality, relevant only for parametric sensitivity analyses that use
#   heterogeneous bias

for ( i in 1:length(dfs) ){
  
  dat = dfs[[i]]
  cat("\n\n***************** ", meta.names[i], " *****************\n")
  
  ### Simple meta-analysis
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
  cat("\n***Naive estimate (HR) and CI:", paste( round( exp(mu), digits ),
                                               " [",
                                               round( exp(mu.lo), digits ),
                                               ", ",
                                               round( exp(mu.hi), digits ),
                                               "]",
                                               sep = "") ) 
  
  #### Normality
  # all seem nicely normal :)
  dat$calib = calib_ests( yi = dat$yi,
                          sei = sqrt(dat$vi) )
  
  cat("\n\nShapiro test: ")
  print( shapiro.test(dat$calib) )
  
  qqnorm(dat$calib)
  qqline(dat$calib, col = "red")
}



# ANALYZE FLEGAL ---------------------------------------------------


dat = dfs[[1]]

# refit naive meta-analysis
( meta = rma.uni( yi = dat$yi, 
                  vi = dat$vi, 
                  method = "REML", 
                  knha = TRUE ) )

# ~ E-value for point estimate -----------------------
# HR for a rare outcome approximates an RR
evalue( est = RR( exp(meta$b) ),
        lo = RR( exp(meta$ci.lb) ),
        hi = RR( exp(meta$ci.ub) ) )

# ~  Homogeneous bias; percentage of meaningfully strong effects -----------------------
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


# ~~ Sensitivity plot (homogeneous bias) -----------------------

# note: y-axis will automatically say "RR"; manually edited PDF to say "HR" but would
#  be easy to do it here
breaks.x1 = seq(1, 2, .1)
p = sens_plot( type = "line",
               method = "calibrated",
               give.CI = TRUE,
               q = log(0.90),
               tail = "below",
               Bmax = log(1.5),
               breaks.x1 = breaks.x1,
               dat = dat,
               yi.name = "yi",
               vi.name = "vi" )

p = p + ylab("Estimated proportion of studies with true HR < 0.9")

# put x-axis on log scale
# from inside sens_plot:
breaks.x2 = round(breaks.x1 + sqrt(breaks.x1^2 - 
                                     breaks.x1), 2)
p = p + scale_x_continuous( breaks = breaks.x1, 
                              sec.axis = sec_axis(~g(.),
                                                  name = "Minimum strength of both confounding RRs", 
                                                  breaks = breaks.x2),
                              trans = "log10")

my_ggsave( name = "flegal_line_plot.pdf",
           width = 4,
           height = 4,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir )

# ~  Heterogeneous bias -----------------------
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




# ANALYZE GBMC ---------------------------------------------------

dat = dfs[[2]]

# refit naive meta-analysis
( meta = rma.uni( yi = dat$yi, 
                  vi = dat$vi, 
                  method = "REML", 
                  knha = TRUE ) )

# ~ E-value for point estimate -----------------------
# HR for a rare outcome approximates an RR
evalue( est = RR( exp(meta$b) ),
        lo = RR( exp(meta$ci.lb) ),
        hi = RR( exp(meta$ci.ub) ) )

# ~  Homogeneous bias; percentage of meaningfully strong effects -----------------------
# uncorrected estimate of effects with RR < 0.90 AND Tmin, Gmin
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

# uncorrected estimate of effects with RR > 1.1 AND Tmin, Gmin
confounded_meta( method = "calibrated",
                 q = log(1.1),
                 r = 0.15,
                 tail = "above",
                 muB = 0,
                 sigB = 0,
                 dat = dat,
                 yi.name = "yi",
                 vi.name = "vi" )


# ~~ Sensitivity plot (homogeneous bias) -----------------------
breaks.x1 = seq(1, 2, .1)
p = sens_plot( type = "line",
               method = "calibrated",
               give.CI = TRUE,
               q = log(1.1),
               tail = "above",
               Bmax = log(1.5),
               breaks.x1 = breaks.x1,
               dat = dat,
               yi.name = "yi",
               vi.name = "vi" )

p = p + ylab("Estimated proportion of studies with true HR > 1.1")

# put x-axis on log scale
# from inside sens_plot:
breaks.x2 = round(breaks.x1 + sqrt(breaks.x1^2 - 
                                     breaks.x1), 2)
p = p + scale_x_continuous( breaks = breaks.x1, 
                            sec.axis = sec_axis(~g(.),
                                                name = "Minimum strength of both confounding RRs", 
                                                breaks = breaks.x2),
                            trans = "log10")


my_ggsave( name = "gbc_line_plot.pdf",
           width = 4,
           height = 4,
           .results.dir = results.dir,
           .overleaf.dir = overleaf.dir )

# ~ Heterogeneous bias (Supplement)  -----------------------
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




