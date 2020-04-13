
rm( list = ls() )

####################### SET UP ####################### 

code.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Code (git)"
private.data.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Data (private)"
results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Results from R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/Metasens tutorial"

setwd(code.dir)
source("helper.R")

setwd(private.data.dir)
dfs = list( read.csv("flegal_prepped.csv"),
                 read.csv("gbc_prepped.csv"),
                 read.csv("kodama_prepped.csv") )


# # choose bias SD by looking at a particular cross-section of the plot (muB)
# B_pct( muB = log(2), sigB = 0.01 )
  
digits = 2
varB.ratio = c(0, 0.8)
boot.reps = 2000  # ~~ increase later!

# start results df from scratch
rm("resE")


library(MetaUtility)
library(ggplot2)
library(dplyr)
library(metafor)
library(testthat)

# use local versions of EValue code
detach(package:EValue, unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/EValue package/evalue_package_git/EValue/R")
source("EValue.R")

####################### NORMALITY ASSUMPTION ####################### 

# investigate normality for each meta-analysis
# all seem nicely normal :)
for ( dat in dfs ){
  
  dat$calib = calib_ests( yi = dat$yi,
                          sei = sqrt(dat$vi) )
  
  print( shapiro.test(dat$calib) )
  
  qqnorm(dat$calib)
  qqline(dat$calib, col = "red")
}


# to do: maybe just put the parametric part on the existing calib plot?


####################### MAIN ANALYSES ####################### 


# maybe use same Bmax throughout so that the plots are more clear?

##### Flegal #####
dat = dfs[[1]]
meta.name = "Flegal"
ql = list( log(1), log(.9) )
tail = c( "above", "below" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 1.3
breaks.x1 = seq(1, Bmax, .1)
est.type = "HR"

setwd(code.dir)
source("general_analysis.R")

##### GBC #####

dat = dfs[[2]]
meta.name = "GBC"
ql = list( log(1), log(1.1) )
tail = c( "below", "above" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 1.5
breaks.x1 = seq(1, Bmax, .1)
est.type = "HR"

setwd(code.dir)
source("general_analysis.R")

##### Kodama #####
dat = dfs[[3]]
meta.name = "Kodama"
ql = list( log(1), log(1.5), log(1.1) )  # different from others because otherwise Phat is too close to 1
tail = c( "below", "above", "above" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 3
breaks.x1 = seq(1, 3, .5)
est.type = "RR"

setwd(code.dir)
source("general_analysis.R")


# # with a different choice of Bmax and 
# dat = dfs[[3]]
# meta.name = "Kodama"
# ql = list( log(1.1) )  # different from others because otherwise Phat is too close to 1
# tail = c( "above" )
# take.exp = TRUE
# r = 0.10
# varB.ratio = c(0, 0.8)
# Bmax = 3
# breaks.x1 = seq(1, 3, .5)
# est.type = "RR"
# 
# setwd(code.dir)
# source("general_analysis.R")
# 
# # choose a good threshold for Kodama
# Phat_causal( .q = log(1.1),
#           .B = 2,
#           .calib = calib_ests( dfs[[3]]$yi, sqrt(dfs[[3]]$vi) ), # assumed on log scale
#           .tail = "above",
#           .causative = TRUE,
#           
#           .give.CI = FALSE,
#           .R = 500,
#           .dat = dfs[[3]] )


##### Look At Aggregated Results #####
resE

setwd(results.dir)
write.csv(resE, "all_results.csv")
  