
# PRELIMINARIES ---------------------------------------------------

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
# renv::restore()

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
                                    pattern = "Linked to OSF \\(MW2\\)/Code \\(git\\)",
                                    replacement = "Results from R" ) 

private.data.dir = str_replace_all( string = here(),
                                    pattern = "Linked to OSF \\(MW2\\)/Code \\(git\\)",
                                    replacement = "Data (private)" )


# get helper fns
setwd(code.dir)
source("helper.R")

# no sci notation
options(scipen=999)



# FLEGAL: PREP AUTHOR-PROVIDED DATA ---------------------------------------------------

setwd(private.data.dir)
setwd("From Flegal")
d = read.csv("flegal_raw_data.csv")


d = cbind(d, 
          MetaUtility::scrape_meta(type = "RR",
                                   est = d$HR,
                                   hi = d$High) )

names(d)[ names(d) == "Author.year" ] = "study"
names(d)[ names(d) == "vyi" ] = "vi"

# remove blank rows at end
d = d[ !is.na(d$EntrezUID), ]

# sanity check vs. reported
( meta = rma.uni( yi = d$yi,
                  vi = d$vi, 
                  method = "DL",
                  knha = FALSE ) )
exp(meta$b)  # inverse because we reversed signs
exp(meta$ci.ub)
# reported: 
# k = 140, 0.94 [0.91, 0.96]
# I^2 = 85%
# for all studies
# quite close

# and reported number of studies (Table 1)
expect_equal( nrow(d), 140)

# reported number of papers doesn't quite match theirs
# their study variable seems a bit messy?
expect_equal( length(unique(d$study)), 97)


# write the prepped data
setwd(private.data.dir)
write.csv(d, "flegal_prepped.csv")


# GLOBAL BMI CONSORTIUM [GBC]: PREP AUTHOR-PROVIDED DATA---------------------------------------------------

setwd(private.data.dir)
setwd("From Emanuele")
d = read.csv("gbmic_rr_acm_by_cohort.csv")


names(d)[ names(d) == "study_name" ] = "study"
names(d)[ names(d) == "lnrr" ] = "yi"
d$vi = d$selnrr^2

# sanity check for dose-response
# as expected, category 2 is the reference
d %>% group_by(bmicat) %>%
  summarize( expmu = exp(mean(yi)) )

# meta-analyze each BMI category
for ( i in c(1,3,4,5,6) ) {
  meta = rma.uni( yi = yi,
                 sei = selnrr,
                 knha = TRUE,
                 data = d %>% filter(bmicat == i) )
  
  cat( paste( "/n/n~~~~~~~~~~~ BMI CATEGORY", i, ":" ) )
  print( paste( round( exp(meta$b), 2 ),
                "[",
                round( exp(meta$ci.lb), 2 ),
                ", ",
                round( exp(meta$ci.ub), 2 ),
                "]" ) )
}
  

# for us, keep only the overweight category
d = d[ d$bmicat == 3, ]

# sanity check
( meta = rma.uni( yi = d$yi,
                  vi = d$vi, 
                  method = "DL",
                  knha = FALSE ) )
exp(meta$b)  # inverse because we reversed signs
exp(meta$ci.lb)
exp(meta$ci.ub)
# this is not equivalent to the analysis they did, but yields
#  HR = 1.07 [1.05, 1.09]
# similar to their 1.11 [1.10, 1.11] in Table 2

# number of studies
nrow(d)

# number of papers vs. reported in Table 2 legend
# also does not quite agree
expect_equal( length(unique(d$study)), 189)


# recode effect sizes as positive since they are reporting a negative effect
#d$flipped = FALSE 

# write the prepped data
setwd(private.data.dir)
write.csv(d, "gbc_prepped.csv")





