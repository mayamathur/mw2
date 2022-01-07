
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

### Sanity check: Reproduce reported pooled estimate ###
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


### Sanity checks: number of studies ###
# number of estimates matches theirs (Table 1)
expect_equal( nrow(d), 140)


### Sanity checks: number of articles ###
# apparent number of papers (94) doesn't quite match theirs (97)
# but the study variable seems a bit messy
expect_equal( length(unique(d$study)), 97)

# try to clean up the study variable
x = str_split(string = d$study, pattern = " ")
firstAuthor = unlist( lapply( x, function(.x) .x[1] ) )
nuni(firstAuthor)
# the years contain letters to distinguish multiple papers with same author and year
yearString = unlist( lapply( x, function(.x) .x[ length(.x) ] ) )
yearNumeric = parse_number(yearString)
# look at parsing failures
# these are studies for which no year was included in "study" variable
d$study[ is.na(yearNumeric) ]
# still doesn't quite match
expect_equal( nuni(myStudy), 97)

# write the prepped data
setwd(private.data.dir)
write.csv(d, "flegal_prepped.csv")


# GLOBAL BMI CONSORTIUM [GBC]: PREP AUTHOR-PROVIDED DATA---------------------------------------------------

setwd(private.data.dir)
setwd("From Emanuele")
d = read.csv("gbmic_rr_acm_by_cohort.csv")

# rename variables
names(d)[ names(d) == "study_name" ] = "study"
names(d)[ names(d) == "lnrr" ] = "yi"
d$vi = d$selnrr^2

### Sanity check: Reproduce reported pooled estimates ###
# but note, per eAppendix, that our analysis is different (e.g., not multivariate)

# sanity check for dose-response
# as expected, category 2 is the reference
d %>% group_by(bmicat) %>%
  summarize( expmu = exp(mean(yi)) )

# meta-analyze each BMI category (standard categories, not the finer-grained ones)
# this is close to eTable 4, "All Studies", final row
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

### Another sanity check: Reproduce reported pooled estimate for overweight only ###
# uses DL rather than REML with KNHA, as in our main analysis, just for the sanity check since
#  former is more common
( meta = rma.uni( yi = d$yi,
                  vi = d$vi, 
                  method = "DL",
                  knha = FALSE ) )
exp(meta$b)  # inverse because we reversed signs
exp(meta$ci.lb)
exp(meta$ci.ub)
# this is not equivalent to the analysis they did, but yields
#  HR = 1.10 [1.07, 1.12]
# reported: HR = 1.11 [1.10, 1.11] in Table 2 and eTable 4 with 189 studies
# similar given differences in analysis


# number of papers vs. reported in Table 2 legend (189)
# below, the "-3" is because, per authors:
# "One study (FINRISK) provided adjusted results for their 4 cohorts combined (FINRISK77, FINRISK82, FINRISK87 and North Karalia study), thus you needed to add 3 to the number of cohorts in the study labelled as "Finrisk Cohort 1977"
expect_equal( length(unique(d$study)), 189-3)


# write the prepped data
setwd(private.data.dir)
write.csv(d, "gbc_prepped.csv")





