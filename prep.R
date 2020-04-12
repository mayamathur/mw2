

####################### SET UP ####################### 

private.data.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Data (private)"
results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Results from R"
  
library(MetaUtility)
library(EValue)
library(ggplot2)
library(dplyr)
library(metafor)
library(testthat)

####################### KODAMA: HAND-ENTER DATA ####################### 

# data-entry checked 2020-4-11
# using the 16 studies in the first meta-analysis shown in Fig 3 (all-cause mortality; low vs. high CRF)
# using the adjusted RRs
study = c("Slattery 1988",
          "Hein 1992",
          "Aijaz 2008",
          "Villeneuve 1998",
          "Stevens 2002a",
          "Farrell 2002",
          "Stevens 2002b",
          "Sandvik 1993",
          "Kampert 1996",
          "Stevens 2004",
          "Laukken 2008",
          "Sawada 1999",
          "Erikksen 1998",
          "Arraiz 1992",
          "Gulati 2003",
          "Myers 2002")
rr = c( 1.47,
        1.53,
        1.25,
        1.22,
        1.92,
        1.25,
        2.23,
        1.38,
        1.42,
        2.58,
        1.83,
        2.38,
        2.09,
        2.17,
        1.89,
        2.96 )
hi = c( 1.55,
        1.85,
        1.47,
        2.57,
        2.59,
        1.8,
        3.3, 
        2.83,
        1.87,
        4.09,
        2.73,
        4.14,
        2.96,
        4.30,
        2.86,
        4.46 )

# already have RRs
library(MetaUtility)
d = MetaUtility::scrape_meta( type = "RR",
                              est = rr, 
                              hi = hi,
                              sqrt = FALSE )
d$study = study
names(d) = c("yi", "vi")

# we did not flip the point estimate signs
d$flipped = FALSE

# sanity check: matches reported 1.70 and 1.92 :)
( meta = rma.uni( yi = d$yi,
                  vi = d$vi,
                  method = "DL",
                  knha = FALSE ) )

exp(meta$b)
exp(meta$ci.ub)
 
# write the prepped data
setwd(private.data.dir)
write.csv(d, "kodama_prepped.csv")



# # for Shapiro test (seems okay :) )
# MetaUtility::prop_stronger(q = q,
#                            M = meta$b,
#                            t2 = meta$tau2,
#                            se.t2 = meta$se.tau2,
#                            tail = "above",
#                            dat = d, 
#                            yi.name = "yi",
#                            vi.name = "vyi" )
# 
# 
# sens_plot( type = "line",
#            q = q, 
#            yr = meta$b,
#            vyr = meta$se^2,
#            t2 = meta$tau2,
#            vt2 = meta$se.tau2^2,
#            Bmax = log(2.5) )
# 
# # y_R^t = meta$b - log(1.5)
# # here q < y_R^t
# # y_R^c > 0
# # so Phat under homogeneity is a *lower* bound and it will only look more robust
# #  with heterogeneity
# confounded_meta( 
#   q = q, 
#   r = 0.1,
#   muB = log(1.5),
#   sigB = 0,
#   yr = meta$b,
#   vyr = meta$se^2,
#   t2 = meta$tau2,
#   vt2 = meta$se.tau2^2,
#   tail = "above" )
# # this bias factor of 1.5 corresponds to RR strengths of 2.4
# # Ghat = 3.46
# 
# 
# # E-value for point estimate: 3.06 and 2.27
# evalues.RR(est = meta$b, hi = meta$ci.ub)



####################### FLEGAL: PREP AUTHOR-PROVIDED DATA ####################### 

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

# and reported number of studies
expect_equal( nrow(d), 140)

# recode effect sizes as positive since they are reporting a negative effect
d$yi = -1 * d$yi
d$flipped = TRUE 

# write the prepped data
setwd(private.data.dir)
write.csv(d, "flegal_prepped.csv")


