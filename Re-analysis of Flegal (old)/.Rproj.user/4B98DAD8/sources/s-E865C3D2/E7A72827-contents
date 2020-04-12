
# look at Flegal's sensitivity to unmeasured confounding
setwd("~/Dropbox/Personal computer/Independent studies/Metawars/Manuscript Metawars/Other metawars/Metawars on obesity/From Flegal")

f = read.csv("flegal.csv")

library(MetaUtility)

f = cbind(f, 
          scrape_meta(type = "RR", est = f$HR, hi = f$High) )

# blank rows at end
f = f[ !is.na(f$EntrezUID), ]

# recode effect sizes as positive since they are reporting a negative effect
f$yi = -1 * f$yi

#temp = f[ f$Type.of.wt.ht.data == "measured", ]
library(testthat)
expect_equal( nrow(f), 140)


##### First Try Regular RE Meta-Analysis ######
library(metafor)
( meta = rma.uni( yi = f$yi,
                  vi = f$vyi, 
                  method = "REML",
                  knha = TRUE ) )
1/exp(meta$b)  # inverse because we reversed signs
1/exp(meta$ci.lb)
# reported: 
# k = 140, 0.94 [0.91, 0.96]
# I^2 = 85%
# for all studies
# quite close



##### Check Normality ######
std = (f$yi- c(meta$b)) / sqrt(c(meta$tau2) + f$vyi)
hist(std, breaks=20)
shapiro.test(std)


##### Significance Funnel ######
library(PublicationBias)
significance_funnel(yi = f$yi,
                    vi = f$vyi)

# DEBUG THE PLOT CODE 


# cannot really use selection models given substantial clustering
res = svalue( yi = f$yi,
              vi = f$vyi,
              q = 0, 
              clustervar = f$EntrezUID,
              model = "robust" )
res$svals

# meta-analyze only the nonaffirmatives
# 2-sided pval
# affirmative vs. non-affirmative
f$sei = sqrt(f$vyi)
f$pval = 2 * ( 1 - pnorm( abs(f$yi) / f$sei ) )
f$affirm = f$pval < 0.05 & f$yi > 0
table(f$affirm)

( meta.worst = rma.uni( yi = f[f$affirm == FALSE,]$yi,
                  vi = f[f$affirm == FALSE,]$vyi, 
                  method = "REML",
                  knha = TRUE ) )
1/exp(meta.worst$b)
1/exp(meta.worst$ci.lb)



##### Prop Stronger #####

# 42% stronger than 1.1 (protective effect of overweight)
# [0.34, 0.50]
prop_stronger( q = log(1.1), 
               M = meta$b, 
               t2 = meta$tau2,
               se.M = meta$se, 
               se.t2 = meta$se.tau2,
               tail = "above",
               R = 50,
               dat = f,
               yi.name = "yi", 
               vi.name = "vyi")

# 9% smaller than .9 (harmful effect of overweight)
# [0.04, 0.18]
prop_stronger( q = log(.9), 
               M = meta$b, 
               t2 = meta$tau2,
               se.M = meta$se, 
               se.t2 = meta$se.tau2,
               tail = "below",
               R = 50,
               dat = f,
               yi.name = "yi", 
               vi.name = "vyi")

# *** proportion as strong as observed estimate in the protective direction
# 17% [.11, .24]
# so not too negligible! 
prop_stronger( q = log(.95), 
               M = meta$b, 
               t2 = meta$tau2,
               se.M = meta$se, 
               se.t2 = meta$se.tau2,
               tail = "below",
               R = 50,
               dat = f,
               yi.name = "yi", 
               vi.name = "vyi")


##### Sensitivity Analysis for Unmeasured Confounding #####

# look at conservatism
meta$b - log(1.5) # estimated true mean if muB = log(1.5)
log(1.1)  # q

# q > estimated true mean
# and est true mean < 0
# so homo bias is a LOWER bound on strong true effects and is conservative as a skeptic


library(EValue)

( res.conf = confounded_meta( q = log(1.1),
                            r = 0.10,
                            muB = log(1.5),
                            sigB = 0,
                            yr = meta$b,
                            vyr = meta$se^2,
                            t2 = meta$tau2,
                            vt2 = meta$se.tau2^2,
                            tail = "above" ) )

# bias required to have >60% effects in the *harmful*, rather than protective, direction
# ** confounding RRs of only 1.45 in each study! 
confounded_meta( q = log(1),
                 r = 0.40,
                 sigB = 0,
                 yr = meta$b,
                 vyr = meta$se^2,
                 t2 = meta$tau2,
                 vt2 = meta$se.tau2^2,
                 tail = "above" )



# look at conservatism
meta$b - log(1.5) # estimated true mean if muB = log(1.5)
log(1.1)  # q

# q > estimated true mean
# and est true mean < 0
# so homo bias is a LOWER bound on strong true effects and is conservative as a skeptic

sens_plot( type = "line", 
           q = log(1.1), 
           muB = log(1.5),
           sigB = 0, 
           Bmax = log(1.8),
           yr = meta$b,
           vyr = meta$se^2,
           breaks.x1 = seq(.5, 5, .1),
           t2 = meta$tau2,
           vt2 = meta$se.tau2^2 )


sens_plot( type = "line", 
           q = log(1), 
           muB = log(1.5),
           sigB = 0, 
           Bmax = log(1.8),
           breaks.x1 = seq(.5, 5, .1),
           yr = meta$b,
           vyr = meta$se^2,
           t2 = meta$tau2,
           vt2 = meta$se.tau2^2 )









# seems borderline for this one