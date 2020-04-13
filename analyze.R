
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
boot.reps = 300  # ~~ increase later!

# start results df from scratch
rm("resE")


# maybe use same Bmax throughout so that the plots are more clear?

# analyze Kodama
analyze_one_meta( dat = dfs[[3]],
                  meta.name = "Kodama",
                  ql = list( log(1), log(1.5) ),  # different from others because otherwise Phat is too close to 1
                  tail = c( "below", "above" ),
                  take.exp = TRUE,
                  boot.reps = 2000,
                  digits = digits,
                  r = 0.10,
                  varB.ratio = c(0, 0.8),
                  Bmax = 3,
                  breaks.x1 = seq(1, 3, .5),
                  est.type = "RR" )


# analyze Flegal
analyze_one_meta( dat = dfs[[1]],
                  meta.name = "Flegal",
                  ql = list( log(1), log(.9) ),
                  tail = c( "above", "below" ),
                  take.exp = TRUE,
                  boot.reps = 2000,
                  digits = digits,
                  r = 0.10,
                  varB.ratio = c(0, 0.8),
                  Bmax = 1.5,
                  breaks.x1 = seq(1, 1.5, .1),
                  est.type = "HR" )

# analyze GBC
analyze_one_meta( dat = dfs[[2]],
                  meta.name = "GBC",
                  ql = list( log(1), log(1.1) ),
                  tail = c( "below", "above" ),
                  take.exp = TRUE,
                  boot.reps = 2000,
                  digits = digits,
                  r = 0.10,
                  varB.ratio = c(0, 0.8),
                  Bmax = 1.5,
                  breaks.x1 = seq(1, 1.5, .1),
                  est.type = "HR" )


resE

setwd(results.dir)
write.csv(resE, "all_results.csv")
  