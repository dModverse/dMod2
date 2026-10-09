# Phytoene conversion p -> pf -> z by the enzyme Fad: fits with and without error model
# Reference: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0187628

## Library dependencies and plot theme --------------------------------------------
library(deSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(cOde)
library(dMod2)

qplot <- function(...) ggplot2::qplot(...) + theme_few() + scale_color_colorblind() + scale_fill_colorblind()

do.fit <- TRUE
do.fiterrors <- TRUE # fit an error model or take sigma from data

do.compile <- TRUE

do.startbest <- TRUE #use bestfit parameters

source("inst/examples/example_phytoene/setup.R")


## Howto proceed -------------------------------------------------
print("now predicting and plotting")
# Data times
timesD <- unique(sort(c(unlist(sapply(data, function(d) c(0, d$time))),31)))
times <- sort(union(timesD,seq(min(timesD), max(timesD), len=101)))

# Fitting
ini <- pouter

if(do.startbest){
  bestfit <- readRDS("inst/examples/example_phytoene/bestfit_pub.Rds")
  common <- intersect(names(ini),names(bestfit))
  ini[common] <- bestfit[common]
}
center <- ini
fixedVar <- NULL
cat("test obj with ini:",obj(ini, fixed = fixedVar)$value)

if(do.fit){
  print("start fit")

  myfit <- trust(obj, ini, fixed = fixedVar, rinit=1, rmax=10, iterlim=500,blather=TRUE)
  print(paste("fit",myfit$converged,"with",myfit$iterations,"iterations to",myfit$value))
  
  cat("data =", normL2(data, g*x*p)(myfit$argument)$value,"\n")

for(con in conditions){
  if(do.fiterrors)
    cat(con," = ",normL2(data, g*x*p,errmodel = err)(myfit$argument, conditions = con)$value,"\n")
  else
    cat(con," = ",normL2(data, g*x*p)(myfit$argument, conditions = con)$value,"\n")
}

prediction <- (g*x*p)(times, myfit$argument, conditions = conditions)

plotCombined( prediction, data) + theme_few()

dataTot <- lapply(data, function(d){
  td <- unique(d$time)
  ptot <- do.call(rbind,lapply(td, function(t){
    tot = sum(subset(d, time==t)$value)
    sig = sqrt(sum(subset(d, time==t)$sigma^2))
    if(TRUE){
      tot = sum(c(mean(subset(d, time ==t & name == "op")$value),mean(subset(d, time ==t & name == "opf")$value),mean(subset(d, time ==t & name == "oz")$value)), na.rm = TRUE)
    }
    data.frame(time=t, name="ptot_obs",value=tot, sigma=sig) #, condition=d$condition[1])
  }))
  rbind(d, ptot)
})

pobs <- plotCombined( prediction, dataTot, name%in%names(observables), facet = "wrap") 
print(pobs)

}


# More data
if(FALSE){
  # Datasets taken at different times differ in activity and setup; mutantB (_B, parB) is a less active enzyme
  # standard_*: p(0) > 0, pf(0) = 0; pfl_*: p(0) = 0, pf(0) > 0
  # Parallel measurements and dose responses (0 and 15 min per dose): see column global
  
  
conditionsTot <- c(
  "standard",

  #  "standard_30C", # other temperature

  "standard_nP",
  "pds_nP", # new enyzme added after ~30min

  "pfl_standard",

  "standard_WT",
  "standard_B",

  "pfl_par",
  "pfl_parB", 
  "standard_par", 
  "standard_parB",
  "standard_parh"
)

datasheet <- read.table(paste0("inst/examples/example_phytoene/datafile_err_sp_badn.csv"),sep=",",header=TRUE) # with columns condition, name, time, value, 
data <- lapply(intersect(conditionsTot,levels(datasheet$condition)), function(mycondition) subset(datasheet, condition == mycondition, select=-condition))
names(data) <- intersect(conditionsTot,levels(datasheet$condition))

dataTot <- lapply(data, function(d){
  td <- unique(d$time)
  ptot <- do.call(rbind,lapply(td, function(t){
    tot = sum(subset(d, time==t)$value)  
    sig = sqrt(sum(subset(d, time==t)$sigma^2))
    data.frame(time=t, name="ptot_obs",value=tot, sigma=sig,global = "ptot",n=NaN) 
  }))
  rbind(d, ptot)
})

}
