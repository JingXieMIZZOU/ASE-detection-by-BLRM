library("lme4")
library("mvtnorm")
library("MCMCpack")
library("nloptr")
library("ggplot2")
library("dplyr")
library("stats4")

# load a sample dataset
mysample<-read.csv(file="example.csv")

# load self-defined functions
source("ASE_paraEst.R")
source("ASE_LogPPs.R")
source("ASE_detect.R")

# detect ASE, SNP variation and GS effect
result<-detection(mysample)
