# ASE detection by BLRM

There are basically 4 steps to do ASE detection using our proposed method.
1. read in data and re-structurize it.
2. Hyper parameter estimation.
3. Posterior probability computation.
4. FDR control and ASE detection.

Inside R file "ASE_paraEst.R":
GDD(i,data) is a function to re-structurize the data to a structure that contains necessary information for analysis, i.e., SNPs, Replicates, counts  from maternal allel (YI) and total counts (NI). This structure is friendly to GLMM fitting in glmer() function in "lme4" package. This function can be easily modified to be applied to other similar real data structures.             
filter(data) is a function to filter the genes with computational problems; para.est(data) is a function to estimate hyper parameters.

Inside R file "ASE_LogPPs.R":
logpost11G(D,par), logpost21G(D,par),logpost31G(D,par),logpost41G(D,par) calculate the log posterior probability of model 1, model 2, model 3, model 4 respectively, for one specifc gene.
PPsFUN<-function(par,index,data) returns log posterior probabilites of all the genes.

Inside R file "ASE_detect.R":
detection(data) returns a list contains  genes which has gene, effect, SNP effect, and gene & SNP effect.

"ASE_example.R" provides an example of ASE detection using proposed BLRM method.
