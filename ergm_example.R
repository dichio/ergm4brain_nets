#################
#
# ERGMs for large-scale brain networks: a minimal implementation
# See Dichio, V & De Vico Fallani, F, "Statistical modeling of complex brain networks: a maximum entropy approach", https://arxiv.org/abs/2209.05829
#
# Created: 12 Mar 2023 (Rstudio version: 2022.12.0.353, ergm version: 4.3.2)
#
# Author: Vito Dichio
#
# Last update: 12 Mar 2023
#
################

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ergm)


##### Step 0: Import brain network 
# Example data: filtered brain networks from functional EEG, subj=1, band=alpha // from Obando & De Vico Fallani, 2017 (http://dx.doi.org/10.1098/rsif.2016.0940)

brainet <- as.matrix(read.table(file="data/subj1_alpha.txt"))


##### Step 1: Build network
# Here simple networks: undurected, unweighted and with no self-loops

brainnet <- network(bnet, directed = FALSE, loops = FALSE, multiple = FALSE)
brainnet 


##### Step 2 (optional) : Compute statistics 
# For a comprehensive list of statistics, type ?ergm.terms  

summary(brainnet ~ edges + twopath + triangle + # some simple stats
                    gwnsp(decay=0.75,T) + gwdegree(decay=0.75,T) + gwesp(decay=0.75,T) + # some weighted stats
                    degree(0:10) # nodes with degree = (a,b)
                    ) 


##### Step 3 : Model fit

# Model of Simpson et al. (2011,2012)
m1 <- 'edges + gwesp(decay=0.75, fixed=T) + gwnsp(decay=0.75, fixed=T), 
        verbose=T,
        control=snctrl(init=c(-1.0,1.0,1.0))'

# Model "M1" of Obando & De Vico Fallani (2017)
m2 <- 'gwdegree(decay=0.75, fixed=T) + gwesp(decay=0.75, fixed=T), 
        verbose=T,
        constraints = ~ edges,
        control=snctrl(init=c(1.0,1.0))'

set.seed(160318)

# ergm - estimation (might take a while)
brainfit <- eval(parse(text = paste('ergm(brainnet ~', m1, ')',sep="")))

# Get results of the ergm fit
summary(bfit)


##### Step 4 (optional) : mcmc diagnostics
# Each statistic measured by its value relative to the corresponding value for the observed network

mcmc.diagnostics(bfit)


##### Step 5 : Goodness of Fit assessment 
# Solid lines: experimental; boxes: interquartile range; light gray lines: 95% confidence interval based on simulated values

bfit.gof <- gof(bfit) 
plot(bfit.gof)


##### Step 6 (optional): Simulate new networks
# Simulating nsim=100 new networks based on the fitted model

new_brains <- simulate(brainfit,nsim=100)
summary(new_brains)






##### Useful resources
# General: https://statnet.org/
# Introductory tutorial: https://statnet.org/workshop-ergm/ergm_tutorial.html
# Less introductory tutorial: https://statnet.org/workshop-advanced-ergm/advanced_ergm_tutorial.html

##### References
# [1] Hunter, D R et al. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. Journal of statistical software, 24(3), nihpa54860.
# [2] Schweinberger, M et al. (2020). Exponential-family models of random graphs: inference in finite, super and infinite population scenarios.
# [3] Simpson, S L et al. (2011). Exponential random graph modeling for complex brain networks. PloS one, 6(5), e20039.
# [4] Dichio, V & De Vico Fallani, F (2022). Statistical modeling of complex brain networks: a maximum entropy approach, arXiv preprint arXiv:2209.05829.