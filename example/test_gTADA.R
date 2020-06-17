###################################################
## Load codes
dirFile <- "../script/"
fileR <- dir(dirFile, ".R$")

for (ii in fileR){
    print(ii)
    source(paste0(dirFile, ii))
}


##################################################
## Get data
data <- read.table("../data/data_DD.txt", header = TRUE, as.is = TRUE)

ntrio = rep(4293, 2)
allDNData <- data[, paste0("dn_", c("damaging_", "lof_"), "DD")]
allMutData <- data[, paste0("mut_", c("damaging", "lof"))]

head(data.frame(allMutData, allDNData))

################################################
## Format data
inputData <- data.frame(Gene = data[, 1], allMutData, allDNData)

#################################################
## Get gene sets
fmrp <- read.table("../data/FMRP_targets.txt")
fmrp <- data.frame(V1 = fmrp[, 1], rep(1, dim(fmrp)[1]))

f1 <- data.frame(Gene = data[, 1], fmrpGene = rep(0, dim(data)[1]))
f1[is.element(f1[, 1], fmrp[, 1]), 2] <- 1

fmrp <- f1 
#write.table(fmrp, paste0("../data/fmrp_formatted.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
#fmrp <- read.table("../data/fmrp_formatted.txt", header = FALSE)

################################################
## Set parameters
nIteration = 1000


mcmcDD <- gTADA(modelName = DNextTADA,
                geneSet = data.frame(fmrp[, 2]), #Don't need gene names
                inputData = inputData, ## Input data should be formatted as above
                Ndn = array(c(ntrio)), #rep(ntrio, 1), ##Two de novo categories
                                        #                                    Ncase = array(ncase), #rep(N$ca, 1), ##Two case categories
                                        #                                   Ncontrol = array(ncontrol), #rep(N$cn, 1), ##Two control categories
                nIteration = nIteration ## Number of iterations: should be upto higher for real data
                )
##############################################
## Take a look at estimated parameters
mcmcDD$pars

### alpha0[2] shows the estimated alpha's value of the tested gene set 

##############################################
## Obtain full results
head(mcmcDD$dataPP)

#############################################
## Draw heatmap to check convergence
plotParHeatmap(pars = c('alpha0[1]', 'alpha0[2]'), mcmcResult = mcmcDD$gTADAmcmc)


##The section below is only for the comparison between extTADA/TADA and gTADA
## if you want to see differences between extTADA and gTADA results
data1 <- merge(data, mcmcDD$dataPP, by.x = 'geneName', by.y = 'Gene')
data1 <- merge(data1, fmrp, by.x = 'geneName', by.y = 'Gene')
data1[, 'fmrpGene'] <- as.factor(ifelse(data1$fmrpGene == 0, "No", "Yes"))

library('ggplot2')
p1 <- ggplot(data1, aes(PP.x, PP.y, col = fmrpGene)) + geom_point() +  xlab("extTADA's PP") + ylab("gTADA's PP") + geom_abline(slope = 1, intercept = 0)
p1
#save(mcmcDD, file = "mcmcDD_result.RData")
