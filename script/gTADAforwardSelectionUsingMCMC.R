##Forward selection using MCMC

gTADAforwardSelection <- function(modelName  , inputData,
                         Ndn = NULL,
                         Ncase = NULL, Ncontrol = NULL,
                    geneSet = NULL, geneSetName = NULL,
                    loglk0 = NULL, #If users have the result of a top gene set

                    thresholdTest = 2,
                    
                    nIteration = NULL, nIteration2 = NULL,
                    
                    nThin = NULL, nThin2 = NULL, nCore = 1, nChain = 1,
                         hyperBetaDN0 = NULL,
                         hyperBetaCC0 = NULL,
                         hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                         hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,

                         upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                         lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                         betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                    adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                    autoAdjustHyperBeta = FALSE,
                    drawHeatMap = FALSE, 
                    writeResult = FALSE,
                  resultDir = NULL, saveResult = FALSE)     {

         geneName <- data.frame(inputData[, "Gene"])
         dataDN <- data.frame(inputData[, grep("dn_", colnames(inputData))])
#         colnames(dataDN) <- paste0("dn_", 1:dim(dataDN)[2])

                mutRate <- data.frame(inputData[, grep("mut_", colnames(inputData))])
 #        colnames(mutRate) <- paste0("mut_", 1:dim(mutRate)[2])
         dataCCcase <- data.frame(inputData[, grep("cc_case", colnames(inputData))])
  #       colnames(dataCCcase) <- paste0("cc_case", 1:dim(dataCCcase)[2])
         dataCCcontrol <- data.frame(inputData[, grep("cc_control", colnames(inputData))])
   #      colnames(dataCCcontrol) <- paste0("cc_control", 1:dim(dataCCcontrol)[2])

          if (is.null(geneSetName))
              geneSetName <- paste0("GeneSet", 1:dim(geneSet)[2])
         if (is.null(resultDir))
             resultDir <- "."
         if (dim(dataCCcontrol)[2] == 0)
             dataCCcontrol = NULL
         if (dim(dataCCcase)[2] == 0)
             dataCCcase = NULL
         if (dim(dataDN)[2] == 0)
             dataDN = NULL
         if (dim(mutRate)[2] == 0)
             mutRate = NULL

         if (is.null(geneSet)){
             stop("\nGene sets have to be input")
         } else {
             geneSet <- data.frame(geneSet)
             if (dim(geneSet)[2] < 2)
                 stop("The forward selection needs at least two gene sets")
         }
         
             
##Start the selection process
    iStart = 1
         nGeneSet = dim(geneSet)[2]
         loglkresult = NULL
         addOption = NULL
         if (!is.null(loglk0)){
             loglkresult <- c(loglkresult, loglk0)
             addOption <- c(addOption, "YES")
             iStart = 2 ##Don't result the first gene set
         }
#######################Run for the first gene set
         cI <- 1
         for (ig in iStart:nGeneSet){
             message("=========\nRUNNING FOR THE GENE SET ", ig)
             mcmcResult <- gTADA(modelName = modelName,
                 geneSet = data.frame(geneSet[,  unique(c(cI, ig))]),
                                    inputData = inputData, ## Input data should be formatted as above
                                    Ndn = Ndn, #rep(ntrio, 1), ##Two de novo categories
                                    Ncase = Ncase, #rep(N$ca, 1), ##Two case categories
                                    Ncontrol = Ncontrol, #rep(N$cn, 1), ##Two control categories
                                    nIteration = nIteration, ## Number of iterations: should be upto higher for real data
                                    nThin = nThin ## Depend on users, but it can be floor(nIteration/1000)
                                   )
         
             loglknew <- mcmcResult$loglk
             loglkresult <- c(loglkresult, loglknew)

             if (ig == 1){

                 loglk0 = loglknew
  #               message("The gene set ", geneSetName[ig], " is added, and llk = ", loglknew)
                outResult = mcmcResult
                addOption <- c(addOption, "YES")
                 
             } else {

                 lCI <- as.numeric(mcmcResult$pars[, 3][-1])

                 if (((loglknew - loglk0) > thresholdTest) & (all(lCI > 0))) {
                     loglk0 = loglknew
                   
                     cI <- unique(c(cI, ig))
                     outResult = mcmcResult
                     addOption <- c(addOption, "YES")
 #                    print(data.frame(loglkresult, addOption))
                 } else {
                     addOption <- c(addOption, "NO")
#                     print(data.frame(geneSetName[1:ig], loglkresult, addOption))
                 }

             } #else

             print(data.frame(geneSetName[1:ig], loglkresult, addOption))


             if (saveResult){
                 write.table(mcmcResult$dataPP, paste0(resultDir, "/TempDataFRDgeneSet.", ig,
                                                       ".Option.", addOption[ig], ".txt"), row.names = FALSE, quote = FALSE)
                 save(mcmcResult, file = paste0(resultDir, "/TempgTADAresult.", ig,  ".Option.", addOption[ig], ".RData"))

                 write.table(data.frame(geneSetName[1:ig], loglkresult, addOption),
                             paste0(resultDir, "/TempGeneSetInfo.", ig, ".txt"),
                             row.names = FALSE, quote = FALSE)
                 
             }
         

         } ##for

    outLLK <- data.frame(geneSetName,                          loglkresult, addOption)
    print(outLLK)
    colnames(outLLK) <- c("GeneSet", "logLK", "Forward_Selection_Information")
    
       return(list(chosenList = cI, outputResult = outResult, outLogLK = outLLK))
     }
