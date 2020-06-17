gTADA <- function(modelName = NULL,
                  inputData,
                         Ndn = NULL,
                         Ncase = NULL, Ncontrol = NULL,
                  geneSet = NULL,
                  useBF = FALSE, # Using directly bayes factors
                  bf0 = NULL, # Bayes Factors from gTADA/extTADA without using gene sets

                  nIteration = NULL,
                  nIteration2 = NULL,
                    
                  nThin = NULL, nThin2 = NULL, nCore = 1, nChain = 1,
                  estimationMethodforUsingBF = c("optimizing", "mcmc"),
                         hyperBetaDN0 = NULL,
                         hyperBetaCC0 = NULL,
                         hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                         hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,
                  alpha01Mean = 0, alpha01SD = 2,
                  rho0 = NULL, nu0 = NULL,

                         upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                         lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                         betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                    adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                    autoAdjustHyperBeta = FALSE,
                    drawHeatMap = FALSE, 
                    writeResult = FALSE,
                  resultDir = NULL,
                  seed0 = sample.int(.Machine$integer.max, 1))
     {


###MCMC process
         geneName <- data.frame(inputData[, "Gene", drop = FALSE])
         dataDN <- data.frame(inputData[, grep("dn_", colnames(inputData)), drop = FALSE])
#         colnames(dataDN) <- paste0("dn_", 1:dim(dataDN)[2])
         mutRate <- data.frame(inputData[, grep("mut_", colnames(inputData)), drop = FALSE])
 #        colnames(mutRate) <- paste0("mut_", 1:dim(mutRate)[2])
         dataCCcase <- data.frame(inputData[, grep("cc_case", colnames(inputData)), drop = FALSE])
  #       colnames(dataCCcase) <- paste0("cc_case", 1:dim(dataCCcase)[2])
         dataCCcontrol <- data.frame(inputData[, grep("cc_control", colnames(inputData)), drop = FALSE])
   #      colnames(dataCCcontrol) <- paste0("cc_control", 1:dim(dataCCcontrol)[2])

         if (!is.null(resultDir))
             resultDir <- "."
         if (dim(dataCCcontrol)[2] == 0)
             dataCCcontrol = NULL
         if (dim(dataCCcase)[2] == 0)
             dataCCcase = NULL
         if (dim(dataDN)[2] == 0)
             dataDN = NULL
         if (dim(mutRate)[2] == 0)
             mutRate = NULL

         if (is.null(nIteration)){
             nIteration <- 5000
             message("No interation input; therefore, 5000 is used\n")
                 
         }
         if (is.null(nThin))
             nThin <- floor(nIteration/1000)
#                  message("print(head(dataDN)) :", head(dataDN))

         ##There are two sections: using BFs (faster) or not using BFs (slow, using MCMC)
         if (useBF){
             message("\n========Bayes Factors are used in the estimation processs==========\n")

             estimationMethodforUsingBF <- match.arg(estimationMethodforUsingBF)
             
             pH0 <- p0.temp.dn <- p0.temp.cc <- rep(1, dim(inputData)[1])
             if (!is.null(dataDN)){
                 p0.temp.dn <- pH0.model.DN(data.dn = dataDN, data.mut = mutRate, n.dn = Ndn)
                 pH0 <- pH0*apply(p0.temp.dn, 1, prod)
             }
             if (!is.null(dataCCcase)){
                 data.cc <- data.frame(dataCCcase, dataCCcontrol)
                 p0.temp.cc <- pH0.model.CC(data.cc = data.cc,
                                        n.cc  = list(ncase = Ncase, ncontrol = Ncontrol),
                                        rho0 = rho0, nu0 = nu0)
                 pH0 <- pH0*apply(p0.temp.cc, 1, prod)
             }
             Ngs = dim(geneSet)[2]
             Ngene = dim(inputData)[1]
             

             modelData <- list(NN = Ngene, #Gene numbers
                   Ngs = Ngs,
                   bf0 = bf0,
                   pH0 = pH0,
                   geneSet = data.frame(geneSet),
                   sigmaPrior = alpha01SD, #sigmaPrior,
                   alpha01Mean = alpha01Mean,
                   alpha01SD = alpha01SD)
             if (estimationMethodforUsingBF == "mcmc"){
                 mcmcData <- stan(#model_code = modelName,
                     model_code = gTADAfromBF,
                     pars = c("alpha0"),
                     data = modelData, ##Model data as described above
                     iter = nIteration, chains = nChain, cores = nCore, thin = nThin)
                 pars0 <- estimatePars(mcmcResult = mcmcData,
                                    pars = paste0("alpha0[", 1:(1+dim(geneSet)[2]), "]"))# pars = c('alpha0[1]', 'alpha0[2]'))
             }
             if (estimationMethodforUsingBF == "optimizing"){
                 #seed0 = sample.int(.Machine$integer.max, 1)
                 tempModel <- stan_model(model_code = gTADAfromBF)
                 mcmcModel2 <- optimizing(object = tempModel,
                                          data = modelData, hessian = TRUE,
                                          seed = seed0,
                                          iter = nIteration)
                 par1 <- mcmcModel2$par
                 par1 <- par1[grep('alpha0', names(par1))]
                 ##Get Hessian
                 tempH <- (mcmcModel2$hessian)
                 se.mle = sqrt(-diag(solve(tempH))) #a1$hessian)))
                 par2 <- cbind(par1, par1 - 1.96*se.mle, par1 + 1.96*se.mle, se.mle,
                               par1/se.mle,
                               pnorm(abs(par1/se.mle), lower.tail = FALSE)) #p values for one side

  
                 
                 mcmcData <- NULL
                 message("\nNo MCMC results will be generated\n")
                 pars0 <- par2
             }
             
             e.alpha0 <- pars0[, 1]
             pi0 = e.alpha0[1]
             for (i in 1:dim(geneSet)[2])
                 pi0 = pi0 + e.alpha0[i + 1] * geneSet[, i]
             pi0 <- exp(pi0)/(1 + exp(pi0))
             dataFDR <- inputData
             dataFDR$BF <- bf0
             dataFDR$PP <- bf0*pi0/(1 - pi0 + pi0 * bf0)

             colnames(pars0) <- c("Mode", "lCI", "uCI", "pSD", "zScore", 
        "pValue")

             pars0 <- cbind(rownames(pars0), pars0)
             colnames(pars0)[1] <- c("Parameter")

             
             loglk = NULL


         } else {
         
         message("MCMC is running")
         mcmcData <- gTADAmcmc(modelName = modelName,
                               geneSet = geneSet,
                                 dataDN = dataDN,
                               mutRate = mutRate, Ndn = Ndn,
                                 dataCCcase = dataCCcase, dataCCcontrol = dataCCcontrol, Ncase = Ncase, Ncontrol = Ncontrol,
                                 nIteration = nIteration, nIteration2 = nIteration2,
                               alpha01Mean = alpha01Mean, alpha01SD = alpha01SD,
                                 nThin = nThin, nThin2 = nThin2, nCore = nCore, nChain = nChain,
                                 hyperBetaDN0 = hyperBetaDN0, hyperBetaCC0 = hyperBetaCC0,
                                 hyper2GammaMeanDN = hyper2GammaMeanDN, hyper2BetaDN = hyper2BetaDN, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                                 hyper2GammaMeanCC = hyper2GammaMeanCC, hyper2BetaCC = hyper2BetaCC,
                                 upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                                 lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                                 betaPars = betaPars, #Adjust beta's values: should be default
                                 adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                    autoAdjustHyperBeta =  autoAdjustHyperBeta)



###############Estimate genetic parameters
         message("\nEstimate genetic parameters from MCMC results")
         ############Add parameters
         parNames <- c("alpha0[1]")
         if (!is.null(geneSet))
             parNames <- paste0('alpha0[', 1:(dim(geneSet)[2]+1), ']')
         if (!is.null(dataCCcontrol)){
             parNames <- c(parNames, paste0('hyperGammaMeanCC[', 1:dim(dataCCcontrol)[2], ']'))
             parNames <- c(parNames, paste0('hyperBetaCC[', 1:dim(dataCCcontrol)[2], ']'))
         }
         if (!is.null(dataDN)){
             parNames <- c(parNames, paste0('hyperGammaMeanDN[', 1:dim(dataDN)[2], ']'))
             parNames <- c(parNames, paste0('hyperBetaDN[', 1:dim(dataDN)[2], ']'))
         }
         

         message("Estimate parameters from MCMC results\n")
         
         pars0 <- estimatePars(pars = parNames,
                     mcmcResult = mcmcData)
         print(pars0)
         pars1 <- as.numeric(pars0[, 1])
         names(pars1) <- rownames(pars0)

          parsFDR <- list(alpha0 = as.numeric(pars1[grep("alpha0", names(pars1))]),
             gammaMeanDN = as.numeric(pars1[grep("hyperGammaMeanDN", names(pars1))]),
                         betaDN = as.numeric(pars1[grep("hyperBetaDN", names(pars1))]),
                         gammaMeanCC = as.numeric(pars1[grep("hyperGammaMeanCC", names(pars1))]),
                         betaCC = as.numeric(pars1[grep("hyperBetaCC", names(pars1))]),

                         nfamily = Ndn,
                         ncase = Ncase,
                         ncontrol = Ncontrol
                         )

         message("\nCalculate posterior probabilities and FDRs")
         colnames(geneName) <- "Gene"
         dataOut <- calculateFDR(pars = parsFDR, geneSet = geneSet,
                                 dnData = dataDN, mutData = mutRate,
                                 caseData = dataCCcase, controlData = dataCCcontrol,
                                 geneName = geneName)

         dataFDR <- dataOut$dataFDR
         loglk <- dataOut$loglk
         message("\nDraw heatmaps")
         outTime <-  format(Sys.time(), "%a_%b_%d_%H_%M_%S_%Y")
         if (drawHeatMap) {
             pdf(paste0(resultDir, "/heatMap", outTime, ".pdf"))
                          allHyperGamma <- rownames(pars0[grep("hyperGammaMean", rownames(pars0)), ])
             for (i1 in 1:length(allHyperGamma))
                 plotParHeatmap(pars = c("alpha0[1]", allHyperGamma[i1]), mcmcResult = mcmcData)
             dev.off()
         }
         if (writeResult){
             write.table(dataFDR, paste0(resultDir, "/Result_extTADA_PosteriorAndFDR", outTime, ".txt"),
                         row.names = FALSE, quote = FALSE)
             write.table(pars0,   paste0("Result_extTADA_estimatedPars", outTime, ".txt"), quote = FALSE)
         }
########################
         message(paste0("\nThe analysis is completed.\nIf you want to analyse steps seperately, please take a look at the example in the manual"))

         pars0 <- cbind(rownames(pars0), pars0)
             colnames(pars0)[1] <- c("Parameter")

         }
         return(list(dataPP = dataFDR, pars = pars0, gTADAmcmc = mcmcData, loglk = loglk))


          }

