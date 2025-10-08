#Install the package:
setwd("/path/modelFits")
Rcpp::compileAttributes()
setwd("/path")
devtools::install("modelFits")

#Simulations
library(Rcpp)
source("modelFit.R")

methods <- "MCP"
#Specify the index of simulation dataset k.
k <- 1

#Specify ratios (lambda/delta)
ratios <- c(10 ,20, 30, 40)/1000

#delta: the parameter for subgroups.
Deltas = c(seq(0.1,2,by=0.1))

#Load data for simulations.
data_simu <- readRDS("/path/data_simu.rds")
data <- data_simu[[k]]
n <- dim(data$Z)[1]

BICs <- rep(1e+10, 16)
finals <- replicate(16, NULL, simplify = FALSE)

#This loop will find model estimates with best BIC.
for(i in 1:length(ratios)){
  
  ratio <- ratios[i]
  Lambda <- Deltas*ratio
  
  for(j in 1:length(Deltas)){
    BIC_new <- BICs
    
    if(j == 1){
      results <- modelFit(data=data,
                          #Tuning parameters:
                          Theta=1,
                          Lambda=Lambda[j],
                          Delta=Deltas[j],
                          Gamma=3,
                          method = methods,
                          warm=T,
                          beta0=NULL, mu0=NULL, epsilon0=NULL,v0=NULL,eta0=NULL,
                          phi = 1,
                          alpha = 1)
    }
    else{
      results <- modelFit(data=data,
                          #Tuning parameters:
                          Theta=1,
                          Lambda=Lambda[j],
                          Delta=Deltas[j],
                          Gamma=3,
                          method = methods,
                          warm=T,
                          phi = 1,
                          alpha = 1,
                          beta0=results$beta.est, mu0=results$mu.est, epsilon0=results$epsilon.est,v0=results$v.est,eta0=results$eta.est
      )
    }   
    print(paste("ratio:",ratio))
    print(paste("Delta:",Deltas[j]))
    
    K <- length(unique(round(results$mu.est, 5)))
    p <- sum(results$beta.est != 0)
    
    all_p <- 4999
    print(paste0("K: ",K,";p: ", p))
    
    #Calculate 
    #results$residual is the RSS of our fitted estimates.
    BIC_new[1:4] <- log(results$residual) + c(0.5, 1,2,3) * log(log(n))*log(all_p)/n*(K+p-1)
    BIC_new[5:8] <- log(results$residual) + c(0.5, 1,2,3) * log(log(n))*log(all_p)/n*(5*K+p-1)
    BIC_new[9:12] <- log(results$residual) + c(0.5, 1,2,3) * log(log(n))*log(all_p)/n*(10*K+p-1)
    BIC_new[13:16] <- log(results$residual) + c(0.5, 1,2,3) * log(log(n))*log(all_p)/n*(20*K+p-1)
    

    for(ii in 1:16){
      if(BIC_new[ii] < BICs[ii]){
        BICs[ii] <- BIC_new[ii]
        finals[[ii]] <- list(results, Delta = Deltas[j], ratio = ratio)
      }
    }
  }
  
  result <- list(BICs = BICs, finals = finals)
  file_name <- paste0("/path/result_",k,".rds")
  saveRDS(result, file_name)
}
