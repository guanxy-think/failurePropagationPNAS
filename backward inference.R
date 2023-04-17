###This code is used to run the backward model for congestion cascade
###Input of this program includes node failure times, node coordinates (for spatial dependence) and network adjacency matrix/null models (for functional depedence).
###This is the version executed with parallel computing.

##Setup the parallel computing environment
.libPaths("/thinklab/R_packages")
library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores,outfile="/thinklab/guan/logfile.txt")
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/thinklab/R_packages"))
clusterEvalQ(cl, setwd("/thinklab/guan"))
setwd("/thinklab/guan")
##

##Calculate node distance matrix
coordinates_node <- read.csv("/thinklab/guan/SiouxFalls_node.txt",header=T,sep="\t")
distance_matrix <- matrix(0,nrow(coordinates_node),nrow(coordinates_node))
for (i in 1:nrow(coordinates_node)) {
  for (j in 1:nrow(coordinates_node)) {
    distance_matrix[i,j] <- sqrt((coordinates_node$X_coordinate[i]-coordinates_node$X_coordinate[j])^2 + (coordinates_node$Y_coordinate[i]-coordinates_node$Y_coordinate[j])^2)
  }
}
##

##Get node failure time
##Node that does not fail has failure time NA 
failed_nodes <- read.csv("/thinklab/guan/failed_nodes.csv",header=T)
merged_node_info <- merge(coordinates_node,failed_nodes,by="node_ID",all.x=T)
node_failure_time <- merged_node_info$failure_time
##

##Custom functions
colProds <- function(ma) {
  return(apply(ma, 2, prod))
}
rowProds <- function(ma) {
  return(apply(ma, 1, prod))
}
##

##Calculate failure propagation probability between two nodes if their failure times have an interval of time_diff
propagationProbOneDay <- function(para,distance_matrix_data,adjacency_matrix_data,time_diff) {
  spatial_dependence <- para[1]*exp(-para[2]*distance_matrix_data)
  functional_dependence <- adjacency_matrix_data*para[3]
  one_over_lambda <- spatial_dependence + functional_dependence - spatial_dependence*functional_dependence    #Equation 11 in the paper
  propagation_prob <- exp(-((time_diff-1)*one_over_lambda)^para[4]) - exp(-(time_diff*one_over_lambda)^para[4])    #Equation 9 in the paper
  diag(propagation_prob) <- 0    #No propagation to oneself
  return(propagation_prob)
}
##

##Calculate the probability that each node fails on the given date_data
failureProb <- function(para,distance_matrix_data,adjacency_matrix_data,date_data) {
  propagationProb <- list()
  for (l in 1:3) {    #Consider 4 time steps (and right censoring), so the largest time interval between two failures is 3.
    propagationProb[[l]] <- propagationProbOneDay(para,distance_matrix_data,adjacency_matrix_data,l)
  }
  hazard <- matrix(NA,76,4)    #There are a total of 76 nodes in the network
  hazard[,1] <- 0
  hazard[which(date_data == 1),1] <- 1    #Seeds of failure propagation
  survival <- matrix(NA,76,4)
  survival[,1] <- 1
  densities <- hazard*survival
  for (i in 2:4) {    #Calculate failure hazard, survival probability and failure probability for every node at every time step.
    Prop_matrix <- matrix(0,76,76)
    for (j in 1:(i-1)) {    #Equation 4 in the paper
      Prop_matrix <- Prop_matrix + densities[,i-j]*propagationProb[[j]]
    }
    #Prop_matrix[which(Prop_matrix > 1)] <- 1
    noFailure_nodes <- colProds(1-Prop_matrix)    #Equation 3 in the paper (without external factors)
    failure_at_i <- 1 - noFailure_nodes
    hazard[,i] <- failure_at_i
    survival[,i] <- survival[,(i-1)]*(1-hazard[,(i-1)])    #Equation 2 in the ppaer
    densities[,i] <- hazard[,i]*survival[,i]
  }
  return(list(densities,survival,hazard))
}
##

##Formulation of the likelihood function as Equation 1 in the paper
likelihood_function <- function(para,distance_matrix_data,adjacency_matrix_data,date_data) {
  failure_info <- failureProb(para,distance_matrix_data,adjacency_matrix_data,date_data)
  failure_density <- failure_info[[1]]
  failure_survival <- failure_info[[2]]
  failure_hazard <- failure_info[[3]]
  
  daily_LL <- rep(NA,4)    #Write the likelihoods for failures that happen at different time steps separately
  daily_LL[1] <- sum(log(failure_density[,1][which(date_data == 1)]))
  for (i in 2:4) {
    node_indexes <- which(date_data == i)
    LL_day_i <- sum(log(failure_density[node_indexes,i]))
    daily_LL[i] <- LL_day_i
  }
  
  RC_node_index <- which(is.na(date_data))    #Right censoring for nodes that did not fail
  RC_LL <- sum(log(failure_survival[RC_node_index,4]*(1-failure_hazard[RC_node_index,4])))
  
  LL <- sum(daily_LL) + RC_LL    #Total likelihood 
  return(-LL)
}
##

##Implement parameter estimation
foreach(pp=1:1000) %dopar% {    #We constructed 1000 null models due to unknown network adjacency matrix.
  tryCatch(
    {
      network_file <- paste("/thinklab/guan/null model ", pp, ".csv",sep="")
      adjacency_matrix <-unname(as.matrix(read.csv(transportation_file,header=T)))
      
      solution <- optim(c(0.1,0.1,0.1,0.1),likelihood_function,
                        distance_matrix_data=distance_matrix,adjacency_matrix_data=adjacency_matrix,date_data=node_failure_time,
                        method="L-BFGS-B",lower=rep(0.000000001,4),control = list(maxit = 30000,trace = TRUE,ndeps=rep(0.0000000001,4),factr=0.45e9,REPORT=1),hessian = TRUE)
      parameter_estimates <- solution$par
      likelihood <- solution$value
      converge <- solution$convergence
      hessianM <- solution$hessian
      parameter_file <- paste(paste("parameters",pp, sep=""),".csv",sep="")
      likelihood_file <- paste(paste("likelihood",pp, sep=""),".csv",sep="")
      hessianM_file <- paste(paste("hessian",pp, sep=""),".csv",sep="")
      converge_file <- paste(paste("converge",pp, sep=""),".csv",sep="")
      write.csv(parameter_estimates,parameter_file,row.names=F,col.names=F)
      write.csv(likelihood,likelihood_file,row.names=F,col.names=F)
      write.csv(hessianM,hessianM_file,row.names=F,col.names=F)
      write.csv(converge,converge_file,row.names=F,col.names=F)
    },
    error=function(e) {mess<-paste(pp,e$message);cat(mess,file="error_log.txt",append=T,fill=T)},
    warning=function(w) {mess<-paste(pp,w$message);cat(mess,file="warning_log.txt",append=T,fill=T)},
    finally={
      write.csv(parameter_estimates,parameter_file,row.names=F,col.names=F)
      write.csv(likelihood,likelihood_file,row.names=F,col.names=F)
      write.csv(hessianM,hessianM_file,row.names=F,col.names=F)
      write.csv(converge,converge_file,row.names=F,col.names=F)
    }
  )
}
##

stopCluster(cl)
stopImplicitCluster()
