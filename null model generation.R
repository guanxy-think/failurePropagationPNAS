rm(list=ls())
library(matrixStats)
library(poweRlaw)
library(numDeriv)
library(igraph)

BPR_variables <- read.csv("/benchmark model/SiouxFallsDataBPR.csv",header=T,colClasses = c("character",rep(NA,5)))
link_set <- unname(sapply(BPR_variables$Link_from_to,sub,pattern="\\.",replacement=','))
##for transportation network
##generating null models
for (null_model_index in 1:100) {
  degrees <- rpldis(length(link_set),1,4)
  while (sum(degrees)%%2 != 0) {
    degrees <- rpldis(length(link_set),1,4)
  }
  rowSum <- degrees
  A.transp <- matrix(0,length(link_set),length(link_set))
  cou <- 0
  while (max(rowSum) > 0) {
    max_indexes <- which(rowSum == max(rowSum))
    for (i in max_indexes) {
      if (rowSum[i]>0) {
        cand_index <- which(rowSum > 0)
        cand_index1 <- cand_index[cand_index != i]
        chosen_index <- sample(cand_index1,rowSum[i]) #problem: if cand_index1 is a number, it will sample [0,cand_index1]
        #print(cand_index)
        #print(cand_index1)
        #print(chosen_index)
        if (length(cand_index1)==1) {chosen_index = cand_index1}
        for (j in chosen_index) {
          #print(j)
          #print(rowSum[j])
          A.transp[i,j] = 1
          A.transp[j,i] = 1
          rowSum[j] <- rowSum[j] - 1
          #print(min(rowSum))
        }
      }
      #print(cou)
      #print(min(rowSum))
      #cou <- cou+1
      rowSum[i] <- 0
    }
  }
  file_name <- paste("C:/Users/user/Google Drive/dissertation/simulation paper/benchmark model/null models/transportation", null_model_index,".csv", sep="")
  colnames(A.transp) <- link_set
  rownames(A.transp) <- link_set
  write.csv(A.transp,file_name)
}
# sum(A.transp)
# sum(degrees)
# sum(rowSum)
# isSymmetric(A.transp)
# which(colSums(A.transp) != rowSums(A.transp))
# which(colSums(A.transp) != degrees)
# which(rowSum == -1)


##Null models generation for interdependent power and transportation networks.
for (m in 1:100) {
#for transportation network
	repeat{
		degrees <- rpldis(2166,1,4)
		if (sum(degrees)%%2==0 & max(degrees) >=10) break
	}
	rowSum <- degrees
	colSum <- degrees
	A.transp.null <- matrix(0,2166,2166)
	cou <- 0
	while (max(rowSum) > 0) {
		max_indexes <- which(rowSum == max(rowSum))
		for (i in max_indexes) {
			if (rowSum[i]>0) {
				cand_index <- which(rowSum > 0)
				cand_index1 <- cand_index[cand_index != i]
				chosen_index <- sample(cand_index1,rowSum[i])
				#print(cand_index)
				#print(cand_index1)
				#print(chosen_index)
				if (length(cand_index1)==1) {chosen_index = cand_index1}
				for (j in chosen_index) {
					#print(j)
					#print(rowSum[j])
					A.transp.null[i,j] = 1
					A.transp.null[j,i] = 1
					rowSum[j] <- rowSum[j] - 1
					#print(min(rowSum))
				}
			}
			#print(cou)
			#print(min(rowSum))
			#cou <- cou+1
			rowSum[i] <- 0
		}
	}
	sum(A.transp.null)
	sum(degrees)
	sum(rowSum)
	isSymmetric(A.transp.null)
	which(colSums(A.transp.null) != rowSums(A.transp.null))
	which(colSums(A.transp.null) != degrees)
	which(rowSum == -1)
	#for power network
	repeat {
		degrees_p <- ceiling(rexp(2166, rate = 0.5)) 
		if (sum(degrees_p)%%2==0 & max(degrees_p) >=10) break
	}
	rowSum_p <- degrees_p
	colSum_p <- degrees_p
	A.power.null <- matrix(0,2166,2166)
	while (max(rowSum_p) > 0) {
		max_indexes <- which(rowSum_p == max(rowSum_p))
		for (i in max_indexes) {
			if (rowSum_p[i]>0) {
				cand_index <- which(rowSum_p > 0)
				cand_index1 <- cand_index[cand_index != i]
				chosen_index <- sample(cand_index1,rowSum_p[i]) #problem: if cand_index1 is a number, it will sample [0,cand_index1]
				#print(cand_index)
				#print(cand_index1)
				#print(chosen_index)
				if (length(cand_index1)==1) {chosen_index = cand_index1}
				for (j in chosen_index) {
					#print(j)
					#print(rowSum[j])
					A.power.null[i,j] = 1
					A.power.null[j,i] = 1
					rowSum_p[j] <- rowSum_p[j] - 1
					#print(min(rowSum))
				}
			}
			#print(cou)
			#print(min(rowSum))
			#cou <- cou+1
			rowSum_p[i] <- 0
		}
	}
	isSymmetric(A.power.null)
	which(colSums(A.power.null) != degrees_p)
	which(rowSum_p == -1)
	#Interdependencies
	A.inter.null <- matrix(0,2166,2166)
	index_transp_10 <- which(degrees >= 10)
	index_power_10 <- which(degrees_p >= 10)
	for (i in index_transp_10) {
		for (j in index_power_10) {
			A.inter.null[i,j] <- 1
		}
	}
	file_transp <- paste(paste("/simulation data/null models/transportation_network",m, sep=""),".csv",sep="")
	file_power <- paste(paste("/simulation data/null models/power_network",m, sep=""),".csv",sep="")
	file_inter <- paste(paste("/simulation data/null models/interdependencies",m, sep=""),".csv",sep="")
	write.csv(A.transp.null,file_transp,row.names=F,col.names=F)
	write.csv(A.power.null,file_power,row.names=F,col.names=F)
	write.csv(A.inter.null,file_inter,row.names=F,col.names=F)
}
