DijkstraAlgor <- function(origin_index, network_matrix) {
	number_of_nodes <- dim(network_matrix)[1]
	shortest_distance <- rep(Inf,number_of_nodes)
	shortest_distance[origin_index] <- 0
	solved_nodes_with_parents <- list()
	unsolved_nodes <- 1:number_of_nodes
	unsolved_nodes <- unsolved_nodes[ unsolved_nodes != origin_index ]
	while(length(unsolved_nodes) > 0) {
		shortest_distance_matrix <- matrix(rep(shortest_distance,number_of_nodes),number_of_nodes,number_of_nodes)
		distance_uptator_matrix <- shortest_distance_matrix + network_matrix
		distance_updator <- apply(distance_uptator_matrix,2,min)
		new_solved_nodes_distance <- min(distance_updator[unsolved_nodes])
		new_solved_nodes <- intersect(which(distance_updator == new_solved_nodes_distance),unsolved_nodes)
		shortest_dist_to_unsolved_nodes_m <- matrix(rep(distance_updator,number_of_nodes),nrow=number_of_nodes,byrow=T)
		parent_nodes <- which(distance_uptator_matrix == shortest_dist_to_unsolved_nodes_m, arr.ind=T)
		shortest_distance[new_solved_nodes] <- new_solved_nodes_distance
		for(node in new_solved_nodes) {
			if(is.infinite(new_solved_nodes_distance)) {          #some nodes do not have a path from the origin
				solved_nodes_with_parents[[node]] <- "unreachable"
			} else {
				solved_nodes_with_parents[[node]] <- parent_nodes[,"row"][parent_nodes[,"col"] == node]
			}
			unsolved_nodes <- unsolved_nodes[unsolved_nodes != node]
		}
	}
	solved_nodes_with_parents[[origin_index]] <- "origin"
	return(solved_nodes_with_parents)
}

linkPathMatrix <- function(network_matrix) {
	number_of_nodes <- dim(network_matrix)[1]
	link_matrix <- which(!is.infinite(network_matrix),arr.ind=T)
	link_set <- apply(link_matrix,1,paste,collapse=",")
	link_path_matrix <- matrix(0,nrow=number_of_nodes^2,ncol=length(link_set))
	colnames(link_path_matrix) <- link_set
	path_set <- c()
	for(i in 1:number_of_nodes) {
		for(j in 1:number_of_nodes) {
			path_set <- c(path_set, paste(c(j,i), collapse=","))
		}
	}
	rownames(link_path_matrix) <- path_set
	for(node in 1:number_of_nodes) {
		one_origin_solution <- DijkstraAlgor(node,network_matrix)
		for(other_node in 1:number_of_nodes) {		
			path <- paste(c(node,other_node),collapse=",")
			link_end <- sample(one_origin_solution[[other_node]],1)         #if there is a tie between two nodes, randomly select one
			if(length(one_origin_solution[[other_node]]) == 1) link_end <- one_origin_solution[[other_node]]
			other_node_temp <- other_node
			while(!is.character(link_end)) {       #back-tracking
				link <- paste(c(link_end,other_node_temp),collapse=",")				
				link_path_matrix[path,link] <- 1
				other_node_temp <- link_end
				link_end <- sample(one_origin_solution[[other_node_temp]],1)
				if(length(one_origin_solution[[other_node_temp]]) == 1) link_end <- one_origin_solution[[other_node_temp]]
			}
		}
	}
	return(link_path_matrix)
}

##get the flow on each link based on link-path matrix and OD
linkFlowVector <- function(LPM,demand_matrix) {
	link_flow <- c()     
	for(link in colnames(LPM)) {
		link_flow[link] <- sum(LPM[,link]*demand_matrix)
	}
	return(link_flow)
}

pathToLinkFlow <- function(path_s,path_f) {
  link_flow <- 0
  for (OD in names(path_s)) {
    OD_paths_m <- path_s[[OD]]
    OD_paths_flow <- path_f[[OD]]
    link_flow_m <- OD_paths_m*OD_paths_flow
    OD_link_flow <- colSums(link_flow_m)
    link_flow <- link_flow + OD_link_flow
  }
  return(link_flow)
}

##Converting a link-cost vector to an adjacency matrix format
vectorToMatrix <- function(vec,node_number) {
  vec_matrix <- matrix(Inf,node_number,node_number)
  for(ele in names(vec)) {
    indexes <- as.integer(strsplit(ele,",")[[1]])
    vec_matrix[indexes[1],indexes[2]] <- vec[ele]
  }
  return(vec_matrix)
}

secondDerivative <- function(Tfree,X,C) {
  return(0.6*Tfree*X^3/C^4)
}

PathBasedAlgorithm <- function(initial_cost,demand_matrix,number_of_nodes,BPR_information) {
  old_criter <- 1
  initial_cost_matrix <- vectorToMatrix(initial_cost,number_of_nodes)
  link_path <- linkPathMatrix(initial_cost_matrix)
  initial_cost_reorder <- initial_cost[match(colnames(link_path),BPR_information$Link_from_to)]
  path_set <- list()
  path_cost <- list()
  path_flow <- list()
  for (OD in rownames(link_path)) {
    OD_path <- matrix(link_path[OD,],nrow=1)
    colnames(OD_path) <- names(link_path[OD,])
    path_set[[OD]] <- OD_path
    path_cost[[OD]] <- sum(initial_cost_reorder*path_set[[OD]][1,])
    indexes <- as.integer(strsplit(OD,",")[[1]])
    path_flow[[OD]] <- demand_matrix[indexes[1],indexes[2]]
  }
  
  flows <- linkFlowVector(link_path,demand_matrix)
  if (!identical(colnames(link_path),names(flows))) {
    stop("the order does not match")
  }
  Ta <- BPR_information$Free_Flow_Travel_Time
  Ca <- BPR_information$Capacity
  Ta <- Ta[match(names(flows),BPR_information$Link_from_to)]
  Ca <- Ca[match(names(flows),BPR_information$Link_from_to)]
  repeat {
    cost <- Ta*(1+0.15*(flows^4)/(Ca^4))
    for (OD in rownames(link_path)) {
      for (i in 1:length(path_cost[[OD]])) {
        path_cost[[OD]][i] <- sum(cost*path_set[[OD]][i,])
      }
    }
    
    cost_matrix <- vectorToMatrix(cost,number_of_nodes)
    link_path_update <- linkPathMatrix(cost_matrix)
    for (OD in rownames(link_path_update)) {
      OD_path_update <- link_path_update[OD,]
      OD_paths_old <- path_set[[OD]]
      OD_path_cost_update <- sum(cost*OD_path_update)
      OD_path_cost_old <- path_cost[[OD]]
      flag <- F
      for (i in 1:dim(OD_paths_old)[1]) {
        if (OD_path_cost_old[i] == OD_path_cost_update) flag <- T     #comparing path cost (instead of the exact path) is better, bacause there may be ties
      }
      if (!flag) {
        path_set[[OD]] <- rbind(OD_paths_old,OD_path_update)
        path_cost[[OD]] <- c(path_cost[[OD]],OD_path_cost_update)   #check whether the added path is the shortest
        path_flow[[OD]] <- c(path_flow[[OD]],0)
      }
    }
    
    S_all <- secondDerivative(Ta,flows,Ca)   ##derivatives for all links
    for (OD in rownames(link_path_update)) {
      path_costs <- path_cost[[OD]]
      path_sets <- path_set[[OD]]
      path_flows <- path_flow[[OD]]
      shortest_path_index <- which(path_costs == min(path_costs))
      shortest_path_cost <- min(path_costs)
      shortest_path <- path_sets[shortest_path_index,]
      non_shortest_paths_index <- c()
      if (length(path_costs)>1) {
        non_shortest_paths_index <- 1:length(path_costs)
        non_shortest_paths_index <- non_shortest_paths_index[-shortest_path_index]
      }
      for (j in non_shortest_paths_index) {
        path_k <- path_sets[j,]
        S_k <- sum(S_all[path_k != shortest_path])
        path_cost_k <- path_costs[j]
        path_flow_k <- path_flows[j]
        path_flow[[OD]][j] <- max(0,path_flow_k-0.01*(path_cost_k - shortest_path_cost)/S_k)
      }
      indexes <- as.integer(strsplit(OD,",")[[1]])
      OD_amount <- demand_matrix[indexes[1],indexes[2]]
      path_flow[[OD]][shortest_path_index] <- OD_amount - sum(path_flow[[OD]][-shortest_path_index])
    }
    
    flows_update <- pathToLinkFlow(path_set,path_flow)
    flows_update <- flows_update[match(names(flows),names(flows_update))]
    #new_criter <- sqrt(sum((flows-flows_update)^2))
    new_criter <- sum(abs(flows-flows_update))
    print(new_criter)
    if(new_criter < 1) {
      break
    }
    flows <- flows_update
    #old_criter <- new_criter
  }
  return(list(path_set,path_flow))
}

linkSetToPath <- function(OD,pathSet) {
  ordered_links <- c()
  origin <- strsplit(OD,",")[[1]][1]
  destination <- strsplit(OD,",")[[1]][2]
  pathSet_temp <- pathSet
  for (a_link in pathSet_temp) {
    if(grepl(origin, a_link)) {
      ordered_links <- c(ordered_links,a_link)
      pathSet_temp <- pathSet_temp[pathSet_temp != a_link]
      break
    }
  }
  
  last_link <- ordered_links[length(ordered_links)]
  current_end <- strsplit(last_link,",")[[1]][2]
  while(current_end != destination) {
    for (a_link in pathSet_temp) {
      if(grepl(current_end, a_link)) {
        ordered_links <- c(ordered_links,a_link)
        pathSet_temp <- pathSet_temp[pathSet_temp != a_link]
        break
      }
    }
    last_link <- ordered_links[length(ordered_links)]
    current_end <- strsplit(last_link,",")[[1]][2]
  }
  return(ordered_links)
}


BPR_variables <- read.csv("/benchmark model/SiouxFallsDataBPR.csv",header=T,colClasses = c("character",rep(NA,5)))
demand <- read.csv("/benchmark model/SiouxFallsDataDemand.csv",header=T)
demand <- as.matrix(unname(demand[,-1]))
BPR_variables$Link_from_to <- unname(sapply(BPR_variables$Link_from_to,sub,pattern="\\.",replacement=','))
##for each rerun, the following 3 lines need to be run
nodes_num <- nrow(demand)
free_flow_cost <- BPR_variables$Free_Flow_Travel_Time
names(free_flow_cost) <- BPR_variables$Link_from_to
outcome <- PathBasedAlgorithm(free_flow_cost,demand,nodes_num,BPR_variables)


free_cost_matrix <- vectorToMatrix(free_flow_cost,nodes_num)
link_path_matrix <- linkPathMatrix(free_cost_matrix)

##incorporate random failures later
##initialization
paths_set <- outcome[[1]]
paths_flow <- outcome[[2]]
link_flows <- pathToLinkFlow(paths_set,paths_flow)
link_set <- names(link_flows)
capacity <- BPR_variables$Capacity
capacity <- capacity[match(names(link_flows),BPR_variables$Link_from_to)]
free_flow_cost <- free_flow_cost[match(names(link_flows),BPR_variables$Link_from_to)]
initial_failure <- which(link_flows > 2.403*capacity)
failed_links <- names(initial_failure)
survival_links <- setdiff(link_set,failed_links)
random_failed_links <- survival_links[which(runif(length(survival_links)) < 0.1)]
failed_links <- c(failed_links,random_failed_links)
timeS <- 1
failure_time <- rep(1,length(failed_links))
propagation <- list()
paths_link_flow <- list()
for (OD in names(paths_set)) {
  paths_link_flow[[OD]] <- paths_flow[[OD]]*paths_set[[OD]]
  propagation[[OD]] <- c()
  #  index_zero_flow <- which(paths_flow[[OD]] == 0)
  #  paths_link_flow[[OD]] <- paths_link_flow[[OD]][-index_zero_flow,]
}


##repeat
link_cost <- free_flow_cost*(1+0.15*(link_flows^4)/(capacity^4))
link_cost[failed_links] <- NA
link_cost <- link_cost[!is.na(link_cost)]
cost_matrix <- vectorToMatrix(link_cost,nodes_num)
new_shortest_paths <- linkPathMatrix(cost_matrix)

new_paths <- list()
##as long as the current path goes through at least one failed link, a new path needs to be find.
for (OD in names(paths_set)) {
  all_paths_for_OD <- paths_link_flow[[OD]]
  new_paths[[OD]] <- matrix(NA,length(paths_flow[[OD]]),dim(new_shortest_paths)[2])
  colnames(new_paths[[OD]]) <- colnames(new_shortest_paths)
  for (i in 1:dim(all_paths_for_OD)[1]) {
    a_path <- all_paths_for_OD[i,]
    a_path_flow <- max(a_path)
    if (sum(a_path[failed_links],na.rm=T) != 0) {
      ##get the path
      traversed_links <- names(a_path)[which(a_path > 0)]
      path_ordered_links <- linkSetToPath(OD,traversed_links)
      
      ##Assuming that rerouting happens after encountering the failed links.
      # first_failed_link_indxes <- min(which(path_ordered_links %in% failed_links))
      # new_origin <- strsplit(path_ordered_links[first_failed_link_indxes],",")[[1]][1]
      # new_destin <- strsplit(OD,",")[[1]][2]
      # new_OD <- paste(new_origin,new_destin,sep=",")
      # new_path <- new_shortest_paths[new_OD,]
      failed_link_indxes_on_path <- which(path_ordered_links %in% failed_links)
      new_path <- new_shortest_paths[OD,]
      
      # j_temp <- 1
      # while(j_temp < first_failed_link_indxes) {
      #   if(!(path_ordered_links[j_temp] %in% names(new_path))) stop("nonexistent link")
      #   new_path[path_ordered_links[j_temp]] <- 1
      #   j_temp <- j_temp + 1
      # }
      new_path <- new_path*a_path_flow
      new_paths[[OD]][i,] <- new_path
      retouted_to <- setdiff(names(new_path)[new_path > 0], names(a_path)[a_path > 0])
      for (index in failed_link_indxes_on_path) {
        propagation[[path_ordered_links[index]]] <- c(propagation[[path_ordered_links[index]]], retouted_to)
      }
    } else {
      new_paths[[OD]][i,] <- 0
      links_modify <- names(a_path)[a_path > 0]
      new_paths[[OD]][i,][links_modify] <- a_path_flow
    }
  }
}

updated_link_flows <- 0
for (OD in names(new_paths)) {
  updated_link_flows <- updated_link_flows +colSums(new_paths[[OD]])  
}
update_capacity <- BPR_variables$Capacity[match(names(updated_link_flows),BPR_variables$Link_from_to)]
update_failure <- names(updated_link_flows)[which(updated_link_flows > 2.403*update_capacity)]
failed_links <- c(failed_links,update_failure)
survival_links <- setdiff(link_set,failed_links)
random_failed_links <- survival_links[which(runif(length(survival_links)) < 0.1)]
failed_links <- c(failed_links,random_failed_links)
timeS <- timeS+1
failure_time <- c(failure_time, rep(timeS,length(update_failure)+length(random_failed_links)))

capacity <- update_capacity
free_flow_cost <- BPR_variables$Free_Flow_Travel_Time[match(names(updated_link_flows),BPR_variables$Link_from_to)]
link_flows <- updated_link_flows
paths_link_flow <- new_paths

# OD_link_flow <- matrix(NA,length(names(paths_set)),length(failed_links))
# rownames(OD_link_flow) <- names(paths_set)
# colnames(OD_link_flow) <- failed_links
# for (OD in rownames(OD_link_flow)) {
#   if(length(paths_flow[[OD]]) > 1) {
#     OD_link_flow[OD,] <- colSums(paths_flow[[OD]]*paths_set[[OD]][,failed_links])
#   } else {
#     OD_link_flow[OD,] <- paths_flow[[OD]]*paths_set[[OD]][,failed_links]
#   }
# }
# for (link in colnames(OD_link_flow)) {
#   redistributed_flows <- OD_link_flow[,link]
#   destination_flows <- c()
#   for (n in names(redistributed_flows)) {
#     if (redistributed_flows[n] > 0) {
#       destin <- strsplit(n,",")[[1]][2]
#       if (destin %in% names(destination_flows)) {
#         destination_flows[destin] <- destination_flows[destin] + redistributed_flows[n]
#       } else {
#         destination_flows[destin] <- redistributed_flows[n]
#       }
#     }
#   }
# }

failed_links
failure_time
propagation
for (aLink in names(propagation)) {
  propagation[[aLink]] <- unique(propagation[[aLink]])
}
write.csv(failure_time,"/benchmark model/failure_time.csv")
write.csv(failed_links,"/benchmark model/failed_links.csv")
propagation_matrix <- matrix(0,76,76)
colnames(propagation_matrix) <- BPR_variables$Link_from_to
rownames(propagation_matrix) <- BPR_variables$Link_from_to
for (aLink in names(propagation)) {
  for (anotherLink in propagation[[aLink]]) {
    propagation_matrix[aLink,anotherLink] <- 1
  }
}
write.csv(propagation_matrix,"/benchmark model/propagation_matrix.csv")



current_flows <- outcome[[4]]
Ta_global <- BPR_variables$Free_Flow_Travel_Time
Ca_global <- BPR_variables$Capacity
Ta_global <- Ta_global[match(names(current_flows),BPR_variables$Link_from_to)]
Ca_global <- Ca_global[match(names(current_flows),BPR_variables$Link_from_to)]
current_cost <- Ta_global*(1+0.15*(current_flows^4)/(Ca_global^4))
current_cost_matrix <- vectorToMatrix(current_cost,nodes_num)
current_LPM <- linkPathMatrix(current_cost_matrix)