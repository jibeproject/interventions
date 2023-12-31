#SET UP
rm(list = ls())

install.packages("sf")
install.packages("igraph")
install.packages("sfnetworks")
install.packages("tidygraph")
install.packages("stplanr")
install.packages("influenceR")

library(sf)
library(dplyr)
library(igraph)
library(sfnetworks)
library(tidygraph)
library(stplanr)
library(influenceR)

####################
#PART 1: read files and clean values
####################
setwd("")
#read City of Manchester boundary
manchester_bound <- st_read(file.path("./GM_Districts.gpkg")) #the city boundaries are found in Teams WP2>Data_WP2>Processed_Data>Greater Manchester>GM_Administrative units>GM_districts
#read network
network <- st_read(file.path("./network2way.gpkg")) #the MATSim generated network is found in Teams WP2>Data_WP2>Processed_Data>Greater Manchester>GM_Network
network <- st_intersection(network, manchester_bound[,1]) %>% st_cast("LINESTRING") #clip to GM
pct <- st_read(file.path("./network_v3.12pct.gpkg")) #read base network with PCT demand attached found in Teams WP2>Data_WP2>Processed_Data>Greater Manchester>GM_Network

#add cycling infra
network <- merge(network, st_drop_geometry(pct[,c("edgeID", "aadb_f", "aadb_b", "cyclesm", "govnearmkt_slc_majority", "quitnss")]), by = "edgeID", all.y=FALSE)

#OFFROAD = 4(best) / PROTECTED = 3 / LANE(painted) = 2 / MIXED = 1(worst)
network$bikeProtectionType <- ifelse(stringr::str_detect(network$bikeProtectionType, "OFFROAD") == TRUE, paste0(as.numeric(4)), network$bikeProtectionType)
network$bikeProtectionType <- ifelse(stringr::str_detect(network$bikeProtectionType, "PROTECTED") == TRUE, paste0(as.numeric(3)), network$bikeProtectionType)
network$bikeProtectionType <- ifelse(stringr::str_detect(network$bikeProtectionType, "LANE") == TRUE, paste0(as.numeric(2)), network$bikeProtectionType)
network$bikeProtectionType <- ifelse(stringr::str_detect(network$bikeProtectionType, "MIXED") == TRUE, paste0(as.numeric(1)), network$bikeProtectionType)

#give NA values lowest demand weights; replace NAs and 0s with 0.00001s (we later invert these values because the algorithm for the shortest path always chooses the least cost path among two vertices, and we want the least cost path to be the one with the highest cycling demand)
network$aadb_f <- ifelse(stringr::str_detect(network$aadb_f, "0") == TRUE, paste0(as.numeric(0.00001)), network$aadb_f)
network$aadb_b <- ifelse(stringr::str_detect(network$aadb_b, "0") == TRUE, paste0(as.numeric(0.00001)), network$aadb_b)
network$aadb_f <- as.numeric(network$aadb_f)
network$aadb_b <- as.numeric(network$aadb_b)
network$aadb_f[is.na(network$aadb_f)] <- 0.00001
network$aadb_b[is.na(network$aadb_b)] <- 0.00001

#give NA values lowest demand weights; replace NAs and 0s with 0.00001s (we later invert these values because the algorithm for the shortest path always chooses the least cost path among two vertices, and we want the least cost path to be the one with the highest cycling demand)
network$govnearmkt_slc_majority[is.na(network$govnearmkt_slc_majority)] <- 0.00001
#convert NaN values to infinite
network$bikeJibeMarginalDisutility[is.na(network$bikeJibeMarginalDisutility)] <- Inf

####################
#PART 2: create Tidy Geospatial Network (Graph) and add weights to edges and nodes
####################

#convert to sf
network <- as_sfnetwork(network)

#add various attributes as links' and nodes' weights
network <- network %>% activate("edges") %>% mutate(weight = st_length(network %>% activate("edges"))) %>%
                                             mutate(weight_1 = E(network)$govnearmkt_slc_majority) %>%
                                             mutate(weight_2 = (E(network)$aadb_f + E(network)$aadb_b)) %>%
                                             mutate(weight_3 = as.numeric(bikeProtectionType)) %>%
                                             mutate(weight_4 = E(network)$f_bikeStressJct)

network <- network %>% activate("nodes") %>% set_vertex_attr("nodeID", value = 1:length(V(network)))

#the PageRank score is a composite measure of incoming links with highets cycling demand as per the Propensity to Cycle Tool (PCT) and cycling stress crossing a junction
network <- network %>% activate("nodes") %>% mutate(pg_rnk = page_rank(network, algo = "prpack", V(network), directed = TRUE, weights = (network %>% activate("edges") %>% pull(weight_1)) + (network %>% activate("edges") %>% pull(weight_4)))$vector)

network <- network %>% activate("nodes") %>% mutate(burt = igraph::constraint((network %>% activate("nodes")), weights = network %>% activate("edges") %>% pull(weight_1)))

network <- network %>% activate("nodes") %>% mutate(bc = centrality_betweenness(weights = network %>% activate("edges") %>% pull(weight), normalized = TRUE, directed = TRUE))

network <- network %>% activate("nodes") %>% mutate(eigen = eigen_centrality(network, weights = (network %>% activate("edges") %>% pull(weight_1)), scale = TRUE, directed = TRUE)$vector)

network <- network %>% activate("nodes") %>% mutate(ens = influenceR::ens(network))

#NB. We have tested Page Rank score, Burt constraint and Burt's Effective Network Size as node weights and the best performing was the PageRank algorithm

####################
#PART 3: Split the network (graph) into sub-graphs
####################
#subgraph g2 containing streets with cycling infra or quiet streets (as defined by CyclingStreets with their "quietness index" >= 75 );
#subgraph g1 containing streets that don't have any or have cycling infra of poor quality (mixed use)

g2 <- subgraph.edges(graph=network, eids=which(E(network)$weight_3==4 | E(network)$weight_3==3 | E(network)$weight_3==2 | E(network)$quitnss>=75), delete.vertices = TRUE) #NB. form and to get changed with sub-setting the graph, they correspond with feature IDs of the nodes rather than nodeIDs
g1 <- subgraph.edges(graph=network, eids=which(E(network)$weight_3==1 & E(network)$quitnss<75), delete.vertices = TRUE)

#decompose subgraphs to find disconnected parts (for now the code makes use only of the 20 largest disconnected subgraphs in g2; future work can extend to include most of the subgraphs)
dc1 <- decompose(g1, mode = "weak") #Decompose a graph into components
c1 <- components(g1) #get the components
dc2 <- decompose(g2, mode = "weak") #Decompose a graph into components
c2 <- components(g2) #get the components

#get top 20 subgraphs nodeIDs of g2
cluster_order <- order(c2$csize, decreasing = TRUE)
top20_nodeIDs <- NULL
rm(i)
for (i in 1:20){
  top20_nodeIDs[[i]] <- V(g2)$nodeID[c2$membership == cluster_order[i]]
}

#get intersections (nodes) that bridge GOOD and BAD infra, aka. boundary nodes
bound.nodes <- vertex_attr(g1)$nodeID[(vertex_attr(g1)$nodeID %in% vertex_attr(g2)$nodeID)==TRUE]

#get the boundary nodes coordinates of the largest 20 sub-graphs with quality infra
boundlist <- NULL

for (j in 1:20){
  nam <- paste0("bound", j)
  assign(nam, bound.nodes[(bound.nodes %in% top20_nodeIDs[[j]])==TRUE]) #top 20 largest sub-graphs
  for(l in 1:length(eval(parse(text=nam)))){
    subgraph_names <- paste("bound_coor", 1:l, sep="")
    subgraph <- vector("list", length(subgraph_names))
    names(subgraph) <- subgraph_names
    subgraph <- vertex_attr(g2)$geometry[(vertex_attr(g2)$nodeID %in% eval(parse(text=nam)))==TRUE][l]
    names(subgraph) <- vertex_attr(g2)$nodeID[(vertex_attr(g2)$nodeID %in% eval(parse(text=nam)))==TRUE][l]
    boundlist <- c(boundlist, subgraph)
  }
  rm(nam)
}

rm(j, l, i, subgraph, subgraph_names)
rm(list = ls()[grep("bound", ls())])

##########################
#PART 4: Construct boundary nodes' pairs
#########################
#check PageRank cut-off values distribution
summary(V(network)$pg_rnk)
#greater cut-off should be selected to create less node pairs and reduced run time: e.g. top 0.1% of nodes with highest PageRanking score
pg_rnk <- min(tail(sort(V(network)$pg_rnk), length(V(network))*0.001))

rm(i,j, nodeID_clust1, nodeID_clust2, pg_rnk_clust1, pg_rnk_clust2, dist, df, euc)

df <- data.frame(nodeID_clust1=as.character(), pg_rnk_clust1=as.numeric(), nodeID_clust2=as.character(), pg_rnk_clust2=as.numeric())
euc <- NULL
nodeIDs <- network %>% activate("nodes") %>% as.data.frame() #make a data.frame to speed up search time

####################
for (i in 1:length(boundlist)){
    nodeID_clust1 <- names(boundlist)[i]
    pg_rnk_clust1 <- nodeIDs[which(nodeIDs$nodeID==nodeID_clust1),"pg_rnk"]
    if (pg_rnk_clust1 > pg_rnk){
      for (j in 1:length(boundlist)){
      nodeID_clust2 <- names(boundlist)[j]
      pg_rnk_clust2 <- nodeIDs[which(nodeIDs$nodeID==nodeID_clust2),"pg_rnk"]
      if ((nodeID_clust1 != nodeID_clust2) & pg_rnk_clust2 > pg_rnk){
            dist <- st_distance(boundlist[[i]], boundlist[[j]], which = "Euclidean")
            df <- data.frame(nodeID_clust1=nodeID_clust1, pg_rnk_clust1=pg_rnk_clust1, nodeID_clust2=nodeID_clust2, pg_rnk_clust2=pg_rnk_clust2, dist=dist)
            euc <- dplyr::bind_rows(df, euc)
      }
    }
    }else{
  }
}

euc <- euc[order(euc$dist),]
euc <- euc[!duplicated(euc$dist),]

###########################
#PART 5: Connect disconnected boundary nodes
#########################
rm(i,paths_euc, inter_euc, interventions)

paths_euc <- NULL
inter_euc <- NULL
interventions <- NULL

#parallelize in future work
for (i in 1:nrow(euc)){                                                   #add least cost path constraint constructed of highest cycling demand and highest AADT
  paths_euc[[i]] <- shortest_paths(network, euc[i,1], euc[i,3], weights =(1/E(network)$govnearmkt_slc_majority + 1/E(network)$aadt), output = "both")
  inter_euc[[i]] <- pct[pct$edgeID %in% E(network)$edgeID[which(E(network) %in% unlist(paths_euc[[i]]$epath))],]
  interventions <- dplyr::bind_rows(interventions, inter_euc[[i]])
}

rm(list=setdiff(ls(), c("interventions", "g1"))) #clear everything except interventions and subgraph g1
.rs.restartR()
rm(list=setdiff(ls()))

###########################
#PART 6: Add cycling interventions to the base network
#########################
#read netwrok v3.13
network <- st_read("D:/JIBE/02_DataOutput/network/gm/network_v3.13.gpkg") #file found in Teams WP2>Data_WP2>Processed_Data>Greater Manchester>GM_Network

interventions <- interventions[!duplicated(interventions$edgeID),]
interventions <- interventions[(interventions$edgeID %in% E(g1)$edgeID),] #keep only edgeIDs found in g1 (graph with BAD cycling infra)

#add PROTECTED cycling lanes to edges without cycling lanes
network$cyclesm <- ifelse((network$edgeID %in% interventions$edgeID)==TRUE, as.character("protected"), network_copy$cyclesm)

#write network with new interventions
sf::write_sf(network,"./network_interventions.gpkg")
