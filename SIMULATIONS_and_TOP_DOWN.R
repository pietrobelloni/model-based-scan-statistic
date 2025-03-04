# Libraries can be installed from CRAN, with the exception of INLA (please refer to https://www.r-inla.org/download-install)
library(INLA)
library(sp)
library(SpatialEpi)
library(ggplot2)
library(patchwork)
library(DClusterm)

# load the map of Padua province with the mortality information about deaths for major causes in males population
# Since the data cannot be shared because of privacy, we simulated the denominator for the risk using the same distribution of the 
# true expected cases of deaths used in our study, adjusted for the population of the municipalities and multiplied by 1000. 
# The resulting variable follows a Gaussian distribution N(24.55, 18.76)
# The map need to be a SpatialPolygonsDataFrame object with the information in the @data part.
# In particular, each row is a specific area of the map, and, 
# for each area, there are the variable "Denominator" (in our study, the expected cases of death) and "Pop" for the population.

setwd("C:/Users/enric/OneDrive/Desktop/Lavoro/ArticoloSubmission/0CODICE")
load("PaduaMap_and_AdjMatrix.Rdata")
map <- PaduaMap

##### SIMULATE_MAPS ################################################################################################################################
# functions:
# - simulate_maps: simulates M maps with a cluster and some areas with excesses using a given map and the associated neighborhood matrix.
# - simulatedMap_Plot: shows the location of the true simulated cluster and the areas with excesses across the map.


simulate_maps <- function(map, adjMatr, M, areaClust, nAreas, beta = 0, rho = 0.999, sigma2u, sigma2v, delta, gamma, showPlot=0){
  # map is the map used to simulate maps and scenarios
  # adjMatr is the adjacency matrix or neighborhood matrix
  # M is the number of desired simulated maps
  # areaClust is the minimum number of areas to form the true cluster 
  # nAreas is the number of "random" areas with mortality excess
  # delta represents the excess value for the areas in the clusters
  # gamma represents the excess value for the random "nAreas" areas 
  # showPlot = 0 does not show plot
  # showPlot = 1 shows the plot of simulated Theta (theoretical excesses) for areas inside and outside the clusters for the M maps
  # showPlot = 2 shows the plot of the distribution of simulated Theta (theoretical excesses) for the M maps
  # showPlot = 3 shows both plots
  
  n <- nrow(map) #number of areas
  g <- inla.read.graph(adjMatr) #read the adjacency matrix
  vi <- rnorm(n*M, mean=0, sd=sqrt(sigma2v)) #generate vector v
  ui <- LaplacesDemon:::rmvnp(M, rep(0,n), (diag(rowSums(adjMatr))-rho*adjMatr)/sigma2u) #generate vector u
  zi <- rep(0, n*M) #initialise the vector for areas in the cluster  
  ri <- rep(0, n*M) #initialise the vector for random areas with mortality excess
  numCl <- round(runif(M, 1, n))

  for(i in 1:M){
    zi[(1+n*(i-1)):(n*i)][numCl[i]] <- 1
    zi[(1+n*(i-1)):(n*i)][which(map$Area %in% names(adjMatr[numCl[i],][which(adjMatr[numCl[i],]==1)]))] <- 1
    ri[(1+n*(i-1)):(n*i)][round(runif(nAreas, 0, n))] <- 1
    while(sum(zi[(1+n*(i-1)):(n*i)])<areaClust){
      zi[(1+n*(i-1)):(n*i)]<-rep(0,n)
      j = round(runif(1,1,n))
      zi[(1+n*(i-1)):(n*i)][j] <- 1
      zi[(1+n*(i-1)):(n*i)][which(map$Area %in% names(adjMatr[j,][which(adjMatr[j,]==1)]))] <- 1
      if(sum(zi[(1+n*(i-1)):(n*i)])>=5) break
    }
  }
  eta <- rep(NA, n*M); Teta <- rep(NA, n*M); SimNum <- rep(NA, n*M) #initialise eta, Theta and the simulated numerator (SimNum)
  for(i in 1:M){
    eta[(1+n*(i-1)):(n*i)] <- beta + ui[i,] + vi[(1+n*(i-1)):(n*i)] + delta*zi[(1+n*(i-1)):(n*i)] + gamma*ri[(1+n*(i-1)):(n*i)] 
    eta[(1+n*(i-1)):(n*i)] <- eta[(1+n*(i-1)):(n*i)] - mean(eta[(1+n*(i-1)):(n*i)])  
    Teta[(1+n*(i-1)):(n*i)] <- exp(eta[(1+n*(i-1)):(n*i)])
    SimNum[(1+n*(i-1)):(n*i)] <- rpois(n, map$Denominator*Teta[(1+n*(i-1)):(n*i)])
  }
  for(i in 1:M){
    assign(paste0("map", i), map) 
    mapName <- paste0("map", i) 
    maps <- get(mapName) 
    maps$NumSim <- SimNum[(1+n*(i-1)):(n*i)]; maps$TetaSim <- Teta[(1+n*(i-1)):(n*i)]
    maps$z <- zi[(1+n*(i-1)):(n*i)]; maps$r <- ri[(1+n*(i-1)):(n*i)] 
    assign(mapName, maps)
  }
  
  # plots 
  dati <- matrix(NA, n*M, 3); dati <- as.data.frame(dati); colnames(dati) <- c("Ind","Teta","z")
  dati[,1] <- seq(1, n*M, 1); dati[,2] <- Teta; dati[,3] <- zi+1
  dati$z <- as.factor(dati$z); v=16
  p1 <- ggplot(aes(x=Ind, y=Teta, colour = z), data=dati)+geom_point()+
    scale_colour_manual(values=c("#B8B8B8","#505050"), labels=c("Not cluster","Cluster"))+
    theme_bw() + xlab("") + ylab(expression("Simulated  "*theta))+
    theme(axis.text = element_text(size=v), axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), axis.title = element_text(size=v), 
          legend.title = element_blank(), legend.text = element_text(size=14),
          strip.text = element_text(size=v), legend.key.size = unit(1, 'cm'))
  p2 <- ggplot(aes(as.factor(rep(1,length(Teta))), Teta), data=dati)+geom_boxplot()+
    theme_bw() + xlab("") + ylab(expression("Simulated  "*theta))+
    theme(axis.text = element_text(size=v), axis.text.x = element_blank(),
          axis.title = element_text(size=v), legend.title = element_text(size=v), 
          legend.text = element_text(size=v), strip.text = element_text(size=v),
          plot.title = element_text(size=v), legend.key.size = unit(1, 'cm'))
  
  if(showPlot == 1) print(p1)
  if(showPlot == 2) print(p2)
  if(showPlot == 3) print(p1+theme(legend.position = "bottom")+p2)
  if(showPlot == 0) {}
    
  # M maps in output 
  map_names <- paste0("map", 1:M)
  return(mget(map_names))
}  

# EXAMPLES:
# SIMULATION 1, SCENARIO 1
SM_simulation1 <- simulate_maps(map, adjMatr, M=50, areaClust=5, nAreas=2, beta=0, rho=0.999, 
                              sigma2u=0.0003, sigma2v=0.00001, delta=0.45, gamma=0.4, showPlot=0)
# SIMULATION 2, SCENARIO 9
# SM_simulation2 <- simulate_maps(map, adjMatr, M=50, areaClust=5, nAreas=2, beta=0, rho=0.999, 
#                                sigma2u=0.03, sigma2v=0.001, delta=0.35, gamma=0.3, showPlot=3)


simulatedMap_Plot <- function(map){
  par(mfrow=c(1,2), mar= margin(1, 0.2, 1, 0.2))
  plot(map, col=map$z, main="True cluster") 
  plot(map, col=map$r, main="Random areas with excess") 
}


# EXAMPLE:
head(SM_simulation1$map1@data)
simulatedMap_Plot(SM_simulation1$map1)



##### TOP_DOWN_COVARIATES ##########################################################################################################################
# functions:
# - top_down: clustering procedure with top-down and bottom-up approaches,
#   using the given map, the associated neighborhood matrix and, eventually, some covariates.
# - clusters_Plot: shows various plot about the clustering procedure. It plots the risk, the adjusted and smoothed risk with the BYM2 model,
#   the areas with mortality excess resulting from the BYM2 model, the true cluster and the detected clusters with top-down and bottom-up approaches.
# - covariatesEffect_Plot: shows the effect of the covariates on the risk using the posteriori distribution of the beta parameters

# generate two covariates 
x1 <- rnorm(nrow(map))
x2 <- runif(nrow(map))
covariates <- cbind(x1,x2)

top_down <- function(map, adjMatr, covariates=NULL, fracPop=1, minAreas=3){
  # map is the map where to study the presence of clusters with mortality excesses
  # adjMatr is the adjacency matrix or neighborhood matrix
  # covariates represents the covariates to put into BYM2 and in DClusterm (matrix or data.frame with named columns)
  # fracPop is the maximum fraction of population to include in the cluster
  # minAreas is the number of minimum contiguous areas to define a cluster 
  
  # BYM2
  g <- inla.read.graph(adjMatr)
  map$re_u <- 1:nrow(map@data)
  prior.try <- list(prec = list(prior = "pc.prec", param = c(0.5 / 0.31, 0.01)),
                    phi = list(prior = "pc", param = c(0.5, 2 / 3)))
  if(is.null(covariates)==T){
    print("The top-down approach is using the bym2 model without covariates")
    formula0 <- NumSim ~ f(re_u, model = "bym2", graph = g, hyper = prior.try)
    res0 <- inla(formula0, family = "poisson", data = map@data, E = Denominator, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE))
    res <- res0
    dicm <- summary(res0)$dic$dic
  } else if(is.null(covariates)==F){
    covar <- paste(colnames(covariates), collapse = " + ")
    formula0 <- NumSim ~ f(re_u, model = "bym2", graph = g, hyper = prior.try)
    formula1 <- as.formula(paste0("NumSim ~ ", covar," + f(re_u, model = \"bym2\", graph = g, hyper = prior.try)"))
    res0 <- inla(formula0, family = "poisson", data = map@data, E = Denominator, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE))
    res1 <- inla(formula1, family = "poisson", data = map@data, E = Denominator, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE))
    dicm <- min(summary(res0)$dic$dic, summary(res1)$dic$dic)
    if(summary(res0)$dic$dic==dicm) {
      res <- res0
      print("The better bym2 model is the null model") } else if(summary(res1)$dic$dic==dicm) {
        res <- res1
        print(paste0("The better bym2 model is the model with ", covar, " covariates"))}
  }  
  bym2_results <- res
  map$RR <- res$summary.fitted.values[, "mean"]
  map$LR <- res$summary.fitted.values[, "0.025quant"]; map$UR <- res$summary.fitted.values[, "0.975quant"]
  map$excess <- round(map$RR,3)
  for(i in 1:nrow(map@data)){
    if(map$LR[i]<=1 & map$UR[i]>=1) map$excess[i]<-0
    if(map$LR[i]<=1 & map$UR[i]<=1) map$excess[i]<-0
  }
  if(summary(res0)$dic$dic==dicm) {
    model0 <- glm(NumSim ~ offset(log(Denominator)) + 1, family = "poisson", data = map)
  } else if(summary(res1)$dic$dic==dicm) {
    model0 <- glm(as.formula(paste0("NumSim ~ offset(log(Denominator)) + 1 +", covar)), family = "poisson", data = map)
  }
  
  sensspec <- matrix(NA, 12, 3)
  sensspec[,2]<-rep(c("TP","FP","TN","FN","Sens","Spec"),2)
  sensspec[,3]<-c(rep("Top-down",6),rep("Bottom-up",6))
  sensspec <- as.data.frame(sensspec)
  colnames(sensspec)<-c("Value","Criteria","Method")
  
  # TOP-DOWN APPROACH 
  centr <- map@data[which(map$excess!=0),]$Area; centroids <- coordinates(map)[centr,]
  gc()
  start <- Sys.time()
  DCl_td <- DetectClustersModel(map, thegrid=centroids, fractpop=fracPop, alpha = 0.05, typeCluster = "S", R = NULL, model0 = model0, ClusterSizeContribution = "Pop")
  time <- difftime(Sys.time(), start, units=c("secs"))
  time_td <- as.numeric(time[1])
  if(is.data.frame(DCl_td)==T){
    clusterFin_td <- slimknclusters(map, DCl_td, minsize=minAreas)
    clusterFin_td <- clusterFin_td[clusterFin_td$pvalue < 0.05/seq(nrow(map@data),nrow(map@data)-nrow(clusterFin_td)+1),]
    elements_td <- get.stclusters(map, clusterFin_td)
    map$cluster <- 0
    for(i in 1:nrow(clusterFin_td)){ map$cluster[elements_td[[i]]] <- paste("cluster",i,sep="") }
    map$num <- rep(0, nrow(map@data))
    for(i in 1:nrow(map@data)){
      for(j in 1:nrow(clusterFin_td)){ if(map$cluster[i]==paste0("cluster",j)) map$num[i] = j }
    }
    ic_td <- matrix(NA, nrow(clusterFin_td), 4); ic_td <- as.data.frame(ic_td)
    colnames(ic_td) <- c("Risk","Lower_0.025","Upper_0.975","Method")
    ic_td[,4] <- rep("Top-down", nrow(clusterFin_td))
    m2 <- glm(NumSim ~ offset(log(Denominator)) + 1, family = "poisson", data = map)
    c2 <- rep(NA, nrow(map@data))
    for(i in 1:nrow(ic_td)){
      c2 <- ifelse(map$cluster==paste0("cluster",i),1,0)
      m3 <- glm(NumSim ~ -1+c2, offset=log(fitted(m2)), family = "poisson", data = map)
      ic_td[i,1:3] <- round(c(exp(summary(m3)$coefficients[1]),
                              exp(summary(m3)$coefficients[1]+c(-1,1)*qnorm(0.975)*summary(m3)$coefficients[2])),3)
    }
    rand_td <- fossil::rand.index(map$z, map$num)
    num2 <- rep(NA, nrow(map@data)); num2 <- ifelse(map$num==0,0,1)
    tab <- table(num2, map$z)
    if(nrow(tab)==1){ tab <- rbind(tab,c(0,0)) }
    sensspec[1,1] <- tab[2,2]; sensspec[2,1] <- tab[2,1]; sensspec[3,1] <- tab[1,1]; sensspec[4,1] <- tab[1,2];
    sensspec[5,1] <- round(tab[2,2]/sum(tab[,2]),3); sensspec[6,1] <- round(tab[1,1]/sum(tab[,1]),3)
  } else{
    map$cluster <- rep(0, nrow(map@data))
    map$num <- rep(0, nrow(map@data))
    ic_td <- matrix(NA, 1, 3); ic_td <- as.data.frame(ic_td)
    colnames(ic_td) <- c("Risk","Lower_0.025","Upper_0.975","Method")
    ic_td[,4] <- rep("Top-down", 1)
    rand_td <- fossil::rand.index(map$z, map$num)
    tab <- table(map$num, map$z)
    if(nrow(tab)==1){ tab <- rbind(tab,c(0,0)) }
    sensspec[1,1] <- tab[2,2]; sensspec[2,1] <- tab[2,1]; sensspec[3,1] <- tab[1,1]; sensspec[4,1] <- tab[1,2];
    sensspec[5,1] <- round(tab[2,2]/sum(tab[,2]),3); sensspec[6,1] <- round(tab[1,1]/sum(tab[,1]),3)
  }
  
  # BOTTOM-UP APPORACH
  centroids_bu <- matrix(NA, nrow(map@data), 2)
  for(i in 1:nrow(map@data)) centroids_bu[i,] <- map@polygons[[i]]@labpt
  gc()
  start <- Sys.time()
  DCl_bu <- DetectClustersModel(map, thegrid=centroids_bu, fractpop=fracPop, alpha = 0.05, typeCluster = "S", R = NULL, model0 = model0, ClusterSizeContribution = "Pop")
  time <- difftime(Sys.time(), start, units=c("secs"))
  time_bu <- as.numeric(time[1])
  if(is.data.frame(DCl_bu)==T){
    clusterFin_bu <- slimknclusters(map, DCl_bu, minsize=minAreas)
    clusterFin_bu <- clusterFin_bu[clusterFin_bu$pvalue < 0.05/seq(nrow(map@data),nrow(map@data)-nrow(clusterFin_bu)+1),]
    elements_bu <- get.stclusters(map, clusterFin_bu)
    map$cluster_bu <- 0
    for(i in 1:nrow(clusterFin_bu)){ map$cluster_bu[elements_bu[[i]]] <- paste("cluster",i,sep="") }
    map$num_bu <- rep(0, nrow(map@data))
    for(i in 1:nrow(map@data)){
      for(j in 1:nrow(clusterFin_bu)){ if(map$cluster_bu[i]==paste0("cluster",j)) map$num_bu[i] = j }
    }
    ic_bu <- matrix(NA, nrow(clusterFin_bu), 4); ic_bu <- as.data.frame(ic_bu)
    colnames(ic_bu) <- c("Risk","Lower_0.025","Upper_0.975","Method")
    ic_bu[,4] <- rep("Bottom-up", nrow(clusterFin_bu))
    m2 <- glm(NumSim ~ offset(log(Denominator)) + 1, family = "poisson", data = map)
    c2 <- rep(NA, nrow(map@data))
    for(i in 1:nrow(ic_bu)){
      c2 <- ifelse(map$cluster_bu==paste0("cluster",i),1,0)
      m3 <- glm(NumSim ~ -1+c2, offset=log(fitted(m2)), family = "poisson", data = map)
      ic_bu[i,1:3] <- round(c(exp(summary(m3)$coefficients[1]),
                              exp(summary(m3)$coefficients[1]+c(-1,1)*qnorm(0.975)*summary(m3)$coefficients[2])),3)
    }
    rand_bu <- fossil::rand.index(map$z, map$num_bu)
    num2_bu <- rep(NA, nrow(map@data)); num2_bu <- ifelse(map$num==0,0,1)
    tab <- table(num2_bu, map$z)
    if(nrow(tab)==1){ tab <- rbind(tab,c(0,0)) }
    sensspec[7,1] <- tab[2,2]; sensspec[8,1] <- tab[2,1]; sensspec[9,1] <- tab[1,1]; sensspec[10,1] <- tab[1,2]; 
    sensspec[11,1] <- round(tab[2,2]/sum(tab[,2]),3); sensspec[12,1] <- round(tab[1,1]/sum(tab[,1]),3)
  } else{
    map$cluster_bu <- rep(0, nrow(map@data))
    map$num_bu <- rep(0, nrow(map@data))
    ic_bu <- matrix(NA, 1, 4); ic_bu <- as.data.frame(ic_bu)
    colnames(ic_bu) <- c("Risk","Lower_0.025","Upper_0.975","Method")
    ic_bu[,4] <- rep("Bottom-up", 1)
    rand_bu <- fossil::rand.index(map$z, map$num_bu)
    tab <- table(map$num_bu, map$z)
    if(nrow(tab)==1){ tab <- rbind(tab,c(0,0)) }
    sensspec[7,1] <- tab[2,2]; sensspec[8,1] <- tab[2,1]; sensspec[9,1] <- tab[1,1]; sensspec[10,1] <- tab[1,2]; 
    sensspec[11,1] <- round(tab[2,2]/sum(tab[,2]),3); sensspec[12,1] <- round(tab[1,1]/sum(tab[,1]),3)
  }

  clusterTime <- as.data.frame(matrix(c(round(time_td,3),round(time_bu,3),"Top-down","Bottom-up"),2,2)); colnames(clusterTime) <- c("Time", "Method")
  clusterRisk <- rbind(ic_td, ic_bu)
  clusterRand <- as.data.frame(matrix(c(round(rand_td,3),round(rand_bu,3),"Top-down","Bottom-up"),2,2)); colnames(clusterRand) <- c("Rand Index", "Method")
  
  # the output is a list with the time of the procedures, the risks of the clusters and the performance ability through 
  # the Rand index, sensitivity and specificity 
  return(list(
    map = map,
    bym2_results = bym2_results,
    clusterTime = clusterTime,
    clusterRisk = clusterRisk,
    clusterRand = clusterRand,
    sensspec = sensspec
  ))
  
}

# EXAMPLE:
# our simulations
TD <- top_down(SM_simulation1$map1, adjMatr, covariates=covariates, fracPop=1, minAreas=3)


clusters_Plot <- function(map, p=1){
  # map is the map where to study the presence of clusters with mortality excesses
  # p is a parameter to choose the graph:
  # p = 0 shows the risk of the map
  # p = 1 shows the adjusted and smoothed risk with the BYM2 model
  # p = 2 shows the areas with mortality excess resulting from the BYM2 model
  # p = 3 shows the true cluster and the detected clusters with top-down and bottom-up approaches
  
  if(p==0) mapvariable(map$NumSim/map$Denominator, map, main="Risk")
  if(p==1) mapvariable(map$RR, map, main="Adjusted and smoothed risk", lower=min(map$NumSim/map$Denominator), upper=max(map$NumSim/map$Denominator)) 
  if(p==2) mapvariable(map$excess, map, main="Areas with mortality excess")
  if(p==3) {par(mfrow=c(1,3), mar=c(5, 0, 5, 0))
    plot(map, col=map$z, main="True cluster")
    plot(map, col=map$num, main="Detected cluster with top-down")
    plot(map, col=map$num_bu, main="Detected cluster with bottom-up")}
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
}

# EXAMPLE:
clusters_Plot(TD$map, p=3)


plot_dens <- list()
covariatesEffect_Plot <- function(bym2_results){
  if(nrow(bym2_results$summary.fixed)==1) {print("The covariates have no effect on the risk")}
  if(nrow(bym2_results$summary.fixed)>1){
    for(i in 1:(nrow(bym2_results$summary.fixed)-1)){
      marginal <- inla.smarginal(res$marginals.fixed[[i+1]]); marginal <- data.frame(marginal)
      plot_dens[[i]]<-ggplot(marginal, aes(x = x, y = y)) + geom_line() +
        labs(x = expression(beta), y = "Density") +
        geom_vline(xintercept = 0, col = "blue") + 
        geom_vline(xintercept = res$summary.fixed$`0.025quant`[i+1], col = "red", "dashed") + 
        geom_vline(xintercept = res$summary.fixed$`0.975quant`[i+1], col = "red") + 
        theme_bw()
    }} 
  return(plot_dens)
}

# ESAMPLE:
listPlot <- covariatesEffect_Plot(TD$bym2_results)
listPlot[1]; listPlot[2] 

