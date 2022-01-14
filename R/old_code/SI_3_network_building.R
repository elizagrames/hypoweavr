# Workspace set up -------------------------------------------------------------

palette <- c("#0776C5", "#65B4E2", "#FFBA49", "#06A9B2", "#E75371", "#F1E0D9")
par(las=1, pty="s")

library(igraph)
library(Matrix)

# Network dissimilarity functions ----------------------------------------------

# Network dissimilarity functions from Schieber et al. (2017)
# https://github.com/tischieber/Quantifying-Network-Structural-Dissimilarities

entropia<-function(a) {
  a<-a[which(a>0)];
  -sum(a*log(a));
}
node_distance<-function(g){
  n<-length(V(g))
  if(n==1){
    retorno=1
  }
  
  if(n>1){
    a<-Matrix(0,nrow=n,ncol=n,sparse=TRUE)
    m<-shortest.paths(g,algorithm=c("unweighted"))
    m[which(m=="Inf")]<-n
    quem<-setdiff(intersect(m,m),0)
    
    for(j in (1:length(quem))){
      l<-which(m==quem[j])/n
      linhas<-floor(l)+1
      posicoesm1<-which(l==floor(l))
      if(length(posicoesm1)>0){
        linhas[posicoesm1]<-linhas[posicoesm1]-1
      }
      a[1:n,quem[j]]<-hist(linhas,plot=FALSE,breaks=(0:n))$counts
    }
    #m<-c()
    retorno=(a/(n-1))
  }
  return(retorno)
}
nnd<-function(g){
  N<-length(V(g))
  nd<-node_distance(g)
  pdfm<-colMeans(nd)
  norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))
  return(c(pdfm,max(c(0,entropia(pdfm)-entropia(nd)/N))/norm))
}
alpha<-function(g){
  N<-length(V(g))
  r<-sort(alpha.centrality(g,exo=degree(g)/(N-1),alpha=1/N))/((N^2))
  return(c(r,max(c(0,1-sum(r)))))
}
gDist <-function(g,h,w1,w2,w3){
  
  first<-0
  second<-0
  third<-0
  g<-g
  h<-h
  N<-length(V(g))
  M<-length(V(h))
  PM<-matrix(0,ncol=max(c(M,N)))
  
  if(w1+w2>0){
    pg=nnd(g)
    PM[1:(N-1)]=pg[1:(N-1)]
    PM[length(PM)]<-pg[N]
    ph=nnd(h)
    PM[1:(M-1)]=PM[1:(M-1)]+ph[1:(M-1)]
    PM[length(PM)]<-PM[length(PM)]+ph[M]
    PM<-PM/2
    
    first<-sqrt(max(c((entropia(PM)-(entropia(pg[1:N])+entropia(ph[1:M]))/2)/log(2),0)))
    second<-abs(sqrt(pg[N+1])-sqrt(ph[M+1]))
  }
  
  if(w3>0){
    pg<-alpha(g)
    ph<-alpha(h)
    m<-max(c(length(pg),length(ph)))
    Pg<-matrix(0,ncol=m)
    Ph<-matrix(0,ncol=m)
    Pg[(m-length(pg)+1):m]<-pg
    Ph[(m-length(ph)+1):m]<-ph
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
    g<-graph.complementer(g)
    h<-graph.complementer(h)
    
    pg<-alpha(g)
    ph<-alpha(h)
    m<-max(c(length(pg),length(ph)))
    Pg<-matrix(0,ncol=m)
    Ph<-matrix(0,ncol=m)
    Pg[(m-length(pg)+1):m]<-pg
    Ph[(m-length(ph)+1):m]<-ph
    third<-third+sqrt((entropia((Pg+Ph)/2)-(entropia(pg)+entropia(ph))/2)/log(2))/2
  }
  
  return(w1*first+w2*second+w3*third)
  
}

# Data cleaning functions ------------------------------------------------------

# Take a list of synonyms and recode data 
recode_tags <- function(newcodes, field){
  tmp <- field
  for(i in 1:length(newcodes)){
    for(j in 1:length(newcodes[[i]])){
      tmp[grep(newcodes[[i]][j], field)] <- names(newcodes)[i]
    }
  }
  return(tmp)
}

# Find valid levels of study characteristics
valid_levels <- function(feature){
  which(levels(feature) %in% feature)
}

# Graph creation functions -----------------------------------------------------

# Create graphs from a list of studies
create_DAGS <- function(input.data, index){
  studies <- list()
  studies <- split(input.data, f=input.data[index])
  
  count.paths <- c()
  all.DAGS <- list()
  length(all.DAGS) <- length(studies)
  
  blank <- as.data.frame(matrix(c(NA, NA), ncol=2, byrow=TRUE))
  colnames(blank) <- c("from", "to")
  blank.DAG <- graph_from_data_frame(blank, directed=TRUE)
  
  for (i in 1:length(studies)) {
    study <- studies[[i]]
    hyp <- c()
    study.DAG <- blank.DAG
    if(nrow(study)>0){
      study <- strsplit(as.character(study$Path), ";")
      study <- unlist(lapply(study, trimws))
      ## For each study, split the paths 
      for (j in 1:length(study)){
        path <- strsplit(as.character(study[[j]]), " > ")
        line <- c()
        ## For each path, split into pairs and merge to hypothesis
        for (k in 1:(length(path[[1]])-1)){
          element <- c(path[[1]][k], path[[1]][k+1])
          line <- append(line, element)
          count.pair <- paste(element[1], element[2], sep = " > ")
          count.paths <- append(count.paths, count.pair)
        }
        
        ## For each hypothesis in a study, create a DAG
        path.edges <- as.data.frame(matrix(c(line), ncol=2,byrow=TRUE))
        colnames(path.edges) <- c("from", "to")
        path.nodes <- unique(line)
        path.DAG <- graph_from_data_frame(path.edges, vertices=path.nodes, directed=TRUE)
        study.DAG <- igraph::union(study.DAG, path.DAG)
      }    
    }
    ## Merge each hypothesis into a study DAG
    ## Create a list of all study DAGs
    all.DAGS[[i]] <- study.DAG
    
  }
  
  all.DAGS <- lapply(all.DAGS, igraph::delete.vertices, "NA")
  return(all.DAGS)
}

# Create cumulative graphs by successively adding studies
create_cumulative <- function(dags, startpoint=2){
  tmp <- list(dags[[1]])
  for(i in startpoint:length(dags)){
    tmp[i] <- list(igraph::union(tmp[[i-1]], dags[[i]]))
  }
  return(tmp)
}

# Create sliding window graphs by successively adding studies
create_windows <- function(dags, windowsize=5, startpoint){
  window_DAG <- list(dags[[1]])
  
  for(i in max(startpoint, windowsize):length(dags)){
    
    tmp <- list(dags[[i]])
    for(j in 2:(windowsize-1)){
      
      tmp[[j]] <- igraph::union(tmp[[j-1]], dags[[i-j]])
      
    }
    
    window_DAG[i] <- list(tmp[[length(tmp)]])
    
  }
  return(window_DAG)
}


# Find changes/dissimilarities between two graphs
DAG_changes <- function(dag, field, startpoint) {
  changes <- rep(NA, length(dag))
  used <- valid_levels(field)
  for (i in startpoint:length(dag)) {
    if(i %in% used){
      previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
      changes[i] <- gDist(dag[[i]], dag[[previous_valid]], 0.45, 0.45, 1)
    }
  }
  return(changes)
}

# Plotting functions -----------------------------------------------------------

# Lookup counts of how many times an edge appears across all graphs
lookup_widths <- function(allDAGs, currentDAG){
  width_lookup <- table(unlist(lapply(allDAGs, function(x){
    attr(igraph::E(x), "vnames")
  })))
  widths <- width_lookup[match(attr(igraph::E(currentDAG), "vnames"), names(width_lookup))]
  return(widths)
}

# Calculate the number of features in a graph
calc_features <- function(graph, startpoint){
  n.nodes <- n.edges <- c()
  for(i in startpoint:length(graph)){
    n.nodes[i] <- length(igraph::V(graph[[i]]))
    n.edges[i] <- length(igraph::E(graph[[i]]))
  }
  return(data.frame(n.nodes=n.nodes, n.edges=n.edges))
}

# Calculate the number of new features in a graph compared to the previous graph
new_features <- function(graph, startpoint){
  n.nodes <- n.edges <- c()
  for(i in startpoint:length(graph)){
    n.nodes[i] <- sum(!(igraph::V(graph[[i]]) %in% igraph::V(graph[[i-1]])))
    n.edges[i] <- sum(!(igraph::E(graph[[i]]) %in% igraph::E(graph[[i-1]])))
  }
  return(data.frame(new.nodes=n.nodes, new.edges=n.edges))
}

# Plot network features
plot_features <- function(features, type, xlab, ymax=NULL){
  if(is.null(ymax)){
    ymax <- ceiling((max(features, na.rm=T)*10)/100)*10
  }
    ylim <- c(0, ymax)
  plot(features[,1] ~ as.numeric(levels(type)), axes=F,
       type="l", lwd=4, ylab="Total count", xlab=xlab, col=paste(palette[2], 88, sep=""), ylim=ylim)
axis(1)
axis(2)
    lines(features[,2] ~ as.numeric(levels(type)), lty=1, lwd=4, col=palette[4])
  legend("topleft", legend=c("Pathways", "Factors and processes"), 
         lty=c(1), lwd=2.5, bty="n", col=palette[c(4,2)])
  
}

plot_changes <- function(x, y, xlab, bw=10){
  graphics::plot(y ~ as.numeric(levels(x)), 
                 ylab="Change in DAG dissimilarity", ylim=c(0,1),lwd=0,
                 xlab=xlab, pch=19, cex=2, col=paste(palette[2], "66", sep=""), axes=F)
  axis(1); axis(2)
  newy <- as.numeric(y)[!is.na(y)]
  newx <- as.numeric(levels(x))[!is.na(y)]
  lines(ksmooth(newx, newy, kernel="normal", bandwidth=bw), col=palette[4], lwd=5)
  
}

# Read in data  ----------------------------------------------------------------

metadata <- read.csv("./data/Study_characteristics_new.csv")[1:145,]

# Note: this is the cleaned version of the spreadsheet with data entry errors 
# removed and tags standardized (i.e. all the forest types reclassified, bird
# codes standardized to four-letter alpha codes, references formatted, etc.)

# For the raw, messy data, please email the corresponding author. The messy data
# includes all articles, including unscreened articles not in the random sample,
# along with the exact text entered in the spreadsheet, typos and all.

# Other datasets we need for tagging metadata and grouping terms

# Define the four letter alpha codes used for bird species
birdcodes <- read.csv("./data/IBP-Alpha-Codes20.csv")

# Read in list of synonyms and terms to group together
synonyms <- read.csv("./data/terms_grouped2.csv")


# Publication year  ------------------------------------------------------------

# Extract just the first year of data collection
metadata$FirstYear <-
  as.numeric(unlist(lapply(strsplit(metadata$Years.of.data.collection, "-"), function(x) {
    x[1]
  })))

# To include years with no studies initiated, we need to make this a factor

metadata$FirstYear <- factor(append(metadata$FirstYear, 
                             min(metadata$FirstYear):2020))[1:nrow(metadata)]

hist(
  as.numeric(as.character(metadata$FirstYear)),
  20,
  border = F,
  col = paste(palette[4], "bb", sep=""),
  xlim = c(1970, 2020),
  xlab = "First year of data collection",
  main = ""
)

# Study locations -------------------------------------------------------------- 

metadata$Latitude <- as.numeric(gsub("\\*", "", unlist(lapply(strsplit(metadata$Coordinates, ", "), function(x){x[1]}))))
metadata$Longitude <- as.numeric(gsub("\\*", "", unlist(lapply(strsplit(metadata$Coordinates, ", "), function(x){x[2]}))))

metadata$Latitude <- factor(append(metadata$Latitude, 
                                   min(metadata$Latitude):max(metadata$Latitude)))[1:nrow(metadata)]
metadata$Longitude <- factor(append(metadata$Longitude, 
                                   min(metadata$Longitude):max(metadata$Longitude)))[1:nrow(metadata)]

# Matrix type ------------------------------------------------------------------

# Because study sites can be surrounded by multiple matrices, we can set up a
# dictionary to define these landscape types based on the authors' original
# descriptions of the study locations
agriculture <- c("agric", "pasture", "crop", "field", "farm")
industry <- c("clearcut", "logg", "mining", "energy", "shale", 
              "pipeline", "seismic")
developed <- c("urban", "develop", "road", "residential", "suburban")
other_matrix <- c("plantation", "bog", "forest", "grass", "sapling",
                  "reservoir", "natural", "abrupt", "shrub", "trail")

# Add these to a named list
matrix_defs <- list(developed=developed, 
                     agriculture=agriculture, 
                     industry=industry, 
                     other=other_matrix)

# Add a wildcard to search word forms
matrix_defs <- lapply(matrix_defs, function(x){
  paste(x, "*", sep="")
})

matrix_types <- topictagger::tag_strictly(metadata$Surrounding.landscape, matrix_defs)
matrix_types[matrix_types>0] <- 1 # presence only for each type

# Bird species -----------------------------------------------------------------

# Tag each study based on taxa present
bird_dictionary <- topictagger::create_dictionary(birdcodes[,c(2,2)])
tagged_birds <- topictagger::tag_strictly(tolower(metadata$Focal.taxa), 
                                          bird_dictionary)

# Check that no species were missed and only community studies are untagged
metadata$Focal.taxa[rowSums(tagged_birds)==0]

tagged_birds <- tagged_birds[,colSums(tagged_birds)>0]
colnames(tagged_birds) <- unlist(lapply(colnames(tagged_birds), function(x){
  strsplit(x, "\\.")[[1]][1]
}))

# Path cleaning  ---------------------------------------------------------------

# Set up a placeholder
metadata$Path <- rep("", nrow(metadata))
metadata$Pathways <- gsub("=", ">", metadata$Pathways)

# Clean up the entered pathways to be interpretable for the network by splitting
# pathways that were entered together in the dataset
for(i in 1:nrow(metadata)){
  path1 <- trimws(strsplit(metadata$Pathways[i], ";")[[1]])
  
  groups <- strsplit(path1, ">")
  groups <- lapply(groups, trimws)
  
  paths <- c()
  for(j in 1:length(groups)){
    
    x <- groups[[j]]
    if(any(grep("\\+", x))){
      x <- strsplit(x, "\\+")
      x <- lapply(x, trimws)
      pathset <- c()
      for(k in 1:(length(x)-1)){
        pair.matrix <- matrix(unlist(lapply(x[[k]], paste, sep=" > ", x[[k+1]])))
        pathset <- append(pathset, pair.matrix)
      }
      x <- paste(pathset, collapse=";")
      
    }else{
      x <- paste(x, collapse=" > ")
    }
    paths <- paste(paths, x, sep=";")
  }
  metadata$Path[i] <- paste(paths, ";", sep="")
}

metadata$Path <- trimws(metadata$Path)

# Prep paths for synonym recognition to consolidate similar variables
tmp.path <- gsub(">", "_>_",metadata$Path)
tmp.path <- gsub(";", "_;_", tmp.path)
processes <- strsplit(trimws(tmp.path), "_")


# Replace synonymous terms to avoid identifying processes twice
# For example, grouping 'nest success' and 'nest failure' as the same factor

for(i in 1:length(processes)){
  y <- processes[[i]]
  y <- trimws(y)
  for(j in 1:length(y)){
    if(y[j]%in%synonyms$processes){
      newterm <- synonyms$X.1[match(y[j], synonyms$processes)]
      if(newterm!=""){
        y[j] <- newterm
      }
    }
  }
  processes[[i]] <- paste(y, collapse=" ")
}

metadata$Path <- unlist(processes)

# Create graphs for studies and years ------------------------------------------

# Graphs for all studies individually by title (index = 21)
study_DAGS <- create_DAGS(metadata, 1)

# Graphs for each year, where all studies from a given year at input at once
yearly_DAGS <- create_DAGS(metadata, which(colnames(metadata)=="FirstYear"))

# Graphs for each latitude band, where all studies from a given latitude band
# are input at once
lat_DAGS <- create_DAGS(metadata, which(colnames(metadata)=="Latitude"))

# Graphs for each longitude band, where all studies from a given longitude band
# are input at once
long_DAGS <- create_DAGS(metadata, which(colnames(metadata)=="Longitude"))

# Cumulative and sliding window graphs -----------------------------------------

cumulative_DAGs <- create_cumulative(yearly_DAGS, 
                                     min(valid_levels(metadata$FirstYear))+1)

yearwindow_DAG <- create_windows(yearly_DAGS, windowsize = 5, 
                                 startpoint = min(valid_levels(metadata$FirstYear))+1)

latband_DAG <- create_windows(lat_DAGS, windowsize = 5, 
                              startpoint = min(valid_levels(metadata$Latitude))+1)

longband_DAG <- create_windows(long_DAGS, windowsize = 5, 
                              startpoint = min(valid_levels(metadata$Longitude)+1))

# # Stopping criteria and estimation of missing information ---------------------

n.sim <- 300

# Placholders
factors <- pathways <- modfact <- modpath <- new.factors <- new.pathways <- poismodfact <- poismodpath <- array(dim=c(n.sim, length(study_DAGS)))

for(i in 1:n.sim){
  # Randomly reorder graphs to create a cumulative network
  tmpdag <- create_cumulative(sample(study_DAGS))

  # New graph features added
  factors[i,] <- as.numeric(new_features(tmpdag, 2)[,1]>0)
  pathways[i,] <- as.numeric(new_features(tmpdag, 2)[,2]>0)

  new.factors[i, ] <- new_features(tmpdag, 2)[,1]
  new.pathways[i,] <- new_features(tmpdag, 2)[,2]

  # Fit a binomial GLM and predict the probability of new graph features
  modfact[i,2:length(study_DAGS)] <- plogis(predict(glm(factors[i,] ~ seq(1, length(study_DAGS), 1),
                                         family = "binomial")))
  modpath[i,2:length(study_DAGS)] <- plogis(predict(glm(pathways[i,] ~ seq(1, length(study_DAGS), 1),
                                         family = "binomial")))

  poismodfact[i,2:length(study_DAGS)] <- predict(glm(new.factors[i,] ~ seq(1, length(study_DAGS), 1),
                                                        family = "poisson"), type="response")
  poismodpath[i,2:length(study_DAGS)] <- predict(glm(new.pathways[i,] ~ seq(1, length(study_DAGS), 1),
                                                        family = "poisson"), type="response")

  }


# Stopping criteria: the number of new factors is less than 1 per every two 
# studies (<0.50) and the probability of any new factors is less than 0.25

plot(apply(poismodfact, 2, mean), type="n", ylim=c(0,5), xlim=c(0,150),
     xlab="Number of included studies",
     ylab="Estimated number of new graph features", axes=F)
axis(1, at = seq(0,150,25)); axis(2, at=seq(0,5,.5))
ci <- apply(poismodfact, 2, quantile, c(0.025, 0.975), na.rm=T)
newx <- seq(1, length(study_DAGS), 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste(palette[2], "66", sep=""), border = F)
lines(apply(poismodfact, 2, mean, na.rm=T), col=palette[2], lwd=4)

ci <- apply(poismodpath, 2, quantile, c(0.025, 0.975), na.rm=T)
newx <- seq(1, length(study_DAGS), 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste(palette[4], "44", sep=""), border = F)
lines(apply(poismodpath, 2, mean, na.rm=T), col=palette[4], lwd=4)

legend("topright",
       legend=c("Mean", "95% confidence interval", "Pathways", "Factors"),
       lwd=c(2, NA, NA, NA, NA), pch=c(NA, 15, 15, 15),
       col=c("black", "grey95", palette[c(4,2)]), bty="n",
       pt.cex = c(NA, 2,2,2))

# Prob of any new features is less than 25% CI
plot(apply(modfact, 2, mean), type="n", ylim=c(0,1), xlim=c(0,150),
     xlab="Number of included studies",
     ylab="Probability of new graph features", axes=F)
axis(1, at = seq(0,150,25)); axis(2, at=seq(0,1,.1))
ci <- apply(modfact, 2, quantile, c(0.025, 0.975), na.rm=T)
newx <- seq(1, length(study_DAGS), 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste(palette[2], "66", sep=""), border = F)
lines(apply(modfact, 2, mean, na.rm=T), col=palette[2], lwd=4)

ci <- apply(modpath, 2, quantile, c(0.025, 0.975), na.rm=T)
newx <- seq(1, length(study_DAGS), 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste(palette[4], "44", sep=""), border = F)
lines(apply(modpath, 2, mean, na.rm=T), col=palette[4], lwd=4)

legend("topright",
       legend=c("Mean", "95% confidence interval", "Pathways", "Factors"),
       lwd=c(2, NA, NA, NA, NA), pch=c(NA, 15, 15, 15),
       col=c("black", "grey95", palette[c(4,2)]), bty="n",
       pt.cex = c(NA, 2,2,2))


# Estimate accumulation of network features

factors2 <- names(igraph::V(cumulative_DAGs[[length(cumulative_DAGs)]]))
pathways2 <- attr(E(cumulative_DAGs[[length(cumulative_DAGs)]]), "vnames")

factor.dat <- array(dim=c(length(study_DAGS), length(factors2)))
pathway.dat <- array(dim=c(length(study_DAGS), length(pathways2)))

for(i in 1:length(study_DAGS)){
  factor.dat[i,] <- as.numeric(factors2 %in%  names(igraph::V(study_DAGS[[i]])))
  pathway.dat[i,] <- as.numeric(pathways2 %in%  attr(E(study_DAGS[[i]]), "vnames"))
}

est.factors <- vegan::poolaccum(factor.dat, permutations = 500)
est.pathways <- vegan::poolaccum(pathway.dat, permutations = 500)

plot(apply(est.factors$jack1, 1, mean), type="n", ylim=c(0,600), xlim=c(0,150), axes=F,
     ylab="Extrapolated total count of graph features",
     xlab="Number of included studies")
axis(1); axis(2)
ci <- apply(est.factors$jack1, 1, quantile, c(0.025, 0.975))
newx <- seq(1, length(study_DAGS)-2, 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste(palette[2], "66", sep=""), border = F)
lines(apply(est.factors$jack1, 1, mean), col=palette[2], lwd=4)

ci <- apply(est.pathways$jack1, 1, quantile, c(0.025, 0.975))
newx <- seq(1, length(study_DAGS)-2, 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])),
        col=paste(palette[4], "44", sep=""), border = F)
lines(apply(est.pathways$jack1, 1, mean), col=palette[4], lwd=4)
legend("topright",
       legend=c("Mean", "95% confidence interval", "Pathways", "Factors"),
       lwd=c(2, NA, NA, NA, NA), pch=c(NA, 15, 15, 15),
       col=c("black", "grey95", palette[c(4,2)]), bty="n",
       pt.cex = c(NA, 2,2,2))

# Global graph from union of all studies ---------------------------------------

metaDAG <- cumulative_DAGs[[length(cumulative_DAGs)]]

metaDAG <- igraph::simplify(metaDAG)

# At this point, there are lots of analyses and visualizations that can be done
# Here, we just show a few examples of plotting network features and structural
# changes across time and space, along with how to subset graphs by moderators.

# Plot network features and changes --------------------------------------------

plot_features(calc_features(cumulative_DAGs, 5), 
              metadata$FirstYear, xlab="Year", ymax=350)

plot_features(calc_features(yearwindow_DAG, 5), 
              metadata$FirstYear, xlab="Year", ymax=120)

yearly_changes <- DAG_changes(cumulative_DAGs, 
                              metadata$FirstYear, 1)
plot_changes(metadata$FirstYear,
             yearly_changes, xlab="First year of data collection")

fiveyr.changes <- DAG_changes(yearwindow_DAG, 
                              metadata$FirstYear, 6)
plot_changes(metadata$FirstYear, 
             fiveyr.changes, xlab="First year of data collection")

plot_features(calc_features(latband_DAG, 10), 
              metadata$Latitude, xlab="Latitude", ymax=200)
plot_features(calc_features(longband_DAG, 19), 
              metadata$Longitude, xlab="Longitude", ymax=100)

latchanges_window <- DAG_changes(latband_DAG, metadata$Latitude, 8)
plot_changes(metadata$Latitude, latchanges_window, xlab="Latitude", bw=10)

longchanges_window <- DAG_changes(longband_DAG, metadata$Longitude, 17)
plot_changes(metadata$Longitude, longchanges_window, xlab="Longitude", bw=10)

# Example subgraphs by study characteristics -----------------------------------

agric <- study_DAGS[which(matrix_types[,'agriculture']>0)]
agric2 <- list(agric[[1]])

for(i in 2:length(agric)){
  agric2[i] <- list(igraph::union(agric2[[i-1]], agric[[i]]))
}

agriculture_network <- agric2[[i]]

# Example graph for a single species -------------------------------------------

birds <- which(tagged_birds[,which(colnames(tagged_birds)=="OVEN")]>0)

bird <- study_DAGS[birds]
bird2 <- list(bird[[1]])

for(i in 2:length(bird)){
  bird2[i] <- list(igraph::union(bird2[[i-1]], bird[[i]]))
}

ovenbirds <- bird2[[i]]
