### Workspace set up -----------------------------------------------------------
setwd("~/Documents/research-projects/ebcmd/")
library(igraph)

palette <- c("#87c6d8", "#d1d6dc", "#f2b5b4", "#fe7b81", "#7b4567")
plot(1:5, rep(0, 5), pch=20, cex=25, col=palette, xpd=T, axes=F, ylab="", xlab="")

palette2 <- c("#ff8ba0", "#ffa5a1", "#ffc2a6", "#f9de87", "#fce9c4", "#d3e5cf")
plot(1:6, rep(0, 6), pch=20, cex=25, col=palette2, xpd=T, axes=F, ylab="", xlab="")

# palette3 <- c("#f89da7", "#fcccd4", "#fbdea2", "#f2e2c6", "#8eb695")
# plot(1:5, rep(0, 5), pch=20, cex=25, col=palette3, xpd=T, axes=F, ylab="", xlab="")




### Functions -----------------------------------------------------------------

# from Schieber2017
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



recode_tags <- function(newcodes, field){
  tmp <- field
  for(i in 1:length(newcodes)){
    for(j in 1:length(newcodes[[i]])){
      tmp[grep(newcodes[[i]][j], field)] <- names(newcodes)[i]
    }
  }
  return(tmp)
}

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


valid_levels <- function(feature){
  which(levels(feature) %in% feature)
}

create_cumulative <- function(dags, startpoint=2){
  tmp <- list(dags[[1]])
  for(i in startpoint:length(dags)){
    tmp[i] <- list(igraph::union(tmp[[i-1]], dags[[i]]))
  }
  return(tmp)
}

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
lookup_widths <- function(allDAGs, currentDAG){
  width_lookup <- table(unlist(lapply(allDAGs, function(x){
    attr(igraph::E(x), "vnames")
  })))
  widths <- width_lookup[match(attr(igraph::E(currentDAG), "vnames"), names(width_lookup))]
  return(widths)
}

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

calc_features <- function(graph, startpoint){
  n.nodes <- n.edges <- c()
  for(i in startpoint:length(graph)){
    n.nodes[i] <- length(igraph::V(graph[[i]]))
    n.edges[i] <- length(igraph::E(graph[[i]]))
  }
  return(data.frame(n.nodes=n.nodes, n.edges=n.edges))
}

new_features <- function(graph, startpoint){
  n.nodes <- n.edges <- c()
  for(i in startpoint:length(graph)){
    n.nodes[i] <- sum(!(igraph::V(graph[[i]]) %in% igraph::V(graph[[i-1]])))
    n.edges[i] <- sum(!(igraph::E(graph[[i]]) %in% igraph::E(graph[[i-1]])))
  }
  return(data.frame(new.nodes=n.nodes, new.edges=n.edges))
}



plot_features <- function(features, type, xlab){
  ylim <- c(0, max(features, na.rm = T))
  plot(features[,1] ~ as.numeric(levels(type)), 
       type="l", lwd=4, ylab="Total count", xlab=xlab, col=palette2[1], ylim=ylim)
  lines(features[,2] ~ as.numeric(levels(type)), lty=2, lwd=4, col=palette2[3])
  legend("topleft", legend=c("Factors and processes", "Pathways"), 
         lty=c(1,2), lwd=2.5, bty="n", col=palette2[3])
  
}

# Read in data  ----------------------------------------------------------------

input.data <- readxl::read_xlsx("~/Downloads/fulltext2 (14).xlsx", sheet=2)
input.data <- input.data[1:157,]
excluded <- input.data[!(input.data$include%in%"yes"),1:40]
input.data <- input.data[input.data$include%in%"yes",1:40]

table(excluded$reasons)

reasons_excl <- excluded$reasons

sort(table(reasons_excl))

no_response <- c()
wrong_pop <- c("North America", "artificial", "focal species", "specialist", "overwinter", "season", "not forest")
wrong_comp <- c("comparing")
existing_data <- c("already included", "data comes", "included as", "existing data")
no_sensitivity <- c("no measure of edge or area", "no measure of distance", "no distance", "no mechanisms", "patch size",
                    "sensitivity", "gradient", "harvest")

reasons <- list(no_sensitivity=no_sensitivity, existing_data=existing_data, wrong_comp=wrong_comp, wrong_pop=wrong_pop)

reasons_excl <- recode_tags(reasons, reasons_excl)
table(reasons_excl)


input.data <- data.frame(input.data)


# Publication year  ------------------------------------------------------------

input.data$firstyear <- as.numeric(unlist(lapply(strsplit(input.data$years, "-"), function(x){x[1]})))

possible_years <- factor(append(input.data$firstyear, seq(1975, 2020, 1)))[1:nrow(input.data)]
input.data$firstyear <- possible_years
input.data <- input.data[order(input.data$firstyear),]

png("./figures/firstyears.png", width=8,height=8,units="in",res=1200)
par(las=1, pty="s")
hist(as.numeric(as.character(input.data$firstyear)), 20, border=F, 
     col=palette2[1], xlim=c(1965,2020), xlab="First year of data collection", main="")
dev.off()
#clean up
rm(possible_years)

# Latitude and longitude  ------------------------------------------------------ 

input.data$latitude <- unlist(lapply(strsplit(input.data$coordinates, ", "), function(x){x[1]}))
input.data$longitude <- unlist(lapply(strsplit(input.data$coordinates, ", "), function(x){x[2]}))
input.data$longitude[which(input.data$longitude==106.666667)] <- -106.666667

#input.data$longitude[56] <- "-66.083333"

input.data$latband <- as.numeric(unlist(lapply(strsplit(input.data$latitude, "\\."), function(x){x[1]})))
latbands <- factor(append(input.data$latband, seq(24, 72, 1)))[1:nrow(input.data)]
input.data$latband <- latbands

input.data$longband <- as.numeric(unlist(lapply(strsplit(input.data$longitude, "\\."), function(x){x[1]})))
longbands <- factor(append(input.data$longband, seq(-140, -52, 1)))[1:nrow(input.data)]
input.data$longband <- longbands

maps::map(regions = c("usa(?!:alaska)", "canada"), col="white")
maps::map("state", add=T, col="grey90")
maps::map(region="canada", add=T, col="grey60")
maps::map(region="usa", add=T, col="grey60")

points(input.data$latitude ~ input.data$longitude, 
       pch=20, cex=2, col="white", lwd=1, lty=1)

points(input.data$latitude ~ input.data$longitude, 
       pch=20, cex=2, col="#ff000040", lwd=1, lty=1)


forests <- raster::raster("~/Documents/research-projects/ebcmd/data/NA_TREEAGE_1096/data/ca04_usak06_1km.tif")

states <- raster::getData(country="USA", level=1)
provinces <- raster::getData(country="Canada", level=1)

canada <- sp::spTransform(provinces, sp::proj4string(forests))
usa <- sp::spTransform(states, sp::proj4string(forests))

new_clip <- raster::extent(raster::extent(forests)[1]+3000000,
                           raster::extent(canada)[2]+10000,
                           raster::extent(usa)[3],
                           raster::extent(canada)[4])

forested <- raster::crop(forests, new_clip)

studysites <- data.frame(as.numeric(input.data$latitude), as.numeric(input.data$longitude))
names(studysites) <- c("latitude", "longitude")
sp::coordinates(studysites) <- ~ longitude + latitude
sp::proj4string(studysites) <- sp::proj4string(states)
studysites <- sp::spTransform(studysites, sp::proj4string(forests))
# 
# 
# 
# png("./figures/study_locations.png", width=12,height=8,units="in",res=1200)
# par(xpd=T)
# sp::plot(forested, col=palette2[6], axes=F, box=F, legend=F)
# sp::plot(canada, add=T, border="grey50")
# sp::plot(usa, add=T, border="grey50")
# 
# points(studysites, col="white", pch=20, lwd=1, cex=1.5)
# points(studysites, col=paste(palette2[1], "99", sep=""), pch=20, cex=1.5)
# points(studysites, col="grey20", pch=1, cex=1.25, lwd=0.5)
# dev.off()

# clean up
rm(latbands, longbands)

# Recode matrix type -----------------------------------------------------------

agriculture <- c("agric", "pasture", "crop", "field", "farm")
industry <- c("clearcut", "logg", "mining", "energy", "shale", "pipeline", "seismic")
developed <- c("urban", "develop", "road", "residential")
other_matrix <- c("plantation", "bog", "forest", "grass", "sapling", "reservoir")



matrix_types <- list(developed=developed, 
                     agriculture=agriculture, 
                     industry=industry, 
                     other=other_matrix)
matrix_types <- lapply(matrix_types, function(x){
  paste(x, "*", sep="")
})

input.data$recode_matrix <- recode_tags(matrix_types, input.data$matrix)

matrix_hits <- topictagger::tag_strictly(input.data$matrix, matrix_types)
matrix_hits[matrix_hits>0] <- 1
prop.table(colSums(matrix_hits))

# Recode forest type -----------------------------------------------------------

forest_types <- input.data$forest_class

prop.table(sort(table(append(forest_types, rep("hammock", 2)))))


oak_hickory <- c("oak-hickory", "hickory-oak")

industry <- c("clearcut", "logg", "mining", "energy", "shale", "pipeline", "seismic")
developed <- c("urban", "develop", "road", "residential")
other_matrix <- c("plantation", "bog", "forest", "grass", "sapling", "reservoir")

matrix_types <- list(agriculture=agriculture, 
                     developed=developed, 
                     industry=industry, 
                     other=other_matrix)

input.data$recode_matrix <- recode_tags(matrix_types, input.data$matrix)


# Path cleaning  --------------------------------------------------------------

input.data$Path <- rep("", nrow(input.data))

for(i in 1:nrow(input.data)){
  path1 <- strsplit(input.data$test.pathways[i], ";")[[1]]
  
  groups <- strsplit(path1, ">")[[1]]
  groups <- unlist(lapply(groups, trimws))
  groups <- strsplit(groups, "\\+")
  groups <- lapply(groups, trimws)
  groups
  
  paths <- c()
  for(j in 1:(length(groups)-1)){
    pair.matrix <- matrix(unlist(lapply(groups[[j]], paste, sep=" > ", groups[[j+1]])))
    pathset <- paste(pair.matrix, sep="; ")
    paths <- paste(paths, pathset, sep="; ")
  }
  input.data$Path[i] <- paths
  
}

input.data$Path <- trimws(input.data$Path)
input.data$Path <- substring(input.data$Path, 3)

input.data$Path
processes <- trimws(unlist(strsplit(trimws(gsub(";", " > ", input.data$Path)), " > ")))
sort(table(processes))

bird_abundance <- c("bird abundane", "bird occurrence", "occupancy", "population density")
patch_size <- c("island size", "patch area", "forest width", "fragmentation")
nest_success <- c("nest survival", "daily nest success")
dispersal <- c("emigration")
bird_community <- c("bird species richness", "bird diversity", "avian community composition")
parental_feeding <- c("nest provisionining")
territory_density <- c("nest density")
habitat_selection <- c("nest placement")
predators <- c("mammalian predator abundance", "avian predator")
age <- c("age structure", "male age")
vegetation <- c("vegetation complexity")


input.data$Path <- gsub("movement", "extraterritorial forays", input.data$Path)

input.data$Path <- gsub("extra-pair", "extra pair", input.data$Path)

synonyms <- list(bird_abundance, patch_size, nest_success, bird_community, territory_density, habitat_selection, age,
                 parental_feeding, dispersal, predators, vegetation)
names(synonyms) <- c("bird abundance", "patch size", "nest success", "bird community", "territory density", 
                     "habitat selection", "age structure", "parental feeding", "dispersal", "predator abundance", "vegetation structure")

for(i in 1:length(synonyms)){
  for(j in 1:length(synonyms[[i]])){
    input.data$Path <- gsub(synonyms[[i]][j], names(synonyms)[i], input.data$Path)
  }
}


# Recode sensitivity metrics  --------------------------------------------------

distance_to_edge <- c("distance")
patch_size <- c("patch", "continuous", "island", "width")
forest_dissection <- c("dissection")

sensitivity_types <- list(distance_to_edge=distance_to_edge, 
                     patch_size=patch_size, 
                     forest_dissection=forest_dissection)

input.data$recode_metrics <- recode_tags(sensitivity_types, input.data$metrics)


metricxmatrix <- table(input.data$recode_matrix, input.data$recode_metrics)
metricxmatrix <- metricxmatrix[!rownames(metricxmatrix)=="N/A",]
colnames(metricxmatrix) <- gsub("_", " ", colnames(metricxmatrix))
metricxmatrix <- data.frame(metricxmatrix)

library(ggplot2)

ggplot(metricxmatrix, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_continuous(high = palette2[1],
                        low = "white",
                        na.value = "white") + # fill with color
  theme(
    line = element_blank(),
    # remove the background, tickmarks, etc
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size = 10)
  ) +
  coord_equal()
ggsave("./figures/metricsxmatrix.png")

# Create mini-DAGs  --------------------------------------------------------

study_DAGS <- create_DAGS(input.data, 19)
yearly_DAGS <- create_DAGS(input.data, which(colnames(input.data)=="firstyear"))
lat_DAGS <- create_DAGS(input.data, which(colnames(input.data)=="latband"))
long_DAGS <- create_DAGS(input.data, which(colnames(input.data)=="longband"))

# Network DAGs   --------------------------------------------------------

cumulative_DAGs <- create_cumulative(yearly_DAGS, min(valid_levels(input.data$firstyear))+1)
cumulative_lat_DAGs <- create_cumulative(lat_DAGS, 2)
cumulative_long_DAGs <- create_cumulative(long_DAGS, 2)

yearwindow_DAG <- create_windows(yearly_DAGS, windowsize = 5, 
                                 startpoint = min(valid_levels(input.data$firstyear))+1)

latband_DAG <- create_windows(lat_DAGS, windowsize = 5, 
                              startpoint = min(valid_levels(input.data$latband))+1)

longband_DAG <- create_windows(long_DAGS, windowsize = 5, 
                              startpoint = min(valid_levels(input.data$longband)+1))


# sensitivity analysis -----

factors <- pathways <- changes <- list()
n.sim <- 1
for(i in 1:n.sim){
  tmpdag <- create_cumulative(sample(study_DAGS))
  factors[[i]] <- calc_features(tmpdag, 1)[,1]
  pathways[[i]] <- calc_features(tmpdag, 1)[,2]
#  changes[[i]] <- DAG_changes(tmpdag, factor(input.data$Title), 1)
  }


factors <- names(igraph::V(cumulative_DAGs[[length(cumulative_DAGs)]]))
specdat <- array(dim=c(length(study_DAGS), length(factors)))
                 
for(i in 1:length(study_DAGS)){
  specdat[i,] <- as.numeric(factors %in%  names(igraph::V(study_DAGS[[i]])))
}

estN <- vegan::poolaccum(specdat, permutations = 200, index="chao2")

plot(apply(estN$jack1, 1, quantile, c(0.025, 0.975))[1,], type="l", ylim=c(0,200))
lines(apply(estN$jack1, 1, quantile, c(0.025, 0.975))[2,], type="l")



par(las=1, pty="s")
plot(factors[[1]], ylim=c(0,50), xlim=c(0,100), type="n", ylab="Number of factors", xlab="Number of studies")
for(i in 1:n.sim){
  lines(factors[[i]], col=paste(palette2[1], "66", sep=""))
}
lines((seq(1,100,1) ~ seq(1,100,1)), lty=2)
lines(rowSums(matrix(unlist(factors), ncol=n.sim, byrow=F))/n.sim, lwd=2)


plot(pathways[[1]], ylim=c(0,90), type="n", ylab="Number of pathways", xlab="Number of studies")
for(i in 1:n.sim){
  lines(pathways[[i]], col=paste(palette2[1], "44", sep=""))
}

lines(rowSums(matrix(unlist(pathways), ncol=n.sim, byrow=F))/100, lwd=2)


# plot(changes[[1]], ylim=c(0,1), type="n", ylab="Change in DAG structure", xlab="Number of studies")
# for(i in 1:n.sim){
#   lines(changes[[i]], col=paste(palette2[1], "44", sep=""))
# }
# length(changes)
# lines(rowSums(matrix(unlist(changes), ncol=17, byrow=F))/17, lwd=2)

library(Matrix)




#### Transitive reduction   ----------------------------------------------------

metaDAG <- cumulative_DAGs[[length(cumulative_DAGs)]]


tkplot(metaDAG, 
            edge.arrow.size=.7, 
            edge.color=palette2[1],
            vertex.color="#ffffff00", 
            vertex.frame.color="#ffffff00", 
            vertex.label.family="Arial",
            vertex.label.dist=0, 
            vertex.label.color="black", 
            vertex.label.cex=.75,vertex.size=sqrt(strength(metaDAG))*5,
            edge.width=as.numeric(lookup_widths(study_DAGS, metaDAG)),
            )


topsort.metaDAG <- topo_sort(metaDAG)
metaDAG.adj <- as.matrix(as_adj(metaDAG))
reduced.adj <- mnem::transitive.reduction(metaDAG.adj)
reducedDAG <- igraph::graph_from_adjacency_matrix(reduced.adj)

tkplot(reducedDAG, 
            edge.arrow.size=.5,
            edge.color=palette2[1],
            vertex.color="#ffffff00", 
            vertex.frame.color="#ffffff00", 
            vertex.label.family="Arial",
            vertex.label.dist=0, 
            vertex.label.color="black", 
            vertex.label.cex=.75, 
            edge.width=.5)


edgedat <- cbind(attr(igraph::E(metaDAG), "vnames"), as.numeric(lookup_widths(study_DAGS, metaDAG)))
edgedat[order(as.numeric(edgedat[,2])),]

#### Dissimilarity  --------------------------------------------------------
library(Matrix)
plot_changes <- function(x, y, xlab, bw=10){
  graphics::plot(y ~ as.numeric(levels(x)), 
       ylab="Change in DAG dissimilarity", 
       xlab=xlab, pch=19, cex=2, col=palette2[2])
  newy <- as.numeric(y)[!is.na(y)]
  newx <- as.numeric(levels(x))[!is.na(y)]
  lines(ksmooth(newx, newy, kernel="normal", bandwidth=bw), col=palette2[1], lwd=5)
  
}


par(las=1, pty="s")

yearly_changes <- DAG_changes(cumulative_DAGs, input.data$firstyear, 1)
plot_changes(input.data$firstyear, yearly_changes, xlab="First year of data collection")

fiveyr.changes <- DAG_changes(yearwindow_DAG, input.data$firstyear, 6)
plot_changes(input.data$firstyear, fiveyr.changes, xlab="First year of data collection")


latchanges_window <- DAG_changes(latband_DAG, input.data$latband, 8)
plot_changes(input.data$latband, latchanges_window, xlab="Latitude", bw=10)

valid_levels(input.data$longband)
longchanges <- DAG_changes(cumulative_long_DAGs, input.data$longband, 17)
plot_changes(input.data$longband, longchanges, xlab="Longitude", bw=10)

longchanges_window <- DAG_changes(longband_DAG, input.data$longband, 17)
plot_changes(input.data$longband, longchanges_window, xlab="Longitude", bw=10)


# Network nodes and edges ------------------------------------------------------


par(las=1, pty="s")

yearwindows <- c()
for(i in 5:length(table(input.data$firstyear))){
  yearwindows[i] <- sum(table(input.data$firstyear)[(i-4):i])
}


latcounts <- c()
for(i in 5:length(table(input.data$latband))){
  latcounts[i] <- sum(table(input.data$latband)[(i-4):i])
}

yearwindows[yearwindows==0] <- 1

plot_features(new_features(yearwindow_DAG, 6)/yearwindows[-1], input.data$firstyear, xlab="Year")


plot(new_features(cumulative_DAGs, 2)[,1]/yearwindows[-1] ~ seq(1975, 2020,1), xlab="Year")

lines(yearwindows ~ seq(1975, 2020, 1))
plot_features(calc_features(latband_DAG, 10), input.data$latband, xlab="Latitude")
plot_features(calc_features(longband_DAG, 19), input.data$longband, xlab="Longitude")

plot_features(calc_features(cumulative_DAGs, 5), input.data$firstyear, xlab="Year")
plot_features(calc_features(cumulative_lat_DAGs, 10), input.data$latband, xlab="Latitude")
plot_features(calc_features(cumulative_long_DAGs, 19), input.data$longband, xlab="Longitude")

### Matrix types ####
agric <- study_DAGS[which(matrix_hits[,'agriculture']>0)]
agric2 <- list(agric[[1]])

for(i in 2:length(agric)){
  agric2[i] <- list(igraph::union(agric2[[i-1]], agric[[i]]))
}


plot.igraph(agric2[[i]], main="Agriculture",
            edge.arrow.size=0.5, edge.color=palette2[1], 
            vertex.color="white", 
            vertex.frame.color="white", vertex.label.family="Arial",
            vertex.label.dist=0, vertex.label.color="black",  vertex.label.cex=.7,
            edge.width=as.numeric(lookup_widths(study_DAGS, agric2[[i]])))


develop <- study_DAGS[which(matrix_hits[,'developed']>0)]
develop2 <- list(develop[[1]])

for(i in 2:length(develop)){
  develop2[i] <- list(igraph::union(develop2[[i-1]], develop[[i]]))
}

plot.igraph(develop2[[i]], main="Development",
            edge.arrow.size=.5, edge.color=palette2[1],
            vertex.color="white", vertex.frame.color="white", 
            vertex.label.family="Arial",
            vertex.label.dist=0, vertex.label.color="black",  vertex.label.cex=.7,
            edge.width=as.numeric(lookup_widths(study_DAGS, develop2[[i]])))

industry <- study_DAGS[which(matrix_hits[,'industry']>0)]
industry2 <- list(industry[[1]])

for(i in 2:length(industry)){
  industry2[i] <- list(igraph::union(industry2[[i-1]], industry[[i]]))
}


plot.igraph(industry2[[i]], main="Industry",
            edge.arrow.size=.5, edge.color=palette2[1],
            vertex.color="white", vertex.frame.color="white", vertex.label.family="Arial",
            vertex.label.dist=0, vertex.label.color="black", vertex.label.cex=.7,
            edge.width=as.numeric(lookup_widths(study_DAGS, industry2[[i]])))



study_distances <- sp::spDists(studysites)
#graph_distances <- array(dim=dim(study_distances))

for(i in 1:length(study_DAGS)){
  for(j in 1:length(study_DAGS)){
    if(i!=j){
      graph_distances[i,j] <- gDist(study_DAGS[[i]], study_DAGS[[j]], 0.45, 0.45, 1)
    }
  }
}

plot(study_distances[upper.tri(study_distances)], graph_distances[upper.tri(graph_distances)])
abline(lm(graph_distances[upper.tri(graph_distances)] ~ study_distances[upper.tri(study_distances)]))


### Metric types #### 
# forestwidth <- study_DAGS[which(input.data$recode_metrics=="forest width")]
# forestwidth2 <- list(forestwidth[[1]])
# 
# for(i in 2:length(forestwidth)){
#   forestwidth2[i] <- list(igraph::union(forestwidth2[[i-1]], forestwidth[[i]]))
# }
# width_lookup <- table(unlist(lapply(forestwidth, function(x){
#   attr(igraph::E(x), "vnames")
# })))
# 
# widths <- width_lookup[match(attr(igraph::E(forestwidth2[[i]]), "vnames"), names(width_lookup))]
# 
# 
# plot.igraph(forestwidth2[[i]], main="Forest Width",
#             edge.arrow.size=0.5, edge.color="black", 
#             vertex.color="white", vertex.frame.color="white", vertex.label.family="Arial",
#             vertex.label.dist=0, vertex.label.color="black", edge.width=as.numeric(widths*widths))
# 


patch <- study_DAGS[which(input.data$recode_metrics=="patch_size")]
patch2 <- list(patch[[1]])

for(i in 2:length(patch)){
  patch2[i] <- list(igraph::union(patch2[[i-1]], patch[[i]]))
}
width_lookup <- table(unlist(lapply(patch, function(x){
  attr(igraph::E(x), "vnames")
})))

widths <- width_lookup[match(attr(igraph::E(patch2[[i]]), "vnames"), names(width_lookup))]


plot.igraph(patch2[[i]], main="Patch Size",
            edge.arrow.size=0.75, edge.color=palette2[1], 
            vertex.color="white", vertex.frame.color="white", vertex.label.family="Arial",
            vertex.label.dist=0, vertex.label.color="black", edge.width=as.numeric(widths))



edgedist <- study_DAGS[which(input.data$recode_metrics=="distance_to_edge")]
edgedist2 <- list(edgedist[[1]])

for(i in 2:length(edgedist)){
  edgedist2[i] <- list(igraph::union(edgedist2[[i-1]], edgedist[[i]]))
}
width_lookup <- table(unlist(lapply(edgedist, function(x){
  attr(igraph::E(x), "vnames")
})))

widths <- width_lookup[match(attr(igraph::E(edgedist2[[i]]), "vnames"), names(width_lookup))]


plot.igraph(edgedist2[[i]], main="Distance to edge",
            edge.arrow.size=0.5, edge.color=palette2[1], 
            vertex.color="white", vertex.frame.color="white", vertex.label.family="Arial",
            vertex.label.dist=0, vertex.label.color="black", edge.width=as.numeric(widths*widths))


birds <- which(bird_dat[,which(fourlettercodes=="oven")]>0)
#birds <- birds[!(birds %in% which(rowSums(bird_dat[,which(fourlettercodes%in%c("oven", "woth"))])>1))]

bird <- study_DAGS[birds]
bird2 <- list(bird[[1]])

for(i in 2:length(bird)){
  bird2[i] <- list(igraph::union(bird2[[i-1]], bird[[i]]))
}
width_lookup <- table(unlist(lapply(bird, function(x){
  attr(igraph::E(x), "vnames")
})))

widths <- width_lookup[match(attr(igraph::E(bird2[[i]]), "vnames"), names(width_lookup))]


plot.igraph(bird2[[i]], main="Ovenbird",
            edge.arrow.size=0.5, edge.color=palette2[1], 
            vertex.color="white", vertex.frame.color="white", vertex.label.family="Arial",
            vertex.label.dist=0, vertex.label.color="black", edge.width=as.numeric(widths*widths))



#### Random comments ####
## need to take count.paths and make it weight the different paths 
## also need a way to mark if paths are hypothesized or data-driven
## maybe this can be done separately?
## like once the final graph is done, then we go back to the data and count paths/hypotheses
## we'll want to plot this in DiagrammeR anyways since it is prettier than igraph

## do we want a separate graph for each study first? that seems silly if they're hypothesized

# 1 for all paths i in X, W = union(x1, x2, ...xi)
# 2 if W contains modules, manually remove modules disconnected to the primary node 
# else, 3 D = W
# if D contains cycles, stop for manual fix
# else, topological sort D
# take the adjacency matrix A for D
# make the matrix S, which is equal to (A^2, OR A^3, OR ... A^k) where k is the number of vertices minus 1 and is the number of hops
# the final A is all the elements of A that are not in S

example1 <- delete.vertices(current.DAG[[1]], "NA")
example2 <- delete.vertices(current.DAG[[2]], "NA")
example3 <- delete.vertices(current.DAG[[3]], "NA")
png("graphs/path-combination.png", width = 8, height = 8, units = 'in', res = 400)
par(mfrow=c(2,2))
set.seed(13)
plot(example3,
     edge.arrow.size=.5, edge.color=c("white", "white", "white", "white", "white", "white", "white", "black"),
     vertex.color="white", vertex.frame.color="white", 
     vertex.label.dist=0, vertex.label.color=c("black", "black", "white", "white", "white", "white", "white", "white", "white"))

set.seed(13)
plot(example3,
     edge.arrow.size=.5, edge.color=c("white", "white", "white", "white", "white", "black", "black", "black"),
     vertex.color="white", vertex.frame.color="white", 
     vertex.label.dist=0, vertex.label.color=c("black", "black", "black", "black", "white", "white", "white", "white", "white"))

set.seed(13)
plot(example3,
     edge.arrow.size=.5, edge.color=c("white", "white", "white", "blue", "blue", "black", "black", "black"),
     vertex.color="white", vertex.frame.color="white", 
     vertex.label.dist=0, vertex.label.color=c("black", "black", "black", "black", "blue", "blue", "white", "white", "white"))
legend("topright", legend=c(paste("D(G', G) = ", round(changes[1],2), sep="")))

set.seed(13)
plot(example3,
     edge.arrow.size=.5, edge.color=c("red", "red", "red", "blue", "blue", "black", "black", "black"),
     vertex.color="white", vertex.frame.color="white", 
     vertex.label.dist=0, vertex.label.color=c("black", "black", "black", "black", "blue", "blue", "red", "red", "red"))
legend("topright", legend=c(paste("D(G', G) = ", round(changes[2],2), sep="")))

dev.off()
