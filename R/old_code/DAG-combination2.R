# Workspace set up -------------------------------------------------------------

palette <- c("#ff8ba0", "#ffa5a1", "#ffc2a6", "#f9de87", "#fce9c4", "#d3e5cf")
plot(1:6, rep(0, 6), pch=20, cex=25, col=palette, xpd=T, axes=F, ylab="", xlab="")

par(las=1, pty="s")

library(igraph)
library(Matrix)
library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)

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
plot_features <- function(features, type, xlab){
  ylim <- c(0, max(features, na.rm = T))
  plot(features[,1] ~ as.numeric(levels(type)), 
       type="l", lwd=4, ylab="Total count", xlab=xlab, col=palette[1], ylim=ylim)
  lines(features[,2] ~ as.numeric(levels(type)), lty=2, lwd=4, col=palette[3])
  legend("topleft", legend=c("Factors and processes", "Pathways"), 
         lty=c(1,2), lwd=2.5, bty="n", col=palette2[3])
  
}

plot_changes <- function(x, y, xlab, bw=10){
  graphics::plot(y ~ as.numeric(levels(x)), 
                 ylab="Change in DAG dissimilarity", ylim=c(0,1),
                 xlab=xlab, pch=19, cex=2, col=paste(palette[1], "66", sep=""), axes=F)
  axis(1); axis(2)
  newy <- as.numeric(y)[!is.na(y)]
  newx <- as.numeric(levels(x))[!is.na(y)]
  lines(ksmooth(newx, newy, kernel="normal", bandwidth=bw), col=palette[1], lwd=5)
  
}

# Read in data  ----------------------------------------------------------------

input.data <- readxl::read_xlsx("./data/screening_data.xlsx", sheet=2)

# Subset to only the first 30% random sample
input.data <- input.data[1:236,]

# Pull aside the articles that were excluded
excluded <- input.data[!(input.data$include%in%"yes"),1:40]

# Retain the included articles for the rest of the script
input.data <- data.frame(input.data[input.data$include%in%"yes",1:40])


# Reasons for exclusion -------------------------------------------------------
table(excluded$reasons)

reasons_excl <- excluded$reasons
reasons_excl[is.na(reasons_excl)] <- excluded$include[is.na(reasons_excl)]
sort(table(reasons_excl))

# Recode the reasons for exclusion
wrong_hab <- c("farm fields", "not forest", "matrix surrounding")
notNA <- c("North America")
nofull <- c("needPDF")
wrong_pop <- c("focal species", "specialist")
wrong_time <- c("overwinter", "season", "transient", "post-fledg")
wrong_comp <- c("comparing", "time for space")
existing_data <-
  c(
    "already included",
    "data comes",
    "included as",
    "existing data",
    "review",
    "version included",
    "any data"
  )
no_sensitivity <-
  c(
    "no measure of edge or area",
    "no measure of distance",
    "no distance",
    "no mechanisms",
    "patch size",
    "sensitivity",
    "gradient",
    "harvest",
    "landscape level",
    "chose continuous",
    "not distance"
  )

reasons <-
  list(
    no_sensitivity = no_sensitivity,
    existing_data = existing_data,
    wrong_comp = wrong_comp,
    wrong_pop = wrong_pop,
    notNA = notNA,
    nofull = nofull,
    wrong_time = wrong_time,
    wrong_hab = wrong_hab
  )

reasons_excl <- recode_tags(reasons, reasons_excl)
table(reasons_excl)


# Publication year  ------------------------------------------------------------

input.data$firstyear <- as.numeric(unlist(lapply(strsplit(input.data$years, "-"), function(x){x[1]})))

possible_years <- factor(append(input.data$firstyear, seq(1975, 2020, 1)))[1:nrow(input.data)]
input.data$firstyear <- possible_years
input.data <- input.data[order(input.data$firstyear),]

hist(
  as.numeric(as.character(input.data$firstyear)),
  20,
  border = F,
  col = palette[1],
  xlim = c(1965, 2020),
  xlab = "First year of data collection",
  main = ""
)

# Latitude and longitude  ------------------------------------------------------ 

input.data$latitude <- unlist(lapply(strsplit(input.data$coordinates, ", "), function(x){x[1]}))
input.data$longitude <- unlist(lapply(strsplit(input.data$coordinates, ", "), function(x){x[2]}))
input.data$longitude[which(input.data$longitude==106.666667)] <- -106.666667

input.data$latband <- as.numeric(unlist(lapply(strsplit(input.data$latitude, "\\."), function(x){x[1]})))
latbands <- factor(append(input.data$latband, seq(24, 72, 1)))[1:nrow(input.data)]
input.data$latband <- latbands

input.data$longband <- as.numeric(unlist(lapply(strsplit(input.data$longitude, "\\."), function(x){x[1]})))
longbands <- factor(append(input.data$longband, seq(-140, -52, 1)))[1:nrow(input.data)]
input.data$longband <- longbands

forests <- raster::raster("./data/NA_TREEAGE_1096/data/ca04_usak06_1km.tif")
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



# png("./figures/study_locations.png", width=12,height=8,units="in",res=1200)
# par(xpd=T)
# sp::plot(forested, col=palette[6], axes=F, box=F, legend=F)
# sp::plot(canada, add=T, border="grey50")
# sp::plot(usa, add=T, border="grey50")
#
# points(studysites, col="white", pch=20, lwd=1, cex=1.5)
# points(studysites, col=paste(palette[1], "99", sep=""), pch=20, cex=1.5)
# points(studysites, col="grey20", pch=1, cex=1.25, lwd=0.5)
# dev.off()

# Matrix type ------------------------------------------------------------------

agriculture <- c("agric", "pasture", "crop", "field", "farm")
industry <-
  c("clearcut",
    "logg",
    "mining",
    "energy",
    "shale",
    "pipeline",
    "seismic")
developed <- c("urban", "develop", "road", "residential")
other_matrix <-
  c("plantation",
    "bog",
    "forest",
    "grass",
    "sapling",
    "reservoir",
    "natural",
    "abrupt")

matrix_types <- list(developed=developed, 
                     agriculture=agriculture, 
                     industry=industry, 
                     other=other_matrix)
matrix_types <- lapply(matrix_types, function(x){
  paste(x, "*", sep="")
})

matrix_hits <- topictagger::tag_strictly(input.data$matrix, matrix_types)
matrix_hits[matrix_hits>0] <- 1
prop.table(colSums(matrix_hits))

# Forest type ------------------------------------------------------------------

forest_types <- input.data$forest_class

prop.table(sort(table(forest_types)))

# Bird species -----------------------------------------------------------------

input.data$taxa[grep("Swainson's warbler", input.data$taxa)] <- "SWWA"

input.data$taxa <- gsub("black-throated blue warbler", "btbw", input.data$taxa)
input.data$taxa[input.data$taxa=="buff-breasted flycatchers"] <- "buff-breasted flycatcher"

birdcodes <- read.csv("./data/IBP-Alpha-Codes20.csv")
birdcodes <- data.frame(apply(birdcodes, 2, tolower))
fourlettercodes <- birdcodes[,'SPEC']

bird_dictionary <- topictagger::create_dictionary(birdcodes[,c(2,4,5)])
tagged_birds <- topictagger::tag_strictly(tolower(input.data$taxa), 
                                          bird_dictionary)

bird_dat <- array(dim=c(nrow(tagged_birds), length(fourlettercodes)))

species <- unlist(lapply(colnames(tagged_birds), function(x){
  strsplit(x, "\\.")[[1]][1]
} ))

for(i in 1:length(fourlettercodes)){
  bird_dat[,i] <- rowSums(tagged_birds[,species==fourlettercodes[i]])
}

# Check that no species were missed
input.data$taxa[rowSums(bird_dat)==0]

studied <- fourlettercodes[which(colSums(bird_dat)>0)]
scinames <- birdcodes$SCINAME[birdcodes$SPEC %in% studied]

# Need to match the scientific and common names to BirdLife to map to phylogeny
birdlife <- readxl::read_xls("./data/BirdLife_Checklist_Version_3/BirdLife_Checklist_Version_3.xls")

commonnames <- birdcodes$COMMONNAME[birdcodes$SPEC %in% studied]

# Read in list of forest birds in North America
lott <- readxl::read_xlsx("./data/lott.xlsx", sheet=2)
lott <- data.frame(lott)[-c(1:40),]

# Save the current common names
original_common <- unique(append(birdcodes$COMMONNAME[birdcodes$SPEC %in% studied], tolower(lott$Common.Name)))[!is.na(unique(append(birdcodes$COMMONNAME[birdcodes$SPEC %in% studied], tolower(lott$Common.Name))))]

commonnames <- original_common
commonnames[commonnames=="gray-headed chickadee"] <- NA
commonnames[is.na(match(commonnames, tolower(birdlife$`Common name`)))]
commonnames[commonnames=="blue-gray gnatcatcher"] <- "blue-grey gnatcatcher"
commonnames[commonnames=="gray catbird"] <- "grey catbird"
commonnames[commonnames=="european starling"] <- "common starling"
commonnames[commonnames=="canada jay"] <- "grey jay"
commonnames[commonnames=="gray-cheeked thrush"] <- "grey-cheeked thrush"
commonnames[commonnames=="brown creeper"] <- "american treecreeper"
commonnames[commonnames=="oregon junco"] <- "dark-eyed junco"

scientificnames <- birdlife$`Scientific name`[match(commonnames, tolower(birdlife$`Common name`))]

# Read in phylogenetic trees
phylo <- ape::read.nexus("./data/tree-pruner-0eecac32-f757-4052-91fc-0cef9612d69f/output.nex")

# Are all the species in the phylogeny in our list of forest birds?
gsub("_", " ", (phylo[[1]]$tip.label)) %in% scientificnames

# Are all the birds in our list of forest birds in the phylogeny?
scientificnames[!(scientificnames %in% gsub("_", " ", (phylo[[1]]$tip.label)))]
# No, but these are non-passerines

bird_lookup <- cbind(original_common, commonnames, scientificnames, birdcodes$SPEC[match(original_common, birdcodes$COMMONNAME)])

# Cut down the bird data to only the forest birds
bird_dat2 <- bird_dat[,which(fourlettercodes %in% birdcodes$SPEC[match(original_common, birdcodes$COMMONNAME)])]

included_species <- fourlettercodes[which(fourlettercodes %in% birdcodes$SPEC[match(original_common, birdcodes$COMMONNAME)])]

bird_lookup <- cbind(bird_lookup, colSums(bird_dat2)[match(bird_lookup[,4], included_species)])


tree <- phylo[[1]]
tree$tip.label<-gsub("_"," ",tree$tip.label)

dat <- data.frame(count=as.numeric(bird_lookup[match(gsub("_", " ", phylo[[1]]$tip.label), bird_lookup[,3]),5]))
dat2 <- data.frame(count=as.numeric(bird_lookup[match(gsub("_", " ", phylo[[1]]$tip.label), bird_lookup[,3]),5]),
                   species=(bird_lookup[match(gsub("_", " ", phylo[[1]]$tip.label), bird_lookup[,3]),1]))

dat2[order(dat2[,1]),]

phylodat <- phylobase::phylo4d(tree, dat)

barplot(
  phylodat,
  bar.col = palette[1],
  center = F,
  scale = F,
  tree.type = "phylogram",
  show.trait = T,
  tip.cex = 1,
  show.data.axis = T,
  trait.labels = c("Number of studies"),
  bar.lwd = 5,
  trait.bg.col = "white",
  data.xlim = c(0, 30),
  grid.vertical = F
)






# Path cleaning  ---------------------------------------------------------------

input.data$Path <- rep("", nrow(input.data))
input.data$test.pathways <- gsub("=", ">", input.data$test.pathways)

for(i in 1:nrow(input.data)){
  path1 <- trimws(strsplit(input.data$test.pathways[i], ";")[[1]])
  
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
  input.data$Path[i] <- paste(paths, ";", sep="")
}

input.data$Path <- trimws(input.data$Path)


tmp.path <- gsub(">", "_>_",input.data$Path)
tmp.path <- gsub(";", "_;_", tmp.path)

processes <- strsplit(trimws(tmp.path), "_")


# Read in list of synonyms and terms to group together
synonyms <- read.csv("./data/terms_grouped.csv")

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

input.data$Path <- unlist(processes)

# Create graphs for studies and years ------------------------------------------

# Graphs for all studies individually by title (index = 21)
study_DAGS <- create_DAGS(input.data, 21)

# Graphs for each year, where all studies from a given year at input at once
yearly_DAGS <- create_DAGS(input.data, which(colnames(input.data)=="firstyear"))

# Graphs for each latitude band, where all studies from a given latitude band
# are input at once
lat_DAGS <- create_DAGS(input.data, which(colnames(input.data)=="latband"))

# Graphs for each longitude band, where all studies from a given longitude band
# are input at once
long_DAGS <- create_DAGS(input.data, which(colnames(input.data)=="longband"))

# Cumulative and sliding window graphs -----------------------------------------

cumulative_DAGs <- create_cumulative(yearly_DAGS, 
                                     min(valid_levels(input.data$firstyear))+1)

yearwindow_DAG <- create_windows(yearly_DAGS, windowsize = 5, 
                                 startpoint = min(valid_levels(input.data$firstyear))+1)

latband_DAG <- create_windows(lat_DAGS, windowsize = 5, 
                              startpoint = min(valid_levels(input.data$latband))+1)

longband_DAG <- create_windows(long_DAGS, windowsize = 5, 
                              startpoint = min(valid_levels(input.data$longband)+1))

# Stopping criteria and estimationg of missing information ---------------------

n.sim <- 300
factors <- pathways <- modfact <- modpath <- array(dim=c(n.sim, length(study_DAGS)))

for(i in 1:n.sim){
  # Randomly reorder graphs to create a cumulative network
  tmpdag <- create_cumulative(sample(study_DAGS))
  
  # New graph features added
  factors[i,] <- as.numeric(new_features(tmpdag, 2)[,1]>0)
  pathways[i,] <- as.numeric(new_features(tmpdag, 2)[,2]>0)
  
  # Fit a binomial GLM and predict the probability of new graph features
  modfact[i,2:122] <- plogis(predict(glm(factors[i,] ~ seq(1, length(study_DAGS), 1), 
                                         family = "binomial")))
  modpath[i,2:122] <- plogis(predict(glm(pathways[i,] ~ seq(1, length(study_DAGS), 1), 
                                         family = "binomial")))
  }

plot(apply(modfact, 2, mean), type="n", ylim=c(0,1), xlab="Number of included studies", 
     ylab="Probability of new graph features", axes=F)
axis(1); axis(2)
ci <- apply(modfact, 2, quantile, c(0.025, 0.975), na.rm=T)
newx <- seq(1, length(study_DAGS), 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])), 
        col=paste(palette[1], "66", sep=""), border = F)
lines(apply(modfact, 2, mean))

ci <- apply(modpath, 2, quantile, c(0.025, 0.975), na.rm=T)
newx <- seq(1, length(study_DAGS), 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])), 
        col=paste(palette[3], "66", sep=""), border = F)
lines(apply(modpath, 2, mean))

legend("topright", 
       legend=c("Mean", "95% confidence interval", "Factors", "Pathways"),
       lwd=c(2, NA, NA, NA, NA), pch=c(NA, 15, 15, 15), 
       col=c("black", "grey95", palette[c(1,3)]), bty="n", 
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

plot(apply(est.factors$jack1, 1, mean), type="n", ylim=c(0,600), axes=F,
     ylab="Extrapolated total count of graph features", 
     xlab="Number of included studies")
axis(1); axis(2)
ci <- apply(est.factors$jack1, 1, quantile, c(0.025, 0.975))
newx <- seq(1, length(study_DAGS)-2, 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])), 
        col=paste(palette[1], "66", sep=""), border = F)
lines(apply(est.factors$jack1, 1, mean))

ci <- apply(est.pathways$jack1, 1, quantile, c(0.025, 0.975))
newx <- seq(1, length(study_DAGS)-2, 1)
polygon(x=c(newx, rev(newx)), y=c(ci[1,], rev(ci[2,])), 
        col=paste(palette[3], "66", sep=""), border = F)
lines(apply(est.pathways$jack1, 1, mean))
legend("topright", 
       legend=c("Mean", "95% confidence interval", "Factors", "Pathways"),
       lwd=c(2, NA, NA, NA, NA), pch=c(NA, 15, 15, 15), 
       col=c("black", "grey95", palette[c(1,3)]), bty="n", 
       pt.cex = c(NA, 2,2,2))

# Global graph from union of all studies ---------------------------------------

metaDAG <- cumulative_DAGs[[length(cumulative_DAGs)]]

# Plot graph and rearrange layout
tkplot(metaDAG, 
            edge.arrow.size=.7, 
            edge.color=palette2[1],
            vertex.color="#ffffff00", 
            vertex.frame.color="#ffffff00", 
            vertex.label.family="Arial",
            vertex.label.dist=0, 
            vertex.label.color="black", 
            vertex.label.cex=.75,vertex.size=sqrt(strength(metaDAG))*3,
            edge.width=as.numeric(lookup_widths(study_DAGS, metaDAG)),
            )


# Convert graph to adjacency matrix and perform a transitive reduction
metaDAG.adj <- as.matrix(as_adj(metaDAG))
reduced.adj <- mnem::transitive.reduction(metaDAG.adj)
reducedDAG <- igraph::graph_from_adjacency_matrix(reduced.adj)

# Note: graph is not acyclic, so transitive reduction has multiple solutions
igraph::is.dag(metaDAG)


edgedat <- cbind(attr(igraph::E(metaDAG), "vnames"), 
                 as.numeric(lookup_widths(study_DAGS, metaDAG)))
edgedat[order(as.numeric(edgedat[,2])),]

# Plot network features and changes --------------------------------------------

plot_features(calc_features(cumulative_DAGs, 5), 
              input.data$firstyear, xlab="Year")
plot_features(calc_features(yearwindow_DAG, 5), 
              input.data$firstyear, xlab="Year")

yearly_changes <- DAG_changes(cumulative_DAGs, 
                              input.data$firstyear, 1)
plot_changes(input.data$firstyear,
             yearly_changes, xlab="First year of data collection")

fiveyr.changes <- DAG_changes(yearwindow_DAG, 
                              input.data$firstyear, 6)
plot_changes(input.data$firstyear, 
             fiveyr.changes, xlab="First year of data collection")

plot_features(calc_features(latband_DAG, 10), 
              input.data$latband, xlab="Latitude")
plot_features(calc_features(longband_DAG, 19), 
              input.data$longband, xlab="Longitude")

latchanges_window <- DAG_changes(latband_DAG, input.data$latband, 8)
plot_changes(input.data$latband, latchanges_window, xlab="Latitude", bw=10)

longchanges_window <- DAG_changes(longband_DAG, input.data$longband, 17)
plot_changes(input.data$longband, longchanges_window, xlab="Longitude", bw=10)

# Subset graphs by study characteristics ---------------------------------------

agric <- study_DAGS[which(matrix_hits[,'agriculture']>0)]
agric2 <- list(agric[[1]])

for(i in 2:length(agric)){
  agric2[i] <- list(igraph::union(agric2[[i-1]], agric[[i]]))
}

tkplot(
  agric2[[i]],
  main = "Agriculture",
  edge.arrow.size = 0.5,
  edge.color = palette2[1],
  vertex.color = "white",
  vertex.frame.color = "white",
  vertex.label.family = "Arial",
  vertex.label.dist = 0,
  vertex.label.color = "black",
  vertex.label.cex = .7,
  vertex.size = sqrt(strength(agric2[[i]])) * 3,
  edge.width = (as.numeric(lookup_widths(study_DAGS, agric2[[i]]))) / 2
)

develop <- study_DAGS[which(matrix_hits[,'developed']>0)]
develop2 <- list(develop[[1]])

for(i in 2:length(develop)){
  develop2[i] <- list(igraph::union(develop2[[i-1]], develop[[i]]))
}

tkplot(
  develop2[[i]],
  main = "Development",
  edge.arrow.size = .5,
  edge.color = palette2[1],
  vertex.color = "white",
  vertex.frame.color = "white",
  vertex.label.family = "Arial",
  vertex.label.dist = 0,
  vertex.label.color = "black",
  vertex.label.cex = .7,
  vertex.size = sqrt(strength(develop2[[i]])) * 3,
  edge.width = (as.numeric(lookup_widths(study_DAGS, develop2[[i]]))) / 2
)


industry <- study_DAGS[which(matrix_hits[,'industry']>0)]
industry2 <- list(industry[[1]])

for(i in 2:length(industry)){
  industry2[i] <- list(igraph::union(industry2[[i-1]], industry[[i]]))
}

tkplot(
  industry2[[i]],
  main = "Industry",
  edge.arrow.size = .5,
  edge.color = palette2[1],
  vertex.color = "white",
  vertex.frame.color = "white",
  vertex.label.family = "Arial",
  vertex.label.dist = 0,
  vertex.label.color = "black",
  vertex.label.cex = .7,
  vertex.size = sqrt(strength(industry2[[i]])) * 3,
  edge.width = as.numeric(lookup_widths(study_DAGS, industry2[[i]])) / 2
)

# Example graph for a single species -------------------------------------------

birds <- which(bird_dat[,which(fourlettercodes=="oven")]>0)

bird <- study_DAGS[birds]
bird2 <- list(bird[[1]])

for(i in 2:length(bird)){
  bird2[i] <- list(igraph::union(bird2[[i-1]], bird[[i]]))
}
width_lookup <- table(unlist(lapply(bird, function(x){
  attr(igraph::E(x), "vnames")
})))

widths <- width_lookup[match(attr(igraph::E(bird2[[i]]), "vnames"), names(width_lookup))]

tkplot(
  bird2[[i]],
  main = "Ovenbird",
  edge.arrow.size = 0.5,
  edge.color = palette[1],
  vertex.color = "white",
  vertex.frame.color = "white",
  vertex.label.family = "Arial",
  vertex.label.dist = 0,
  vertex.label.color = "black",
  vertex.size = sqrt(strength(bird2[[i]])) * 3,
  edge.width = as.numeric(lookup_widths(study_DAGS, bird2[[i]])) / 2
)
