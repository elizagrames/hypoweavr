#' Make machine-readable ontologies
#' @description Given a human-readable ontology or hierarchical dictionary, fills in blank rows
#' @param x a data.frame
#' @return the same data.frame with blank rows filled in
fill_rows <- function (x){
  # Note: this function is pulled from topictagger
  if (any(x == "")) {
    x[x == ""] <- NA
  }
  for (j in ncol(x):2) {
    for (i in nrow(x):1) {
      if (!is.na(x[i, j])) {
        if (is.na(x[i, (j - 1)])) {
          parents <- x[1:i, (j - 1)]
          parent <- utils::tail(parents[!is.na(parents)],
                                1)
          x[i, (j - 1)] <- parent
        }
      }
    }
  }
  return(x)
}

#' Tidy up punctuation
#' @description Standardizes punctuation in a string
#' @param x a string
#' @param punc a list of punctuation marks
#' @param default a punctuation mark that should be standard
#' @return the input string with punctuation marks replaced by the default
#' @examples clean_punctuation("a;b/c", c("/", ";"), ";")
clean_punctuation <- function(x, punc, default){
  punc <- unique(append(default, punc))
  if(length(punc)>1){
    for(i in 1:length(punc)){
      x <- gsub(punc[i], default, x)
    }
  }
  return(x)
}

#' Standardizes pathways stored as text
#' @description Given pathways stored as text, cleans and preps them to be converted to a network
#' @param path a character vector describing relationships between concepts
#' @param join a character vector that indicates two concepts are both affected similarly
#' @param cause a character vector that indicates implied causal relationships
#' @param sep a character vector that separates different pathways in the same text
#' @return a character vector the same length as the input path with standardized format
clean_path <- function(path, join="+", cause=c(">", "="), sep=c(";")){
  # use standard symbols for distinguishing paths
  subfunc <- function(z){
    z <- clean_punctuation(z, cause, ">")
    z <- clean_punctuation(z, sep, ";")
    z <- clean_punctuation(z, join, "+")

    pathways <- trimws(strsplit(z, ";")[[1]])
    groups <- lapply(strsplit(pathways, ">"), trimws)

    # collapse joined terms with causal terms

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

    output1 <- trimws(paste(paths, ";", sep=""))

    return(output1)
  }

  output <- unlist(lapply(path, subfunc))

  return(output)
}


#' Merges synonymous terms
#' @description Replace synonyms with standard terms in a pathway created by clean_path()
#' @param path a character vector describing relationships between concepts
#' @param terms a character vector listing replacement terms
#' @param synonyms a character vector listing synonyms to be replaced
#' @return a character vector the same length as the input path with synonyms replaced
replace_terms <- function(path, terms, synonyms){

  subfunc <- function(x){
    tmp.path <- gsub(">", "_>_", x)
    tmp.path <- gsub(";", "_;_", tmp.path)
    processes <- strsplit(trimws(tmp.path), "_")

    for(i in 1:length(processes)){
      y <- trimws(processes[[i]])
      for(j in 1:length(y)){
        if(y[j]%in%synonyms){
          newterm <- terms[match(y[j], synonyms)]
          if(newterm!=""){
            y[j] <- newterm
          }
        }
      }
      processes[[i]] <- paste(y, collapse=" ")
    }
    output1 <- unlist(processes)

    return(output1)
  }

  output <- unlist(lapply(path, subfunc))
  return(output)

}

#' Create a network from a path description
#' @description Given an implied causal pathway described in text, creates a network
#' @param paths a character vector describing relationships between concepts
#' @return a list of graphs the same length as input paths
generate_graph <- function(paths){

  subfunc <- function(x){

    blank <- as.data.frame(matrix(c(NA, NA), ncol=2, byrow=TRUE))
    colnames(blank) <- c("from", "to")
    study.graph <- suppressWarnings(igraph::graph_from_data_frame(blank, directed=TRUE))
    count.paths <- c()

    study <- unlist(lapply(strsplit(as.character(x), ";"), trimws))

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
      path.edges <- suppressWarnings(as.data.frame(matrix(c(line), ncol=2,byrow=TRUE)))
      colnames(path.edges) <- c("from", "to")
      path.nodes <- unique(line)
      path.DAG <- suppressWarnings(igraph::graph_from_data_frame(path.edges, vertices=path.nodes, directed=TRUE))
      study.graph <- igraph::union(study.graph, path.DAG)
    }

    study.graph <- igraph::delete.vertices(study.graph, "NA")

    return(study.graph)

  }
  all.graphs <- lapply(paths, subfunc)
  return(all.graphs)
}

#' Unify graphs by shared characteristics
#' @description Takes the union of graphs
#' @param graphs a list of graph objects
#' @param by a vector indicating which graphs should be merged
#' @return a graph
merge_graphs <- function(graphs, by){
  level.graphs <- graphs[by]
  for(j in 1:length(level.graphs)){
    if(j==1){
      new.graphs <- level.graphs[[j]]
    }else{
      new.graphs <- igraph::union(new.graphs, level.graphs[[j]])
    }
  }

  return(new.graphs)
}

#' Arrange graphs into a series
#' @description Order graphs by metadata
#' @param graphs a list of graph objects
#' @param order.by a character vector of metadata by which to order and group graphs
#' @param sort.order a vector indicating the order in which to sort order.by
#' @return a list of graph objects
create_series <- function(graphs, order.by, sort.order=NULL){
  if(is.null(sort.order)){
    by.list <- sort(unique(order.by))
  }else{
    by.list <- sort.order
  }
  all.graphs <- list()
  length(all.graphs) <- length(by.list)

  for(i in 1:length(by.list)){
    if(by.list[i]%in%order.by){
      tmp.graph <- merge_graphs(graphs, order.by==by.list[i])
      all.graphs[[i]] <- tmp.graph
    }
  }
  names(all.graphs) <- by.list
  return(all.graphs)
}

#' Create a cumulative series of graphs
#' @description Takes the union of all previous graphs in a series at each time step to create a cumulative series of graphs
#' @param graphs a list of graph objects
#' @param order.by a character vector of metadata by which to order and group graphs
#' @param sort.order a vector indicating the order in which to sort order.by
#' @return a list of graph objects
create_cumulative <- function(graphs, order.by, sort.order=NULL){
  if(is.null(sort.order)){
    by.list <- sort(unique(order.by))
  }else{
    by.list <- sort.order
  }

  used <- which(by.list %in% order.by)



  all.graphs <- list()
  length(all.graphs) <- length(by.list)
  all.graphs[[used[1]]] <- merge_graphs(graphs, order.by==by.list[used[1]])

  for(i in (used[1]+1):length(by.list)){
    if(by.list[i]%in%order.by){
      tmp.graph <- merge_graphs(graphs, order.by==by.list[i])
      all.graphs[[i]] <- igraph::union(all.graphs[[i-1]], tmp.graph)
    }else{
      all.graphs[[i]] <- all.graphs[[i-1]]
    }
  }
  names(all.graphs) <- by.list
  return(all.graphs)
}


#' Create a series of graphs based on sliding windows
#' @description Takes the union of previous graphs within a sliding window to create a series of graphs
#' @param graphs a list of graph objects
#' @param window.size an integer indicating how large windows should be
#' @param startpoint where in the list to start the first window
#' @return a list of graph objects
create_windows <- function(graphs, window.size=5, startpoint=1){
  window_DAG <- list()
  length(window_DAG) <- length(graphs)
  for(i in max(startpoint, window.size):length(graphs)){

    tmp <- list(graphs[[i]])

    if(!is.null(tmp[[1]])){
      for(j in 2:(window.size-1)){

        tmp[[j]] <- igraph::union(tmp[[j-1]], graphs[[i-j]])

      }

      window_DAG[i] <- list(tmp[[length(tmp)]])
    }

  }

  dagnames <- paste(round(as.numeric(names(graphs))-5), round(as.numeric(names(graphs))), sep="-")

  dagnames[1:max(startpoint, window.size)] <- NA

  names(window_DAG) <- dagnames

  return(window_DAG)

}

#' Calculate graph dissimilarity
#' @description Calculates dissimilarity between two graphs based on topological distance
#' @param x a graph object or list of graph objects
#' @param y a graph object (if x is a graph object)
#' @return the dissimilarity metric from Scheiber et al. 2017
calc_dissimilarity <- function(x, y=NULL){
  if(!is.null(y)){
    x <- list(x, y)
  }
  changes <- rep(NA, length(x))
  used <- which(!unlist(lapply(x, is.null)))
  for (i in used[2]:length(x)) {
    if(i %in% used){
      previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
      changes[i] <- gDist(x[[i]], x[[previous_valid]], 0.45, 0.45, 1)
    }
  }
  if(length(x)==2){
    changes <- changes[2]
  }

  return(changes)
}







#' Graph distance
#' @description See Schieber et al. (2017) for details
#' @param g a graph object to compare to h
#' @param h a graph object to compare to g
#' @param w1 weight of term1
#' @param w2 weight of term2
#' @param w3 weight of term3
#' @return graph distance
gDist <-function(g,h,w1,w2,w3){

  # Network dissimilarity functions from Schieber et al. (2017)
  # https://github.com/tischieber/Quantifying-Network-Structural-Dissimilarities

  # Shannon entropy
  entropia<-function(a) {
    a<-a[which(a>0)];
    -sum(a*log(a));
  }

  # Node distance
  node_distance<-function(g){
    n<-length(igraph::V(g))
    if(n==1){
      retorno=1
    }

    if(n>1){
      a<-Matrix::Matrix(0,nrow=n,ncol=n,sparse=TRUE)
      m<-igraph::shortest.paths(g,algorithm=c("unweighted"))
      m[which(m=="Inf")]<-n
      quem<-setdiff(intersect(m,m),0)

      for(j in (1:length(quem))){
        l<-which(m==quem[j])/n
        linhas<-floor(l)+1
        posicoesm1<-which(l==floor(l))
        if(length(posicoesm1)>0){
          linhas[posicoesm1]<-linhas[posicoesm1]-1
        }
        a[1:n,quem[j]]<-graphics::hist(linhas,plot=FALSE,breaks=(0:n))$counts
      }
      #m<-c()
      retorno=(a/(n-1))
      retorno <- as.matrix(retorno) # was returning a sparse matrix
    }
    return(retorno)
  }

  # Network node dispersion
  nnd<-function(g){
    N<-length(igraph::V(g))
    nd<-node_distance(g)
    pdfm<-colMeans(nd)
    norm<-log(max(c(2,length(which(pdfm[1:(N-1)]>0))+1)))
    return(c(pdfm,max(c(0,entropia(pdfm)-entropia(nd)/N))/norm))
  }

  # Alpha centrality
  alpha<-function(g){
    N<-length(igraph::V(g))
    r<-sort(igraph::alpha.centrality(g,exo=igraph::degree(g)/(N-1),alpha=1/N))/((N^2))
    return(c(r,max(c(0,1-sum(r)))))
  }




  first<-0
  second<-0
  third<-0
  g<-g
  h<-h
  N<-length(igraph::V(g))
  M<-length(igraph::V(h))
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
    g<-igraph::graph.complementer(g)
    h<-igraph::graph.complementer(h)

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

#' Calculate the number and identity of new features in a graph
#' @description Returns the number and/or identity of new nodes and edges between two graphs
#' @param x a graph object or list of graph objects
#' @param y a graph object (if x is a graph object)
#' @param return a string indicating if 'counts' or 'features' should be returned
#' @return either a list of numbers or list of feature names
new_features <- function(x, y=NULL, return="counts"){
  if(!is.null(y)){
    x <- list(x, y)
  }

  used <- which(!unlist(lapply(x, is.null)))

  if(return=="counts"){
    n.nodes <- n.edges <- rep(NA, length(x))
    for (i in used[2]:length(x)) {
      if(i %in% used){
        previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
        n.nodes[i] <- sum(!(igraph::V(x[[i]]) %in% igraph::V(x[[previous_valid]])))
        n.edges[i] <- sum(!(igraph::E(x[[i]]) %in% igraph::E(x[[previous_valid]])))
      }
    }
    if(length(x)==2){
      n.nodes <- n.nodes[2]
      n.edges <- n.edges[2]
    }
    output <- list(nodes=n.nodes, edges=n.edges)
  }


  if(return=="features"){
    n.nodes <- n.edges <- list()
    length(n.nodes) <- length(n.edges) <- length(x)

    for (i in used[2]:length(x)) {
      if(i %in% used){
        previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
        n.nodes[[i]] <- names(igraph::V(x[[i]])[!(igraph::V(x[[i]]) %in% igraph::V(x[[previous_valid]]))])
        n.edges[[i]] <- attr(igraph::E(x[[i]]), "vnames")[!(igraph::E(x[[i]]) %in% igraph::E(x[[previous_valid]]))]
      }
    }
    if(length(x)==2){
      n.nodes <- n.nodes[[2]]
      n.edges <- n.edges[[2]]
    }
    output <- list(nodes=n.nodes, edges=n.edges)
  }

  return(output)
}

#' Calculate the number and identity of shared features between graphs
#' @description Returns the number and/or identity of shared nodes and edges between two graphs
#' @param x a graph object or list of graph objects
#' @param y a graph object (if x is a graph object)
#' @param return a string indicating if 'counts' or 'features' should be returned
#' @return either a list of numbers or list of feature names
shared_features <- function(x, y=NULL, return="counts"){
  if(!is.null(y)){
    x <- list(x, y)
  }

  used <- which(!unlist(lapply(x, is.null)))

  if(return=="counts"){
    n.nodes <- n.edges <- rep(NA, length(x))

    for (i in used[2]:length(x)) {
      if(i %in% used){
        previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
        n.nodes[i] <- length(igraph::V(igraph::intersection(x[[i]], x[[previous_valid]])))
        n.edges[i] <- length(igraph::E(igraph::intersection(x[[i]], x[[previous_valid]])))
      }
    }
    if(length(x)==2){
      n.nodes <- n.nodes[2]
      n.edges <- n.edges[2]
    }
    output <- list(nodes=n.nodes, edges=n.edges)
  }

  if(return=="features"){
    n.nodes <- n.edges <- list()
    length(n.nodes) <- length(n.edges) <- length(x)

    for (i in used[2]:length(x)) {
      if(i %in% used){
        previous_valid <- max((1:(i-1))[1:(i-1) %in% used])
        n.nodes[[i]] <- names(igraph::V(igraph::intersection(x[[i]], x[[previous_valid]])))
        n.edges[[i]] <- attr(igraph::E(igraph::intersection(x[[i]], x[[previous_valid]])), "vnames")
      }
    }
    if(length(x)==2){
      n.nodes <- n.nodes[[2]]
      n.edges <- n.edges[[2]]
    }
    output <- list(nodes=n.nodes, edges=n.edges)
  }
  return(output)

}

#' Calculate node and edge metrics
#' @description Given a graph, calculates node and edge metrics using functions from the igraph package
#' @param x a graph object or list of graphs
#' @param metric a string indicating which metric to use; see details for options
#' @param return.df logical; should a list of metrics be returned or a unified data frame
#' @details optional metrics: page_rank, edge_betweenness, strength, alpha_centrality, FILL IN MORE OF THESE LATER. See igraph help pages for details on each function.
graph_metrics <- function(x, metric="page_rank", return.df=FALSE){
  subfunc <- function(z){
    if(class(z)=="igraph"){
      if(metric%in%c("page_rank", "eigen_centrality", "page.rank")){
        tmp <- eval(as.name(metric))(z)$vector
      }else{
        tmp <- eval(as.name(metric))(z)
      }

      if(metric%in%c("edge.betweenness", "edge_betweenness")){
        names(tmp) <- attr(igraph::E(z), "vnames")
      }
      return(tmp)

    }else{
      NA
    }

  }

  if(class(x)!="igraph"){
    output <- lapply(x, subfunc)
    names(output) <- names(x)
    if(return.df){
      full_set <- unique(unlist(lapply(output, names)))
      dat <- array(dim=c(length(full_set), length(output)))
      for(i in 1:length(output)){
        dat[,i] <- output[[i]][match(full_set, names(output[[i]]))]
      }
      rownames(dat) <- full_set
      colnames(dat) <- names(output)

      output <- dat

    }
  }else{
    output <- subfunc(x)
  }

  return(output)

}
