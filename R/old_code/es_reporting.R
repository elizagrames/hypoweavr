pied <- read.csv("~/Downloads/insect_MASEM_effects - 636 (1).csv")

# Function to report effect size data for continuous measures
es_cor <- function(xvar, yvar){
  retained <- (!is.na(xvar) & !is.na(yvar))
  N <- sum(retained)
  tmp <- cor.test(xvar[retained], yvar[retained])
  round(cbind(df=tmp$parameter, t=tmp$statistic, p=tmp$p.value, r=tmp$estimate, N=N), 5)
}

# Given a list of continuous, numeric variables, reports effects in a table
report_effects <- function(variables){
  lookups <- array(dim=c(length(variables), length(variables)))
  rownames(lookups) <- colnames(lookups) <- names(variables)
  lookups[upper.tri(lookups)] <- seq(1, sum(upper.tri(lookups)), )
  
  cor_dat <- array(dim=c(sum(upper.tri(lookups)), 5))
  
  for(i in 1:sum(upper.tri(lookups))){
    x <- names(sort(rowSums(lookups==i, na.rm = T), decreasing = T))[1]
    y <- names(sort(colSums(lookups==i, na.rm = T), decreasing = T))[1]
    yvar <- variables[names(variables)==y][[1]]
    xvar <- variables[names(variables)==x][[1]]
    cor_dat[i,] <- es_cor(xvar, yvar)
  }
  colnames(cor_dat) <- c("df", "t", "p", "r", "N")
  return(cor_dat)
}

# List of named variables
# All variables need to be continuous and numeric
# Missing and NA values are fine
variables <- list(chicks=pied$Maximum.number.of.hatched.young.averaged.over.first.clutches, 
                  fledglings=pied$Average.number.of.young.fledged.from.first.clutches, 
                  nestsuccess=(pied$N-pied$Number.of.failed.first.clutches)/pied$N, 
                  moths=pied$Density.of.winter.moths.in.Wytham.Wood)

# Print out a table that can be passed to metafor
report_effects(variables)