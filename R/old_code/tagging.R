# Definitions of search components ---------------------------------------------
time <-
  c(
    "reproductive success",
    "fledging success",
    "reproductive performance",
    "pairing success",
    "parasitism rates",
    "provisioning rates",
    "reduced fecundity",
    "singing males",
    "territorial males",
    "territory density",
    "territory occupancy",
    "ground nests",
    "young fledged",
    "breeding",
    "reproduction",
    "nesting",
    "fledgling"
  )
outcome <-
  c(
    "imperfect detection",
    "population density",
    "population persistence",
    "detection probability",
    "local extinction",
    "total number",
    "population status",
    "remaining populations",
    "adult mortality",
    "area-sensitive species",
    "population dynamics",
    "population growth",
    "population estimate",
    "breeding population",
    "total population",
    "demographic parameters",
    "density increased",
    "population sizes",
    "wildlife population",
    "higher density",
    "highest density",
    "effective population",
    "population trend",
    "lower density",
    "occupied patch",
    "small population",
    "population sinks",
    "population viability",
    "relative density",
    "population decline",
    "local population",
    "population change",
    "species increased",
    "species occurrence",
    "stable population",
    "species present",
    "population model",
    "survival rates",
    "abundance",
    "area sensitive",
    "area sensitivity",
    "occupancy",
    "point count",
    "survivorship",
    "species composition"
  )

population <-
  c(
    "forest interior species",
    "avian responses",
    "avian species",
    "dendroica cerulea",
    "dumetella carolinensis",
    "hylocichla mustelina",
    "oporornis formosus",
    "piranga olivacea",
    "seiurus aurocapilla",
    "seiurus aurocapillus",
    "sensitive species",
    "setophaga chrysoparia",
    "singing males",
    "wilsonia citrina",
    "young fledged",
    "passerine",
    "passeriformes",
    "bird",
    "warbler",
    "thrush",
    "vireo",
    "nuthatch",
    "kinglet",
    "flycatcher",
    "tanager",
    "avian community",
    "chickadee",
    'wren',
    "veery"
  )

intervention <-
  c(
    "landscape management",
    "vegetation structure",
    "island biogeography",
    "agricultural landscape",
    "continuous forest",
    "landscape metrics",
    "larger patches",
    "smaller patches",
    "forest patch",
    "habitat amount",
    "landscape composition",
    "total forest",
    "forest cover",
    "forest structure",
    "forested areas",
    "forest opening",
    "large forest",
    "habitat alteration",
    "habitat change",
    "habitat configuration",
    "habitat degradation",
    "habitat disturbance",
    "habitat heterogeneity",
    "habitat patch",
    "habitat quality",
    "habitat requirement",
    "heterogeneous landscape",
    "forest edges",
    "landscape attributes",
    "landscape characteristics",
    "landscape context",
    "landscape ecology",
    "landscape effects",
    "landscape features",
    "landscape level",
    "landscape pattern",
    "landscape scale",
    "landscape structure",
    "landscape variable",
    "large areas",
    "large patches",
    "large tracts",
    "natural disturbance",
    "occupied patch",
    "intact forest",
    "patch isolation",
    "forest stand",
    "forest tract",
    "patch occupancy",
    "patch scale",
    "quality habitat",
    "sensitive species",
    "small patch",
    "small woodlots",
    "spatial configuration",
    "surrounding landscape",
    "urbanizing landscape",
    "urban landscape",
    "woodland patch",
    "fragment",
    "remnant",
    "patch size",
    "area sensitive",
    "forest corridor",
    "forest area",
    "woodland area",
    "woodlot area",
    "grove size",
    "habitat island",
    "habitat area",
    "area requirement",
    "forest island",
    "forest size", "small forest")

space <- c("canopy cover", "forest", "woodland", "woodlot", "broadleaf", "coniferous", 'deciduous')



should_stem <- function(word){
  splitup <- strsplit(word, " ")[[1]]
  for(i in 1:length(splitup)){
    wordcut <- SnowballC::wordStem(splitup[i], language="en")
    stem_length <- nchar(wordcut)
    
    if(i==1){
      if(stem_length > 3){
        words <- paste(wordcut, "* ", sep="")
      }
      if(stem_length <= 3){
        words <- paste(splitup[i], " ", sep="")
      }
    }
    if(i > 1){
      if(stem_length > 3){
        words <- paste(words, wordcut, "* ", sep="")
      }
      if(stem_length <= 3){
        words <- paste(words, splitup[i], " ", sep="")
      }
    }
  }
  
  words <- trimws(words)
  
  return(words)
}


load("./output/search_results_deduped.rda")

texts <- tolower(paste(search_results$title, search_results$abstract, search_results$keywords, search_results$taxa_notes, search_results$concepts))

# re-run search terms
source("./scripts/definitions.R")

time <- gsub("i\\*", "\\*", unlist(lapply(time, should_stem)))
outcome <- gsub("i\\*", "\\*", unlist(lapply(outcome, should_stem)))
population <- gsub("i\\*", "\\*", unlist(lapply(population, should_stem)))
intervention <- gsub("i\\*", "\\*", unlist(lapply(intervention, should_stem)))
space <- gsub("i\\*", "\\*", unlist(lapply(space, should_stem)))

search_dictionary <- list(outcome=append(time, outcome), population=population, intervention=intervention, space=space)


search_info <- topictagger::tag_strictly(texts,
                                         scheme = search_dictionary,
                                         allow_multiple = TRUE)


search_info[search_info>0] <- 1

all_cats <- rowSums(search_info)

missing_terms <- all_cats!=4

# priorities

# geotag 
setwd("../gap-mapping/")
geodat <- read.csv("./data/ontologies/worldcities.csv")

geodat <- geodat[-which(tolower(geodat$city) %in% litsearchr::custom_stopwords), ]
geodat <- geodat[-which(nchar(geodat$city) < 4), ]
geodat <- geodat[-which(geodat$city == "Male"), ]
geodat <- geodat[-which(tolower(geodat$admin_name) == "pest"), ]

# We also need to correct for any duplicated region names, which are often
# common ways to refer to geographic areas (e.g. "Center").
checkgeo <-  as.data.frame.matrix(table(geodat$admin_name, geodat$country))
twocountries <- checkgeo
twocountries[twocountries > 0] <- 1
twocountries <- rowSums(twocountries)
twocountries <- twocountries[twocountries > 1]

geodat$admin_name[geodat$admin_name %in% names(twocountries)] <-
  paste(geodat$admin_name[geodat$admin_name %in% names(twocountries)],
        geodat$country[geodat$admin_name %in% names(twocountries)])

geodat$admin_name[which(
  tolower(geodat$admin_name) %in% c("central", "southern", 
                                    "northern", "western", "eastern"))] <-
  paste(geodat$admin_name[which(
    tolower(geodat$admin_name) %in% c("central", "southern", "northern", 
                                      "western", "eastern"))],
    geodat$country[which(
      tolower(geodat$admin_name) %in% c("central", "southern", "northern", 
                                        "western", "eastern")
    )])

# For good measure, we drop all duplicated cities anywhere
geodat <- geodat[-which(duplicated(geodat$city)), ]


# We may want to include other types of information associated with the 
# geographic tags we are able to pull, such as UN classification. 
undata <- read.csv("./data/ontologies/UN_countries.csv")
geodat$continents <-
  undata$Continent[match(geodat$iso3, undata$ISO3.Alpha.code)]
geodat$region <-
  undata$Region[match(geodat$iso3, undata$ISO3.Alpha.code)]
geodat$subregion <-
  undata$Subregion[match(geodat$iso3, undata$ISO3.Alpha.code)]
geodat$admin_name[geodat$country == "Sri Lanka"] <-
  paste(geodat$admin_name[geodat$country == "Sri Lanka"], "Sri Lanka")

setwd("../ebcmd/")

# We need to put our geotagging references in the proper hierarchy 
geodat <- geodat[,colnames(geodat) %in% c('continents', 'region', 'subregion', 'country', 'admin_name', 'city')]

l2 <- topictagger::make_entry(geodat$country, geodat$admin_name)
l3 <- topictagger::make_entry(geodat$region, geodat$country, l2)
l4 <- topictagger::make_entry(geodat$continents, geodat$region, l3)

geotags <- quanteda::dictionary(l4)


geo_info <- topictagger::tag_strictly(texts,
                                      scheme = geotags,
                                      allow_multiple = FALSE)


rm(l2,l3,l4, time, space, outcome, population, twocountries,intervention, checkgeo, geodat, undata)

continents <- unlist(lapply(geo_info, function(x) {
  strsplit(x, "\\.")[[1]][1]
}))
countries <- unlist(lapply(geo_info, function(x) {
  strsplit(x, "\\.")[[1]][4]
}))

# should also tag taxa?

aou <- read.csv("./data/NACC_list_species.csv")
aou$common <- rep("", nrow(aou))
n.fam <- length(unique(aou$family))
for(i in 1:n.fam){
  tmp <- aou[aou$family==unique(aou$family)[i],]
  common_name <- sort(table(strsplit(paste(tmp$common_name, collapse = " "), " ")[[1]]), decreasing = T)
  if(common_name[1]>1){
    aou$common[aou$family==unique(aou$family)[i]] <- names(common_name)[1]
  }
}


taxa <- aou[,c("order", "family", "common", "genus")]

taxa <- unique(taxa)
taxa[taxa==""] <- NA
taxa$common <- paste(taxa$common, "*", sep="")
taxa <- apply(taxa, 2, tolower)
taxa <- topictagger::create_dictionary(taxa)
taxa

  
taxa_tags <- topictagger::tag_strictly(tolower(paste("x", search_results$title, "x", sep=" ")),
                                      scheme = taxa,
                                      allow_multiple = TRUE)


passerines <- rowSums(taxa_tags[,grep("passeriformes", colnames(taxa_tags))])
length(which(passerines>0))
which(rowSums(taxa_tags)==0)
table(rowSums(taxa_tags))

rm(aou, taxa, n.fam, i)

retain <- which(missing_terms==FALSE & continents%in%c("Northern America", NA) & (passerines >0 | rowSums(taxa_tags)==0))

drop <- which(missing_terms==TRUE)

single_screen <- which(missing_terms==FALSE)[!which(missing_terms==FALSE)%in%retain]

length(drop) + length(retain)+length(single_screen) # cool, checks out

head(search_results$title[retain], 20)
grep("goshawk", tolower(search_results$title))

rowSums(taxa_tags)[45]
sort(taxa_tags[45,])

x <- synthesisr::write_ris(synthesisr::as.bibliography(search_results[retain,]))

x <- x[-grep("TY  - ", x)]

x[1] <- paste("TY  - JOUR", x[1], sep="\n")
end_rows <- grep("ER  -", x)
x[end_rows[-4747]+1] <- paste(x[end_rows[-4747]+1], "\n", "TY  - JOUR", sep="")

x1 <- x[1:end_rows[2500]]
x2 <- x[(end_rows[2500]+1):end_rows[4747]]

x2[1] <- "TY  - JOUR"

writeLines(x1, "./output/retained_results1.ris")
writeLines(x2, "./output/retained_results2.ris")



table(x2[grep("TY  -", x2)])
table(x2[grep("ER  -", x2)])

x2[length(x)-1]

