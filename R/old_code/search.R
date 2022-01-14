source("./scripts/definitions.R")
search <- litsearchr::write_search(list(append(time, outcome), population, intervention, space), languages = "English", closure = "none", stemming = T, exactphrase = T)
writeLines(search, "./output/stemmed_search.txt")
