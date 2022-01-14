zoorec <- list.files("./data/EBCMD_searches-20210119T023018Z-001/EBCMD_searches/ZooRec/", recursive = T, full.names = T)
zoorec <- synthesisr::read_refs(zoorec)

biosis <- list.files("./data/EBCMD_searches-20210119T023018Z-001/EBCMD_searches/BIOSIS", recursive = T, full.names = T)
biosis <- synthesisr::read_refs(biosis)

pq <- list.files("./data/EBCMD_searches-20210119T023018Z-001/EBCMD_searches/PQDiss", recursive = T, full.names = T)
pq <- synthesisr::read_refs(pq)

biorxiv <- list.files("./data/EBCMD_searches-20210119T023018Z-001/EBCMD_searches/bioRxiv", recursive = T, full.names = T)
biorxiv <- synthesisr::read_refs(biorxiv)

Scopus <- list.files("./data/EBCMD_searches-20210119T023018Z-001/EBCMD_searches/Scopus", recursive = T, full.names = T)
Scopus <- synthesisr::read_refs(Scopus)

ebsco <- list.files("./data/EBCMD_searches-20210119T023018Z-001/EBCMD_searches", recursive = F, include.dirs = F, pattern = ".ris", full.names = T)
ebsco <- synthesisr::read_refs(ebsco)

full_results <- synthesisr::merge_columns(list(zoorec, biosis, pq, biorxiv, Scopus, ebsco))
save(full_results, file="./data/full_results.rda")

exact_titles <- which(duplicated(tm::removePunctuation(tolower(full_results$title))) & !is.na(full_results$title))
exact_dois <- which(duplicated(tm::removePunctuation(tolower(full_results$doi))) & !is.na(full_results$doi))

partial_dedupe <- full_results[-unique(append(exact_dois, exact_titles)),]

# now abstracts
exact_abstracts <- which(duplicated(tm::removePunctuation(tolower(partial_dedupe$abstract))) & !is.na(partial_dedupe$abstract))

mostly_deduped <- partial_dedupe[-exact_abstracts,]

nospaces <- which(duplicated(tm::removePunctuation(gsub(" ", "", tolower(mostly_deduped$title)))) & !is.na(mostly_deduped$title))
search_results <- mostly_deduped[-nospaces,]

save(search_results, file="./output/search_results_deduped.rda")






 
