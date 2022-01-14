articles <- read.csv("~/Downloads/Articles_P46561_0222_A4747.csv")
answers <- read.csv("~/Downloads/Answers_P46561_0222_A4747.csv")
users <- read.csv("~/Downloads/UserAnswers_P46561_0222_A4747.csv")

resolutions <- users[users$User.Name=="chris.elphick",c(1,5,14)]
colnames(resolutions)[2] <- "resolution"

documents1 <- synthesisr::read_ref("~/Documents/research-projects/ebcmd/output/retained_results1.ris")
documents2 <- synthesisr::read_ref("~/Documents/research-projects/ebcmd/output/retained_results2.ris")

documents <- synthesisr::merge_columns(documents1, documents2)

documents$title[which(is.na(match(documents$title, articles$Title)))]
answers$Title[grep("Neotropical migratory", answers$Title)]
documents$title[4726] <- paste(documents$title[4726], ".", sep="")
colnames(answers)[16] <- "title"

documents <- dplyr::left_join(documents, answers, by="title")

documents <- dplyr::left_join(documents, resolutions, by="Article.ID")

documents$resolution[documents$Status!="conflict"] <- documents$Include[documents$Status!="conflict"]


full_text <- documents[documents$resolution=="true",]

table(full_text$Tested.mechanisms)


unclear <- grep("unclear", full_text$North.America)
full_text$keywords[unclear]

not_NA <- unclear[c(1,6,9,16,17,25,26,28,33,36,37,38,39,41,43,53,54,55,79,91)]

full_text$resolution[not_NA] <- "false"

full_text <- full_text[full_text$resolution=="true",]


aou <- read.csv("./data/NACC_list_species.csv")
songbirds <- aou[aou$order=="Passeriformes",]

birds <- list(birds=tolower(append(songbirds$common_name, songbirds$species)))

birds <- quanteda::dictionary(birds)

single_species <- as.numeric(quanteda::dfm(tolower(paste(full_text$title, full_text$abstract, full_text$keywords)), dictionary=birds))

priorities <- full_text[single_species!=0,]

table(is.na(full_text$doi))

writeLines(litsearchr::write_title_search(gsub("\\)", "", gsub("\\(", "", full_text$title[is.na(full_text$doi)]))), "./output/need_doi_all_titles.txt")

for(i in 1:nrow(full_text)){
  if(is.na(full_text$doi[i])){
    newdoi <- sort(table(full_results$doi[which(tolower(full_results$title) %in% tolower(full_text$title[i]))]), decreasing = T)[1]
    if(!is.na(newdoi)){
      full_text$doi[i] <- names(newdoi)
    }
  }
}

for(i in 1:nrow(full_text)){
  if(is.na(full_text$url[i])){
    newurl <- sort(table(full_results$url[which(tolower(full_results$title) %in% tolower(full_text$title[i]))]), decreasing = T)[1]
    if(!is.na(newurl)){
      full_text$url[i] <- names(newurl)
    }
  }
}



full_results$url[which(tolower(full_results$title) %in% tolower(full_text$title[i]))]








refs <- synthesisr::write_ris(synthesisr::as.bibliography(full_text[,1:25]))

writeLines(refs, "./output/full_text_screening.ris")


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




need_doi <- synthesisr::read_ref("./output/priority_need_dois.bib")

writeLines(litsearchr::write_title_search(gsub("\\)", "", gsub("\\(", "", need_doi$title))), "./output/need_doi.txt")




