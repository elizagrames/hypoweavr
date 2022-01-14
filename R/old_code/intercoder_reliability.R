answers <- read.csv("~/Downloads/Articles_P46561_0308_A4747.csv")
users <- read.csv("~/Downloads/UserAnswers_P46561_0308_A4747.csv")


eliza <- users[users$User.Name=="eliza.grames",c(1,5)]
danielle <- users[users$User.Name=="danielle.schwartz",c(1,5)]
chris <- users[users$User.Name=="chris.elphick",c(1,5)]


colnames(eliza)[2] <- "eliza"
colnames(danielle)[2] <- "danielle"
colnames(chris)[2] <- "chris"

answers <- dplyr::left_join(answers, eliza, "Article.ID")
answers <- dplyr::left_join(answers, danielle, "Article.ID")
answers <- dplyr::left_join(answers, chris, "Article.ID")

duplicates <- which(!is.na(answers$eliza) & !is.na(answers$danielle))

answers <- answers[duplicates,]

answers$resolution <- rep(NA, nrow(answers))

answers$resolution[which(answers$eliza=="false" & answers$danielle=="false")] <- "false"
answers$resolution[which(answers$eliza=="true" & answers$danielle=="true")] <- "true"
answers$resolution[answers$eliza!=answers$danielle] <- answers$chris[answers$eliza!=answers$danielle]


tp.e <- sum(answers$resolution=="true" & answers$eliza=="true")
tn.e <- sum(answers$resolution=="false" & answers$eliza=="false")
fp.e <- sum(answers$resolution=="false" & answers$eliza=="true")
fn.e <- sum(answers$resolution=="true" & answers$eliza=="false")


tp.d <- sum(answers$resolution=="true" & answers$danielle=="true")
tn.d <- sum(answers$resolution=="false" & answers$danielle=="false")
fp.d <- sum(answers$resolution=="false" & answers$danielle=="true")
fn.d <- sum(answers$resolution=="true" & answers$danielle=="false")









