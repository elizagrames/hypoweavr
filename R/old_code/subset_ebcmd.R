input.data <- readxl::read_xlsx("~/Downloads/fulltext2.xlsx")

samplesize <- ceiling(nrow(input.data)*.2)

seedno <- sample(1:100) #55
seedno


set.seed(seedno)
random_subset <- sample(1:nrow(input.data))
use_these <- input.data[random_subset,]

done <- yes <- c()
for(i in 1:100){
  set.seed(i)
  random_subset <- sample(1:nrow(input.data))
  use_these <- input.data[random_subset,]
  
  done[i] <- sum(!is.na(use_these$include[1:79]))
  yes[i] <- sum(use_these$include[1:79]=="yes", na.rm = T)
}
which.max(done)
which.max(yes)

plot(done, type="l", ylim=c(0,.5))
lines(yes)

colnames(use_these)
write.csv(use_these, file="~/Documents/research-projects/ebcmd/data/random_subset_seed55.csv")


which(use_these$include=="yes")

which(!is.na(use_these$include))



sum(!is.na(use_these$include[1:200]))
