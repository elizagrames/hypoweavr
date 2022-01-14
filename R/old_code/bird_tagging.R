input.data$taxa[grep("Swainson's warbler", input.data$taxa)] <- "SWWA"

input.data$taxa[103] <- gsub("black-throated blue warbler", "btbw", input.data$taxa[103])
input.data$taxa[input.data$taxa=="buff-breasted flycatchers"] <- "buff-breasted flycatcher"

birdcodes <- read.csv("~/Documents/research-projects/ebcmd/data/IBP-Alpha-Codes20.csv")

birdcodes <- apply(birdcodes, 2, tolower)
birdcodes <- data.frame(birdcodes)

bird_dictionary <- topictagger::create_dictionary(birdcodes[,c(2,4,5)])
tagged_birds <- topictagger::tag_strictly(tolower(input.data$taxa), bird_dictionary)

fourlettercodes <- birdcodes[,'SPEC']
bird_dat <- array(dim=c(nrow(tagged_birds), length(fourlettercodes)))

species <- unlist(lapply(colnames(tagged_birds), function(x){
  strsplit(x, "\\.")[[1]][1]
} ))


for(i in 1:length(fourlettercodes)){
  bird_dat[,i] <- rowSums(tagged_birds[,species==fourlettercodes[i]])
}



input.data$taxa[rowSums(bird_dat)==0]

studied <- fourlettercodes[which(colSums(bird_dat)>0)]
scinames <- birdcodes$SCINAME[birdcodes$SPEC %in% studied]


phylo <- ape::read.nexus("~/Documents/research-projects/ebcmd/data/tree-pruner-0eecac32-f757-4052-91fc-0cef9612d69f/output.nex")

birdlife <- readxl::read_xls("~/Documents/research-projects/ebcmd/data/BirdLife_Checklist_Version_3/BirdLife_Checklist_Version_3.xls")


commonnames <- birdcodes$COMMONNAME[birdcodes$SPEC %in% studied]


lott <- readxl::read_xlsx("~/Downloads/13750_2019_170_MOESM3_ESM.xlsx", sheet=2)
lott <- data.frame(lott)[-c(1:40),]

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
#birdlife <- birdlife[birdlife$]


scientificnames <- birdlife$`Scientific name`[match(commonnames, tolower(birdlife$`Common name`))]
#writeLines(scientificnames, "~/Documents/research-projects/ebcmd/data/scientific_names_Apr1.txt")



phylo <- ape::read.nexus("~/Documents/research-projects/ebcmd/data/tree-pruner-0eecac32-f757-4052-91fc-0cef9612d69f/output.nex")

gsub("_", " ", (phylo[[1]]$tip.label)) %in% scientificnames

scientificnames[!(scientificnames %in% gsub("_", " ", (phylo[[1]]$tip.label)))]


# BEWICK

bird_lookup <- cbind(original_common, commonnames, scientificnames, birdcodes$SPEC[match(original_common, birdcodes$COMMONNAME)])


bird_dat <- bird_dat[,which(fourlettercodes %in% birdcodes$SPEC[match(original_common, birdcodes$COMMONNAME)])]

included_species <- fourlettercodes[which(fourlettercodes %in% birdcodes$SPEC[match(original_common, birdcodes$COMMONNAME)])]

bird_lookup <- cbind(bird_lookup, colSums(bird_dat)[match(bird_lookup[,4], included_species)])


library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)

tre <- phylo[[1]]
tre$tip.label<-gsub("_"," ",tre$tip.label)

dat <- data.frame(count=as.numeric(bird_lookup[match(gsub("_", " ", phylo[[1]]$tip.label), bird_lookup[,3]),5]))
dat2 <- data.frame(count=as.numeric(bird_lookup[match(gsub("_", " ", phylo[[1]]$tip.label), bird_lookup[,3]),5]),
                  species=(bird_lookup[match(gsub("_", " ", phylo[[1]]$tip.label), bird_lookup[,3]),1]))

dat2[order(dat2[,1]),]

p4d <- phylo4d(tre, dat)

png("./figures/phylogeny.png", width=6,height=6,units="in",res=1200)

gridplot(p4d, col="black", center=F, scale=F, tree.type="fan", show.color.scale=T,
         show.trait=F, tip.cex=.7, 
         bar.lwd=7, trait.bg.col="white", data.xlim=c(0,30), grid.vertical=F,
         cell.col =  colorRampPalette(c("white", palette2[c(1)]))(5))

barplot(p4d, bar.col=palette2[1], center=F, scale=F, tree.type="phylogram",
        show.trait=T, tip.cex=1,  show.data.axis = T, trait.labels = c("Number of studies"),
        bar.lwd=5, trait.bg.col="white", data.xlim=c(0,30), grid.vertical=F)



dev.off()
