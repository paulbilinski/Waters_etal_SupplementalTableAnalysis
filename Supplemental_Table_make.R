#paralog analysis

#load lib
library(plyr)
#get data
x<-read.table("newblasts",header=T)
#make subset where query=subject
y=subset(x,as.character(x$query_id)==as.character(x$subject_id))
#get max length for each gene
long_y=ddply(y,.(query_id),summarise,full_length=max(length))
#summarize length and %id for each pair
sum_x=ddply(x,.(query_id,subject_id),summarise,length=sum(length),id=mean(identity))
#merge and make table of all hits with summed length
long_x<-merge(sum_x,long_y)
#calculate percentage
percent_length <- as.data.frame(long_x$length/long_x$full_length)
long_xall <- cbind(long_x,percent_length)
colnames(long_xall)[6] <- "perclength"
long_xall <- as.data.frame(long_xall)
#minimum length 75%
long_x75=long_xall[which(long_xall$perclength>=0.75),]
#long_x75 <- subset(long_xall, long_xall$%length>=0.75)
#minimum ID 90%
long_x75_90=long_x75[which(long_x75$id>=90),]
#minimum 95%
long_x75_95=long_x75[which(long_x75$id>=95),]

#Above90

query_parent75_90 <- as.data.frame(long_x75_90$query_id)
query_parent75_90 <- as.data.frame(substr(long_x75_90$query_id, 1, 13))
colnames(query_parent75_90)[1]<-"query_parent"
subject_parent75_90 <- as.data.frame(long_x75_90$subject_id)
subject_parent75_90 <- as.data.frame(substr(long_x75_90$subject_id, 1, 13))
colnames(subject_parent75_90)[1]<-"subject_parent"

full_x75_90 <- cbind(long_x75_90,query_parent75_90,subject_parent75_90)
final_x75_90 <- subset(full_x75_90, as.character(full_x75_90$query_parent) != full_x75_90$subject_parent)
above90 <- as.data.frame(unique(final_x75_90$query_id))
colnames(above90)[1]<-"query_id"

pegs <- read.csv("List2_Pegs_renamed.csv")
colnames(pegs)[1]<-"query_id"
above90pegs <- merge(pegs,above90,by="query_id")

megs <- read.csv("List1_Megs_renamed.csv")
colnames(megs)[1]<-"query_id"
above90megs <- merge(megs,above90,by="query_id")

#90megs=22/69, 90pegs=33/112

##Above 95

query_parent75_95 <- as.data.frame(long_x75_95$query_id)
query_parent75_95 <- as.data.frame(substr(long_x75_95$query_id, 1, 13))
colnames(query_parent75_95)[1]<-"query_parent"
subject_parent75_95 <- as.data.frame(long_x75_95$subject_id)
subject_parent75_95 <- as.data.frame(substr(long_x75_95$subject_id, 1, 13))
colnames(subject_parent75_95)[1]<-"subject_parent"

full_x75_95 <- cbind(long_x75_95,query_parent75_95,subject_parent75_95)
final_x75_95 <- subset(full_x75_95, as.character(full_x75_95$query_parent) != full_x75_95$subject_parent)
above95 <- as.data.frame(unique(final_x75_95$query_id))
colnames(above95)[1]<-"query_id"

pegs <- read.csv("List2_Pegs_renamed.csv")
colnames(pegs)[1]<-"query_id"
above95pegs <- merge(pegs,above95,by="query_id")

megs <- read.csv("List1_Megs_renamed.csv")
colnames(megs)[1]<-"query_id"
above95megs <- merge(megs,above95,by="query_id")

#above95megs = 12/69, 95pegs=6/112

##Generation of Table S4
Schnablemaize1 <- read.csv("MaizeSubgenome1.csv")
Schnablemaize2 <- read.csv("MaizeSubgenome2.csv")
Waters_subset <- read.csv("Waters_Genesubset.csv")
MEGS <- read.csv("List1_Megs.csv")
colnames(MEGS) <- c("locus")
PEGS <- read.csv("List2_Pegs.csv")
colnames(PEGS) <- c("locus")

colnames(Waters_subset)[1] <- "locus"
Not <- Waters_subset$locus
NotMegs <- setdiff(Waters_subset$locus, MEGS$name2)
NotMegsPegs <- setdiff(NotMegs, PEGS$name2)
ExcludedMP <- as.data.frame(x=NotMegsPegs)
colnames(ExcludedMP)[1] <- "locus"
subgenome1 <- as.data.frame(Schnablemaize1)
colnames(subgenome1)[1] <- "locus"
colnames(ExcludedMP)[1] <- "locus"
watersinsubgenome1 <- merge(subgenome1, ExcludedMP, by="locus")
subgenome2 <- as.data.frame(Schnablemaize2)
colnames(subgenome2)[1] <- "locus"
watersinsubgenome2 <- merge(subgenome2, ExcludedMP, by="locus")
colnames(MEGS)[1] <- "locus"
colnames(PEGS)[1] <- "locus"
megsinsubgenome1 <- merge(subgenome1, MEGS, by="locus")
megsinsubgenome2 <- merge(subgenome2, MEGS, by="locus")
pegsinsubgenome1 <- merge(subgenome1, PEGS, by="locus")
pegsinsubgenome2 <- merge(subgenome2, PEGS, by="locus")

zvz <- read.csv("CDS_zeavzea.csv")
zvzname1 <- data.frame(DNDS=zvz$DNDS, name=zvz$name1)
zvzname2 <- data.frame(DNDS=zvz$DNDS, name=zvz$strand2)
zvznames <- rbind(zvzname1, zvzname2)
duplicated(zvznames$name)
unique(zvznames$name)

colnames(zvz)[8] <- "locus"
temp1 <- merge(watersinsubgenome1, zvz, by="locus")
colnames(zvz)[8] <- "name1"
colnames(zvz)[20] <- "locus"
temp2 <- merge(watersinsubgenome1, zvz, by="locus")
watersinsubgenome1_dnds <- merge(temp1, temp2, by="locus") #if exists in column 1 remove from column 2
zvz <- read.csv("CDS_zeavzea.csv")
colnames(zvz)[8] <- "locus"
temp1 <- merge(watersinsubgenome2, zvz, by="locus")
colnames(zvz)[8] <- "name1"
colnames(zvz)[20] <- "locus"
temp2 <- merge(watersinsubgenome2, zvz, by="locus")
watersinsubgenome2_dnds <- merge(temp1, temp2, by="locus")
#Megs in subgenomes
zvz <- read.csv("CDS_zeavzea.csv")
colnames(zvz)[8] <- "locus"
temp1 <- merge(megsinsubgenome1, zvz, by="locus")
colnames(zvz)[8] <- "name1"
colnames(zvz)[20] <- "locus"
temp2 <- merge(megsinsubgenome1, zvz, by="locus")
megsinsubgenome1_dnds <- merge(temp1, temp2, by="locus")
zvz <- read.csv("CDS_zeavzea.csv")
colnames(zvz)[8] <- "locus"
temp1 <- merge(megsinsubgenome2, zvz, by="locus")
colnames(zvz)[8] <- "name1"
colnames(zvz)[20] <- "locus"
temp2 <- merge(megsinsubgenome2, zvz, by="locus")
megsinsubgenome2_dnds <- merge(temp1, temp2, by="locus")
#Pegs in subgenomes
zvz <- read.csv("CDS_zeavzea.csv")
colnames(zvz)[8] <- "locus"
temp1 <- merge(pegsinsubgenome1, zvz, by="locus")
colnames(zvz)[8] <- "name1"
colnames(zvz)[20] <- "locus"
temp2 <- merge(pegsinsubgenome1, zvz, by="locus")
pegsinsubgenome1_dnds <- merge(temp1, temp2, by="locus",all=TRUE)
zvz <- read.csv("CDS_zeavzea.csv")
colnames(zvz)[8] <- "locus"
temp1 <- merge(pegsinsubgenome2, zvz, by="locus")
colnames(zvz)[8] <- "name1"
colnames(zvz)[20] <- "locus"
temp2 <- merge(pegsinsubgenome2, zvz, by="locus")
pegsinsubgenome2_dnds <- merge(temp1, temp2, by="locus", all=TRUE)

##Generation of Table S6, Domestication Candidates
fgs_1 <- read.csv("FGS_v1_genes.csv")
fgs_1_nodp <- subset(fgs_1, !duplicated(fgs_1[,1]))
colnames(fgs_1_nodp)[1] <- "locus"
refgen1waters <- merge(fgs_1_nodp,Waters_subset,by="locus")
colnames(fgs_1_nodp)[1] <- "locus"
refgen1megs <- merge(fgs_1_nodp,MEGS)
refgen1pegs <- merge(fgs_1_nodp,PEGS)
length(refgen1waters$locus)
length(refgen1megs$locus)
length(refgen1pegs$locus)

improv <- read.csv("improvement_candidates.csv",header=FALSE)
domest <- read.csv("domestication_candidates.csv",header=FALSE)
colnames(improv)[1] <- "locus"
colnames(domest)[1] <- "locus"
colnames(MEGS)[1] <- "locus"
colnames(PEGS)[1] <- "locus"
colnames(MEGS)[1] <- "locus"
length(merge(refgen1waters,domest,by="locus")$locus)
length(merge(PEGS,domest,by="locus")$locus)
length(merge(MEGS,domest,by="locus")$locus)
length(merge(refgen1waters,improv,by="locus")$locus)
length(merge(PEGS,improv,by="locus")$locus)
length(merge(MEGS,improv,by="locus")$locus)

##Generation of Table S7, Population Genetic Statistics for Imprinted Genes
#Analysis of Popgen stats
full <- read.csv("Sorghum_Zearefgen2_knks.csv")
data <- subset(full, full$Ks<1)
MEGS <- read.csv("List1_Megs.csv")
MEGS_dnds <- merge(data,MEGS,by="name2")
PEGS <- read.csv("List2_Pegs.csv")
PEGS_dnds <- merge(data,PEGS,by="name2")
KnKs6 <- subset(data,data$KnKs<6)
computeall <- read.csv("Fullcompute.csv")
colnames(KnKs6)[20] <- "locus"
MEGS <- read.csv("List1_Megs.csv")
colnames(MEGS) <- c("locus")
PEGS <- read.csv("List2_Pegs.csv")
colnames(PEGS) <- c("locus")
Not <- Waters_subset$locus
NotMegs <- setdiff(Waters_subset$locus, MEGS$name2)
NotMegsPegs <- setdiff(NotMegs, PEGS$name2)
ExcludedMP <- as.data.frame(x=NotMegsPegs)
colnames(ExcludedMP)[1] <- "locus"
computesyn <- subset(computeall, computeall$sitetype=="SYN")
computesynMZ <- subset(computesyn, computesyn$taxa=="mz")
dadcomputesynMZ <- merge(computesynMZ,PEGS,by="locus")
momcomputesynMZ <- merge(computesynMZ,MEGS,by="locus")
stats <- read.delim("Maize_Genestats.txt", header=TRUE, sep="\t", na.strings="NA")
stats <- cbind(stats, stats$ThetaPi / stats$seqbp)
colnames(stats)[8] <- "ThetaPibp"
MEGS <- read.csv("List1_Megs.csv")
momstats <- merge(stats,MEGS,by="name2")
PEGS <- read.csv("List2_Pegs.csv")
dadstats <- merge(stats,PEGS,by="name2")
Waters_subset <- read.csv("Waters_Genesubset.csv")
colnames(Waters_subset)[1] <- "locus"
colnames(stats)[1] <- "locus"
Waters_synMZ <- merge(Waters_subset,computesynMZ,by="locus")
Waters_stats <- merge(Waters_subset,stats,by="locus")
Waters_synMZ_noMP <- merge(Waters_synMZ, ExcludedMP, by="locus")
Waters_stats_noMP <- merge(Waters_stats, ExcludedMP, by="locus")
Waters_KnKs <- merge(Waters_subset,KnKs6,by="locus")
Waters_KnKs_noMP <- merge(Waters_KnKs, ExcludedMP, by="locus")

#Hprime
summary(Waters_synMZ_noMP$Hprime)
summary(dadcomputesynMZ$Hprime)
summary(momcomputesynMZ$Hprime)
wilcox.test(dadcomputesynMZ$Hprime, Waters_synMZ_noMP$Hprime)
wilcox.test(momcomputesynMZ$Hprime, Waters_synMZ_noMP$Hprime)
#TajD
summary(Waters_synMZ_noMP$TajD)
summary(dadcomputesynMZ$TajD)
summary(momcomputesynMZ$TajD)
wilcox.test(dadcomputesynMZ$TajD, Waters_synMZ_noMP$TajD)
wilcox.test(momcomputesynMZ$TajD, Waters_synMZ_noMP$TajD)
#Hapdiv
summary(Waters_synMZ_noMP$hapdiv)
summary(dadcomputesynMZ$hapdiv)
summary(momcomputesynMZ$hapdiv)
wilcox.test(dadcomputesynMZ$hapdiv, Waters_synMZ_noMP$hapdiv)
wilcox.test(momcomputesynMZ$hapdiv, Waters_synMZ_noMP$hapdiv)
#ThetaPibp
summary(Waters_stats_noMP$ThetaPibp)
summary(dadstats$ThetaPibp)
summary(momstats$ThetaPibp)
wilcox.test(dadstats$ThetaPibp, Waters_stats_noMP$ThetaPibp)
wilcox.test(momstats$ThetaPibp, Waters_stats_noMP$ThetaPibp)
#DnDs
summary(Waters_KnKs_noMP$KnKs)
summary(PEGS_dnds$KnKs)
summary(MEGS_dnds$KnKs)
wilcox.test(PEGS_dnds$KnKs, Waters_KnKs_noMP$KnKs)
wilcox.test(MEGS_dnds$KnKs, Waters_KnKs_noMP$KnKs)


