##Generation of Table S3


##Generation of Table S4


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


