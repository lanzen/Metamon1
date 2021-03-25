setwd("/home/alanzen/projects/Metamon/WP2/SWARM_WP2_all_20200708")

require(vegan)
load(".RData")

table(row.names(otusR.COI.nc) == row.names(md.COI.nc))

selectPCR = (md.COI.nc$Experiment=="GrabEXandPCRRep")

md.pcr.COI = droplevels(md.COI.nc[selectPCR ,])
md.pcr.COI$Grab=as.factor(md.pcr.COI$Grab)
summary(md.pcr.COI) 
# 20 samples: 5 extract groups each with (5 and 3 replicates) x (all and 1 PCR).

otus.pcr.COI = otusR.COI.nc[selectPCR ,]
summary(rowSums(otus.pcr.COI)) # 1286 -- 60k reads.
writeDivStats("PCRReps_diversity_COI.csv", otudistro=otus.pcr.COI)
divR.pcr.COI = read.csv("PCRReps_diversity_COI.csv",row.names=1)
divR.pcr.COI = divR.pcr.COI[,c("Reads","Richness","Rarefied.richness","H","J",
                       "Chao1")]

ffDiv = divR.pcr.COI[rep(c(T,F,F,F),5),]
foDiv = divR.pcr.COI[rep(c(F,T,F,F),5),]
ttDiv = divR.pcr.COI[rep(c(F,F,T,F),5),]
toDiv = divR.pcr.COI[rep(c(F,F,F,T),5),]

## Five vs 1 PCR reps with 5 extraction reps:
fiveToOne = ffDiv / foDiv
summary(fiveToOne)
# Reads           Richness     Rarefied.richness       H                J              Chao1      
# Min.   :0.5640   Min.   :1.294   Min.   :1.004     Min.   :0.9159   Min.   :0.8427   Min.   :1.370  
# Median :1.0754   Median :1.533   Median :1.227     Median :1.0594   Median :0.9343   Median :1.798  
# Mean   :1.4141   Mean   :1.645   Mean   :1.223     Mean   :1.1078   Mean   :1.0194   Mean   :1.840  
# Max.   :3.2041   Max.   :2.092   Max.   :1.426     Max.   :1.4611   Max.   :1.3903   Max.   :2.312 
# Actually looks nicer than 18S

wilcox.test(fiveToOne$Rarefied.richness-1, alternative="greater")
# V = 15, p-value = 0.03
# alternative hypothesis: true location is greater than 0

## Three vs 1 PCR reps with 3 extraction reps:
threeToOne = ttDiv / toDiv
summary(threeToOne)
# Reads           Richness     Rarefied.richness       H                J              Chao1      
# Min.   :0.7741   Min.   :1.306   Min.   :1.028     Min.   :0.9784   Min.   :0.9278   Min.   :1.366  
# Median :0.9731   Median :1.478   Median :1.254     Median :1.0365   Median :0.9422   Median :1.596  
# Mean   :1.3922   Mean   :1.540   Mean   :1.225     Mean   :1.1281   Mean   :1.0453   Mean   :1.606  
# Max.   :3.3896   Max.   :1.858   Max.   :1.398     Max.   :1.5546   Max.   :1.4577   Max.   :1.815  
# Increases more than 18S

wilcox.test(threeToOne$Rarefied.richness-1, alternative="greater")
# V = 15, p-value = 0.03
# alternative hypothesis: true location is greater than 0


## Five vs 3 extraction reps with all PCRs:
fiveToThree = ffDiv / ttDiv
summary(fiveToThree)



