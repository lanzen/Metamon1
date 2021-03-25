setwd("/home/alanzen/projects/Metamon/WP2/SWARM_WP2_all_20200708")

require(vegan)
load(".RData")

table(row.names(otusR.18S) == row.names(mdR.18S))

selectPCR = (mdR.18S$Experiment=="GrabEXandPCRRep")

md.pcr.18S = droplevels(mdR.18S[selectPCR ,])
md.pcr.18S$Grab=as.factor(md.pcr.18S$Grab)
summary(md.pcr.18S) 
# 20 samples: 5 extract groups each with (5 and 3 replicates) x (all and 1 PCR).

otus.pcr.18S = otusR.18S[selectPCR ,]
summary(rowSums(otus.pcr.18S)) # 91 -- 139k reads.
writeDivStats("PCRReps_diversity_18S.csv", otudistro=otus.pcr.18S)
divR.pcr.18S = read.csv("PCRReps_diversity_18S.csv",row.names=1)
divR.pcr.18S = divR.pcr.18S[,c("Reads","Richness","Rarefied.richness","H","J",
                       "Chao1")]

# Make one dataset for each in vitro pooling type (5-5, 5-1, 3-3, 3-1)
ffDiv = divR.pcr.18S[rep(c(T,F,F,F),5),]
foDiv = divR.pcr.18S[rep(c(F,T,F,F),5),]
ttDiv = divR.pcr.18S[rep(c(F,F,T,F),5),]
toDiv = divR.pcr.18S[rep(c(F,F,F,T),5),]

## Five vs 1 PCR reps with 5 extraction reps, relative differences:
fiveToOne = ffDiv / foDiv
summary(fiveToOne)
# Reads           Richness      Rarefied.richness       H               J              Chao1      
# Min.   :0.9923   Min.   :0.9984   Min.   :0.983     Min.   :0.954   Min.   :0.9542   Min.   :1.158  
# Median :1.0673   Median :1.3200   Median :1.258     Median :1.017   Median :0.9764   Median :1.541  
# Mean   :1.0738   Mean   :1.2624   Mean   :1.199     Mean   :1.002   Mean   :0.9715   Mean   :1.470  
# Max.   :1.1634   Max.   :1.4409   Max.   :1.373     Max.   :1.027   Max.   :0.9813   Max.   :1.697  
# -> 20% higher rarefied richness v consistent compared to richness per repliacte for PCR (goes from 1000 to 1200) <-

wilcox.test(fiveToOne$Rarefied.richness-1, alternative="greater")
# V = 14, p-value = 0.0625
# alternative hypothesis: true location is greater than 0

ano(log(fiveToOne$Rarefied.richness))

## Three vs 1 PCR reps with 3 extraction reps:
threeToOne = ttDiv / toDiv
summary(threeToOne)
# Reads           Richness     Rarefied.richness       H                J              Chao1      
# Min.   :0.9172   Min.   :1.074   Min.   :1.056     Min.   :0.9491   Min.   :0.9403   Min.   :1.213  
# Median :1.0513   Median :1.264   Median :1.235     Median :1.0013   Median :0.9663   Median :1.384  
# Mean   :1.0377   Mean   :1.227   Mean   :1.188     Mean   :0.9965   Mean   :0.9693   Mean   :1.372  
# Max.   :1.1350   Max.   :1.370   Max.   :1.309     Max.   :1.0363   Max.   :0.9930   Max.   :1.527  
# -> 18% higher rarefied richness also consistent compared to richness per repliacte for PCR (goes from 1000 to 1150-ish) <-

wilcox.test(threeToOne$Rarefied.richness-1, alternative="greater")
# V = 15, p-value = 0.03
# alternative hypothesis: true location is greater than 0


## Five vs 3 extraction reps with all PCRs:
fiveToThree = ffDiv / ttDiv
summary(fiveToThree)
# Reads          Richness     Rarefied.richness       H                J              Chao1      
# Min.   :0.987   Min.   :1.064   Min.   :1.062     Min.   :0.9811   Min.   :0.9691   Min.   :1.070  
# Median :1.013   Median :1.109   Median :1.070     Median :1.0035   Median :0.9880   Median :1.200  
# Mean   :1.054   Mean   :1.124   Mean   :1.095     Mean   :1.0021   Mean   :0.9867   Mean   :1.185  
# Max.   :1.139   Max.   :1.180   Max.   :1.156     Max.   :1.0146   Max.   :1.0014   Max.   :1.278  



