setwd("/home/alanzen/projects/Metamon/WP2/SWARM_WP2_spatio_202101/")

require(vegan)
require(decontam)
source('../R/utils/filtering.R')
source('../R/utils/diversity.r')
source('../R/utils/taxaplot.R')
source('../R/utils/mergeOTUTable.R')
source('../R/utils/correlationTests.r')


## Read and filter data

md = read.csv(file="metadata.csv", header=T, row.names=1)

#order metadata alphabetically
md = md[order(row.names(md)),]

md$Plate18S = as.factor(md$Plate18S)
md$PlateCOI = as.factor(md$PlateCOI)
md$SSU_ID  = as.character(md$SSU_ID)
summary(md)

otus.18S.all = read.delim("18S/CREST_LULU/SWARM_table_curated.tsv",
                          row.names=1,header=T,sep="\t")

names(otus.18S.all) = gsub("\\.","-",names(otus.18S.all))
table(md$SSU_ID %in% names(otus.18S.all))  # all 36

otus.18S = as.data.frame(t(otus.18S.all[,md$SSU_ID]))
table(row.names(otus.18S) == md$SSU_ID)
row.names(otus.18S) = row.names(md)
table(row.names(otus.18S) == row.names(md))

# Classification df
bestTx = array(dim=dim(otus.18S.all)[1])
for (i in 1:length(bestTx)) bestTx[i] <- tail(unlist(strsplit(as.character(otus.18S.all$classification[i]), 
                                                              split=";", 
                                                              fixed=TRUE)), 1)
tax.18S = data.frame(row.names=row.names(otus.18S.all), 
                     classification=otus.18S.all$classification, bestTx=bestTx)


rowSums(otus.18S)
# Only blanks have few reads
writeDivStats("QCDiversity_stats_intermediate_18S.csv", otudistro=otus.18S)
writeDivStats("QCDiversity_stats_intermediate_Met_18S.csv", otudistro=otus.18S[,grep("Metazoa",tax.18S$classification)])

# ------ Contamination control ------

dim(otus.18S) # we start with 6371 OTUs
sum(otus.18S) # and 3,783,689 reads

## Cross contaminants
otus.18S = filterCrossContaminants2(otus.18S,100)
sum(otus.18S) #3,783,268 (419 reads less)

##  Make barchart and check abundant OTUs for contamination (just for control)

otus.named = decostand(otus.18S,method="total")
names(otus.named) = paste(names(otus.named),tax.18S[names(otus.named),]$bestTx)

grouping_info<-data.frame(row.names=row.names(md), md$Type)
pdf("img/18S/SV_chart_18S_w_controls.pdf",width=17, height = 10)
taxaplot(30,grouping_info,otus.named,grouping_info)
dev.off()

# SWARM_7727.Anthoathecata looks overrepresented, and also SWARM_182.Calanoida
# but both are marine

## Prevalence based filtering based on negative controls (decontam)

om = as.matrix(otus.18S)
predicted_contaminants = isContaminant(om, neg=(md$Type=="C"))#,
                                       #batch=md$Plate18S)

summary(predicted_contaminants) #-> 1 predicted contaminant
table(row.names(predicted_contaminants) == names(otus.18S)) # Yes
sc = (!is.na(predicted_contaminants$p) & predicted_contaminants$p<=.05) # 1
sigCont = predicted_contaminants[sc,]
row.names(sigCont)
sigCont$Tax = tax.18S$bestTx[row.names(tax.18S) %in% row.names(sigCont)]
summary(as.factor(sigCont$Tax))
# Unclassified SAR, SWARM_2089

## Remove the suspected contaminant, which however is a marine diatom

cont = (names(otus.18S) %in% row.names(sigCont))
otus.18S = otus.18S[,-cont]
sum(otus.18S) #3,783,266: removed only 2 reads :D

# Remove the blank and mock samples from further analysis
otusR.18S = otus.18S[md$Type == "S",]
sum(otusR.18S) #3,783,025 (241 reads removed with the 6 control samples)
otusR.18S = otusR.18S[,colSums(otusR.18S)>0] 
dim(otusR.18S) #Retains 6358 OTUs

mdR = droplevels(md[md$Type == "S",])
table(row.names(otusR.18S) == row.names(mdR))

## Taxon based filter (replace with plankton based??)

for (taxonDel in c("No hits","Insecta","Mammalia","Aves","Myxini",
                   "Arachnida","Actinopterygii","Chondrichthyes","Collembola")) {
  toDel =  row.names(tax.18S)[grep(taxonDel,tax.18S$classification)]
  if (length(toDel)>0){
    otusR.18S = otusR.18S[,!(names(otusR.18S) %in% toDel)]
    print(paste("Removing ",taxonDel))
    print(dim(otusR.18S)[2])
    print(sum(otusR.18S))
  }
}

# starting with 3,783,025 reads in 6358 OTUs
# 1] "Removing  No hits"
# [1] 6271
# [1] 3778965
# [1] "Removing  Mammalia"
# [1] 6270
# [1] 3778960
# [1] "Removing  Myxini"
# [1] 6269
# [1] 3778958
# [1] "Removing  Arachnida"
# [1] 6266
# [1] 3778943
# [1] "Removing  Actinopterygii"
# [1] 6260
# [1] 3778901



# ---- Relative abundance ------

otus.18S.ra.all = decostand(otusR.18S, method="total")
otus.18S.H.all = decostand(otusR.18S, method="hell")

summary(rowSums(otusR.18S)) 
# 103 -- 153 thousand reads per sample -> no filtering needed

# For dissimilarity, we require a minimum of 3 expected reads assuming average constant relative OTU abundance
minAb = 3 / min(rowSums(otusR.18S)) #min is 103k so min. is 3E-5
otusR18S.ab = dropRareByAvgAbundance(otusR.18S, minAb)
dim(otusR18S.ab)
dim(otusR18S.ab)/dim(otusR.18S) #1030 OTUs remain at this cutoff, or 16%

otusR18S.ab.ra = decostand(otusR18S.ab,method="total")
otusR18S.ab.H = decostand(otusR18S.ab,method="hellinger")

tax.18SR = tax.18S[names(otusR.18S),]
tax.18S.ab = tax.18S[names(otusR18S.ab),]
table(row.names(tax.18S.ab) == names(otusR18S.ab))


otusR18S.m = otusR.18S[,grep("Metazoa",tax.18SR$classification)]
otusR18S.ab.m = otusR18S.ab[,grep("Metazoa",tax.18S.ab$classification)]
dim(otusR18S.ab.m)
dim(otusR18S.ab.m)/dim(otusR18S.ab) # 330 (32%) of OTUs are metazoan
rowSums(otusR18S.ab.m) / rowSums(otusR18S.ab)
# G1-I-1A    G1-I-1B    G1-I-1C    G1-I-1D    G1-I-1E     G1-I-5    G1-II-5   G1-III-5 
# 0.51865306 0.48336771 0.55579501 0.50068780 0.54011384 0.52531455 0.50176200 0.47249607 
# G1-IV-5     G1-V-5    G2-I-1A    G2-I-1B    G2-I-1C    G2-I-1D    G2-I-1E     G2-I-5 
# 0.05675120 0.37479629 0.19581893 0.46045407 0.35453928 0.39061202 0.20639173 0.32430633 
# G2-II-5   G2-III-5    G2-IV-5     G2-V-5    G3-I-1A    G3-I-1B    G3-I-1C    G3-I-1D 
# 0.33196260 0.37187478 0.10490844 0.09799080 0.09706375 0.11779707 0.16986783 0.21711127 
# G3-I-1E     G3-I-5    G3-II-5   G3-III-5    G3-IV-5     G3-V-5 
# 0.18070052 0.15612275 0.09873560 0.22401985 0.14608221 0.16606075 


otusR18S.ab.ra.m = decostand(otusR18S.ab.m,method="total")
otusR18S.ab.H.m = decostand(otusR18S.ab.m,method="hellinger")

## OTU plot with abundance filtered  
otus.named = decostand(otusR18S.ab,method="total")
names(otus.named) = paste(names(otus.named),tax.18S[names(otus.named),]$bestTx)
grouping_info<-data.frame(row.names=row.names(mdR), paste("Grab",mdR$Grab))
pdf("img/18S/SV_chart_18S_abFiltered.pdf",width=10, height = 5)
taxaplot(30,grouping_info,otus.named,grouping_info)
dev.off()

## OTU plot with abundance filtered metazoa
otus.named = decostand(otusR18S.ab.m,method="total")
names(otus.named) = paste(names(otus.named),tax.18S[names(otus.named),]$bestTx)
pdf("img/18S/SV_chart_18S_abFiltered_metazoa.pdf",width=10, height = 5)
taxaplot(30,grouping_info,otus.named,grouping_info)
dev.off()
  
#-------- Diversity stats --------

writeDivStats("Diversity_filtered_18S_spatial_heterogeneity.csv", otusR.18S)
writeDivStats("Diversity_filtered_18S_spatial_heterogeneity_met.csv", otusR18S.m)
writeDivStats("Diversity_filtered_18S_spatial_heterogeneity_ab.csv", otusR18S.ab)
writeDivStats("Diversity_filtered_18S_spatial_heterogeneity_ab_met.csv", otusR18S.ab.m)


divR.18S = read.csv("Diversity_filtered_18S_spatial_heterogeneity_ab.csv",row.names=1)
divR.18S = divR.18S[,c("Reads","Richness","Rarefied.richness","H","J",
                       "Chao1")]

divR.18S.met = read.csv("Diversity_filtered_18S_spatial_heterogeneity_ab_met.csv",row.names=1)
divR.18S.met = divR.18S.met[,c("Reads","Richness","Rarefied.richness","H","J",
                               "Chao1")]

cor.test(divR.18S$Reads,divR.18S$Richness)
cor.test(divR.18S$Reads,divR.18S$Rarefied.richness)

# ---- NMDS --------

nmds = metaMDS(otusR18S.ab.ra) # converges after 3 iterations, Stress 0.13, resid 3E-4
ordiplot(nmds,type="none",xlim=c(-.7,.5),ylim=c(-.5,.5))
points(nmds,pch=as.numeric(mdR$PooledExtracts),col=as.numeric(mdR$Grab),lwd=2,cex=1.2)
text(nmds,labels=paste("G",mdR$Grab,":",mdR$Sediment,sep=""),pos=1,cex=.7,col=as.numeric(mdR$Grab))
legend("bottomleft",pch=c(1,5,1,1,1), col=c("grey","grey","black","red","green"),
       legend=c("Individual extract", "5 pooled extracts", "Grab 1", "Grab 2", "Grab 3"),
       cex=.8)

nmds.m = metaMDS(otusR18S.ab.ra.m) # Does not converge, stress 0.22 after 20 iterations
ordiplot(nmds.m,type="none")#,xlim=c(-.3,.3),ylim=c(-.2,.2))
points(nmds.m,pch=as.numeric(mdR$PooledExtracts),col=as.numeric(mdR$Grab),lwd=2,cex=1.2)
text(nmds.m,labels=paste("G",mdR$Grab,":",mdR$Sediment,sep=""),pos=1,cex=.7,col=as.numeric(mdR$Grab))
legend("bottomleft",pch=c(1,5,1,1,1), col=c("grey","grey","black","red","green"),
       legend=c("Individual extract", "5 pooled extracts", "Grab 1", "Grab 2", "Grab 3"),
       cex=.8)
# ---- Dissimilarity comparisons ------

ssu.intraSub1 = otusR18S.ab.ra[mdR$Grab==1 & mdR$PooledExtracts==1,]
ssu.intraSub2 = otusR18S.ab.ra[mdR$Grab==2 & mdR$PooledExtracts==1,]
ssu.intraSub3 = otusR18S.ab.ra[mdR$Grab==3 & mdR$PooledExtracts==1,]

is1.bc.ssu = vegdist(ssu.intraSub1)
is2.bc.ssu = vegdist(ssu.intraSub2)
is3.bc.ssu = vegdist(ssu.intraSub3)
ssu.intraSub.bc = c(is1.bc.ssu,is2.bc.ssu,is3.bc.ssu)
summary(ssu.intraSub.bc) # .16 -- .46, median .24

## Samples with 5 pooled replicats from grabs 1--3
ssu.intraGrab1 = otusR18S.ab.ra[mdR$Grab==1 & mdR$PooledExtracts==5,]
ig1.bc.ssu = vegdist(ssu.intraGrab1)
summary(ig1.bc.ssu) #.36 -- .56 med .53 

ssu.intraGrab2 = otusR18S.ab.ra[mdR$Grab==2 & mdR$PooledExtracts==5,]
ig2.bc.ssu = vegdist(ssu.intraGrab2)
summary(ig2.bc.ssu) #.22 -- .44 med .36

ssu.intraGrab3 = otusR18S.ab.ra[mdR$Grab==3 & mdR$PooledExtracts==5,]
ig3.bc.ssu = vegdist(ssu.intraGrab3)
summary(ig3.bc.ssu) #.25 -- .43 med .33

igAll.bc.ssu = c(unlist(ig1.bc.ssu),unlist(ig2.bc.ssu),unlist(ig3.bc.ssu)) 
summary(igAll.bc.ssu)
# .23 -- .56, med .38

## Is BC significantly higher between pooled sub-grabs compared to between extractions?
wilcox.test(igAll.bc.ssu, ssu.intraSub.bc, alternative="greater")
# -> Yes (p=1.7E-6)

## Samples with 5 pooled replicats from all grabs sorted by grab
intraGrabAll.ssu = rbind(intraGrab1, intraGrab2, intraGrab3)
intraGrabAll.bc.ssu = as.matrix(vegdist(intraGrabAll.ssu))
betweenGrabs.bc.ssu = c(intraGrabAll.bc.ssu[c(6:15),c(1:5)],
                        intraGrabAll.bc.ssu[c(11:15),c(6:10)])
summary(betweenGrabs.bc.ssu) # .22 -- .55, med .44

## Is BC significantly higher between grabs compared to between intra-grab pooled
wilcox.test(betweenGrabs.bc.ssu, igAll.bc.ssu, alternative="greater")
# -> Yes but close call (p=0.02)


pdf("img/BCBoxplot_for_Fig3_18S.pdf",width=2.5,height=3)
ggplot(bcDist, aes(x=type, y=bc)) + geom_boxplot(notch=T, outlier.colour="red", outlier.shape=16,
                                                 outlier.size=3) + geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()


# # ---- Dissimilarity comparisons Metazoa ------
# 
# ssu.intraSub1 = otusR18S.ab.ra.m[mdR$Grab==1 & mdR$PooledExtracts==1,]
# ssu.intraSub2 = otusR18S.ab.ra.m[mdR$Grab==2 & mdR$PooledExtracts==1,]
# ssu.intraSub3 = otusR18S.ab.ra.m[mdR$Grab==3 & mdR$PooledExtracts==1,]
# 
# is1.bc = vegdist(ssu.intraSub1)
# is2.bc = vegdist(ssu.intraSub2)
# is3.bc = vegdist(ssu.intraSub3)
# ssu.intraSub.bc = c(is1.bc,is2.bc,is3.bc)
# summary(ssu.intraSub.bc) # .18 -- .90, median .46
# 
# ## Samples with 5 pooled replicats from grabs 1--3
# intraGrab1 = otusR18S.ab.ra.m[mdR$Grab==1 & mdR$PooledExtracts==5,]
# ig1.bc = vegdist(intraGrab1)
# summary(ig1.bc) #.55 -- .94 med .89
# 
# intraGrab2 = otusR18S.ab.ra.m[mdR$Grab==2 & mdR$PooledExtracts==5,]
# ig2.bc = vegdist(intraGrab2)
# summary(ig2.bc) #.75 -- .86 med .80
# 
# intraGrab3 = otusR18S.ab.ra.m[mdR$Grab==3 & mdR$PooledExtracts==5,]
# ig3.bc = vegdist(intraGrab3)
# summary(ig3.bc) #.74 -- .87 med .78
# 
# igAll.bc = c(unlist(ig1.bc),unlist(ig2.bc),unlist(ig3.bc)) 
# summary(igAll.bc)
# # .55 -- .94, med .81
# 
# ## Is BC significantly higher between pooled sub-grabs compared to between extractions?
# wilcox.test(igAll.bc, ssu.intraSub.bc, alternative="greater")
# # -> Yes (p=4E-7)
# 
# ## Samples with 5 pooled replicats from all grabs sorted by grab
# intraGrabAll = rbind(intraGrab1, intraGrab2, intraGrab3)
# intraGrabAll.bc = as.matrix(vegdist(intraGrabAll))
# betweenGrabs.bc = c(intraGrabAll.bc[c(6:15),c(1:5)],intraGrabAll.bc[c(11:15),c(6:10)])
# summary(betweenGrabs.bc) # .55 -- .94, med .87
# 
# ## Is BC significantly higher between grabs compared to between intra-grab pooled
# wilcox.test(betweenGrabs.bc, igAll.bc, alternative="greater")
# # -> Yes but close call (p=0.013)
# 
# require(vioplot)
# 
# boxplot(ssu.intraSub.bc, igAll.bc, betweenGrabs.bc, names=c("Btw. extracts same sub-sample",
#                                                         "Btw. pooled intra-grabs",
#                                                         "Btw. grabs"),las=2, notch=T,col="grey")
# 
# require(vioplot)
# vioplot(ssu.intraSub.bc, igAll.bc, betweenGrabs.bc, names=c("Btw. extracts same sub-sample",
#                                                         "Btw. pooled intra-grabs",
#                                                         "Btw. grabs"),las=2)

# --- Alpha div. comparison ----

## Sample rarefaction curve

sac_ssu.intraSub1=specaccum(ssu.intraSub1)
sac_ssu.intraSub2=specaccum(ssu.intraSub2)
sac_ssu.intraSub3=specaccum(ssu.intraSub3)

col1=rgb(255, 255, 0, max = 255, alpha = 125, names = "y1")
col2=rgb(0, 255, 0, max = 255, alpha = 125, names = "g1")
col3=rgb(0, 0, 255, max = 255, alpha = 125, names = "b1")

plot(sac_ssu.intraSub1, ci.type="polygon", ci.col=col1,ylim=range(600,1000))
lines(sac_ssu.intraSub2, ci.type="polygon", ci.col=col2)
lines(sac_ssu.intraSub3, ci.type="polygon", ci.col=col3)
legend("bottomright",col=c("yellow","green","blue"),pch=16,
       legend=c("Grab 1","Grab 2","Grab 3"))


sac_ig1=specaccum(intraGrab1)
sac_ig2=specaccum(intraGrab2)
sac_ig3=specaccum(intraGrab3)

plot(sac_ig1, ci.type="polygon", ci.col=col1,ylim=range(600,1000))
lines(sac_ig2, ci.type="polygon", ci.col=col2)
lines(sac_ig3, ci.type="polygon", ci.col=col3)
legend("bottomright",col=c("yellow","green","blue"),pch=16,
       legend=c("Grab 1","Grab 2","Grab 3"))

## ---- Pooled grab diversity, 5 pooled extracts only ----

otusR18S.ab.temp = otusR18S.ab[mdR$PooledExtracts==5,]
md.temp = mdR[mdR$PooledExtracts==5,]
otusR18S.ab.pooled = mergeOTUTable(otusR18S.ab.temp, md.temp, "Grab")
specnumber(otusR18S.ab.pooled)
# 1   2   3 
# 959 933 933 

dim(otusR18S.ab.pooled)
#1030

#Exclusive to G1: 37
sum(otusR18S.ab.pooled[1,]>0 & otusR18S.ab.pooled[2,]==0 & otusR18S.ab.pooled[3,]==0)

#Exclusive to G2: 25
sum(otusR18S.ab.pooled[1,]==0 & otusR18S.ab.pooled[2,]>0 & otusR18S.ab.pooled[3,]==0)

#Exclusive to G3: 18
sum(otusR18S.ab.pooled[1,]==0 & otusR18S.ab.pooled[2,]==0 & otusR18S.ab.pooled[3,]>0)

# Shared G1_G2: 35
sum(otusR18S.ab.pooled[1,]>0 & otusR18S.ab.pooled[2,]>0 & otusR18S.ab.pooled[3,]==0)

# Shared G1_G3: 42
sum(otusR18S.ab.pooled[1,]>0 & otusR18S.ab.pooled[2,]==0 & otusR18S.ab.pooled[3,]>0)

# Shared G2_G3: 28
sum(otusR18S.ab.pooled[1,]==0 & otusR18S.ab.pooled[2,]>0 & otusR18S.ab.pooled[3,]>0)

# Shared by all: 845
sum(otusR18S.ab.pooled[1,]>0 & otusR18S.ab.pooled[2,]>0 & otusR18S.ab.pooled[3,]>0)

diversity(otusR18S.ab.pooled)
# 1        2        3 
# 4.508055 4.840315 4.637109 

## --- AMBI distribution ------- 

require(BBI)

table(names(otusR18S.ab) == row.names(tax.18S.ab))

mt = data.frame(Taxonpath = tax.18S.ab$classification, t(otusR18S.ab))

bbi.18S = BBI(mt) # 117 found, 913 not found
write.csv(bbi.18S$BBI,"BIs_18S.csv")
bis = read.csv("BIs_18S.csv",row.names=1)
plot(bis$Shannon~divR.18S$H)

plot (bis$AMBI~bis.COI$AMBI) # no trend
