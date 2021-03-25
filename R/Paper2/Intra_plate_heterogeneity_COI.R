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
md$COI_ID  = as.character(md$COI_ID)
summary(md)

otus.COI.all = read.delim("COI/CREST_LULU/SWARM_table_curated.tsv",
                          row.names=1,header=T,sep="\t")

names(otus.COI.all) = gsub("\\.","-",names(otus.COI.all))
table(md$COI_ID %in% names(otus.COI.all))  # all 36

otus.COI = as.data.frame(t(otus.COI.all[,md$COI_ID]))
table(row.names(otus.COI) == md$COI_ID)
row.names(otus.COI) = row.names(md)
table(row.names(otus.COI) == row.names(md))

# Classification df
bestTx = array(dim=dim(otus.COI.all)[1])
for (i in 1:length(bestTx)) bestTx[i] <- tail(unlist(strsplit(as.character(otus.COI.all$classification[i]), 
                                                              split=";", 
                                                              fixed=TRUE)), 1)
tax.COI = data.frame(row.names=row.names(otus.COI.all), 
                     classification=otus.COI.all$classification, bestTx=bestTx)


rowSums(otus.COI)
# Only blanks have few reads
writeDivStats("QCDiversity_stats_intermediate_COI.csv", otudistro=otus.COI)
#writeDivStats("QCDiversity_stats_intermediate_Met_COI.csv", otudistro=otus.COI[,grep("Metazoa",tax.COI$classification)])

# ------ Contamination control ------

dim(otus.COI) # we start with 19,870 OTUs
sum(otus.COI) # and 3,314,749 reads

## Cross contaminants
otus.COI = filterCrossContaminants2(otus.COI,100)
sum(otus.COI) #3,314,331 (400 reads less)

## Filter data by abundance before removing unclassified (test purpose only)
# We require a minimum of 3 expected reads assuming average constant relative OTU abundance
# 
# minAb = 3 / 82175 #(3.6E-5)
# otusCOI.ab = dropRareByAvgAbundance(otus.COI, minAb)
# dim(otusCOI.ab)
# dim(otusCOI.ab)/dim(otus.COI) #1653 OTUs remain at this cutoff, or 8%
# otus.COI = otusCOI.ab
#tax.COI = tax.COI[names(otusCOI.ab),]


## Remove all taxa not classified to at least phylum level

phylumOrBetter = rep(NA,dim(tax.COI)[1])
for (i in c(1:length(phylumOrBetter))){
  tt = strsplit(as.character(tax.COI[i,"classification"]),";")
  phylumOrBetter[i] = (length(tt[[1]]) > 4)
}
sum(phylumOrBetter) #6052 (569 for pre-abundance filtered)
otus.COI = otus.COI[,phylumOrBetter]
tax.COI = tax.COI[phylumOrBetter,]
table(row.names(tax.COI) == names(otus.COI))
rowSums(otus.COI)

# Limit to metazoa
otus.COI = otus.COI[,grep("Metazoa",tax.COI$classification)]
tax.COI = tax.COI[grep("Metazoa",tax.COI$classification),]

dim(otus.COI) #3806 (379 for pre-abundance filtered)
sum(otus.COI) #1,192,688 (ca 2/3 removed) (1,112,437 for pre-ab. filtered)
rowSums(otus.COI)


##  Make barchart and check abundant OTUs for contamination (just for control)

otus.named = decostand(otus.COI,method="total")
names(otus.named) = paste(names(otus.named),tax.COI[names(otus.named),]$bestTx)

grouping_info<-data.frame(row.names=row.names(md), md$Type)
pdf("img/COI/SV_chart_COI_w_controls.pdf",width=17, height = 10)
taxaplot(30,grouping_info,otus.named,grouping_info)
dev.off()

## Prevalence based filtering based on negative controls (decontam)

om = as.matrix(otus.COI)
predicted_contaminants = isContaminant(om, neg=(md$Type=="C"))#,
#batch=md$PlateCOI)

summary(predicted_contaminants) #-> 2 predicted contaminants
table(row.names(predicted_contaminants) == names(otus.COI)) # Yes
sc = (!is.na(predicted_contaminants$p) & predicted_contaminants$p<=.05) # 2
sigCont = predicted_contaminants[sc,]
row.names(sigCont)
sigCont$Tax = tax.COI$bestTx[row.names(tax.COI) %in% row.names(sigCont)]
summary(as.factor(sigCont$Tax))
# Polychaeta and Spionida - do not remove becaues cannot be real contaminants

# cont = (names(otus.COI) %in% row.names(sigCont))
# otus.COI = otus.COI[,-cont]

# Remove the blank and mock samples from further analysis
otusR.COI = otus.COI[md$Type == "S",]
sum(otusR.COI) #119,277 (ca 300 reads removed with the 6 control samples)
otusR.COI = otusR.COI[,colSums(otusR.COI)>0] 
dim(otusR.COI) #Retains 3800 OTUs (only 5 removed) (373 for pre-ab filtered)

mdR = droplevels(md[md$Type == "S",])
table(row.names(otusR.COI) == row.names(mdR))

## Taxon based filter (replace with plankton based??)

for (taxonDel in c("No hits","Insecta","Mammalia","Aves","Myxini",
                   "Arachnida","Actinopterygii","Chondrichthyes","Collembola",
                   "Onychophorida")) {
  toDel =  row.names(tax.COI)[grep(taxonDel,tax.COI$classification)]
  if (length(toDel)>0){
    otusR.COI = otusR.COI[,!(names(otusR.COI) %in% toDel)]
    print(paste("Removing ",taxonDel))
    print(dim(otusR.COI)[2])
    print(sum(otusR.COI))
  }
}
# 1] "Removing  Insecta"
# [1] 2197
# [1] 823994
# [1] "Removing  Mammalia"
# [1] 2163
# [1] 822935
# [1] "Removing  Myxini"
# [1] 2161
# [1] 822850
# [1] "Removing  Arachnida"
# [1] 1802
# [1] 796708
# [1] "Removing  Actinopterygii"
# [1] 1769
# [1] 792010
# [1] "Removing  Collembola"
# [1] 1768
# [1] 792000
# [1] "Removing  Onychophorida"
# [1] 1759
# [1] 788501

# (191 for pre-ab filtered)

# ---- Relative abundance ------

otus.COI.ra.all = decostand(otusR.COI, method="total")

summary(rowSums(otusR.COI)) 
# 3 -- 112 thousand reads per sample (2.4k -- 111 for pre-ab)

## Testing purpose - pre-abundance filtering instead of here:
#otusRCOI.ab = otusR.COI (191 OTUs)

# For dissimilarity, we require a minimum of 3 expected reads assuming average constant relative OTU abundance
minAb = 3 / min(rowSums(otusR.COI)) #ca 1E-3
otusRCOI.ab = dropRareByAvgAbundance(otusR.COI, minAb)
dim(otusRCOI.ab)
dim(otusRCOI.ab)/dim(otusR.COI) #108 OTUs remain at this cutoff, or 6%


otusRCOI.ab.ra = decostand(otusRCOI.ab,method="total")
#otusRCOI.ab.H = decostand(otusRCOI.ab,method="hellinger")

tax.COIR = tax.COI[names(otusR.COI),]
tax.COI.ab = tax.COI[names(otusRCOI.ab),]
table(row.names(tax.COI.ab) == names(otusRCOI.ab))

## OTU plot with abundance filtered  
otus.named = decostand(otusRCOI.ab,method="total")
names(otus.named) = paste(names(otus.named),tax.COI[names(otus.named),]$bestTx)
grouping_info<-data.frame(row.names=row.names(mdR), paste("Grab",mdR$Grab))
pdf("img/COI/SV_chart_COI_abFiltered.pdf",width=10, height = 5)
taxaplot(30,grouping_info,otus.named,grouping_info)
dev.off()

#-------- Diversity stats --------

writeDivStats("Diversity_filtered_COI_spatial_heterogeneity.csv", otusR.COI)
writeDivStats("Diversity_filtered_COI_spatial_heterogeneity_ab.csv", otusRCOI.ab)


divR.COI = read.csv("Diversity_filtered_COI_spatial_heterogeneity_ab.csv",row.names=1)
divR.COI = divR.COI[,c("Reads","Richness","Rarefied.richness","H","J",
                       "Chao1")]

cor.test(divR.COI$Reads,divR.COI$Richness)
cor.test(divR.COI$Reads,divR.COI$Rarefied.richness)
plot(divR.COI$Reads,divR.COI$Rarefied.richness) 
#Strong negative trend - richness would be better

# ---- NMDS --------

nmds = metaMDS(otusRCOI.ab.ra) #Stress .17, converges w max resid 6.5E-5
ordiplot(nmds,type="none")#,xlim=c(-.7,.5),ylim=c(-.5,.5))
points(nmds,pch=as.numeric(mdR$PooledExtracts),col=as.numeric(mdR$Grab),lwd=2,cex=1.2)
text(nmds,labels=paste("G",mdR$Grab,":",mdR$Sediment,sep=""),pos=1,cex=.7,col=as.numeric(mdR$Grab))
legend("bottomleft",pch=c(1,5,1,1,1), col=c("grey","grey","black","red","green"),
       legend=c("Individual extract", "5 pooled extracts", "Grab 1", "Grab 2", "Grab 3"),
       cex=.8)

# ---- Dissimilarity comparisons ------

intraSub1 = otusRCOI.ab.ra[mdR$Grab==1 & mdR$PooledExtracts==1,]
intraSub2 = otusRCOI.ab.ra[mdR$Grab==2 & mdR$PooledExtracts==1,]
intraSub3 = otusRCOI.ab.ra[mdR$Grab==3 & mdR$PooledExtracts==1,]

is1.bc = vegdist(intraSub1)
is2.bc = vegdist(intraSub2)
is3.bc = vegdist(intraSub3)
intraSub.bc = c(is1.bc,is2.bc,is3.bc)
summary(intraSub.bc) # .03 -- .76, median .49

## Samples with 5 pooled replicats from grabs 1--3
intraGrab1 = otusRCOI.ab.ra[mdR$Grab==1 & mdR$PooledExtracts==5,]
ig1.bc = vegdist(intraGrab1)
summary(ig1.bc) #.26 -- .99 med .93

intraGrab2 = otusRCOI.ab.ra[mdR$Grab==2 & mdR$PooledExtracts==5,]
ig2.bc = vegdist(intraGrab2)
summary(ig2.bc) #.68 -- .98 med .86

intraGrab3 = otusRCOI.ab.ra[mdR$Grab==3 & mdR$PooledExtracts==5,]
ig3.bc = vegdist(intraGrab3)
summary(ig3.bc) #.71 -- .90 med .87

igAll.bc = c(unlist(ig1.bc),unlist(ig2.bc),unlist(ig3.bc)) 
summary(igAll.bc)
# .26 -- .99, med .89

## Is BC significantly higher between pooled sub-grabs compared to between extractions?
wilcox.test(igAll.bc, intraSub.bc, alternative="greater")
# -> Yes (p=5E-13)

## Samples with 5 pooled replicats from all grabs sorted by grab
intraGrabAll = rbind(intraGrab1, intraGrab2, intraGrab3)
intraGrabAll.bc = as.matrix(vegdist(intraGrabAll))
betweenGrabs.bc = c(intraGrabAll.bc[c(6:15),c(1:5)],intraGrabAll.bc[c(11:15),c(6:10)])
summary(betweenGrabs.bc) # .26 -- .995, med .92

## Is BC significantly higher between grabs compared to between intra-grab pooled
wilcox.test(betweenGrabs.bc, igAll.bc, alternative="greater")
# -> No (p=.3)

boxplot(intraSub.bc, igAll.bc, betweenGrabs.bc, names=c("Btw. extracts same sub-sample",
                                                        "Btw. pooled intra-grabs",
                                                        "Btw. grabs"),las=2, notch=T,col="grey")

require(vioplot)
vioplot(intraSub.bc, igAll.bc, betweenGrabs.bc, names=c("Btw. extracts same sub-sample",
                                                        "Btw. pooled intra-grabs",
                                                        "Btw. grabs"),las=2)

# --- Alpha div. comparison ----

## Sample rarefaction curve
sac_intrasub1=specaccum(intraSub1)
sac_intrasub2=specaccum(intraSub2)
sac_intrasub3=specaccum(intraSub3)

col1=rgb(255, 255, 0, max = 255, alpha = 125, names = "y1")
col2=rgb(0, 255, 0, max = 255, alpha = 125, names = "g1")
col3=rgb(0, 0, 255, max = 255, alpha = 125, names = "b1")

plot(sac_intrasub1, ci.type="polygon", ci.col=col1,ylim=range(30,70))
lines(sac_intrasub2, ci.type="polygon", ci.col=col2)
lines(sac_intrasub3, ci.type="polygon", ci.col=col3)
legend("bottomright",col=c("yellow","green","blue"),pch=16,
       legend=c("Grab 1","Grab 2","Grab 3"))


sac_ig1=specaccum(intraGrab1)
sac_ig2=specaccum(intraGrab2)
sac_ig3=specaccum(intraGrab3)

plot(sac_ig1, ci.type="polygon", ci.col=col1,ylim=range(30,90))
lines(sac_ig2, ci.type="polygon", ci.col=col2)
lines(sac_ig3, ci.type="polygon", ci.col=col3)
legend("bottomright",col=c("yellow","green","blue"),pch=16,
       legend=c("Grab 1","Grab 2","Grab 3"))



## ---- Pooled grab diversity, 5 pooled extracts only ----

otusRCOI.ab.temp = otusRCOI.ab[mdR$PooledExtracts==5,]
md.temp = mdR[mdR$PooledExtracts==5,]
otusRCOI.ab.pooled = mergeOTUTable(otusRCOI.ab.temp, md.temp, "Grab")
specnumber(otusRCOI.ab.pooled)
# 1   2   3 
# 83 84   87


dim(otusRCOI.ab.pooled)
# 108

#Exclusive to G1: 7
sum(otusRCOI.ab.pooled[1,]>0 & otusRCOI.ab.pooled[2,]==0 & otusRCOI.ab.pooled[3,]==0)

#Exclusive to G2: 10
sum(otusRCOI.ab.pooled[1,]==0 & otusRCOI.ab.pooled[2,]>0 & otusRCOI.ab.pooled[3,]==0)

#Exclusive to G3: 8
sum(otusRCOI.ab.pooled[1,]==0 & otusRCOI.ab.pooled[2,]==0 & otusRCOI.ab.pooled[3,]>0)

# Shared G1_G2: 4
sum(otusRCOI.ab.pooled[1,]>0 & otusRCOI.ab.pooled[2,]>0 & otusRCOI.ab.pooled[3,]==0)

# Shared G1_G3: 9
sum(otusRCOI.ab.pooled[1,]>0 & otusRCOI.ab.pooled[2,]==0 & otusRCOI.ab.pooled[3,]>0)

# Shared G2_G3: 7
sum(otusRCOI.ab.pooled[1,]==0 & otusRCOI.ab.pooled[2,]>0 & otusRCOI.ab.pooled[3,]>0)

# Shared by all: 63
sum(otusRCOI.ab.pooled[1,]>0 & otusRCOI.ab.pooled[2,]>0 & otusRCOI.ab.pooled[3,]>0)


require(BBI)
# ass.COI = read.table("COI/CREST_LULU/All_Assignments.tsv",
#                      sep="\t", header=T,row.names=3, check.names = F)
# 
# names(ass.COI) = gsub("\\.","-",names(ass.COI))
# 
# asra.COI =  decostand(as.data.frame(t(ass.COI[,as.character(mdR$COI_ID)])),
#                        method="total")
# mt = ass.COI[,c("Taxonpath",as.character(mdR$COI_ID))]
#ass.COI[,c("Taxonpath",as.character(mdR$COI_ID))]
#names(mt) = c("Taxonpath",row.names(mdR))
#mt$Taxonpath = row.names(mt)

table(names(otusRCOI.ab) == row.names(tax.COI.ab))

mt = data.frame(Taxonpath = tax.COI.ab$classification, t(otusRCOI.ab))

bbi.COI = BBI(mt) # 71 found, 37 not found
write.csv(bbi.COI$BBI,"BIs_COI.csv")
bis = read.csv("BIs_COI.csv",row.names=1)
plot(bis$Shannon~divR.COI$H) # Different!

mtPA =  data.frame(Taxonpath = tax.COI.ab$classification, 
                   t(decostand(otusRCOI.ab,method="pa")))
bbi.COIpa = BBI(mtPA) # 71 found, 37 not found
write.csv(bbi.COIpa$BBI,"BIs_COI_pa.csv")

