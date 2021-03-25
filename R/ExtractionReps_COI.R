setwd("/home/alanzen/projects/Metamon/WP2/SWARM_WP2_all_20200708")

require(vegan)
require(VennDiagram)

source('../R/utils/heterogeneity_rarefaction_functions.R')
source('../R/utils/filtering.R')
source('../R/utils/diversity.r')
source('../R/utils/correlationTests.r')
source('../R/utils/taxaplot.R')
source('../R/utils/mergeOTUTable.R')

load(".RData")

# Check agreement metadata vs OTU samples
table(row.names(otusR.COIf) == row.names(md.COIf))

# ---------- Make subset of data for only "experiment 1" (H,T,P) --------------

## We are only concerned with grab 4, including in silico and in vitro reps
# these will be named mdX.COI, otusX.COI etc.

selection = (md.COIf$Experiment=="WP2HomXTest" | md.COIf$Experiment=="WP2PCRRep")

mdX.COI = droplevels(md.COIf[selection,])
summary(mdX.COI)
# We have all ten / four reps for x reps. 2 PCR reps missing

otusX.COI = otusR.COIf[selection,]
otusX.COI.ab = otusRCOI.ab[selection,]


## Make diversity stats and subsample using 40k for consistency with Fig 4
writeDivStats("Diversity_COI_3300r_xreps.csv", otusX.COI, rr_cutoff=3300, rarefyForShannon = TRUE)
divX.COI = read.csv("Diversity_COI_3300r_xreps.csv",row.names=1)
divX.COI = divX.COI[,c(1,3,5,6)]



## Make subset of only individual extracion replicates (not including PCR reps)
divReps.COI = divX.COI[mdX.COI$Group=="ExtractionRep",]
summary(divReps.COI) # _> min is 3319 reads
mdReps.COI = droplevels(mdX.COI[mdX.COI$Group=="ExtractionRep",])

# ---------- Compare alpha diversity between treatment  groups -------

printANOVA1Factor(mdReps.COI$KitHom, divReps.COI, .05)
# Precellys2 has significantly higher rarefied richness vs all others

## Add points for pools
divISPool = divX.COI[mdX.COI$Group=="InSilicoPool",][c(1,4,2,3),]
mdISPool =mdX.COI[mdX.COI$Group=="InSilicoPool",][c(1,4,2,3),]

pdf("../img/COI/Fig2_RarefiedRichnessA_allMet.pdf",width=3.5,height=5)
kolore=c("red", "green", "blue", "grey")
boxplot(divReps.COI$Rarefied.richness ~ mdReps.COI$KitHom,
        notch=T, beside=T,las=2, ylim=c(50,350),
        col=kolore)

points(divISPool$Rarefied.richness, 
       col=c("red","green","blue","black"), pch=3, lw=2)
dev.off()


pdf("../img/COI/Fig2_ShannonB_allMet.pdf",width=3.5,height=5)
boxplot(divReps.COI$H ~ mdReps.COI$KitHom,
        notch=T, beside=T,las=2, #ylim=c(600,1700),
        col=kolore)

points(factor(mdISPool$KitHom), divISPool$H, col=c("white","green","blue","white"), pch=3, lw=2)
dev.off()


## Test statistical significance:

wilcox.test(divReps.COI$Rarefied.richness[mdReps.COI$KitHom=="Precellys1"],
            divReps.COI$Rarefied.richness[mdReps.COI$KitHom=="Precellys2"])
# RS PC2 > PC1: p=2E-4

wilcox.test(divReps.COI$H[mdReps.COI$KitHom=="Precellys1"],
            divReps.COI$H[mdReps.COI$KitHom=="Precellys2"])
# NS for Shannon


wilcox.test(divReps.COI$H[mdReps.COI$KitHom=="Vortex_0.5g"],
            divReps.COI$H[mdReps.COI$KitHom=="Vortex_5g"])
#NS

# -------- Compare dissimilarity and make NMDS  --------

otusX.COI.abH = decostand(otusX.COI.ab, method="hell")

## Hierarchical clustering
hc = hclust(vegdist(otusX.COI.abH))
plot(hc, labels=paste(row.names(mdX.COI),mdX.COI$KitHom), cex=.7)
# Interstingly has in vitro and in silico pools very close to each other but w longer branch length

nmds = metaMDS(otusX.COI.abH)
ordiplot(nmds,type="none")#,xlim=c(-1.2,1),ylim=c(-1,.5))
points(nmds,pch=as.numeric(droplevels(mdX.COI$Group))+1,col=mdX.COI$Colour,lwd=2,cex=1)
legend("topright",pch=as.numeric(unique(mdX.COI$Group)),legend=unique(mdX.COI$Group),cex=.8)
legend("topleft",pch=1,col=unique(mdX.COI$Colour)[c(1,4,2,3,5)],
       legend=unique(mdX.COI$KitHom)[c(1,4,2,3,5)],cex=.8)
#text(nmds,pos=1,cex=.7,col=as.numeric(md.COIf$KitHom))


for (i in c(1:5)){
  
  kh = unique(mdX.COI$KitHom)[i]
  ordispider(nmds,groups=mdX.COI$KitHom,label=F,cex=.7,
             col=unique(mdX.COI$Colour)[i],
             show.groups=kh)
}

# ------ Degree of heterogeneity ------

pspcl1.hBC = unlist(vegdist(otusX.COI.abH[mdX.COI$KitHom=="Precellys1" & mdX.COI$Group=="ExtractionRep",]))
pspcl2.hBC = unlist(vegdist(otusX.COI.abH[mdX.COI$KitHom=="Precellys2" & mdX.COI$Group=="ExtractionRep",]))
psv.hBC = unlist(vegdist(otusX.COI.abH[mdX.COI$KitHom=="Vortex_0.5g" & mdX.COI$Group=="ExtractionRep",]))
pmv.hBC = unlist(vegdist(otusX.COI.abH[mdX.COI$KitHom=="Vortex_5g" & mdX.COI$Group=="ExtractionRep",]))
pcrrep.hBC= unlist(vegdist(otusX.COI.abH[mdX.COI$Group=="PCRRepX58",]))


boxplot(pspcl1.hBC,pspcl2.hBC,psv.hBC,pmv.hBC,pcrrep.hBC,
        names=c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCR Rep"),
        col=c("red","green","blue","darkgrey","cyan"),notch=T,
        las=2)

wilcox.test(pspcl1.hBC,pspcl2.hBC) #p=2E-6
wilcox.test(psv.hBC,pspcl2.hBC) #p=2E-15
wilcox.test(pcrrep.hBC,pspcl2.hBC) #p=0.0009

require(vioplot)
vioplot(pspcl1.hBC,pspcl2.hBC,psv.hBC,pmv.hBC,pcrrep.hBC,
        names=c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCR Rep"),las=2)


# ------ taxonomy --------
  
taxa.all.COI = read.table("COI/CREST_Filtered/Relative_Abundance.tsv",sep="\t",
                          header=T,row.names=3, check.names = F)


taxa.ass.COI = read.table("COI/CREST_Filtered/All_Assignments.tsv",sep="\t",
                          header=T,row.names=3, check.names = F)

tra = taxa.all.COI[,c("Rank","Taxonpath",row.names(mdX.COI))] #or md.COIf
tass = taxa.ass.COI[,c("Rank","Taxonpath",row.names(mdX.COI))] 

groups=paste(mdX.COI$Group,mdX.COI$KitHom)
groups[c(43:47)] = "In silico pools"
grouping_info<-data.frame(row.names=row.names(mdX.COI), groups)

ranks = data.frame(rank=c("phylum","class","order","family","genus","species"),
                   levels=c(18,25,25,25,25,25))

summary(tra$Rank)
# 18 phyla, 43 classes, 75 orders, 88 families, 93 genera

pdf("../img/COI/taxon_barplots_xheterogeneity/Assignmets.pdf",height=9,width=19)
assignedTaxaRA = decostand(as.data.frame(t(tass[,-c(1:2)])),method="total")
taxaplot(25,grouping_info,assignedTaxaRA)
dev.off()

#93% at class, 72 at order and 41 at spp


for (i in c(1:6)){
  r=as.character(ranks$rank[i])
  leveltaxa = as.data.frame(t(tra[tra$Rank==r,-c(1:2)]))
  print(paste(mean(rowSums(leveltaxa))*100,"% classified at rank",r))
  pdf(paste("../img/COI/taxon_barplots_xheterogeneity/",r,".pdf",sep=""),height=9,width=19)
  taxaplot(ranks$levels[i],grouping_info,leveltaxa)
  dev.off()
}
# [1] "100 % classified at rank phylum"
# [1] "93.6222814839119 % classified at rank class"
# [1] "72.5167796396647 % classified at rank order"
# [1] "56.3228131511792 % classified at rank family"
# [1] "56.260922989907 % classified at rank genus"
# [1] "40.8210150570109 % classified at rank species"

# -------- Resampling by combining N replicates  -------

## 
md.n.COIR=mdX.COI
grp = c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCRRep")
colour = c("red","green","blue","grey","cyan")
otusS = rrarefy(otusX.COI, 3319)

noSamples = 5
noReps = 10
iterations = 100
totalRichness = 2366
metRichness = 2366
##
nonMetRichness = totalRichness - metRichness

res = array(dim=c(noSamples,noReps,iterations,totalRichness))
reads = 3300  

resAll = repeatedSubsamplingPlate1(otus=otusS, res=res, 
                          noSample=noSamples, 
                          noReps=noReps, iterations=iterations, reads=reads, 
                          classification=tax.COIR, metadata=md.n.COIR ,PCRRepSamples=8)

ppPlot("../img/COI/Pooling_effect_COI_3300reads.svg", resAll@total, 
       c(50,400), colour, grp)#)

ppHPlot ("../img/COI/Pooling_effect_COI_3300reads_H_allMet.svg", resAll@total, 
        c(1,4), colour, grp)# customOrder=c(1,5,2:4))


