# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon1/")
# 

require(vegan)

source('R/utils/heterogeneity_rarefaction_functions.R')
source('R/utils/filtering.R')
source('R/utils/diversity.r')
source('R/utils/correlationTests.r')
source('R/utils/taxaplot.R')
source('R/utils/mergeOTUTable.R')


# Check agreement metadata vs OTU samples
table(row.names(otusR.18S) == row.names(mdR.18S))

# ---------- Make subset of data for only "experiment 1" (H,T,P) --------------

## We are only concerned with grab 4, including in silico and in vitro reps
# these will be named mdX.18S, otusX.18S etc.

selection = (mdR.18S$Experiment=="WP2HomXTest" | mdR.18S$Experiment=="WP2PCRRep")

mdX.18S = mdR.18S[selection,]
otusX.18S = otusR.18S[selection,]
otus.18S.pmX = otus.18S.pm[selection,]
otusX.18S.ab = otusR18S.ab[selection,]

## Make diversity stats and subsample using 40k for consistency with Fig 4
writeDivStats("Diversity_18S_40k_xreps.csv", otusX.18S, rr_cutoff=40000, rarefyForShannon = TRUE)
divX.18S.40k = read.csv("Diversity_18S_40k_xreps.csv",row.names=1)
divX.18S.40k = divX.18S.40k[,c(1,3,5,6)]

divX.18S.met = divR.18S.met[selection,]

## Make subset of only individual extracion replicates (not including PCR reps)
divReps.18S = divX.18S.40k[mdX.18S$Group=="ExtractionRep",]
mdReps.18S = droplevels(mdX.18S[mdX.18S$Group=="ExtractionRep",])

# ---------- Compare alpha diversity between treatment  groups -------

## Data for table 3
pc1xOut = divX.18S[c(2:10),]
pc1xOut.met = divX.18S.met[c(2:10),]
dx = rbind(pc1xOut,divX.18S)
dx.met = rbind(pc1xOut.met,divX.18S.met)
xo =mdX.18S[c(2:10),]
xo$KitHom="PC1 x outlier"
mdT3 =rbind(xo,mdX.18S)
for (kh in c("Precellys1","PC1 x outlier","Precellys2","Vortex_5g","Vortex_0.5g")){
  dt = dx[mdT3$KitHom==kh & mdT3$Group=="ExtractionRep",]
  dm = dx.met[mdT3$KitHom==kh & mdT3$Group=="ExtractionRep",]
  print(kh)
  print(summary(dm$Richness))
  print(summary(dm/dt))
  print(paste("Met read share 2SD:", 2*sd(dm$Reads/dt$Reads)))
  print(paste("Met S 2SD:", 2*sd(dm$Richness)))
  print(paste("Met S share 2SD:",2*sd(dm$Richness/dt$Richness)))
}

divX.18S.met[mdX.18S$Group=="PooledReps",]
divX.18S.met[mdX.18S$Group=="PooledReps",] / divX.18S[mdX.18S$Group=="PooledReps",]

## Add points for pools
divIVPool = divX.18S.40k[mdX.18S$Group=="PooledReps",][c(1,4,2,3),]
mdIVPool = droplevels(mdX.18S[mdX.18S$Group=="PooledReps",])[c(1,4,2,3),]

divISPool = divX.18S.40k[mdX.18S$Group=="InSilicoPool",][c(1,4,2,3),]
mdISPool = droplevels(mdX.18S[mdX.18S$Group=="InSilicoPool",])[c(1,4,2,3),]

pdf("img/18S/Fig2_RarefiedRichnessA.pdf",width=3.5,height=5)
kolore=c("red", "green", "blue", "grey")
boxplot(divReps.18S$Rarefied.richness ~ mdReps.18S$KitHom,
        notch=T, beside=T,las=2, ylim=c(550,1600),
        col=kolore)

points(factor(mdISPool$KitHom), divISPool$Rarefied.richness, col=mdISPool$Colour, pch=3)
points(factor(mdIVPool$KitHom), divIVPool$Rarefied.richness, col=mdIVPool$Colour, pch=5)
dev.off()


pdf("img/18S/Fig2_ShannonB.pdf",width=3.5,height=5)
kolore=c("red", "green", "blue", "grey")
boxplot(divReps.18S$H ~ mdReps.18S$KitHom,
        notch=T, beside=T,las=2, #ylim=c(600,1700),
        col=kolore)

points(factor(mdISPool$KitHom), divISPool$H, col=mdISPool$Colour, pch=3)
points(factor(mdIVPool$KitHom), divIVPool$H, col="black", pch=5)
dev.off()

## Test statistical significance:

wilcox.test(divReps.18S$Rarefied.richness[mdReps.18S$KitHom=="Precellys1"],
            divReps.18S$Rarefied.richness[mdReps.18S$KitHom=="Precellys2"])
# RS PC2 > PC1: p=3E-4

wilcox.test(divReps.18S$H[mdReps.18S$KitHom=="Precellys1"],
            divReps.18S$H[mdReps.18S$KitHom=="Precellys2"])
# NS for Shannon


wilcox.test(divReps.18S$H[mdReps.18S$KitHom=="Vortex_0.5g"],
            divReps.18S$H[mdReps.18S$KitHom=="Vortex_5g"])

wilcox.test(divReps.18S$Rarefied.richness[mdReps.18S$KitHom=="Vortex_0.5g"],
            divReps.18S$Rarefied.richness[mdReps.18S$KitHom=="Vortex_5g"])

# Vortex 5g> 0.5g (p=0.002 for both H' and rarefied richness)

# -------- Compare dissimilarity and make NMDS  --------

otusX.18S.abH = decostand(otusX.18S.ab, method="hell")

## Hierarchical clustering
hc = hclust(vegdist(otusX.18S.abH))
plot(hc, labels=paste(row.names(mdX.18S),mdX.18S$KitHom), cex=.7)
# Interstingly has in vitro and in silico pools very close to each other but w longer branch length

## Non-metric dimensional scaling
nmds = metaMDS(otusX.18S.abH)
ordiplot(nmds,type="none",xlim=c(-.5,.5),ylim=c(-.2,.6))
points(nmds,pch=as.numeric(mdX.18S$Group),col=mdX.18S$Colour,lwd=2,cex=1.1)


for (i in c(1:4)){
  
  kh = unique(mdX.18S$KitHom)[i]
  ordispider(nmds,groups=mdX.18S$KitHom,label=F,cex=.7,
             col=unique(mdX.18S$Colour)[i],
             show.groups=kh)
}

# ------ Degree of heterogeneity ------

pspcl1.hBC = unlist(vegdist(otusX.18S.abH[mdX.18S$KitHom=="Precellys1" & mdX.18S$Group=="ExtractionRep",]))
pspcl2.hBC = unlist(vegdist(otusX.18S.abH[mdX.18S$KitHom=="Precellys2" & mdX.18S$Group=="ExtractionRep",]))
psv.hBC = unlist(vegdist(otusX.18S.abH[mdX.18S$KitHom=="Vortex_0.5g" & mdX.18S$Group=="ExtractionRep",]))
pmv.hBC = unlist(vegdist(otusX.18S.abH[mdX.18S$KitHom=="Vortex_5g" & mdX.18S$Group=="ExtractionRep",]))
pcrrep.hBC= unlist(vegdist(otusX.18S.abH[mdX.18S$Group=="PCRRepX58",]))


boxplot(pspcl1.hBC,pspcl2.hBC,psv.hBC,pmv.hBC,pcrrep.hBC,
        names=c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCR Rep"),
        col=c("red","green","blue","darkgrey","cyan"),notch=T,
        las=2)

wilcox.test(pspcl1.hBC,pspcl2.hBC) #p=.16
wilcox.test(pspcl1.hBC,psv.hBC) #p=6E-7
wilcox.test(pcrrep.hBC,pmv.hBC) #p=1E-7

require(vioplot)
vioplot(pspcl1.hBC,pspcl2.hBC,psv.hBC,pmv.hBC,pcrrep.hBC,
        names=c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCR Rep"),las=2)

# -------- Compare dissimilarity Metazoa only  --------

otus.18S.pmX.ab = dropRareByAvgAbundance(otus.18S.pmX,3 / min(rowSums(otusR.18S)))
otus.18S.pmX.ab.H = decostand(otus.18S.pmX.ab, method="hell")

nmds = metaMDS(otus.18S.pmX.ab.H)
ordiplot(nmds,type="none",xlim=c(-1.5,1.5),ylim=c(-.5,.5))
points(nmds,pch=as.numeric(mdX.18S$Group),col=mdX.18S$Colour,lwd=2,cex=1.1)


for (i in c(1:4)){
  kh = unique(mdX.18S$KitHom)[i]
  ordispider(nmds,groups=mdX.18S$KitHom,label=F,cex=.7,
             col=unique(mdX.18S$Colour)[i],
             show.groups=kh)
  }

pspcl1.m.hBC = unlist(vegdist(otus.18S.pmX.ab.H[mdX.18S$KitHom=="Precellys1" & mdX.18S$Group=="ExtractionRep",]))
pspcl2.m.hBC = unlist(vegdist(otus.18S.pmX.ab.H[mdX.18S$KitHom=="Precellys2" & mdX.18S$Group=="ExtractionRep",]))
psv.m.hBC = unlist(vegdist(otus.18S.pmX.ab.H[mdX.18S$KitHom=="Vortex_0.5g" & mdX.18S$Group=="ExtractionRep",]))
pmv.m.hBC = unlist(vegdist(otus.18S.pmX.ab.H[mdX.18S$KitHom=="Vortex_5g" & mdX.18S$Group=="ExtractionRep",]))
pcrrep.m.hBC= unlist(vegdist(otus.18S.pmX.ab.H[mdX.18S$Group=="PCRRepX58",]))


boxplot(pspcl1.m.hBC,pspcl2.m.hBC,psv.m.hBC,pmv.m.hBC,pcrrep.m.hBC,
        names=c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCR Rep"),
        col=c("red","green","blue","darkgrey","cyan"),notch=T,
        las=2)

# Compare met to total euk BC
wilcox.test(c(pspcl1.m.hBC,pspcl2.m.hBC,psv.m.hBC),c(pspcl1.hBC,pspcl2.hBC,psv.hBC)) #<2E-16
wilcox.test(pspcl1.m.hBC,pspcl1.hBC) #<2E-16
wilcox.test(pspcl2.m.hBC,pspcl2.hBC) #<2E-16
wilcox.test(psv.m.hBC,psv.hBC) #<2E-16

# Compare treatment groups
wilcox.test(pspcl1.hBC,pspcl2.hBC) #p=1E-11
wilcox.test(pspcl1.hBC,psv.hBC) #p1E-10
wilcox.test(pspcl2.hBC,pmv.hBC) #p2E-5
wilcox.test(pcrrep.hBC,pmv.hBC) #p=1E-7


# ------ taxonomy --------

taxa.all.18S = read.table("18S/CREST_Filtered/Relative_Abundance.tsv",sep="\t",
                          header=T,row.names=3, check.names = F)

taxa.ass.18S = read.table("18S/CREST_Filtered/All_Assignments.tsv",sep="\t",
                          header=T,row.names=3, check.names = F)

tra = taxa.all.18S[,c("Rank","Taxonpath",row.names(mdX.18S))] #or mdR.18S
tass = taxa.ass.18S[,c("Rank","Taxonpath",row.names(mdX.18S))] 

groups=paste(mdX.18S$Group,mdX.18S$KitHom)
groups[c(49:53)] = "In silico pools"
groups[c(25:27,48)] = "In vitro pools"
grouping_info<-data.frame(row.names=row.names(mdX.18S), groups)

ranks = data.frame(rank=c("superkingdom","kingdom","phylum","class","order","family","genus"),
                   levels=c(8,11,25,25,25,25,25))

summary(tra$Rank)
# 13 superkingdoms, 26 kingdoms, 102 phyla, 204 classes, 259 orders, 164 families, 
# 150 genera

pdf("img/18S/taxon_barplots_xheterogeneity/Assignmets.pdf",height=9,width=19)
assignedTaxaRA = decostand(as.data.frame(t(tass[,-c(1:2)])),method="total")
taxaplot(25,grouping_info,assignedTaxaRA)
dev.off()

for (i in c(1:7)){
  r=as.character(ranks$rank[i])
  leveltaxa = as.data.frame(t(tra[tra$Rank==r,-c(1:2)]))
  print(paste(mean(rowSums(leveltaxa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots_xheterogeneity/",r,".pdf",sep=""),height=9,width=19)
  taxaplot(ranks$levels[i],grouping_info,leveltaxa)
  dev.off()
}
# [1] "69.1676110799428 % classified at rank phylum"
# [1] "65.3210395377336 % classified at rank class"
# [1] "48.1136103973555 % classified at rank order"
# [1] "19.9565541019454 % classified at rank family"
#18 % at genus

metKingdom = tra["Metazoa (Animalia)",c(3:55)]

tra.met = tra[grep("Metazoa",tra$Taxonpath),] 
summary(tra.met$Rank)
# 22 phyla, 50 classes, 76 orders, 12 fams, 1 genus, no spp.
ranks = data.frame(rank=c("phylum","class","order","family"),
                   levels=c(18,18,18,12))

for (i in c(1:4)){
  r=as.character(ranks$rank[i])
  leveltax.18Sa = as.data.frame(t(tra.met[tra.met$Rank==r,-c(1:2)])) / t(metKingdom)
  print(paste(mean(rowSums(leveltax.18Sa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots_xheterogeneity/",r,"_Metazoa.pdf",sep=""),height=7,width=19)
  taxaplot(ranks$levels[i],grouping_info,leveltax.18Sa)
  dev.off()
}

# ------ Comparing taxa between methods --------

tax.18Sa.comp = tra[,c(T,T,md.comp$Experiment=="WP2HomXTest" & md.comp$Xreps==1)]
dim(tax.18Sa.comp)
tax.18Sa.comp = tax.18Sa.comp[rowMeans(tax.18Sa.comp[,-c(1:2)])>1E-4,]
dim(tax.18Sa.comp) #417 tax.18Sa
mdc = droplevels(md.comp[md.comp$Experiment=="WP2HomXTest"  & md.comp$Xreps==1,])

for (i in c(1:7)){
  r=as.character(ranks$rank[i])
  leveltax.18Sa = as.data.frame(t(tax.18Sa.comp[tax.18Sa.comp$Rank==r,-c(1:2)]))
  printANOVA(mdc[,c("Xkit","KitHom")], 
             leveltax.18Sa, a = .05/dim(leveltax.18Sa)[2],
             imgDir="img/18S/tax.18Sa_systematic_diffs")
}


# -------- Resampling by combining N replicates  -------


md.n.18SR=mdX.18S
grp = c("Precellys1","Precellys2","Vortex_0.5g","Vortex_5g","PCRRep")
colour = c("red","green","blue","grey","cyan")

noSamples = 5
noReps = 10
iterations = 100
totalRichness = 9818
metRichness = 1481
##
nonMetRichness = totalRichness - metRichness

res = array(dim=c(noSamples,noReps,iterations,totalRichness))
res_p = array(dim=c(noSamples,noReps,iterations,nonMetRichness))
res_m = array(dim=c(noSamples,noReps,iterations,metRichness))
reads = 40000

resAll = repeatedSubsamplingPlate1(otus.18S=otusX.18S, res=res, res_p=res_p, res_m=res_m, 
                          noSample=noSamples, 
                          noReps=noReps, iterations=iterations, reads=reads, 
                          classification=tax.18SR)

ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads.svg", resAll@total, 
       c(600,1600), colour, grp)#)
ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads_non-metazoa.svg", 
       resAll@nonmet, c(500,1400), colour, grp)#)
ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads_metazoa.svg", 
       resAll@met, c(70,300), colour, grp)#)


ppHPlot("img/18S/Pooling_effect_no_Reps_H_40kreads.svg", resAll@total, 
        c(3.5,5), colour, grp)# customOrder=c(1,5,2:4))
ppHPlot("img/18S/Pooling_effect_no_Reps_H_40kreads_non-metazoa.svg", 
        resAll@nonmet, c(3.9,4.6), colour, grp)#)
ppHPlot("img/18S/Pooling_effect_no_Reps_H_40kreads_metazoa.svg", 
        resAll@met, c(1,4), colour, grp)#)


# --------- No. of samples effect removing rare -------

avgTaxonAbundances = colMeans(decostand(otusX.18S,method="total"))

minAb=0.001

otusAb = otusX.18S[,avgtax.18SonAbundances >= minAb]
classificationAb = tax.18SR[avgTaxonAbundances >= minAb]

dim(otusAb) #Retains 127 of 9818
length(grep("Metazoa",classificationAb))

totalRichness = 127
metRichness = 35
nonMetRichness = totalRichness - metRichness

res = array(dim=c(noSamples,noReps,iterations,totalRichness))
res_p = array(dim=c(noSamples,noReps,iterations,nonMetRichness))
res_m = array(dim=c(noSamples,noReps,iterations,metRichness))

resAb = repeatedSubsamplingPlate1(otus.18S=otusAb, res=res, res_p=res_p, res_m=res_m, 
                            noSample=noSamples, 
                            noReps=noReps, iterations=iterations, reads=39750, 
                            classification=classificationAb)

ppPlot("img/18S/Pooling_effect_no_Reps_S_39kreads_min01p.svg", resAb@total, 
       c(95,125), colour, grp)

ppPlot("img/18S/Pooling_effect_no_Reps_S_39kreads_min01p_met.svg", resAb@met, 
       c(13,33), colour, grp)

ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads_min01p_nonmet.svg", resAb@nonmet,
       c(80,93), colour, grp)

quantile(avgTaxonAbundances, probs=seq(0,1,0.1))
# 90% percentile is 3.1E-5
# The 10% most abundant
minAb = 3.1E-5
otusAb = otusX.18S[,avgTaxonAbundances >= minAb]
classificationAb = classification[avgTaxonAbundances >= minAb]

dim(otusAb) #Retains 977 (so yes, 10%, but ~20% of metazoa i.e. 257)
length(grep("Metazoa",classificationAb))

totalRichness = 977
metRichness = 257
nonMetRichness = totalRichness - metRichness

res = array(dim=c(noSamples,noReps,iterations,totalRichness))
res_p = array(dim=c(noSamples,noReps,iterations,nonMetRichness))
res_m = array(dim=c(noSamples,noReps,iterations,metRichness))

resAb = repeatedSubsamplingPlate1(otus.18S=otusAb, res=res, res_p=res_p, res_m=res_m, 
                            noSample=noSamples, 
                            noReps=noReps, iterations=iterations, reads=40000, 
                            classification=classificationAb)

ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads_10pMostAb.svg", resAb@total, 
       c(400,850), colour, grp)

ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads_10pMostAb_met.svg", resAb@met, 
       c(50,170), colour, grp)

ppPlot("img/18S/Pooling_effect_no_Reps_S_40kreads_10pMostAb_nonmet.svg", resAb@nonmet,
       c(340,690), colour, grp)

# ----- Beyond 40k rarefaction -----

calcRarefaction = function(otus.18S, steplength=5000, max=50000){
  steps = max/steplength
  rarefied = matrix(nrow=steps,ncol=15)
  # Not doing std. error as was done in Lanzen et al 2016 on rarefied as not comparable
  stdev = matrix(nrow=steps,ncol=6) 
  mdRep = md.n.18SR[(md.n.18SR$Group=="ExtractionRep" | md.n.18SR$Group=="PCRRepX58"),]
  otusRep = otus.18S[(md.n.18SR$Group=="ExtractionRep" | md.n.18SR$Group=="PCRRepX58"),]
  otusRealPooled = otus.18S[(md.n.18SR$Group=="PooledReps"),]
  otusISPooled = otus.18S[md.n.18SR$Group=="InSilicoPool",]
  
  reads = c(rowSums(otusRep),rowSums(otusRealPooled),rowSums(otusISPooled))
  
  for (i in 1:steps){
    print(paste(i,"out of",steps))
    rarefied[i,1] = stdev[i,1] = i*steplength
    rarefied[i,1] = stdev[i,1] = i*steplength 
    
    rsrep = rarefy(otusRep,sample=i*steplength)
    rspooled = rarefy(otusRealPooled,sample=i*steplength)
    rsis = rarefy(otusISPooled,sample=i*steplength)
    
    for (j in 1:44){
      if (reads[j] < i*steplength)  rsrep[j] = NaN
    }
    for (j in 45:48){
      if (reads[j] < i*steplength) rspooled[j-44] = NaN
    }
    for (j in 49:53){
      if (reads[j] < i*steplength) rsis[j-48] = NaN
    }
    
    for (s in c(1:noSamples)) {
      rarefied[i,s+1] = mean(rsrep[mdRep$KitHom==grp[s]],na.rm=T)
      stdev[i,s+1] = sd(rsrep[mdRep$KitHom==grp[s]],na.rm=T)
    }
    rarefied[i,7:10] = rspooled
    
    rarefied[i,11:15] = rsis
  }
  return(new("Rarefied", rarefied=rarefied, stdev=stdev))
}


# ----- Plot beyond 40k color -------
max=250000
sl=2000
cr = calcRarefaction(otus.18S.pmX,steplength=sl, max=max)

plotRare(rarefied=cr@rarefied, stdev=cr@stdev,
         file="img/18S/rarefaction_beyond40k.svg",height=6,width=7,
         ylims=range(300,3000),xlims=c(sl,max),diff=100,colour=colour)

plotRare(rarefied=cr@rarefied, stdev=cr@stdev,
         file="img/18S/rarefaction_beyond40k_zoom.svg",height=6,width=7,
         ylims=c(300,1800),xlims=c(steplength,50000),diff=50,colour=colour)

# ---- Plots (b/w one for each dataset) ----

plotBW = function(raf, rpPos=c(7,10,8,9,NA), rsPos=c(11,14,12,13,15), 
                  lims=range(0,300)){
  for (s in c(1:5)){
    plot(raf@rarefied[c(1:20),1],raf@rarefied[c(1:20),s+1],type="l",ylim=lims,
         xlab="Reads per replicate",ylab="OTU richness",lwd=2,
         main=toString(grp[s]))
    arrows(raf@rarefied[c(1:20),1], raf@rarefied[c(1:20),s+1]-
             raf@stdev[c(1:20),s+1], raf@rarefied[c(1:20),1], 
           raf@rarefied[c(1:20),s+1]+raf@stdev[c(1:20),s+1], 
           length=0.05, angle=90, code=0)
    lines(raf@rarefied[c(1:20),1],raf@rarefied[c(1:20),rpPos[s]],lwd=2,col="grey")
    lines(raf@rarefied[c(1:20),1],raf@rarefied[c(1:20),rsPos[s]],lwd=2,col="grey",lty=2)
   
    pause = readline()
  }
}

plotBW(cr)

