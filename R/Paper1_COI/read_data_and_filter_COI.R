setwd("/home/alanzen/projects/Metamon/WP2/SWARM_WP2_all_20200708/")

require(vegan)
require(VennDiagram)

source('../R/heterogeneity_rarefaction_functions.R')
source('~/kode/R/filtering.R')
source('~/kode/R/diversity.r')
source('~/kode/R/correlationTests.r')
source('~/kode/R/taxaplot.R')
source('~/kode/R/mergeOTUTable.R')


# Reads metadata that is common for COI and COI
md.all = read.csv(file="metadata.csv",header=T,row.names=1)

# Limit to those samples that have COI data
md.coiR = md.all[md.all$COI,]

# order metadata alphabetically
md.coiR = md.coiR[order(md.coiR$COI_ID),]

# Read OTUs
otus.all.COI = read.delim("COI/CREST_LULU/SWARM_table_curated.tsv",
                          row.names=1,header=T,sep="\t")
# Fix names
names(otus.all.COI) = gsub("\\.","-",names(otus.all.COI))

# order OTU samples alphabetically as well
otus.all.COI = otus.all.COI[,order(names(otus.all.COI))]

# Put taxonomic classification in a separate vector and remove
tax.COI=otus.all.COI$classification
otus.all.COI = otus.all.COI[,-1] # Remove classification column

# ------ Make in silico pooled grp for comparison and get proper OTU table -------
otus.all.COI$X01_10_isPool = rowSums(otus.all.COI[,md.coiR$Group=="ExtractionRep" & md.coiR$KitHom=="Precellys1"])
otus.all.COI$X11_20_isPool = rowSums(otus.all.COI[,md.coiR$Group=="ExtractionRep" & md.coiR$KitHom=="Vortex_0.5g"])
otus.all.COI$X21_24_isPool = rowSums(otus.all.COI[,md.coiR$Group=="ExtractionRep" & md.coiR$KitHom=="Vortex_5g"])
otus.all.COI$X55_64_isPool = rowSums(otus.all.COI[,md.coiR$Group=="ExtractionRep" & md.coiR$KitHom=="Precellys2"])
otus.all.COI$XPCR = rowSums(otus.all.COI[,!is.na(md.coiR$Group) & md.coiR$Group=="PCRRepX58"])

table(names(otus.all.COI) == md.coiR$COI_ID)

names(otus.all.COI) = row.names(md.coiR)

otus.t.COI = as.data.frame(t(otus.all.COI))
rowSums(otus.t.COI)

rowSums(otus.t.COI[md.coiR$Group=="Controls",])
rowSums(otus.t.COI[md.coiR$Group=="InSilicoPool",])
sum(otus.t.COI) #10,311,906 reads
dim(otus.t.COI) #29,217 OTUs


# Write diversity statistics before taxonomic filtering and contaminant removal

writeDivStats("QCDiversity_stats_intermediate_COI.csv", otus.t.COI)
writeDivStats("QCDiversity_stats_intermediate_Met_COI.csv", otus.t.COI[,grep("Metazoa",tax.COI)])

# ------ Contamination control ------

# Remove all taxa not classified to at least phylum level
# 
phylumOrBetter = rep(NA,29217)
for (i in c(1:length(phylumOrBetter))){
  tt = strsplit(as.character(tax.COI[i]),";")
  phylumOrBetter[i] = (length(tt[[1]]) > 4)
}
sum(phylumOrBetter) #8248
tax.COI = tax.COI[phylumOrBetter]

otus.COI = otus.t.COI[,phylumOrBetter]

# Limit to metazoa
otus.COI = otus.COI[,grep("Metazoa",tax.COI)]
tax.COI = tax.COI[grep("Metazoa",tax.COI)]
dim(otus.COI) #5095 (679 order or better)
sum(otus.COI) #2864040 (1,437,043  order or better)
sum(otus.COI)/sum(otus.t.COI) # 28%
rowSums(otus.COI)

for (taxonDel in c("Mammalia","Aves","Myxini",
                   "Actinopterygii", "Insecta","Arachnida")){#, "No hits"
  toDel = grep(taxonDel,tax.COI)
  if (length(toDel)>0){
    otus.COI = otus.COI[,-toDel]
    tax.COI = tax.COI[-toDel]
    print(paste("Removing ",taxonDel))
    print(dim(otus.COI)[2])
    print(sum(otus.COI))
  }
}

# 35 Mammal OTUs (~3k reads), 1 bird (60 r), 2 Myxini, 56 fishes (16k reads)
# 2148 insect OTUs (733k reads!) and 476 arachnida left in - removing them leaves 
# 2043060 reads (order filtering: 1,408,134 reads)

otus.COI = filterCrossContaminants2(otus.COI,100)
sum(otus.COI) #2042626

otus.COI.all.ra = decostand(otus.COI, method="total")

otus.COI.blank.ra = otus.COI.all.ra[md.coiR$Group!="InSilicoPool",
                            colSums(otus.COI.all.ra[md.coiR$Group=="Controls",]) > 0]
otus.COI.blank = otus.COI[md.coiR$Group!="InSilicoPool",
                  colSums(otus.COI.all.ra[md.coiR$Group=="Controls",]) > 0]
md.noIS = md.coiR[md.coiR$Group!="InSilicoPool",]
dim(otus.COI.blank)
# 35 OTUs

tx = tax.COI[colSums(otus.COI.all.ra[md.coiR$Group=="Controls",]) > 0]
# get the last rank predicted
bestTx = array(dim=length(tx))
for (i in 1:length(tx)) bestTx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", 
                                                          fixed=TRUE)), 1)
names(otus.COI.blank.ra) = make.names(paste(names(otus.COI.blank),bestTx,sep="_"))
names(otus.COI.blank) = make.names(paste(names(otus.COI.blank),bestTx,sep="_"))

pdf("../img/controls/COI/OTUs_in_controls_all_barplot_RA.pdf", height=9, width=30)
taxaplot(35, data.frame(row.names=row.names(otus.COI.blank), md.noIS$Experiment),otus.COI.blank.ra)
dev.off()

pdf("../img/controls/COI/OTUs_in_controls_all_barplot_AA.pdf", height=9, width=30)
taxaplot(35, data.frame(row.names=row.names(otus.COI.blank), md.noIS$Experiment),otus.COI.blank)
dev.off()

md.noIS$Ctrl = ifelse(md.noIS$Group=="Controls", as.vector(md.noIS$Experiment)
                      , "Sample")
otus.COI.blank.pool = mergeOTUTable(otus.COI.blank, md.noIS, by="Ctrl")


otus.COI.blank.pool.ra = decostand(otus.COI.blank.pool, method="total")
pdf("../img/controls/COI/OTUs_in_controls_pooled_barplot_RA.pdf", height=9, width=20)
taxaplot(35, data.frame(row.names=row.names(otus.COI.blank.pool), row.names(otus.COI.blank.pool)),
         otus.COI.blank.pool.ra)
dev.off()

onlyInBlank = otus.COI.blank.pool.ra[,otus.COI.blank.pool["Sample",] == 0]# & otus.COI.blank.pool["Mock",] == 0]
dim(onlyInBlank) #5 otus incl. 2 marine gastropods and a shrimp. Where from?

otus.COI.blank.pool.ra.f = otus.COI.blank.pool.ra[,!(names(otus.COI.blank.pool.ra) %in% names(onlyInBlank))]
dim(otus.COI.blank.pool.ra.f) #30

xr = as.data.frame(t(otus.COI.blank.pool.ra.f))

xr$XBlankRatio = xr$Xblank / xr$Sample
xr$PCRBlankRatio = xr$PCRBlank / xr$Sample

sum(xr$XBlankRatio>1)  
summary(xr$XBlankRatio[xr$XBlankRatio>1]) #Median 1526, 3Q 4492, in total 12 of 30

xTop <- xr[order(xr$XBlankRatio, decreasing = T),]

sum(xr$PCRBlankRatio>1)
summary(xr$PCRBlankRatio[xr$PCRBlankRatio>1]) #Median 135, 3Q 9106, in total 7 of 30

pTop <- xr[order(xr$PCRBlankRatio, decreasing = T),]

pdf("../img/controls/COI/BLEX_ratio_order_only.pdf",width=5,height=4)
bp <- barplot(xTop$XBlankRatio[c(12:1)], horiz = T, xlab="Abundance in Extracion controls / real samples", 
              main=paste("Contamination candidates (Extraction)"), las=2, cex.names = 0.6,
              log="x", xlim=c(50,20000))
text(45, bp, labels=row.names(xTop)[c(12:1)], pos=4, cex=0.6, col="red")
dev.off()

pdf("../img/controls/COI/PCRBlank_ratio.pdf",width=5,height=3)
bp <- barplot(pTop$PCRBlankRatio[c(7:1)], horiz = T, xlab="Abundance in PCR blanks / real samples", 
              main=paste("Contamination candidates (PCR blanks)"), las=2, cex.names = 0.6,
              log="x", xlim=c(1,15000))
text(1, bp, labels=row.names(pTop)[c(7:1)], pos=4, cex=0.6, col="red")
dev.off()

# No suspected contaminants!

md.COI.nc = md.coiR[md.coiR$Include,]

otusR.COI.nc = otus.COI[md.coiR$Include,]
## Remove mislabeled in vitro pools
otusR.COI.nc = otusR.COI.nc[md.COI.nc$Group!="PooledReps" | md.COI.nc$Experiment!="WP2HomXTest",]
md.COI.nc = md.COI.nc[md.COI.nc$Group!="PooledReps" | md.COI.nc$Experiment!="WP2HomXTest",]

taxR.COI = tax.COI[colSums(otusR.COI.nc)>0]
otusR.COI.nc = otusR.COI.nc[,colSums(otusR.COI.nc)>0] 
otusR.COIf = otusR.COI.nc[rowSums(otusR.COI.nc)>2000,]
md.COIf = md.COI.nc[rowSums(otusR.COI.nc)>2000,]
dim(otusR.COIf) #Retains 2366 OTUs and 82 samples of 86, missing PCR 9 and 10
row.names(otusR.COIf)

table(row.names(otusR.COI.nc) == row.names(md.COI.nc))
write.table(t(otusR.COI.nc), "COI/SWARM_table_filtered_met.csv", sep="\t", quote=F)

rownames(otusR.COIf)

# -------- Diversity calc ----------

plot(specnumber(otusR.COIf)~rowSums(otusR.COIf),xlab="Reads",ylab="Richness", log="xy")
writeDivStats("Diversity_filtered_COI_WP2_order_only_min2000reads_metazoa_phylum.csv", otusR.COIf,
              rr_cutoff=3300, rarefyForShannon = TRUE)

divR.COIf = read.csv("Diversity_filtered_COI_WP2_order_only_min2000reads_metazoa_phylum.csv",row.names=1)
divR.COIf = divR.COIf[,c("Reads","Richness","Rarefied.richness","H","J",
               "Chao1")]

otus.COI.ra.all = decostand(otusR.COIf, method="total")
otus.COI.H.all = decostand(otusR.COIf, method="hell")

minAb = 3 / min(rowSums(otusR.COIf)) #min is 3319 so min. is 0.0009
otusRCOI.ab = dropRareByAvgAbundance(otusR.COIf, minAb)
dim(otusRCOI.ab)/dim(otusR.COIf) #89 of 2366 (4%)


# ------- Taxonomy --------

mdAllButH4 = md.COI.nc[md.COI.nc$Experiment!="IntraGrabAndXHom",]

taxa.all.COI = read.table("COI/CREST_Filtered/Relative_Abundance.tsv",sep="\t",
                          header=T,row.names=3, check.names = F)

ass.all.COI = read.table("COI/CREST_Filtered/All_Assignments.tsv",sep="\t",
                         header=T,row.names=3, check.names = F)

ass.allbutH4 = decostand(as.data.frame(t(ass.all.COI[,row.names(mdAllButH4)])),
                         method="total")

tra_all = taxa.all.COI[,c("Rank","Taxonpath",row.names(mdAllButH4))]

grouping_info<-data.frame(row.names=row.names(mdAllButH4), mdAllButH4$Group)

ranks = data.frame(rank=c("phylum","class","order","family","genus","species"),
                   levels=c(15,25,25,25,25,25))

summary(tra_all$Rank)
# 18 phyla, 43 classes, 75 orders, 88 families, 
# 93 genera, 77 spp.

for (i in c(1:6)){
  r=as.character(ranks$rank[i])
  leveltax.COIa = as.data.frame(t(tra_all[tra_all$Rank==r,-c(1:2)]))
  print(paste(mean(rowSums(leveltax.COIa))*100,"% classified at rank",r))
  pdf(paste("../img/COI/taxon_barplots_all/",r,".pdf",sep=""),height=9,width=20)
  taxaplot(ranks$levels[i],grouping_info,leveltax.COIa)
  dev.off()
}

pdf("../img/COI/taxon_barplots_all/Assignments_all_butH4.pdf",height=9,width=20)
taxaplot(30,grouping_info,ass.allbutH4)
dev.off()

## Grabs 1-3 and 5

mdExtra = mdAllButH4[mdAllButH4$Experiment=="GrabEXandPCRRep",]

tra_og = taxa.all.COI[,c("Rank","Taxonpath",row.names(mdExtra))]

grouping_info<-data.frame(row.names=row.names(mdExtra), paste("Grab",mdExtra$Grab))


for (i in c(1:6)){
  r=as.character(ranks$rank[i])
  leveltax.COIa = as.data.frame(t(tra_og[tra_og$Rank==r,-c(1:2)]))
  print(paste(mean(rowSums(leveltax.COIa))*100,"% classified at rank",r))
  pdf(paste("../img/COI/taxon_barplots_all/",r,"_pre_v_post_PCR_pooling.pdf",sep=""),
      height=7,width=10)
  taxaplot(ranks$levels[i],grouping_info,leveltax.COIa)
  dev.off()
}

ass.extra = decostand(as.data.frame(t(ass.all.COI[,row.names(mdExtra)])),
                      method="total")
pdf("../img/COI/taxon_barplots_all/Assignments_pre_v_post_PCR_pooling.pdf",
    height=6,width=9)
taxaplot(30,grouping_info,ass.extra)
dev.off()

