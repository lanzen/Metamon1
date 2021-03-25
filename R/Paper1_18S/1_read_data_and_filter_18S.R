# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon1/")
# 

source('R/utils/filtering.R')
source('R/utils/diversity.r')
source('R/utils/taxaplot.R')
source('R/utils/mergeOTUTable.R')
source('R/utils/correlationTests.r')

require(vegan)

# ------ Read data --------

# Read metadata that is common for COI and 18S
md.all = read.csv(file="SWARM_WP2_all_20200708/metadata.csv", header=T, row.names=1)

# Limit analysis to those samples that have 18S data
md.ssuR = md.all[md.all$X18S,] 

# order metadata alphabetically
md.ssuR = md.ssuR[order(md.ssuR$SSU_ID),]

# Read 18S OTUs
otus.18S.all = read.delim("SWARM_WP2_all_20200708/18S/CREST_LULU/SWARM_table_curated.tsv",
                          row.names=1,header=T,sep="\t")
# Fix names
names(otus.18S.all) = gsub("\\.","-",names(otus.18S.all))

# order OTU samples alphabetically as well
otus.18S.all = otus.18S.all[,order(names(otus.18S.all))]

# Put taxonomic classification in a separate vector and remove
tax.18S=otus.18S.all$classification
otus.18S.all = otus.18S.all[,-1]

# ------ Make in silico pooled groups for comparison and get proper OTU table -------
otus.18S.all$X01_10_isPool = rowSums(otus.18S.all[,md.ssuR$Group=="ExtractionRep" & md.ssuR$KitHom=="Precellys1"])
otus.18S.all$X11_20_isPool = rowSums(otus.18S.all[,md.ssuR$Group=="ExtractionRep" & md.ssuR$KitHom=="Vortex_0.5g"])
otus.18S.all$X21_24_isPool = rowSums(otus.18S.all[,md.ssuR$Group=="ExtractionRep" & md.ssuR$KitHom=="Vortex_5g"])
otus.18S.all$X55_64_isPool = rowSums(otus.18S.all[,md.ssuR$Group=="ExtractionRep" & md.ssuR$KitHom=="Precellys2"])
otus.18S.all$XPCR = rowSums(otus.18S.all[,!is.na(md.ssuR$Group) & md.ssuR$Group=="PCRRepX58"])

# Check name agreement metadata SSU ID vs OTU sample names (all should be True)
table(names(otus.18S.all) == md.ssuR$SSU_ID)

# If all agree, proceed, and change sample names to metadata ones
names(otus.18S.all) = row.names(md.ssuR)

# Reformat
otus.18S.t = as.data.frame(t(otus.18S.all))

# Check number of reads per sample, sample group  and total
summary(rowSums(otus.18S.t))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9   90500  104600  131900  117900 1094000 


# Write diversity statistics before taxonomic filtering and contaminant removal

writeDivStats("SWARM_WP2_all_20200708/QCDiversity_stats_intermediate_18S.csv", otudistro=otus.18S.t)
writeDivStats("SWARM_WP2_all_20200708/QCDiversity_stats_intermediate_Met_18S.csv", otudistro=otus.18S.t[,grep("Metazoa",tax.18S)])

# ------ Contamination control ------

dim(otus.18S.t) #we start with 10191 OTUs
otus.18S = otus.18S.t

for (taxonDel in c("No hits","Insecta","Mammalia","Aves","Myxini",
                   "Arachnida","Actinopterygii")) {
  toDel = grep(taxonDel,tax.18S)
  if (length(toDel)>0){
    otus.18S = otus.18S[,-toDel]
    tax.18S = tax.18S[-toDel]
    print(paste("Removing ",taxonDel))
    print(dim(otus.18S)[2])
    print(sum(otus.18S))
  }
}

# Filter cross contamination in a manner analogous to UNCROSS (R Edgar)
# removing all occurences <1% of average occurence
otus.18S = filterCrossContaminants2(otus.18S,100)
sum(otus.18S) #13042977
sum(otus.18S[md.ssuR$Group!="InSilicoPool",]) #8,928,662 

# Transform into relative abundance
otus.18S.all.ra = decostand(otus.18S, method="total")

# Investigate composition of blanks
otus.18S.blank.ra = otus.18S.all.ra[md.ssuR$Group!="InSilicoPool",
                            colSums(otus.18S.all.ra[md.ssuR$Group=="Controls",]) > 0]
otus.18S.blank = otus.18S[md.ssuR$Group!="InSilicoPool",
                  colSums(otus.18S.all.ra[md.ssuR$Group=="Controls",]) > 0]
md.noIS = md.ssuR[md.ssuR$Group!="InSilicoPool",]
dim(otus.18S.blank)
# 254 OTUs in blanks

tx = tax.18S[colSums(otus.18S.all.ra[md.ssuR$Group=="Controls",]) > 0]
# get the last rank predicted
bestTx = array(dim=length(tx))
for (i in 1:length(tx)) bestTx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", 
                                                          fixed=TRUE)), 1)
names(otus.18S.blank.ra) = make.names(paste(names(otus.18S.blank),bestTx,sep="_"))
names(otus.18S.blank) = make.names(paste(names(otus.18S.blank),bestTx,sep="_"))

md.noIS$Ctrl = ifelse(md.noIS$Group=="Controls", as.vector(md.noIS$Experiment)
                      , "Sample")

otus.18S.blank.pool = mergeOTUTable(otus.18S.blank, md.noIS, by="Ctrl")
otus.18S.blank.pool.ra = decostand(otus.18S.blank.pool, method="total")

onlyInBlank = otus.18S.blank.pool.ra[,otus.18S.blank.pool["Sample",] == 0]# & otus.18S.blank.pool["Mock",] == 0]
dim(onlyInBlank) #7 otus

otus.18S.blank.pool.ra.f = otus.18S.blank.pool.ra[,!(names(otus.18S.blank.pool.ra) %in% names(onlyInBlank))]
dim(otus.18S.blank.pool.ra.f) #247

xr = as.data.frame(t(otus.18S.blank.pool.ra.f))

xr$XBlankRatio = xr$Xblank / xr$Sample
xr$PCRBlankRatio = xr$PCRBlank / xr$Sample

sum(xr$XBlankRatio>1) #77 (31%)
summary(xr$XBlankRatio[xr$XBlankRatio>1]) #Median 13, 3Q 69

sum(xr$PCRBlankRatio>1) #99
summary(xr$PCRBlankRatio[xr$PCRBlankRatio>1]) #Median 8, 3Q 24

xTop <- xr[order(xr$XBlankRatio, decreasing = T),]
pTop <- xr[order(xr$PCRBlankRatio, decreasing = T),]


# Top candidate (57466) is Terebellida polychaete and the 2nd and 4th are also polychate (Spionida)
# 3rd is Paranoidae Annelid thus top ones are unlikely to be true contaminants.
# However, 2513 (possibly Lobulomyces sp.) is suspicious as 95% similar to soil / freshwater fungi
# and 1944 and 1070 (freshwater / rumen fungus), and 998 (soil)

pdf("img/controls/18S/PCRBlank_ratio.pdf",width=7,height=8)
bp <- barplot(pTop$PCRBlankRatio[c(40:1)], horiz = T, xlab="Abundance in PCR blanks / real samples", 
              main=paste("Contamination candidates (PCR blanks)"), las=2, cex.names = 0.6,
              log="x", xlim=c(10,20000))
text(10, bp, labels=row.names(pTop)[c(40:1)], pos=4, cex=0.6, col="red")
dev.off()

# SWARM_2, top fungi, is also similar to unclass. eukaryotes from marine sediments (Svalbard) and 
# look like possible Chytridiomycota, that may be marine (https://www.nature.com/articles/srep30120)
# SWARM_31 looks more like Choanoflagellida. SWARM_163 again like Svalbard sediment

# Remove suspected contaminants
true_contaminants = c("SWARM_2513","SWARM_1944","SWARM_1070","SWARM_998") 
cont = (names(otus.18S) %in% true_contaminants)
otus.18S = otus.18S[,-cont]
tax.18S = tax.18S[-cont]
sum(otus.18S) #13,042,974 ca 20k removed

# Remove the blank and mock samples from further analysis
otusR.18S = otus.18S[md.ssuR$Include,]
dim(otusR.18S) #10029

tax.18SR = tax.18S[colSums(otusR.18S)>0]
otusR.18S = otusR.18S[,colSums(otusR.18S)>0] 
dim(otusR.18S) #Retains 10,022
mdR.18S = droplevels(md.ssuR[md.ssuR$Include,])

# Control that cleaned data coincides with metadata names
table(row.names(otusR.18S) == row.names(mdR.18S))

write.table(t(otusR.18S), "SWARM_WP2_all_20200708/18S/SWARM_table_filtered.csv", quote=F, sep="\t")

otus.18S.pm = otusR.18S[,grep("Metazoa",tax.18SR)]
otus.18S.pnm = otusR.18S[,-grep("Metazoa",tax.18SR)]
dim(otus.18S.pm) # only 1538 are metazoan

writeDivStats("SWARM_WP2_all_20200708/Diversity_filtered_18S_WP2_wReps.csv", otusR.18S, rr_cutoff=40000, rarefyForShannon = TRUE)
writeDivStats("SWARM_WP2_all_20200708/Diversity_filtered_18S_WP2_met.csv", otus.18S.pm)

divR.18S = read.csv("SWARM_WP2_all_20200708/Diversity_filtered_18S_WP2_wReps.csv",row.names=1)
divR.18S = divR.18S[,c("Reads","Richness","Rarefied.richness","H","J",
               "Chao1")]

divR.18S.met = read.csv("SWARM_WP2_all_20200708/Diversity_filtered_18S_WP2_met.csv",row.names=1)
divR.18S.met = divR.18S.met[,c("Reads","Richness","Rarefied.richness","H","J",
                       "Chao1")]

# ---- Filter rare sequences below "detection limit" ------

otus.18S.ra.all = decostand(otusR.18S, method="total")

# For dissimilarity, we require a minimum of 3 expected reads assuming average constant relative OTU abundance
# i.e. an average relative abundance of 6.4E-5 across samples
minAb = 3 / min(rowSums(otusR.18S)) #min is 46826 so min. is 6.4E-5
otusR18S.ab = dropRareByAvgAbundance(otusR.18S, minAb)
dim(otusR18S.ab)/dim(otusR.18S) #770 of 10022 OTUs remain at this cutoff

# ------- Taxonomy --------

selected = mdR.18S[mdR.18S$Experiment!="IntraGrabAndXHom",]

taxa.all.18S = read.table("SWARM_WP2_all_20200708/18S/CREST_Filtered/Relative_Abundance.tsv",sep="\t",
                          header=T,row.names=3, check.names = F)

ass.all.18S = read.table("SWARM_WP2_all_20200708/18S/CREST_Filtered/All_Assignments.tsv",sep="\t",
                         header=T,row.names=3, check.names = F)

assignments.selected = decostand(as.data.frame(t(ass.all.18S[,row.names(selected)])),
                      method="total")

tra.selected = taxa.all.18S[,c("Rank","Taxonpath",row.names(selected))]

grouping_info<-data.frame(row.names=row.names(selected), selected$Group)

ranks = data.frame(rank=c("superkingdom","kingdom","phylum","class","order","family","genus"),
                   levels=c(6,9,25,25,25,25,25))

summary(tra.selected$Rank)
# 13 superkingdoms, 26 kingdoms, 102 phyla, 207 classes, 261 orders, 164 families, 
# 148 genera A lot of abundance is fungi (!)


for (i in c(1:7)){
  r=as.character(ranks$rank[i])
  leveltax.18Sa = as.data.frame(t(tra.selected[tra.selected$Rank==r,-c(1:2)]))
  print(paste(mean(rowSums(leveltax.18Sa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots_all/",r,".pdf",sep=""),height=9,width=20)
  taxaplot(ranks$levels[i],grouping_info,leveltax.18Sa)
  dev.off()
}

metKingdom = tra.selected["Metazoa (Animalia)",c(3:92)]

tra.selected.met = tra.selected[grep("Metazoa",tra.selected$Taxonpath),] 

ranks = data.frame(rank=c("phylum","class","order","family"),
                   levels=c(22,24,24,12))

for (i in c(1:4)){
  r=as.character(ranks$rank[i])
  leveltax.18Sa = as.data.frame(t(tra.selected.met[tra.selected.met$Rank==r,-c(1:2)])) / t(metKingdom)
  print(paste(mean(rowSums(leveltax.18Sa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots_all/",r,"_Metazoa.pdf",sep=""),height=9,width=20)
  taxaplot(ranks$levels[i],grouping_info,leveltax.18Sa)
  dev.off()
}

pdf("img/18S/taxon_barplots_all/Assignments.pdf",height=9,width=20)
taxaplot(30,grouping_info,assignments.selected)
dev.off()

## Grabs 1-3 and 5

mdExtra = selected[selected$Experiment=="GrabEXandPCRRep",]

tra_og = taxa.all.18S[,c("Rank","Taxonpath",row.names(mdExtra))]

grouping_info<-data.frame(row.names=row.names(mdExtra), paste("Grab",mdExtra$Grab))

ranks = data.frame(rank=c("superkingdom","kingdom","phylum","class","order","family","genus"),
                   levels=c(6,9,25,25,25,25,25))

for (i in c(1:7)){
  r=as.character(ranks$rank[i])
  leveltax.18Sa = as.data.frame(t(tra_og[tra_og$Rank==r,-c(1:2)]))
  print(paste(mean(rowSums(leveltax.18Sa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots_all/",r,"_pre_v_post_PCR_pooling.pdf",sep=""),
      height=7,width=10)
  taxaplot(ranks$levels[i],grouping_info,leveltax.18Sa)
  dev.off()
}


assignments.extra = decostand(as.data.frame(t(ass.all.18S[,row.names(mdExtra)])),
                      method="total")
pdf("img/18S/taxon_barplots_all/Assignments_pre_v_post_PCR_pooling.pdf",
    height=6,width=9)
taxaplot(30,grouping_info,assignments.extra)
dev.off()
