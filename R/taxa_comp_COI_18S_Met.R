# ---- CLASS RANK -----

tra.18S = taxa.all.18S[,c("Rank","Taxonpath",row.names(mdX.18S))] #or mdR.18S
metKingdom = tra.18S["Metazoa (Animalia)",row.names(mdX.18S)]
tra.met.18S = tra.18S[grep("Metazoa",tra.18S$Taxonpath),]

tra.COI = taxa.all.COI[,c("Rank","Taxonpath",row.names(mdX.COI))] #or md.COIf

classes.18S = as.data.frame(t(tra.met.18S[tra.met.18S$Rank=="class",-c(1:2)]))
classes.COI = as.data.frame(t(tra.COI[tra.COI$Rank=="class",-c(1:2)]))

classes.18S = classes.18S[,colSums(classes.18S)>0]
classes.COI = classes.COI[,colSums(classes.COI)>0]

#Normalize to metazoa
classes.18S.ra = classes.18S/t(metKingdom)

groups.18S=paste("18S",mdX.18S$Group,mdX.18S$KitHom)
groups.18S[c(49:53)] = "18S In silico pools"
groups.18S[c(25:27,48)] = "18S In vitro pools"

groups.COI=paste("COI",mdX.COI$Group,mdX.COI$KitHom)
groups.COI[c(43:47)] = "COI In silico pools"

dim(classes.COI) #38 COI classes
dim(classes.18S.ra) #40 18S
sum(names(classes.COI) %in% names(classes.18S.ra)) #21 out of 57 shared, 19 unique to 18S, 17 to COI

names(classes.COI)[!names(classes.COI) %in% names(classes.18S)]
# [1] "Malacostraca"   "Collembola"     "Hexanauplia"    "Diplopoda"      "Chilopoda"      "Pycnogonida"   
# [7] "Clitellata"     "Gymnolaemata"   "Crinoidea"      "Cephalopoda"    "Polyplacophora" "Adenophorea"   
# [13] "Anopla"         "Onychophorida"  "Rhabditophora"  "Demospongiae"   "Priapulimorpha"

names(classes.18S)[!names(classes.18S) %in% names(classes.COI)]
# [1] "Ostracoda"                           "Enoplea"                            
# [3] "Thaliacea"                           "Appendicularia class incertae sedis"
# [5] "Catenulida"                          "Seriata"                            
# [7] "Macrostomida"                        "Rhabdocoela class incertae sedis"   
# [9] "Phoroniformea"                       "Calcarea"                           
# [11] "Hexactinellida"                      "Typhlocoela"                        
# [13] "Priapulidae class incertae sedis"    "Homalorhagida class incertae sedis" 
# [15] "Cyclorhagida class incertae sedis"   "Heteronemertea"                     
# [17] "Palaeonemertea"                      "Macrodasyida class incertae sedis"  
# [19] "Chaetonotida class incertae sedis"  

row.names(classes.18S.ra) = paste("18S", row.names(classes.18S.ra))
row.names(classes.COI) = paste("COI", row.names(classes.COI))
classes.18S.ra$class = row.names(classes.18S.ra)
classes.COI$class = row.names(classes.COI)

classes.both = merge(classes.18S.ra, classes.COI, all=T)
<<<<<<< HEAD
D2_LUR20/D2_Keystones.xlsxrow.names(classes.both) = classes.both$class
=======
classes.both[is.na(classes.both)] = 0
row.names(classes.both) = classes.both$class
>>>>>>> e8b9ab1f5cd2232eab18be119f56b56f47086952
classes.both = classes.both[c(row.names(classes.18S.ra),row.names(classes.COI)),]
names(classes.both)
classes.both = classes.both[,-c(1:199)[names(classes.both) == "class"]]

grouping_info<-data.frame(row.names=row.names(classes.both),c(groups.18S,groups.COI))


pdf("../img/Classes_18Smet_v_COI.pdf",height=8,width=30)
taxaplot(20,grouping_info,classes.both)
dev.off()

# ---- PHYLUM RANK -----

phyla.18S = as.data.frame(t(tra.met.18S[tra.met.18S$Rank=="phylum",-c(1:2)]))
phyla.COI = as.data.frame(t(tra.COI[tra.COI$Rank=="phylum",-c(1:2)]))

phyla.18S = phyla.18S[,colSums(phyla.18S)>0]
phyla.COI = phyla.COI[,colSums(phyla.COI)>0]

#Normalize to metazoa
phyla.18S.ra = phyla.18S/t(metKingdom)

groups.18S=paste("18S",mdX.18S$Group,mdX.18S$KitHom)
groups.18S[c(49:53)] = "18S In silico pools"
groups.18S[c(25:27,48)] = "18S In vitro pools"

groups.COI=paste("COI",mdX.COI$Group,mdX.COI$KitHom)
groups.COI[c(43:47)] = "COI In silico pools"

dim(phyla.COI) #17 COI phyla
dim(phyla.18S.ra) #20 18S phyla
sum(names(phyla.COI) %in% names(phyla.18S.ra)) #15 out of 22 shared, 5 unique to 18S, 2 to COI

names(phyla.COI)[!names(phyla.COI) %in% names(phyla.18S)]
# [1] "Bryozoa", "Onychophora"

names(phyla.18S)[!names(phyla.18S) %in% names(phyla.COI)]
# [1] "Brachiopoda"     "Ctenophora"      "Xenacoelomorpha" "Kinorhyncha"     "Gastrotricha"   

row.names(phyla.18S.ra) = paste("18S", row.names(phyla.18S.ra))
row.names(phyla.COI) = paste("COI", row.names(phyla.COI))
phyla.18S.ra$class = row.names(phyla.18S.ra)
phyla.COI$class = row.names(phyla.COI)

<<<<<<< HEAD
phyla.both = merge(phyla.18S.ra, phyla.COI, all=T)phyla.both[is.na(phyla.both)] = 0
=======
phyla.both = merge(phyla.18S.ra, phyla.COI, all=T)
phyla.both[is.na(phyla.both)] = 0
>>>>>>> e8b9ab1f5cd2232eab18be119f56b56f47086952
row.names(phyla.both) = phyla.both$class
phyla.both = phyla.both[c(row.names(phyla.18S.ra),row.names(phyla.COI)),]
names(phyla.both)
phyla.both = phyla.both[,-c(1:199)[names(phyla.both) == "class"]]

grouping_info<-data.frame(row.names=row.names(phyla.both),c(groups.18S,groups.COI))


pdf("../img/Phyla_18Smet_v_COI.pdf",height=8,width=30)
taxaplot(20,grouping_info,phyla.both)
dev.off()

