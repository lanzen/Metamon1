source("~/kode/R/dropRareTaxa.R")

writeDivStats = function(outfile, otudistro, rr_cutoff=min(rowSums(otudistro)), rarefyForShannon=F) {

  write(c("Dataset","Reads","Reads_wo_S","Richness",
          "Richness_wo_S","Rarefied.richness","H","J","Rarefied J","Simpson", 
          "Chao1","ACE"),
        file=outfile,ncolumns=12,sep=",")
  
  noS=dropSingletons(otudistro)
  reads = rowSums(otudistro)
  readsNoS=rowSums(noS) 
  
  ds = row.names(otudistro)
  sn=specnumber(otudistro)
  rr=rarefy(otudistro,sample=rr_cutoff)
  
  if (rarefyForShannon){
    shDistro = rrarefy(otudistro, sample=rr_cutoff)
  }
  else shDistro = otudistro
  
  e=estimateR(otudistro)
  for (i in c(1:length(ds))){
    sh= diversity(shDistro[i,],index="shannon")
    ace = e[4,i]
    chao1 = e[2,i]
    write(c(ds[i],reads[i],readsNoS[i],
            sn[i],specnumber(noS[i,]),rr[i],
            sh,sh/log(sn[i]), sh/log(rr[i]),
            diversity(shDistro[i,],index="simpson"),
            ace,chao1),
          file=outfile,ncolumns=12,sep=",",append=T)    
  }
}