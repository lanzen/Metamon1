dropRareTaxa <- function(df, minAbundance){
  require(vegan)
  radist=decostand(df,method="total")
  m=as.matrix(radist)
  fdist=as.data.frame(ifelse(m>minAbundance,m,0))
  names(fdist)=names(radist)
  row.names(fdist)=row.names(radist)
  nonzero = fdist[,colSums(fdist)>0]
  return(nonzero)
}

dropAbundant <- function(df, maxAbundance){
  require(vegan)
  radist=decostand(df,method="total")
  m=as.matrix(radist)
  fdist=as.data.frame(ifelse(m<maxAbundance,m,0))
  names(fdist)=names(radist)
  row.names(fdist)=row.names(radist)
  nonzero = fdist[,colSums(fdist)>0]
  return(nonzero)
}

dropSingletons <- function(df){
  require(vegan)
  good = df[,colSums(df)>1]
  return(good)
}
