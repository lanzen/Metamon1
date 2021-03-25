library(multcomp)
library(lsmeans)

printAllvsAll = function(df,a,method="kendall"){
  totalCols = dim(df)[2]
  for (i in c(1:(totalCols-1))) {
    pn = colnames(df)[i]
    n = df[,i]
    print ("**************************************************")
    print (pn)
    print ("**************************************************")
    
    for (j in c((i+1):totalCols)){
      pm = colnames(df)[j]
      print (pm)
      m = df[,j]
      print(m)
      print(n)
      print(method)
      ct = cor.test(n,m,
                    method=method,use="pairwise.complete.obs")
      
      if (!is.na(ct$p.value) &&  ct$p.value < a){
        print(pm)
        print(ct)
        print(summary(lm(m~n)))
        plot(n,m,xlab=pn,ylab=pm)
        x=max(n,na.rm=T)*.2+min(n,na.rm=T)
        y=max(m)*.9
        text(x,y,paste(paste("p =",signif(ct$p.value,2)),paste(" tau =",signif(ct$estimate,2)))
             ,col="blue")
        lr = lm(m~n)
        abline(lr,col="grey")
        line <- readline()      
      }
    }
  }
}

printVS = function(df1,df2,a,method="kendall"){
  cols1 = dim(df1)[2]
  cols2 = dim(df2)[2]
  for (i in c(1:(cols1))) {
    pn = colnames(df1)[i]
    n = df1[,i]
    print ("**************************************************")
    print (pn)
    print ("**************************************************")
    
    for (j in c(1:cols2)){
      pm = colnames(df2)[j]
      print (pm)
      m = df2[,j]
      ct = cor.test(n,m,
                    method=method,use="pairwise.complete.obs")
      
      if (!is.na(ct$p.value) &&  ct$p.value < a){
        print(pm)
        print(ct)
        print(summary(lm(n~m)))
        plot(m,n,xlab=pm,ylab=pn)
        x=min(m,na.rm=T)
        y=max(n)*.9
        text(x,y,paste(paste("p =",signif(ct$p.value,2)),paste(" tau =",signif(ct$estimate,2)))
             ,col="blue")
        lr = lm(n~m)
        abline(lr,col="grey")
        line <- readline()      
      }
    }
  }
}


printANOVA = function(factors,para.df,a,notch=T,imgDir=""){
  totalp = dim(para.df)[2]
  totalf = dim(factors)[2]
  for (i in c(1:totalf)) {
    fn = colnames(factors)[i]
    n = factors[,i]
   
    
    for (j in c(1:totalp)){      
      pm = colnames(para.df)[j]
      m = para.df[,j]
      an=anova(lm(m~n))
      p=an$Pr[1]
     
      if (!is.na(p) && p<a){
        print ("**************************************************")
        print (paste(fn, pm,sep = " - "))
        print ("**************************************************")
        an2 = aov(m~as.factor(n))
        posthoc=TukeyHSD(an2)
        print(posthoc)
        #print(paste(fn,pm,p,sep=","))
        
        if (imgDir!=""){
          pdf(paste(imgDir,"/",fn,"_",pm,".pdf",sep=""),width=3, height = 6)
          boxplot(m~n,notch=notch,main=paste(fn,"-",pm,", p=",signif(p,2)),cex.main=0.8,las=2)
          dev.off()
        }
        boxplot(m~n,notch=notch,main=paste(fn,"-",pm,", p=",signif(p,2)),cex.main=0.8,las=2)
        
        #cld(ref.grid(an2))
        
        #print(HSD.test(aov(m~as.factor(n)), "n", group=TRUE))
        
        line <- readline()      
      }
    }
  }
}

printANOVA1Factor = function(f,para.df,a,notch=T){
  totalp = dim(para.df)[2]

  for (j in c(1:totalp)){      
      pm = colnames(para.df)[j]
      m = para.df[,j]
      an=anova(lm(m~f))
      p=an$Pr[1]
      print(paste(pm, p))
      if (!is.na(p) && p<a){
        posthoc=TukeyHSD(aov(m~as.factor(f)))
        print(posthoc)
        boxplot(m~f,notch=notch,main=paste(pm,", p=",signif(p,2)),las=2)
        line <- readline()      
      }
    }
  
}

printANOVASimpler = function(factors,para.df,a,notch=T){
  totalp = dim(para.df)[2]
  totalf = dim(factors)[2]
  for (i in c(1:totalf)) {
    fn = colnames(factors)[i]
    n = factors[,i]
    print ("**************************************************")
    print (fn)
    print ("**************************************************")
    
    for (j in c(1:totalp)){      
      pm = colnames(para.df)[j]
      m = para.df[,j]
      an=anova(lm(m~n))
      p=an$Pr[1]
      if (!is.na(p) && p<a){
        print(paste(fn,",",pm,",",p))
      }
    }
  }
}

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$n=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$n) , ]
  return(Tukey.labels)
}


# ------------------------------------
# Partial Effect Plot function (peplot)
# From Miles D. Williams, Rpubs 
# (https://rpubs.com/milesdwilliams15/328471)
# ------------------------------------
peplot <- function(mod,v,ci=.95, plot_points = "n",
                   xlab=v,ylab=names(mod[12]$model)[1],
                   main="Partial Effect Plot",
                   pe_lty=1,pe_lwd=3,pe_col="black",
                   ci_lty=1,ci_lwd=1,ci_col="black",
                   pch_col="black",pch_ty=19,
                   ylim=c(min(pred[,"lwr"]),max(pred[,"upr"]))){
  modDat <- mod[12]$model
  modDat1 <- modDat[,-1]
  modDat2 <- modDat[,which(names(modDat)!=v)]
  x <- resid(lm(modDat1[,v] ~., data=modDat1[,which(names(modDat1)!=v)]))
  y <- resid(lm(modDat2[,1] ~ ., modDat2[,-1]))
  plot(x,y,type=plot_points,xlab=xlab,ylab=ylab,
       ylim=ylim,col=pch_col,pch=pch_ty,
       main=main)
  part <- lm(y~x)
  wx <- par("usr")[1:2]
  new.x <- seq(wx[1],wx[2],len=100)
  pred <- predict(part, new=data.frame(x=new.x), interval="conf",
                  level = ci)
  lines(new.x,pred[,"fit"],lwd=pe_lwd,lty=pe_lty,col=pe_col)
  lines(new.x,pred[,"lwr"],lwd=ci_lwd,lty=ci_lty,col=ci_col)
  lines(new.x,pred[,"upr"],lwd=ci_lwd,lty=ci_lty,col=ci_col)
}