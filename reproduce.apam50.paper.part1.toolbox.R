##
## Necessary functions to train AIMS models.
##
## Author : Eric R. Paquet (eric.r.paquet@gmail.com)
##
## Copyright McGill University, 2015
##

## get the specific cross-validation results
## for k number of rules
get.cv.pam50 <- function(cv,k){
  to.ret <- c()
  probs <- c()
  all.probs <- c()
  for (i in cv){
    to.ret <- rbind(to.ret,cbind(i$test.pred$cl[,k,drop=F],i$to.leave))
    probs <- c(probs,i$test.pred$prob[,k,drop=F])
    all.probs <- rbind(all.probs,i$test.pred$all.probs[[k]])
  }
  to.o <- order(as.numeric(to.ret[,2]))
  to.ret <- to.ret[to.o,]
  probs <- probs[to.o]
  all.probs <- all.probs[to.o,]
  list(idx=as.numeric(to.ret[,2]),pam50=to.ret[,1],probs=probs,all.probs=all.probs)
}

## This will plot a matrix + number in the matrix
## and using a colorcode (red == high number, white = low number) 
plotPam50.mat <- function(to.plot.mismatch.matrix,N,percentage=T){
  require(gplots)
  image(x=1:5,y=1:5,z=t(to.plot.mismatch.matrix)[,5:1],col=colorpanel(50, "white", "red"),
        ylab="",xlab="",axes=F)
  for(i in 1:5){for(j in 1:5){
    cur.num <- t(to.plot.mismatch.matrix)[,5:1][i,j]
    if (percentage){
      text(i,j,sprintf("%d (%.2f %%)",cur.num,cur.num*100/N))
    }
    else{
      text(i,j,sprintf("%d",cur.num))
    }
  }}
  
  abline(v=0:5+.5,lwd=3)
  abline(h=0:5+.5,lwd=3)
  pam.name=rownames(to.plot.mismatch.matrix)
  axis(2, at=5:1, labels=pam.name,las=2,cex.axis=.75,lwd=2,tick=F)
  pam.name <- colnames(to.plot.mismatch.matrix)[5:1]
  axis(3, at=5:1, labels=pam.name,las=2,cex.axis=.75,tick=F,lwd=2)
}

## This the the function used to assign subtypes
## on a data D using EntrezID and the model sel.nbc (eg a-pam50)
apply.nbc <- function(D,EntrezID,sel.nbc){
  D <- apply(D,2,as.numeric)
  EntrezID <- as.character(EntrezID)
  sel.ids.nb <- get.all.pairs.genes(sel.nbc$all.pairs)
  sel.gn <- EntrezID %in% sel.ids.nb
  D <- D[sel.gn,]
  Entrez <- EntrezID[sel.gn]
  col.D.test <- removeDuplicatedEntrezPerPatients(D,Entrez)
  pred.test <- predict.one.vs.all.tsp(D=col.D.test$dataset,
                                      GeneName=col.D.test$EntrezID,
                                      sel.nbc)

  ## Add more information to the output variable
  pred.test$data.used <- col.D.test$dataset
  pred.test$EntrezID.used <- col.D.test$EntrezID
  
  invisible(pred.test)
}

## split the simples rules 
## This function will output a heatmap of the conditional
## probability table for a NB classifier
plotNBClassifier <- function(nbc,mappingEntrez2Gene){
  to.mat <- nbc$apriori/sum(nbc$apriori)
  cur.tags <- c("Prior")
  for (ni in names(nbc$tables)){
    spl <- strsplit(gsub("^X","",ni),"\\.")[[1]]
    cur.tags <- c(cur.tags,paste(mappingEntrez2Gene[spl[1],2],mappingEntrez2Gene[spl[2],2],sep="<"))
    to.mat <- rbind(to.mat,nbc$tables[[ni]][,1])
  }
  rownames(to.mat) <- cur.tags
  colnames(to.mat) <- rownames(nbc$tables[[ni]])
  to.mat
}

## Given a cross-validation reporte the misclassified samples
misClassify <- function(cv,goldstandard,k){
  probs.mis <- c()
  max.cor <- rowMax(goldstandard$cor)
  for (cvi in cv){
    probs.mis <- rbind(probs.mis,cbind(max.cor[cvi$to.leave],
                                       goldstandard$subtype[cvi$to.leave] == cvi$test.pred$cl[,k],
                                       goldstandard$subtype[cvi$to.leave],cvi$test.pred$cl[,k]))
  }
  
  probs.mis <- probs.mis[order(as.numeric(probs.mis[,1])),]
}

## given cross-validation cv and k number of rules
## return the misclassified samples
rankMisclassify <- function(cv,goldstandard,k){
  to.ret <- c()
  for (cvi in cv){
    cur.cor <- goldstandard$cor[cvi$to.leave,]
    cur.cl <- goldstandard$subtype[cvi$to.leave]
    cur.diff <- cvi$test.pred$cl[,k] != cur.cl
    diff.cl <- cvi$test.pred$cl[cur.diff,k]
    cur.cor.diff <- cur.cor[cur.diff,]
    to.ret <- c(to.ret,as.numeric(lapply(1:length(diff.cl),
                                         function(x){temp=6-rank(cur.cor.diff[x,],ties.method = "first");
                                                     temp[grep(diff.cl[x],colnames(cur.cor.diff))]})))
   
  }
  to.ret
}

## This function split the rules eg (a < b, b < c) and
## return a vector with all the genes c(a,b,c)
get.all.pairs.genes <- function(all.pairs){
  genes <- c()
  for (cp in strsplit(all.pairs,"<")){
    genes <- c(genes,cp)
  }
  unique(genes)
}


#
confu.k.cv <- function(cv,cl){
  to.ret <- list()
  for (cv.i in 1:length(cv)){
    for (ki in colnames(cv[[cv.i]]$test.pred$cl)){
      to.add <- cbind(pred=cv[[cv.i]]$test.pred$cl[,ki],
                      gold=cl[cv[[cv.i]]$to.leave])
      if (ki %in% names(to.ret)){
        to.ret[[ki]] <- rbind(to.ret[[ki]],to.add)
      }
      else{
        to.ret[[ki]] <- to.add
      }
    }
  }
  to.ret
}

## From an Entrez Id pair eg 1 < 2 it will return a gene symbol
## version of the pair A < B
## If pam50 is supplied it will mark PAM50 gene symbols with a *
pairs2symbols.only <- function(pairs,probe.info,pam50=c()){
  to.ret <- c()
  for (p in pairs){
    spl <- strsplit(sort.pairs(p),"<")[[1]]
    gi1 <- probe.info$GeneName[which(probe.info$EntrezID == spl[1])[1]]
    gi2 <- probe.info$GeneName[which(probe.info$EntrezID == spl[2])[1]]
    if (spl[1] %in% pam50){
      gi1 <- sprintf("%s*",gi1)
    }
    if (spl[2] %in% pam50){
      gi2 <- sprintf("%s*",gi2)
    }
    to.ret <- c(to.ret,sprintf("%s < %s",gi1,gi2))
  }
  to.ret
}

pairs2symbols <- function(pairs,probe.info,nbc,pam,pam50=c()){
  to.ret <- c()
  for (p in pairs){
    spl <- strsplit(sort.pairs(p),"<")[[1]]
    gi1 <- probe.info$GeneName[which(probe.info$EntrezID == spl[1])[1]]
    if (spl[1] %in% pam50){
      gi1 <- sprintf("%s*",gi1)
    }
    
    gi2 <- probe.info$GeneName[which(probe.info$EntrezID == spl[2])[1]]
    if (spl[2] %in% pam50){
      gi2 <- sprintf("%s*",gi2)
    }
    if (nbc$table[[sprintf("X%s.%s",spl[1],spl[2])]][pam,1] > nbc$table[[sprintf("X%s.%s",spl[1],spl[2])]][pam,2]){
      to.ret <- c(to.ret,paste(gi1,gi2,sep="<"))
    }
    else{
      to.ret <- c(to.ret,paste(gi2,gi1,sep="<"))
    }
  }
  to.ret
}
