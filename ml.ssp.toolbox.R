##
## Necessary functions to train AIMS models.
##
## Author : Eric R. Paquet (eric.r.paquet@gmail.com)
##
## Copyright McGill University, 2015
##

require("ROCR")

color2pam50 <- function(color){
  c("Basal","Normal","Her2","LumB","LumA","ER+","negative","positive","her2+","2")[match(color,c("firebrick2","green4","hotpink2","deepskyblue","blue4","green","white","red","orange","pink"))]
}
pam502color <- function(color){
  c("firebrick2","green4","hotpink2","deepskyblue","blue4","green","white","red","orange","pink")[match(color,c("Basal","Normal","Her2","LumB","LumA","ER+","negative","positive","her2+","2"))]
}

getGeneFu.scale.info <- function(dataset){
  require(genefu)
  data(pam50)
  probes <- as.character(dataset$probe.info[,grep("ProbeID",names(dataset$probe.info))])
  dataset$probe.info$EntrezID[dataset$probe.info$EntrezID == "NULL"]=NA
  rownames(dataset$exprs) <- probes
  dataset.pam50.genefu.scale <- intrinsic.cluster.predict(sbt.model=pam50.scale,
                                                   data=t(dataset$exprs), annot=data.frame(probe=probes,EntrezGene.ID=as.character(dataset$probe.info$EntrezID),stringsAsFactors=F),
                                                   do.mapping=TRUE,
                                                   do.prediction.strength=FALSE, verbose=TRUE)
  dataset.pam50.genefu.scale
}

smaller <- function(x,y){
  x < y
}

## Expect a matrix of log2(exp) values
## rows = genes
## columns = samples
## Output a matrix with computed pairs of features
compute.pairs <- function(d,func=smaller,n.cores=20){
  require(multicore)

  list.pairs <- mclapply(1:(nrow(d)-1),function(i){
    message(sprintf("Computing %d/%d",i,nrow(d)))
    i.j <- (i+1):nrow(d)
    to.ret <-  matrix(0,nrow=length(i.j),ncol=ncol(d))
    to.ret.diff <- matrix(0,nrow=length(i.j),ncol=ncol(d))
    to.ret.idxs <- matrix(0,nrow=length(i.j),ncol=2)
    rendu <- 1
    for (j in i.j){
      to.ret[rendu,] <- func(d[i,],d[j,])
      to.ret.diff[rendu,] <- d[i,] - d[j,]
      to.ret.idxs[rendu,] <- c(i,j)
      rendu <- rendu + 1
    }
    list(to.ret=to.ret,to.ret.diff=to.ret.diff,to.ret.idxs=to.ret.idxs)
  },mc.cores=n.cores)
  
  num.features <- nrow(d)*(nrow(d) + 1)/2 - nrow(d)
  to.ret <- matrix(0,nrow=num.features,ncol=ncol(d))
  to.ret.diff <- matrix(0,nrow=num.features,ncol=ncol(d))
  to.ret.idxs <- matrix(0,nrow=num.features,ncol=2)
  rendu <- 1
  for(lpi in list.pairs){
    message(sprintf("Rendu = %d",rendu))
    cur.ranges <- rendu:(rendu+nrow(lpi$to.ret)-1)
    to.ret[cur.ranges,] <- lpi$to.ret
    to.ret.diff[cur.ranges,] <- lpi$to.ret.diff
    to.ret.idxs[cur.ranges,] <- lpi$to.ret.idxs
    rendu <- rendu+nrow(lpi$to.ret)
  }
  stopifnot(rendu == (num.features+1))
  
  invisible(list(features=to.ret,indices=to.ret.idxs,diff=to.ret.diff))
}

## Expect a matrix of log2(exp) values
## rows = genes
## columns = samples
## Just generate random p-values (runif) to estimate the background
## Output a matrix with computed pairs of features
compute.pairs.pv <- function(d,cl,w,func=smaller,random=F,n.cores=10){
  require(multicore)
  if (random){
    num.features <- nrow(d)*(nrow(d) + 1)/2 - nrow(d)
    to.ret <- matrix(0,nrow=num.features,ncol=3)
    rendu = 1
    for (i in 1:(nrow(d)-1)){
      for (j in (i+1):nrow(d)){
        to.ret[rendu,] <- c(i,j,runif(1))
        rendu = rendu + 1
      }
    }
    to.ret
  }
  else{
    list.pairs <- mclapply(1:(nrow(d)-1),function(i){
      if (i %% 100 == 1){
        message(sprintf("Computing %d/%d",i,nrow(d)))
      }
      i.j <- (i+1):nrow(d)
      ##to.ret <-  matrix(0,nrow=length(i.j),ncol=ncol(d))
      ##to.ret.diff <- matrix(0,nrow=length(i.j),ncol=ncol(d))
      to.ret.idxs <- matrix(0,nrow=length(i.j),ncol=3)
      rendu <- 1
      for (j in i.j){
        curp <- func(d[i,],d[j,])
        cnts00 <- sum((!curp & cl == 0) * w)
        cnts10 <- sum((curp & cl == 0) * w)
        cnts11 <- sum((curp & cl == 1) * w)
        cnts01 <- sum((!curp & cl == 1) * w)
        to.ret.idxs[rendu,] <- c(i,j,chisq.test(cbind(c(cnts00,cnts10),c(cnts01,cnts11)))$p.value)
        
        rendu <- rendu + 1
      }
      to.ret.idxs
    },mc.cores=n.cores)
    
    num.features <- nrow(d)*(nrow(d) + 1)/2 - nrow(d)
    to.ret <- matrix(0,nrow=num.features,ncol=3)
    
    rendu <- 1
    for(lpi in list.pairs){
      ##message(sprintf("Rendu = %d",rendu))
      cur.ranges <- rendu:(rendu+nrow(lpi)-1)
      to.ret[cur.ranges,] <- lpi
      rendu <- rendu+nrow(lpi)
    }
    stopifnot(rendu == (num.features+1))
    to.ret
  }
}

## compute.pairs.pv <- function(d,cl,w,func=smaller,random=F,n.cores=10){
##   require(Rgtsp)
##   if (random){
##     num.features <- nrow(d)*(nrow(d) + 1)/2 - nrow(d)
##     to.ret <- matrix(0,nrow=num.features,ncol=3)
##     rendu = 1
##     for (i in 1:(nrow(d)-1)){
##       for (j in (i+1):nrow(d)){
##         to.ret[rendu,] <- c(i,j,runif(1))
##         rendu = rendu + 1
##       }
##     }
##     to.ret
##   }
##   else{
##     X <- t(d)
##     n <- nrow(d)*(nrow(d) + 1)/2 - nrow(d)
##     n = min(max(n, 1), 0.5 * ncol(X) * (ncol(X) - 1))
##     i = rep(0, n)
##     j = rep(0, n)
##     s = rep(0, n)
##     z = list(I = i, J = j, S = s) 
##     message("Rgtsp starts computing")
##     z = .C("RGTSP_tsp_NW", as.double(X), as.integer(cl), 
##       as.double(w), as.integer(ncol(X)), as.integer(nrow(X)), 
##       i = as.integer(i), j = as.integer(j), s = as.double(s), 
##       n = as.integer(n))
##     message("Rgtsp ends computing")
##     cbind(z$i,z$j,1-z$s)
##   }
## }

compute.pairs.full <- function(d,func=smaller){
  num.features <- nrow(d)^2 - nrow(d)
  to.ret <- matrix(0,nrow=num.features,ncol=ncol(d))
  to.ret.idxs <- matrix(0,nrow=num.features,ncol=2)
  rendu <- 1
  for (i in 1:nrow(d)){
    for (j in 1:nrow(d)){
      if (i != j){
        to.ret[rendu,] <- func(d[i,],d[j,])
        to.ret.idxs[rendu,] <- c(i,j)
        rendu <- rendu + 1
      }
    }
  }
  invisible(list(features=to.ret,indices=to.ret.idxs))
}

# Need list(exprs,GeneName)
comp.sel.pairs <- function(dataset,sel.pairs,func=smaller){
  to.ret <- matrix(0,nrow=length(sel.pairs),ncol(dataset$exprs))
  for (spi in 1:length(sel.pairs)){
    ss.cur <- strsplit(sel.pairs[spi],"<")[[1]]
    gene1 = which(dataset$GeneName == ss.cur[1])
    gene2 = which(dataset$GeneName == ss.cur[2])
    #stopifnot(length(gene1) == 1 & length(gene2) == 1)
    if (length(gene1) == 1 & length(gene2) == 1){
      to.ret[spi,] <- func(dataset$exprs[gene1,],dataset$exprs[gene2,])
    }
    else{
      message(paste("You are missing the pair or have more than one",sel.pairs[spi],"in",dataset$name))
    }
  }
  rownames(to.ret) <- sel.pairs
  ##to.ret <- apply(to.ret,2,as.factor) ## convert to factor
  ##to.ret <- apply(to.ret,2,as.numeric)
  to.ret
}


# Cross validation
# 10-fold CV
cv.nb <- function(D,cl,K.FOLD,feature.names,seed=1234,n.cores=10){
  require(ROCR)
  require(e1071)
  
  set.seed(seed)
  
  k.samples.out <- list()
  all.samples <- 1:ncol(D)

  ## select the samples to left out during cross-validation
  for (i in 1:K.FOLD){
    cur.sel <- sample(all.samples,floor(ncol(D)/K.FOLD))
    all.samples <- setdiff(all.samples,cur.sel)
    k.samples.out[[i]] <- cur.sel
  }

  ## Cross-validate
  CV.stats <- list()
  rendu <- 1
  for(ksi in k.samples.out){
    cat(sprintf("Rendu = %d/%d\n",rendu,length(k.samples.out)))
    rendu <- rendu + 1
    sel.train <- setdiff(1:ncol(D),ksi)
    sel.test <- ksi
    
    train.d <- D[,sel.train]
    test.d <- D[,sel.test]
    
    train.cl <- cl[sel.train]
    test.cl <- cl[sel.test]
    
    class.cur <- factor(train.cl)
    pvalues <- as.numeric(mclapply(1:nrow(train.d),
                                   function(i){
                                     x=train.d[i,];
                                     chisq.test(table(factor(x,levels=c(0,1)),class.cur))$p.value
                                   }))
    
    cur.ordering <- order(pvalues)
    cur.ordering <- cur.ordering[!duplicated(feature.names[cur.ordering])]
    
    current.ml.stats <- mclapply(1:20,function(sel.N){
      sel.ids <- cur.ordering[1:sel.N]
      
      ## train
      nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.d[sel.ids,,drop=F])))
      prob.train <- predict(nb,data.frame(t(train.d[sel.ids,,drop=F])),type="raw")
      prob.test <- predict(nb,data.frame(t(test.d[sel.ids,,drop=F])),type="raw")

      list(sel.N=sel.N,
           prob.train=prob.train,
           class.train=train.cl,
           prob.test=prob.test,
           class.test=test.cl,
           sel.ids=sel.ids)
    })
    CV.stats[[length(CV.stats)+1]] <- list(sel.train=sel.train,
                                           sel.test=sel.test,
                                           current.ml.stats)
  }
  
  invisible(CV.stats)
}

cv.nb.ind <- function(D,cl,K.FOLD,feature.names,ind.data,seed=1234,n.cores=10){
  require(ROCR)
  require(e1071)
  
  set.seed(seed)
  
  k.samples.out <- list()
  all.samples <- 1:ncol(D)

  ## select the samples to left out during cross-validation
  for (i in 1:K.FOLD){
    cur.sel <- sample(all.samples,floor(ncol(D)/K.FOLD))
    all.samples <- setdiff(all.samples,cur.sel)
    k.samples.out[[i]] <- cur.sel
  }

  ## Cross-validate
  CV.stats <- list()
  rendu <- 1
  for(ksi in k.samples.out){
    cat(sprintf("Rendu = %d/%d\n",rendu,length(k.samples.out)))
    rendu <- rendu + 1
    sel.train <- setdiff(1:ncol(D),ksi)
    sel.test <- ksi
    
    train.d <- D[,sel.train]
    test.d <- D[,sel.test]
    
    train.cl <- cl[sel.train]
    test.cl <- cl[sel.test]
    
    class.cur <- factor(train.cl)
    pvalues <- as.numeric(mclapply(1:nrow(train.d),
                                   function(i){
                                     x=train.d[i,];
                                     chisq.test(table(factor(x,levels=c(0,1)),class.cur))$p.value
                                   }))
    
    cur.ordering <- order(pvalues)
    cur.ordering <- cur.ordering[!duplicated(feature.names[cur.ordering])]
    
    current.ml.stats <- mclapply(1:20,function(sel.N){
      sel.ids <- cur.ordering[1:sel.N]
      
      ## train
      nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.d[sel.ids,,drop=F])))
      prob.train <- predict(nb,data.frame(t(train.d[sel.ids,,drop=F])),type="raw")
      prob.test <- predict(nb,data.frame(t(test.d[sel.ids,,drop=F])),type="raw")

      prob.gq.test <- predict(nb,data.frame(t(comp.sel.pairs(ind.data,rownames(train.d)[sel.ids]))),type="raw")

      list(sel.N=sel.N,
           prob.train=prob.train,
           class.train=train.cl,
           prob.test=prob.test,
           class.test=test.cl,
           sel.ids=sel.ids,
           prob.gq.test=prob.gq.test)
    })
    CV.stats[[length(CV.stats)+1]] <- list(sel.train=sel.train,
                                           sel.test=sel.test,
                                           current.ml.stats)
  }
  
  invisible(CV.stats)
}


getAuc <- function(pam50.er.kstats,measure="auc"){
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]][[1]]
      cur.pred=prediction(cur.stats$prob.test[,2],cur.stats$class.test)
      performance(cur.pred,measure=measure)@y.values[[1]]
    }
                                     ))
  }
  k.auc
}

# Simple accuracy
# a = pred, b = real class
# w = weigth given to the samples. Useful when
# estimating acc on a combination of datasets. We might
# want the dataset to weigth equaly.
simple.acc <- function(a,b,w=rep(1,length(a))){
  stopifnot(length(a)==length(b))
  if (round(sum(w)) != length(a)){
    message("Adjusting weigth")
    w <- w*length(w)/sum(w)
  }
  sum((a == b) * w)/length(a)
}

# Get Accuracy instead of AUC
getAcc <- function(pam50.er.kstats){
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]][[1]]
      cur.pred=prediction(cur.stats$prob.test[,2],cur.stats$class.test)
      pred.class <- as.character(apply(cur.stats$prob.test,1,function(x)colnames(cur.stats$prob.test)[which.max(x)]))
      real.class <- as.character(cur.stats$class.test)
      simple.acc(pred.class,real.class)
    }
                                     ))
  }
  k.auc
}

# return the percentage of good classification per class
get.perc.pred <- function(pam50.er.kstats,query.class){
  query.class <- as.character(query.class)
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]][[1]]
      cur.pred=prediction(cur.stats$prob.test[,2],cur.stats$class.test)
      pred.class <- as.character(apply(cur.stats$prob.test,1,function(x)colnames(cur.stats$prob.test)[which.max(x)]))
      real.class <- as.character(cur.stats$class.test)
      stopifnot(length(pred.class) == length(real.class))
      if (!(query.class %in% real.class)){
        message(sprintf("!!! Weird %s not present in real.class !!!",query.class))
      }
      sel.cur <- real.class %in% query.class
      #message(sum(sel.cur))
      sum((pred.class == real.class)[sel.cur])/sum(sel.cur)
    }
                                     ))
  }
  k.auc
}

getAuc.gq <- function(pam50.er.kstats,cl,measure="auc"){
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]][[1]]
      sel.not.na=!is.na(cl)
      cur.pred=prediction(cur.stats$prob.gq.test[sel.not.na,2],cl[sel.not.na])
      performance(cur.pred,measure=measure)@y.values[[1]]
    }
                                     ))
  }
  k.auc
}

getAcc.train <- function(pam50.er.kstats,cl){
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]][[1]]
      sel.not.na=!is.na(cl)

      pred.class <- as.character(apply(cur.stats$prob.gq.test[sel.not.na,],1,function(x)colnames(cur.stats$prob.gq.test)[which.max(x)]))
      real.class <- as.character(cl[sel.not.na])
      stopifnot(length(pred.class) == length(real.class))
      simple.acc(pred.class,real.class)
    }
                                     ))
  }
  k.auc
}

get.perc.pred.train <- function(pam50.er.kstats,cl,query.class){
  query.class <- as.character(query.class)
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]][[1]]
      sel.not.na=!is.na(cl)

      pred.class <- as.character(apply(cur.stats$prob.gq.test[sel.not.na,],1,function(x)colnames(cur.stats$prob.gq.test)[which.max(x)]))
      real.class <- as.character(cl[sel.not.na])
      stopifnot(length(pred.class) == length(real.class))
      if (!(query.class %in% real.class)){
        message(sprintf("!!! Weird %s not present in real.class !!!",query.class))
      }
      sel.cur <- real.class %in% query.class
      #message(sum(sel.cur))
      sum((pred.class == real.class)[sel.cur])/sum(sel.cur)
    }
                                     ))
  }
  k.auc
}

getAuc.train <- function(pam50.er.kstats,measure="auc"){
  k.auc=list()
  for (ki in 1:length(pam50.er.kstats[[1]][[3]])){
    k.auc[[ki]] <- as.numeric(lapply(pam50.er.kstats,function(cur.stats){
      cur.stats=cur.stats[[3]][[ki]]
      cur.pred=prediction(cur.stats$prob.train[,2],cur.stats$class.train)
      performance(cur.pred,measure=measure)@y.values[[1]]
    }
                                     ))
  }
  k.auc
}

getGenes <- function(feature.name){
  strsplit(feature.name,"<")[[1]]
}

selIdsWithMaxConstraint <- function(cur.ordering,feature.names,sel.N,max.repeat){
  toRet <- c()
  cur.genes <- c()
  to.add <- 1
  while(length(toRet) < sel.N){
    # lets say we had this feature
    #stopifnot(to.add <= length(feature.names)) # reach the end without being able to select all the necessary ids
    if (!(to.add <= length(feature.names))){
      message("!!!!!!!!!!!!! Not ABLE to select all features !!!!!!!!!!!!!!!!!!")
      break
    }
    idx <- cur.ordering[to.add]
    cur.genes.i <- getGenes(feature.names[to.add])
    if (all(table(c(cur.genes,cur.genes.i)) <= max.repeat)){ # it seems you could add it
        toRet <- c(toRet,idx)
        cur.genes <- c(cur.genes,cur.genes.i)
    }
    to.add = to.add + 1
  }
  toRet
}

test.selIdsWithMaxConstraint <- function(){
  selIdsWithMaxConstraint(1:10,c("a<b","b<c","e<f"),2,1)
}

## Use compute pairs as first argument
cv.nb.cp <- function(D,cl,K.FOLD,ind.data,seed=1234,n.cores=10,random=F,w=rep(1,ncol(D))){
  require(ROCR)
  require(e1071)
  
  set.seed(seed)
  
  k.samples.out <- list()
  all.samples <- 1:ncol(D)

  ## select the samples to leave out during cross-validation
  for (i in 1:K.FOLD){
    ok <- F
    count.ok <- 0
    while(!ok){ # need both classes in both groups
      cur.sel <- sample(all.samples,floor(ncol(D)/K.FOLD))
      if (all(table(factor(cl[cur.sel],levels=c(0,1))) > 0) &
          all(table(factor(cl[setdiff(1:length(cl),cur.sel)],levels=c(0,1))) > 0)){
        break
      }
      count.ok = count.ok + 1
      message("TRYING KFOLD ONE MORE TIME")
      if (count.ok > 10){
        message("YOU SHOULD THINK ABOUT STOPPING ME")
      }
    }
    all.samples <- setdiff(all.samples,cur.sel)
    k.samples.out[[i]] <- cur.sel
  }

  ## Cross-validate
  CV.stats <- list()
  rendu <- 1
  ## Cross validate
  for(ksi in k.samples.out){
    cat(sprintf("Rendu = %d/%d\n",rendu,length(k.samples.out)))
    rendu <- rendu + 1
    sel.train <- setdiff(1:ncol(D),ksi)
    sel.test <- ksi
     
    train.d <- D[,sel.train]
    test.d <- D[,sel.test]
    
    test.cl <- cl[sel.test]
    train.cl <- cl[sel.train]
    class.cur <- factor(train.cl)
    cur.ordering <- sample(1:nrow(train.d)) # random by default
    if (!random){
      from.idx=seq(1,nrow(train.d),ceiling(nrow(train.d)/10))
      to.idx=c((from.idx-1)[2:length(from.idx)],nrow(train.d))
      pvalues.tmp <- mclapply(1:length(from.idx),
                          function(i){
                            as.numeric(apply(train.d[from.idx[i]:to.idx[i],],
                                             1,function(x){
                                               cnts00 <- (x == 0 & class.cur == 0) * w
                                               cnts10 <- (x == 1 & class.cur == 0) * w
                                               cnts11 <- (x == 1 & class.cur == 1) * w
                                               cnts01 <- (x == 0 & class.cur == 1) * w
                                               
                                               chisq.test(cbind(c(cnts00,cnts10),c(cnts01,cnts11)))$p.value
                                             }))
                          },mc.cores=n.cores)
      pvalues <- rep(0,nrow(train.d))
      for (i in 1:length(from.idx)){
        pvalues[from.idx[i]:to.idx[i]] <- pvalues.tmp[[i]]
      }
      
      cur.ordering <- order(pvalues)
    }
    
    current.ml.stats <- mclapply(seq(1,20,1),function(sel.N){
      toRet <- list()
      feature.names <- rownames(cp.pam50$features)
      for (max.repeat in c(1)){
        
        sel.ids <- selIdsWithMaxConstraint(cur.ordering,feature.names[cur.ordering],sel.N,max.repeat)
        
        ## train
        nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.d[sel.ids,,drop=F])))
        prob.train <- predict(nb,data.frame(t(train.d[sel.ids,,drop=F])),type="raw")
        prob.test <- predict(nb,data.frame(t(test.d[sel.ids,,drop=F])),type="raw")
        prob.gq.test <- NA
        if (sum(is.na(ind.data)) == 0){
          prob.gq.test <- predict(nb,data.frame(t(comp.sel.pairs(ind.data,rownames(cp.pam50$features)[sel.ids]))),type="raw")
        }
        
        toRet[[length(toRet) + 1]] <- list(sel.N=sel.N,
                                           prob.train=prob.train,
                                           class.train=train.cl,
                                           prob.test=prob.test,
                                           class.test=test.cl,
                                           sel.ids=sel.ids,
                                           sel.pairs=rownames(cp.pam50$features)[sel.ids],
                                           prob.gq.test=prob.gq.test,
                                           max.repeat=max.repeat)
      }
      toRet
    })
    CV.stats[[length(CV.stats)+1]] <- list(sel.train=sel.train,
                                           sel.test=sel.test,
                                           current.ml.stats)
  }
  
  invisible(CV.stats)
}

## D represents exprs not pairs
cv.nb.cp.v2 <- function(D,
                        cl,K.FOLD,GeneName,
                        ind.data,seed=1234,n.cores=10,random=F,
                        w=rep(1,ncol(D))){
  require(ROCR)
  require(e1071)
  
  set.seed(seed)
  
  k.samples.out <- list()
  all.samples <- 1:ncol(D)

  ## select the samples to leave out during cross-validation
  if (K.FOLD == 1){
     k.samples.out[[1]] <- c()
  }
  else{
    for (i in 1:K.FOLD){
      ok <- F
      count.ok <- 0
      while(!ok){ # need both classes in both groups
        cur.sel <- sample(all.samples,floor(ncol(D)/K.FOLD))
        if (all(table(factor(cl[cur.sel],levels=c(0,1))) > 0) &
            all(table(factor(cl[setdiff(1:length(cl),cur.sel)],levels=c(0,1))) > 0)){
          break
        }
        count.ok = count.ok + 1
        message("TRYING KFOLD ONE MORE TIME")
        if (count.ok > 10){
          message("YOU SHOULD THINK ABOUT STOPPING ME")
        }
      }
      all.samples <- setdiff(all.samples,cur.sel)
      k.samples.out[[i]] <- cur.sel
    }
  }

  set.seed(as.numeric(Sys.time())) # We want the remaining to be really random

  ## Cross-validate
  CV.stats <- list()
  rendu <- 1
  ## Cross validate
  for(ksi in k.samples.out){
    cat(sprintf("Rendu = %d/%d\n",rendu,length(k.samples.out)))
    rendu <- rendu + 1
    sel.train <- setdiff(1:ncol(D),ksi)
    sel.test <- 1:ncol(D)
    sel.test = ksi
    if (K.FOLD == 1){
      sel.test = 1:ncol(D)
    }
     
    train.d <- D[,sel.train]
    test.d <- D[,sel.test]
    
    test.cl <- cl[sel.test]
    train.cl <- cl[sel.train]
    class.cur <- factor(train.cl)

    ## Only compute the p-values
    ## If random it will random generate p-values using runif
    ## NOTE this is clearly not going to follow a usual p-value like distribution
    pv.pairs <- compute.pairs.pv(train.d,train.cl,w[sel.train],random=random)
    cur.ordering <- order(pv.pairs[,3])
    feature.names <- paste(GeneName[pv.pairs[,1]],GeneName[pv.pairs[,2]],sep="<")
    feature.names <- feature.names[cur.ordering]
    
    current.ml.stats <- mclapply(seq(1,20,1),function(sel.N){
      toRet <- list()
      
      for (max.repeat in c(1)){
        
        sel.ids <- selIdsWithMaxConstraint(1:nrow(pv.pairs),feature.names,sel.N,max.repeat)
        sel.pairs <-  feature.names[sel.ids]
        train.pairs <- comp.sel.pairs(list(exprs=train.d,GeneName=GeneName),sel.pairs)
        test.pairs <- comp.sel.pairs(list(exprs=test.d,GeneName=GeneName),sel.pairs)
        
        ## train
        nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.pairs)))
        prob.train <- predict(nb,data.frame(t(train.pairs)),type="raw")
        prob.test <- predict(nb,data.frame(t(test.pairs)),type="raw")
        
        prob.gq.test <- NA
        if (sum(is.na(ind.data)) == 0){
          ind.pairs <- comp.sel.pairs(ind.data,feature.names[sel.ids])
          prob.gq.test <- predict(nb,data.frame(t(ind.pairs)),type="raw")
        }
        
        toRet[[length(toRet) + 1]] <- list(sel.N=sel.N,
                                           prob.train=prob.train,
                                           class.train=train.cl,
                                           prob.test=prob.test,
                                           class.test=test.cl,
                                           sel.ids=sel.ids,
                                           sel.pairs=sel.pairs,
                                           prob.gq.test=prob.gq.test,
                                           max.repeat=max.repeat)
      }
      toRet
    })
    CV.stats[[length(CV.stats)+1]] <- list(sel.train=sel.train,
                                           sel.test=sel.test,
                                           current.ml.stats)
  }
  
  invisible(CV.stats)
}

cv.nb.rank.ind <- function(D,cl,K.FOLD,feature.names,ind.data,N=100,seed=1234,n.cores=10,random=F){
  require(ROCR)
  require(e1071)
  
  set.seed(seed)
  
  k.samples.out <- list()
  all.samples <- 1:ncol(D)
  D.rank <- apply(D,2,rank)

  ## select the samples to leave out during cross-validation
  for (i in 1:K.FOLD){
    ok <- F
    count.ok <- 0
    while(!ok){ # need both classes in both groups
      cur.sel <- sample(all.samples,floor(ncol(D)/K.FOLD))
      if (all(table(factor(cl[cur.sel],levels=c(0,1))) > 0) &
          all(table(factor(cl[setdiff(1:length(cl),cur.sel)],levels=c(0,1))) > 0)){
        break
      }
      count.ok = count.ok + 1
      message("TRYING KFOLD ONE MORE TIME")
      if (count.ok > 10){
        message("YOU SHOULD THINK ABOUT STOPPING ME")
      }
    }
    all.samples <- setdiff(all.samples,cur.sel)
    k.samples.out[[i]] <- cur.sel
  }

  ## Cross-validate
  CV.stats <- list()
  rendu <- 1
  ## Cross validate
  for(ksi in k.samples.out){
    cat(sprintf("Rendu = %d/%d\n",rendu,length(k.samples.out)))
    rendu <- rendu + 1
    sel.train <- setdiff(1:ncol(D),ksi)
    sel.test <- ksi
    
    # sel genes
    train.cl <- cl[sel.train]
    train.d.rank <- D.rank[,sel.train]
    rank.diff <- rowMeans(train.d.rank[,train.cl==0]) - rowMeans(train.d.rank[,train.cl==1])
    rank.diff.idx <- order(rank.diff)
    symbols=feature.names
    selected.genes <- c(rank.diff.idx[1:(N/2)],rank.diff.idx[(-(N/2 - 1):0)+length(rank.diff.idx)])
    cp.pam50 <- compute.pairs(D[selected.genes,],n.cores=n.cores)
    sel.symbols.func <- symbols[selected.genes]
    rownames(cp.pam50$features) <- paste(sel.symbols.func[cp.pam50$indices[,1]],sel.symbols.func[cp.pam50$indices[,2]],sep="<")
    sorted.names <- as.character(lapply(strsplit(rownames(cp.pam50$features),"<"),function(x){x=sort(x);paste(x[1],x[2],sep=",")}))
    
    train.d <- cp.pam50$features[,sel.train]
    test.d <- cp.pam50$features[,sel.test]
    
    test.cl <- cl[sel.test]

    class.cur <- factor(train.cl)
    cur.ordering <- sample(1:nrow(train.d))
    if (!random){
      pvalues <- as.numeric(mclapply(1:nrow(train.d),
                                     function(i){
                                       x=train.d[i,];
                                       chisq.test(table(factor(x,levels=c(0,1)),class.cur))$p.value
                                     },mc.cores=n.cores))
      
      cur.ordering <- order(pvalues)
    }
    # remove duplicated pairs
    cur.ordering <- cur.ordering[!duplicated(sorted.names[cur.ordering])]
    
    current.ml.stats <- mclapply(seq(1,20,1),function(sel.N){
      toRet <- list()
      feature.names <- rownames(cp.pam50$features)
      #for (max.repeat in 1:sel.N){
      for (max.repeat in c(1)){
        
        sel.ids <- selIdsWithMaxConstraint(cur.ordering,feature.names[cur.ordering],sel.N,max.repeat)
        
        ## train
        nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.d[sel.ids,,drop=F])))
        prob.train <- predict(nb,data.frame(t(train.d[sel.ids,,drop=F])),type="raw")
        prob.test <- predict(nb,data.frame(t(test.d[sel.ids,,drop=F])),type="raw")
        prob.gq.test <- NA
        if (sum(is.na(ind.data)) == 0){
          prob.gq.test <- predict(nb,data.frame(t(comp.sel.pairs(ind.data,rownames(cp.pam50$features)[sel.ids]))),type="raw")
        }
        
        toRet[[length(toRet) + 1]] <- list(sel.N=sel.N,
                                           prob.train=prob.train,
                                           class.train=train.cl,
                                           prob.test=prob.test,
                                           class.test=test.cl,
                                           sel.ids=sel.ids,
                                           sel.pairs=rownames(cp.pam50$features)[sel.ids],
                                           prob.gq.test=prob.gq.test,
                                           max.repeat=max.repeat)
      }
      toRet
    })
    CV.stats[[length(CV.stats)+1]] <- list(sel.train=sel.train,
                                           sel.test=sel.test,
                                           current.ml.stats)
  }
  
  invisible(CV.stats)
}


cv.nb.rank <- function(D,cl,K.FOLD,feature.names,N=100,seed=1234,n.cores=10,random=F){
  require(ROCR)
  require(e1071)
  
  set.seed(seed)
  
  k.samples.out <- list()
  all.samples <- 1:ncol(D)
  D.rank <- apply(D,2,rank)

  ## select the samples to left out during cross-validation
  for (i in 1:K.FOLD){
    cur.sel <- sample(all.samples,floor(ncol(D)/K.FOLD))
    all.samples <- setdiff(all.samples,cur.sel)
    k.samples.out[[i]] <- cur.sel
  }

  ## Cross-validate
  CV.stats <- list()
  rendu <- 1
  for(ksi in k.samples.out){
    cat(sprintf("Rendu = %d/%d\n",rendu,length(k.samples.out)))
    rendu <- rendu + 1
    sel.train <- setdiff(1:ncol(D),ksi)
    sel.test <- ksi

    ## select the most variable genes
    train.cl <- cl[sel.train]
    train.d.rank <- D.rank[,sel.train]
    rank.diff <- rowMeans(train.d.rank[,train.cl==0]) - rowMeans(train.d.rank[,train.cl==1])
    rank.diff.idx <- order(rank.diff)
    symbols=feature.names
    selected.genes <- c(rank.diff.idx[1:(N/2)],rank.diff.idx[(-(N/2 - 1):0)+length(rank.diff.idx)])
    ## compute pairs on selected genes 
    cp.pam50 <- compute.pairs(D[selected.genes,],n.cores=n.cores)
    sel.symbols.func <- symbols[selected.genes]
    rownames(cp.pam50$features) <- paste(sel.symbols.func[cp.pam50$indices[,1]],sel.symbols.func[cp.pam50$indices[,2]],sep="<")
    sorted.names <- as.character(lapply(strsplit(rownames(cp.pam50$features),"<"),function(x){x=sort(x);paste(x[1],x[2],sep=",")}))
    
    train.d <- cp.pam50$features[,sel.train]
    test.d <- cp.pam50$features[,sel.test]
    
    test.cl <- cl[sel.test]

    class.cur <- factor(train.cl)
    cur.ordering <- sample(1:nrow(train.d)) # random ordering
    if (!random){
      pvalues <- as.numeric(mclapply(1:nrow(train.d),
                                     function(i){
                                       x=train.d[i,];
                                       chisq.test(table(factor(x,levels=c(0,1)),class.cur))$p.value
                                     }))
      
      cur.ordering <- order(pvalues)
    }
    
    cur.ordering <- cur.ordering[!duplicated(sorted.names[cur.ordering])]    
    current.ml.stats <- mclapply(1:20,function(sel.N){
      sel.ids <- cur.ordering[1:sel.N]
      
      ## train
      nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.d[sel.ids,,drop=F])))
      prob.train <- predict(nb,data.frame(t(train.d[sel.ids,,drop=F])),type="raw")
      prob.test <- predict(nb,data.frame(t(test.d[sel.ids,,drop=F])),type="raw")
      
      list(sel.N=sel.N,
           prob.train=prob.train,
           class.train=train.cl,
           prob.test=prob.test,
           class.test=test.cl,
           sel.ids=sel.ids,
           sel.pairs=rownames(cp.pam50$features)[sel.ids])
    })
    CV.stats[[length(CV.stats)+1]] <- list(sel.train=sel.train,
                                           sel.test=sel.test,
                                           current.ml.stats)
  }
  
  invisible(CV.stats)
}

keep.rank.most.variable <- function(dataset){
  to.keep <- rep(F,nrow(dataset$exprs))
  symbols <- as.character(dataset$GeneName)
  sd.rank <- apply(apply(dataset$exprs,2,rank),1,sd)
  sel.out <- order(-sd.rank)
  to.keep[sel.out[!duplicated(symbols[sel.out])]] <- T
  dataset$exprs <- dataset$exprs[to.keep,]
  dataset$GeneName <- dataset$GeneName[to.keep]
  dataset$probe.info <- dataset$probe.info[to.keep,]
  dataset$exprs.norm <- dataset$exprs.norm[to.keep,]
  dataset
}


getPairsTable <- function(d){
  stats.pairs <- list()
  stats.pairs.ids <- list()
  # Iterate k-fold
  for (i in 1:length(d)){
    #iterate the number of features
    for (j in 1:length(d[[i]][[3]])){
      if (length(stats.pairs) < j){
        stats.pairs[[j]] <- d[[i]][[3]][[j]]$sel.pairs
        stats.pairs.ids[[j]] <- d[[i]][[3]][[j]]$sel.ids
      }
      else{
        stats.pairs[[j]] <- c(stats.pairs[[j]],d[[i]][[3]][[j]]$sel.pairs)
        stats.pairs.ids[[j]] <- c(stats.pairs.ids[[j]],d[[i]][[3]][[j]]$sel.ids)
      }
    }
  }
  final.stats <- list()
  final.stats.ids <- list()
  for (i in 1:length(stats.pairs)){
    final.stats[[i]] <- sort(table(stats.pairs[[i]]),decreasing=T)/length(d)
  }
}


entrez2symbols <- function(pairs.chr,data){
  to.ret <- c()
  for (pi in strsplit(pairs.chr,"<")){
    eid1 <- pi[1]
    symbol1 <- data$probe.info$GeneName[match(eid1,data$probe.info$EntrezID)[1]]
    eid2 <- pi[2]
    symbol2 <- data$probe.info$GeneName[match(eid2,data$probe.info$EntrezID)[1]]
    
    to.ret <- c(to.ret,sprintf("%s<%s",symbol1,symbol2))
  }
  to.ret
}
# class, (name,ml.stats), CV
getPairsUsageStats <- function(pam50.k.stats.sel.rank,num.pairs=10,repeated=1,data){
  pairs.usage <- list()
  for (ci in pam50.k.stats.sel.rank){ ## iterate classes ie Basal, Her2, LumA, LumB
    pairs.counter.entrez <- c()
    for (cvi in ci[[2]]){ ## iterate CV
      pairs.counter.entrez <- c(pairs.counter.entrez,cvi[[3]][[num.pairs]][[repeated]]$sel.pairs)
    }
    pairs.counter.entrez <- sort(table(pairs.counter.entrez),decreasing=T)
    pairs.counter.symbol <- pairs.counter.entrez
    names(pairs.counter.symbol) <- entrez2symbols(names(pairs.counter.entrez),data)
    pairs.usage[[ ci[[1]] ]] <- list(pairs.counter.symbol=pairs.counter.symbol,
                                     pairs.counter.entrez=pairs.counter.entrez)
  }
  pairs.usage
}

getPairsUsageStats <- function(pam50.k.stats.sel.rank,num.pairs=10,repeated=1,data){
  pairs.usage <- list()
  for (ci in pam50.k.stats.sel.rank){ ## iterate classes ie Basal, Her2, LumA, LumB
    pairs.counter.entrez <- c()
    for (cvi in ci[[2]]){ ## iterate cross-validation
      pairs.counter.entrez <- c(pairs.counter.entrez,cvi[[3]][[num.pairs]][[repeated]]$sel.pairs)
    }
    pairs.counter.entrez <- sort(table(pairs.counter.entrez),decreasing=T)
    pairs.counter.symbol <- pairs.counter.entrez
    names(pairs.counter.symbol) <- entrez2symbols(names(pairs.counter.entrez),data)
    pairs.usage[[ ci[[1]] ]] <- list(pairs.counter.symbol=pairs.counter.symbol,
                                     pairs.counter.entrez=pairs.counter.entrez)
  }
  pairs.usage
}

getPairsCV <- function(pam50.k.stats.sel.rank,num.pairs,repeated,data){
  pairs.usage <- list()
  for (ci in pam50.k.stats.sel.rank){ ## iterate classes ie Basal, Her2, LumA, LumB
    pairs.usage.class <- c()
    for (cvi in ci[[2]]){ ## iterate cross-validation
      pairs.counter.entrez <- cvi[[3]][[num.pairs]][[repeated]]$sel.pairs
      pairs.counter.symbol <- entrez2symbols(pairs.counter.entrez,data)
      
      pairs.usage.class[[ length(pairs.usage.class) + 1 ]] <- list(pairs.counter.entrez=pairs.counter.entrez,
                                                                   pairs.counter.symbol=pairs.counter.symbol)
    }
    pairs.usage[[ ci$pamp50 ]] <- pairs.usage.class
  }
  pairs.usage
}

## Given a dataset with the format
## used in this code, select only some of the genes
focus.dataset <- function(dataset,sel){
  dataset$exprs <- dataset$exprs[sel,]
  dataset$exprs.norm <- dataset$exprs.norm[sel,]
  dataset$GeneName <- dataset$GeneName[sel]
  dataset$probe.info <- dataset$probe.info[sel,]
  invisible(dataset)
}

## Remove the duplicated Entrez by keeping the most higly expressed
## This is closer to the single sample selection
## D is a raw gene expression matrix rows == genes and columns patients
removeDuplicatedEntrezPerPatients <- function(D,EntrezID,probes){
  ## Maybe we have nothing to do already
  if (all(!duplicated(EntrezID))){
    return(list(dataset=D,EntrezID=EntrezID))
  }
  else{
    uniqEntrez <- sort(unique(EntrezID))
    newD <- matrix(0,nrow=length(uniqEntrez),ncol=ncol(D))
    if (!missing(probes)){
      sel.probes <- matrix("",nrow=length(uniqEntrez),ncol=ncol(D))
    }
    for (i in 1:ncol(D)){
      curD <- D[,i]
      curEntrez <- EntrezID
      oi <- order(curD,decreasing=T) ## order by raw expression
      curD <- curD[oi]
      curEntrez <- curEntrez[oi]
      cur.sel <- !duplicated(curEntrez) ## remove duplicated
      curD <- curD[cur.sel]
      curEntrez <- curEntrez[cur.sel]
      reorder <- match(uniqEntrez,curEntrez)
      newD[,i] <- curD[reorder]
      if (!missing(probes)){
        sel.probes[,i] <- probes[oi][cur.sel][reorder]
      }
    }
    colnames(newD) <- colnames(D)
    if (!missing(probes)){
      colnames(sel.probes) <- colnames(D)
      return(list(dataset=newD,EntrezID=uniqEntrez,probes=sel.probes))
    }
    else{
      return(list(dataset=newD,EntrezID=uniqEntrez))
    }
  }
}

## Merge datasets
## What's important :
## exprs, GeneName, probe.info, demo, name, chip, pam50.genefu.scale (cor, subtype)
merge.datasets <- function(list.of.datasets,add.dataset=T,clinical="clinical",use.pam50="pam50.parker"){
  ## Merge names
  name <- ""
  chip <- ""
  demo <- c()
  demo.tag <- c()
  inter.gene.id <- c()
  for (di in list.of.datasets){
    if (length(demo.tag) == 0){
      demo.tag <- colnames(di[[clinical]])
      inter.gene.id <- unique(di$probe.info$EntrezID)
    }
    else{
      demo.tag <- intersect(demo.tag,colnames(di[[clinical]]))
      inter.gene.id <- intersect(inter.gene.id,di$probe.info$EntrezID)
    }
  }
  message(sprintf("Keeping %d Entrez genes and %d clinicals",length(inter.gene.id),length(demo.tag)))
  
  ## Take probe.info from the first dataset
  first.d.idxs <- match(inter.gene.id,list.of.datasets[[1]]$probe.info$EntrezID)
  probe.info <- list.of.datasets[[1]]$probe.info[first.d.idxs,]
  #GeneName <- list.of.datasets[[1]]$GeneName[first.d.idxs]
  pam50.genefu.scale <- list(cor=c(),subtype=c())
  exprs <- c()
  time <- c()
  event <- c()
  for (di in list.of.datasets){
    clean.d <- removeDuplicatedEntrezPerPatients(di$exprs,di$probe.info$EntrezID)
    exprs <- cbind(exprs,clean.d$dataset[match(inter.gene.id,clean.d$EntrezID),])
    if (name == ""){
      name <- di$name
      chip <- di$chip
    }
    else{
      name <- paste(name,di$name,sep="_")
      chip <- paste(chip,di$chip,sep="_")
    }
    if (add.dataset){
      demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],dataset=as.character(di$name)))
    }
    else{
       demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],dataset="NO.WEIGHT"))
    }
    ## demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],
    ##                         dataset=paste("ch",gsub(" ","",as.character(di$chip)),sep=".")))
    ## demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],
    ##                         dataset="ALL.THE.SAME"))
    ## demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],dataset="NO.WEIGHT")
    pam50.genefu.scale[["cor"]] <- rbind(pam50.genefu.scale[["cor"]],di[[use.pam50]]$cor)
    pam50.genefu.scale[["subtype"]] <- c(pam50.genefu.scale[["subtype"]],di[[use.pam50]]$subtype)
    cur.time <- rep(NA,ncol(di$exprs))
    cur.event <- rep(NA,ncol(di$exprs))
    if ("event" %in% names(di) & "time" %in% names(di)){
      cur.time <- di$time
      cur.event <- di$event
    }
    time <- c(time,cur.time)
    event <- c(event,cur.event)
  }
  
  invisible(list(name=name,
                 chip=chip,
                 GeneName=probe.info[,"EntrezID"],
                 probe.info=probe.info,
                 demo=demo,
                 time=time,
                 event=event,
                 exprs=exprs,
                 pam50.genefu.scale=pam50.genefu.scale))
}

getWeights <- function(tags){
  tags=as.character(tags)
  dataset.cnts <- table(tags)
  equal.w <- length(tags)/length(dataset.cnts)
  as.numeric((equal.w/dataset.cnts)[tags])
}

## Merge datasets
## What's important :
## exprs, GeneName, probe.info, demo, name, chip, pam50.genefu.scale (cor, subtype)
merge.datasets.2 <- function(list.of.datasets,add.dataset=T,clinical="clinical",use.pam50="pam50.parker",inter.gene.id=c()){
  ## Merge names
  name <- ""
  chip <- ""
  demo <- c()
  demo.tag <- c()
  
  for (di in list.of.datasets){
    if (length(inter.gene.id) == 0){
      inter.gene.id <- unique(di$probe.info$EntrezID)
    }
    else{
      inter.gene.id <- intersect(inter.gene.id,di$probe.info$EntrezID)
    }
    
    if (length(demo.tag) == 0){
      demo.tag <- colnames(di[[clinical]]) 
    }
    else{
      demo.tag <- intersect(demo.tag,colnames(di[[clinical]]))
      
    }
  }
  message(sprintf("Keeping %d Entrez genes and %d clinicals",length(inter.gene.id),length(demo.tag)))
  
  ## Take probe.info from the first dataset
  first.d.idxs <- match(inter.gene.id,list.of.datasets[[1]]$probe.info$EntrezID)
  probe.info <- list.of.datasets[[1]]$probe.info[first.d.idxs,]
  #GeneName <- list.of.datasets[[1]]$GeneName[first.d.idxs]
  pam50.genefu.scale <- list(cor=c(),subtype=c())

  raw.exprs <- c()
  norm.exprs <- c()
  
  time <- c()
  event <- c()
  for (di in list.of.datasets){
    raw.clean.d <- removeDuplicatedEntrezPerPatients(di$exprs,di$probe.info$EntrezID)
    raw.exprs <- cbind(raw.exprs,raw.clean.d$dataset[match(inter.gene.id,raw.clean.d$EntrezID),])

    norm.clean.d <- removeDuplicatedEntrezPerPatients(di$exprs.norm,di$probe.info$EntrezID)
    norm.exprs <- cbind(norm.exprs,norm.clean.d$dataset[match(inter.gene.id,norm.clean.d$EntrezID),])
    
    if (name == ""){
      name <- di$name
      chip <- di$chip
    }
    else{
      name <- paste(name,di$name,sep="_")
      chip <- paste(chip,di$chip,sep="_")
    }
    if (add.dataset){
      demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],dataset=as.character(di$name)))
    }
    else{
       demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],dataset="NO.WEIGHT"))
    }
    ## demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],
    ##                         dataset=paste("ch",gsub(" ","",as.character(di$chip)),sep=".")))
    ## demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],
    ##                         dataset="ALL.THE.SAME"))
    ## demo <- rbind(demo,cbind(di[[clinical]][,demo.tag],dataset="NO.WEIGHT")
    pam50.genefu.scale[["cor"]] <- rbind(pam50.genefu.scale[["cor"]],di[[use.pam50]]$cor)
    pam50.genefu.scale[["subtype"]] <- c(pam50.genefu.scale[["subtype"]],di[[use.pam50]]$subtype)
    cur.time <- rep(NA,ncol(di$exprs))
    cur.event <- rep(NA,ncol(di$exprs))
    if ("event" %in% names(di) & "time" %in% names(di)){
      cur.time <- di$time
      cur.event <- di$event
    }
    time <- c(time,cur.time)
    event <- c(event,cur.event)
  }
  
  invisible(list(name=name,
                 chip=chip,
                 GeneName=probe.info[,"EntrezID"],
                 probe.info=probe.info,
                 demo=demo,
                 time=time,
                 event=event,
                 exprs=norm.exprs,
                 raw.exprs=raw.exprs,
                 norm.exprs=norm.exprs,
                 pam50.genefu.scale=pam50.genefu.scale))
}
