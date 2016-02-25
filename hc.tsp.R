##
## Necessary functions to train AIMS models.
##
## Author : Eric R. Paquet (eric.r.paquet@gmail.com)
##
## Copyright McGill University, 2015
##

source("ml.ssp.toolbox.R")
require(e1071)

## This function will sort pairs
## "2 < 1" -> "1 < 2"
sort.pairs <- function(pairs){
  as.character(lapply(strsplit(pairs,"<"),function(x){x <- sort(x);paste(x[1],"<",x[2],sep="")}))
}

## this function will compute a TSP model
## similar to what was proposed by Geman (ref bioinformatics paper)
## except at every nodes we are using a naive bayes classifier.
##
## D = a dataset rows == genes and columns == samples
## cl = a character vector of classes
## w = the weights of every samples (used when D is composed of multiple datasets)
## k = how many features tested default 1:20
hc.tsp <- function(D,cl,GeneName,cl.order=NULL,w=rep(1,ncol(D)),k=1:20){
  
  if (ncol(D) != length(cl)){
    stop(sprintf("Need the same number of rows in D and element in cl currently : %d and %d",nrow(D),length(cl)))
  }
  if (!all(!duplicated(GeneName))){
    stop("Provide me with aggregated GeneName now you have duplication (try ?aggregate)")
  }
  
  cl <- as.character(cl)
  message(sprintf("hc.tsp working on %d classes",length(unique(cl))))
  ## Get the most abundant class
  cl.cnt.sort <- sort(table(cl),decreasing=T)
  sorted.classes <- names(cl.cnt.sort)
  if (!is.null(cl.order)){
    sorted.classes <- cl.order
  }
  ## generate a list of NB classifiers
  processed.classes <- c()
  cur.hctsp <- list()
  all.uniq.sel.pairs <- c()
  selected.pairs.list <- list()
  for (cur.class in sorted.classes[1:(length(sorted.classes)-1)]){
    ## current class versus remaining
    processed.classes <- c(processed.classes,cur.class)
    other.classes <- setdiff(names(cl.cnt.sort),processed.classes)
    message(sprintf("%s versus ALL remaining",cur.class))
    c1 <- which(cl == cur.class)
    c2 <- which(cl %in% other.classes)

    train.d <- cbind(D[,c1],D[,c2])
    train.cl <- c(rep(1,length(c1)),
                  rep(0,length(c2)))
    
    pv.pairs <- compute.pairs.pv(train.d,
                                 train.cl,
                                 c(w[c1],w[c2]),
                                 random=F)

    cur.ordering <- order(pv.pairs[,3])
    feature.names <- paste(GeneName[pv.pairs[,1]],GeneName[pv.pairs[,2]],sep="<")
    feature.names <- feature.names[cur.ordering]

    sel.ids <- selIdsWithMaxConstraint(1:nrow(pv.pairs),feature.names,max(k),1)
    sel.pairs <- feature.names[sel.ids]
    
    selected.pairs.list[[cur.class]] <- sel.pairs[1:max(k)]
    
    all.uniq.sel.pairs <- unique(c(all.uniq.sel.pairs,sort.pairs(sel.pairs)))
    ## compute the most significant pairs here
    ## generate one NB for every k
    train.cl <- c(rep(cur.class,length(c1)),
                  rep("ELSE",length(c2)))
    list.of.nb <- list()
    for (ki in k){
      ## train
      train.pairs <- comp.sel.pairs(list(exprs=train.d,GeneName=GeneName),sel.pairs[1:ki])
      nb <- naiveBayes(train.cl ~ ., data=data.frame(t(train.pairs)))
      list.of.nb[[ki]] <- list(k=ki,nbc=nb)
    }
    cur.hctsp[[length(cur.hctsp) + 1]] <- list(class=cur.class,lnb=list.of.nb)
  }
  cur.hctsp[[length(cur.hctsp) + 1]] <- list(class=sorted.classes[length(sorted.classes)],lnb=NULL)
  to.ret <- list(hctsp=cur.hctsp,all.pairs=all.uniq.sel.pairs,k=k,sorted.classes=sorted.classes,selected.pairs.list=selected.pairs.list)
  invisible(to.ret)
}

get.pair <- function(tag){
  gsub("\\.","<",gsub("^X","",tag))
}

## This will reweigth the NBC using pairs + weights
correctNB <- function(nb,train.pairs,cl,w){
  ## Need to correct the nb
  
  for(ri in names(nb$apriori)){
    sel.samples <- cl == ri
    cur.w <- w[sel.samples]
    sum.ri <- sum(cur.w)
    ## Correct the prior
    nb$apriori[ri] = sum.ri

    ## Correct the tables
    for(ti in names(nb$tables)){
      cur.pairs <- get.pair(ti)
      cur.pairs.values <- train.pairs[cur.pairs,sel.samples]
      for (ci in colnames(nb$tables[[ti]])){
        nb$tables[[ti]][ri,ci] <- sum(cur.w[cur.pairs.values == ci])/sum.ri
      }
    }
  }
  nb
}

##apam50.2 <- correctNB(apam50$one.vs.all.tsp[[18]],
##                      yo,merged.dataset$pam50.genefu.scale$subtype,getWeights(merged.dataset$demo$dataset))

## this function will compute a TSP model
## similar to what was proposed by Goeman (ref bioinformatics paper)
## except at every nodes we are using a naive bayes classifier.
##
## D = a dataset rows == genes and columns == samples
## cl = a character vector of classes
## w = the weights of every samples (used when D is composed of multiple datasets)
## k = how many features tested default 1:20
one.vs.all.tsp <- function(D,cl,GeneName,w=rep(1,ncol(D)),k=1:20,random=F){
  
  if (ncol(D) != length(cl)){
    stop(sprintf("Need the same number of rows in D and element in cl currently : %d and %d",nrow(D),length(cl)))
  }
  if (!all(!duplicated(GeneName))){
    stop("Provide me with aggregated GeneName now you have duplication (try ?aggregate)")
  }
  
  cl <- as.character(cl)
  message(sprintf("hc.tsp working on %d classes",length(unique(cl))))

  ## first find all the interesting pairs
  good.pairs <- list()
  all.pairs <- c()
  selected.pairs.list <- list()
  ## 
  for (cur.class in unique(cl)){
    ## current class versus remaining
    message(sprintf("%s versus ALL",cur.class))
    c1 <- which(cl == cur.class)
    c2 <- setdiff(1:length(cl),c1)

    train.d <- cbind(D[,c1],D[,c2])
    train.cl <- c(rep(1,length(c1)),
                  rep(0,length(c2)))
    
    pv.pairs <- compute.pairs.pv(train.d,
                                 train.cl,
                                 c(w[c1],w[c2]),
                                 random=random)

    cur.ordering <- order(pv.pairs[,3])
    feature.names <- paste(GeneName[pv.pairs[,1]],GeneName[pv.pairs[,2]],sep="<")
    feature.names <- feature.names[cur.ordering]

    sel.ids <- selIdsWithMaxConstraint(1:nrow(pv.pairs),feature.names,max(k),1)
    sel.pairs <- feature.names[sel.ids]
    good.pairs[[cur.class]] <- sel.pairs
    selected.pairs.list[[cur.class]] <- sel.pairs[1:max(k)]
    all.pairs <- unique(c(all.pairs,sort.pairs(sel.pairs)))
  }
  ## then combine
  tsps <- list()
  for (ki in k){
    message("Generating NBC for k =",ki)
    ## first get all the unique pairs for this k
    cur.pairs <- c()
    for (gpi in good.pairs){
      cur.pairs <- unique(c(cur.pairs,sort.pairs(gpi[1:ki])))
    }
    train.pairs <- comp.sel.pairs(list(exprs=D,GeneName=GeneName),cur.pairs)
    nb <- naiveBayes(cl ~ ., data=data.frame(t(train.pairs)))
    ## Correct NB
    ## nb <- correctNB(nb,train.pairs,cl,w)
    tsps[[ki]] <- nb
  }
  
  to.ret <- list(all.pairs=all.pairs,
                 k=k,
                 one.vs.all.tsp=tsps,
                 selected.pairs.list=selected.pairs.list)
  invisible(to.ret)
}

predict.one.vs.all.tsp <- function(D,GeneName,one.vs.all.tsp){
  ## Need to add some QC
  ## First compute
  train.pairs <- comp.sel.pairs(list(exprs=D,GeneName=GeneName),one.vs.all.tsp$all.pairs)

  classes <- matrix("",ncol=length(one.vs.all.tsp$k),nrow=ncol(D))
  prob <- matrix(0,ncol=length(one.vs.all.tsp$k),nrow=ncol(D))
  colnames(classes) <- colnames(prob) <- as.character(one.vs.all.tsp$k)
  rownames(classes) <- rownames(prob) <-colnames(D)
  nb.d <- data.frame(t(train.pairs))
  all.probs <- list()
  for (ki in one.vs.all.tsp$k){
    message(sprintf("Current k = %d",ki))
    prob.train <- predict(one.vs.all.tsp$one.vs.all.tsp[[ki]],nb.d,type="raw")
    cur.cl <- apply(prob.train,1,function(prob.cur){colnames(prob.train)[which.max(prob.cur)]})
    cur.prob <- apply(prob.train,1,function(prob.cur){max(prob.cur)})
    prob[,as.character(ki)] <- cur.prob
    classes[,as.character(ki)] <- cur.cl
    all.probs[[as.character(ki)]] <- prob.train
  }
  invisible(list(cl = classes,prob = prob,all.probs = all.probs))
}

predict.hctsp <- function(D,GeneName,hctsp){
  ## Need to add some QC
  ## First compute
  train.pairs <- comp.sel.pairs(list(exprs=D,GeneName=GeneName),hctsp$all.pairs)

  classes <- matrix("",ncol=length(hctsp$k),nrow=ncol(D))
  prob <- matrix(0,ncol=length(hctsp$k),nrow=ncol(D))
  colnames(classes) <- colnames(prob) <- as.character(hctsp$k)
  rownames(classes) <- rownames(prob) <-colnames(D)
  
  for (ki in hctsp$k){
    message(sprintf("Current k = %d",ki))
    remaining <- rep(T,ncol(D))
    for (i in 1:length(hctsp$hctsp)){
      message(sprintf("Remaining %d / %d",sum(remaining),ncol(D)))
      if (!is.null(hctsp$hctsp[[i]]$lnb)){
        prob.train <- predict(hctsp$hctsp[[i]]$lnb[[ki]]$nbc,data.frame(t(train.pairs[,remaining])),type="raw")
        greater.prob <- prob.train[,hctsp$hctsp[[i]]$class] > prob.train[,"ELSE"]
        classes[remaining,as.character(ki)][greater.prob] <- hctsp$hctsp[[i]]$class
        prob[remaining,as.character(ki)][greater.prob] <- prob.train[greater.prob,hctsp$hctsp[[i]]$class]
        prob[remaining,as.character(ki)][!greater.prob] <- prob.train[!greater.prob,"ELSE"] ## mainly just to set up the prob on the last class in the tree
        remaining <- classes[,as.character(ki)] == ""
      }
      else{
        classes[classes[,as.character(ki)] == "",as.character(ki)] <- hctsp$hctsp[[i]]$class
      }
    }
  }
  list(cl=classes,prob=prob)
}

get.acc <- function(pred.cl,cl,w=rep(1,length(cl))){
  to.ret <- list()
  to.ret[["ALL"]] <- as.numeric(apply(pred.cl,2,function(cur.cl)simple.acc(cur.cl,cl,w)))
  for (i in sort(unique(cl))){
    to.ret[[i]] <- as.numeric(apply(pred.cl,2,function(cur.cl){
      cur.sel <- cl == i;
      simple.acc(cur.cl[cur.sel],cl[cur.sel],w[cur.sel])}
                                    ))
  }
  to.ret
}

cv.hctsp <- function(D,cl,GeneName,k.fold,k=1:20,w=rep(1,ncol(D)),train.func=hc.tsp,pred.func=predict.hctsp){
  require(multicore)
  k.split <- list() ## will contain the samples to leave out
  k.size <- floor(ncol(D)/k.fold)
  remaining <- 1:ncol(D)
  set.seed(1234)
  for (i in 1:k.fold){
    cur.sel <- sample(remaining,k.size)
    if (i == 1 & length(cur.sel) >= 6){
      if(!all(cur.sel[1:6] == c(560,3064,2999,3068,4236,3150))){
        message(sprintf("##############\n#############\n\n RANDOM SEEDING : [[[ NOT ]]] OK\n\n ncol = %d, k.size = %d #############\n###########",ncol(D),k.size))
      }else{
        message("##############\n#############\n\n RANDOM SEEDING : OK \n\n#############\n###########")
      }
    }
    remaining <- setdiff(remaining,cur.sel)
    k.split[[length(k.split) + 1]] <- cur.sel
  }
  k.fold.stats <- mclapply(1:length(k.split),function(ki){
    message(paste("Start K-fold ",ki,sep=""))
    ## train
    to.keep <- setdiff(1:ncol(D),k.split[[ki]])
    mm <- train.func(D=D[,to.keep,drop=F],
                     cl=cl[to.keep],
                     GeneName=GeneName,
                     w=w[to.keep],
                     k=k)
    ## test
    train.pred <- pred.func(D[,to.keep,drop=F],GeneName,mm)
    test.pred <- pred.func(D[,k.split[[ki]],drop=F],GeneName,mm)
    message(paste("End K-fold ",ki,sep=""))
    list(train.pred=train.pred,
         train.acc=get.acc(train.pred$cl,cl[to.keep],w=w[to.keep]),
         test.pred=test.pred,
         test.acc=get.acc(test.pred$cl,cl[k.split[[ki]]],w=w[k.split[[ki]]]),
         to.keep=to.keep,
         to.leave=k.split[[ki]],
         nb.model=mm)
  },mc.cores=10)
  k.fold.stats
}

compute.pairs.pv <- function(d,cl,w,func=smaller,random=F,n.cores=10){
  require(Rgtsp)
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
    X <- t(d)
    n <- min(nrow(d)*(nrow(d) + 1)/2 - nrow(d),500000)
    n = min(max(n, 1), 0.5 * ncol(X) * (ncol(X) - 1))
    i = rep(0, n)
    j = rep(0, n)
    s = rep(0, n)
    z = list(I = i, J = j, S = s) 
    message("Rgtsp starts computing")
    z = .C("RGTSP_tsp_NW", as.double(X), as.integer(cl), 
      as.double(w), as.integer(ncol(X)), as.integer(nrow(X)), 
      i = as.integer(i), j = as.integer(j), s = as.double(s), 
      n = as.integer(n))
    message("Rgtsp ends computing")
    cbind(z$i,z$j,1-z$s)
  }
}

cv.plot <- function(cv.out,train.plot=T){
  stats.train <- list()
  stats.test <- list()
  for (ki in cv.out){
    for (ki.cl in names(ki$train.acc)){
      if (!(ki.cl %in% names(stats.train))){
        stats.train[[ki.cl]] <- c()
      }
      stats.train[[ki.cl]] <- rbind(stats.train[[ki.cl]],ki$train.acc[[ki.cl]])
    }
    for (ki.cl in names(ki$test.acc)){
      if (!(ki.cl %in% names(stats.test))){
        stats.test[[ki.cl]] <- c()
      }
      stats.test[[ki.cl]] <- rbind(stats.test[[ki.cl]],ki$test.acc[[ki.cl]])
    }
  }

  ## for(i in names(stats.test)){
  ##   boxplot(stats.test[[i]],ylab=sprintf("Accuracy %d-fold CV",nrow(stats.test[[i]])),
  ##           xlab="K (Number of pairs per individual subtypes)",
  ##           main=i,ylim=c(0,1))
  ## }
  
  ## for(i in names(stats.train)){
  ##   boxplot(stats.test[[i]],ylab="Accuracy training sets",
  ##           xlab="K (Number of pairs per individual subtypes)",
  ##           main=i,ylim=c(0,1))
  ## }
  cv.stats.train.mean <- lapply(stats.train,function(ci){apply(ci,2,mean)})
  cv.stats.train.sd <- lapply(stats.train,function(ci){apply(ci,2,sd)})
  cv.stats.test.mean <- lapply(stats.test,function(ci){apply(ci,2,mean)})
  cv.stats.test.sd <- lapply(stats.test,function(ci){apply(ci,2,sd)})
  for (ni in names(cv.stats.train.mean)){
    matplot(cbind(cv.stats.train.mean[[ni]],cv.stats.test.mean[[ni]],
                  cv.stats.train.mean[[ni]]+1.96*cv.stats.train.sd[[ni]]/sqrt(length(cv.out)),cv.stats.test.mean[[ni]]+1.96*cv.stats.test.sd[[ni]]/sqrt(length(cv.out)),
                  cv.stats.train.mean[[ni]]-1.96*cv.stats.train.sd[[ni]]/sqrt(length(cv.out)),cv.stats.test.mean[[ni]]-1.96*cv.stats.test.sd[[ni]]/sqrt(length(cv.out))),
            ylim=c(0,1),lty=c(1,1,2,2,2,2),col=c("green","red"),lwd=c(3,3,1,1,1,1),main=ni,ylab="Accuracy",xlab="k (number of pairs)",type="l")
    abline(v=1:length(cv.stats.train.mean[[ni]]),lty=2,col="gray",lwd=.5)
    legend("bottomright",c("train",sprintf("%d-fold-CV",length(cv.out))),col=c("green","red"),lwd=3)
  }
  invisible(list(train.mean=cv.stats.train.mean,test.mean=cv.stats.test.mean,
                 train.sd=cv.stats.train.sd,test.sd=cv.stats.test.sd))
}


## test.hc.tsp <- function(){
## require(genefu)
## data(pam50)

## load("/data/epaquet/ngssp/in/TRAIN.datasets.GUEDJ.TCGA.CURTIS.TCGA.EXPO.McGill.GQ.RData.RData")

## sel.gn <- data$probe.info$EntrezID %in% as.character(pam50.scale$centroids.map$EntrezGene.ID)
## D <- data$exprs[sel.gn,]
## Entrez <- data$probe.info$EntrezID[sel.gn]
## col.D <- removeDuplicatedEntrezPerPatients(D,Entrez)

## ## cv <- cv.hctsp(D=col.D$dataset,cl=data$demo$PAM50,GeneName=col.D$EntrezID,k.fold=10,w=getWeights(data$demo$dataset),k=seq(1,20,2))

## cv <- cv.hctsp(D=col.D$dataset,cl=data$demo$PAM50,GeneName=col.D$EntrezID,k.fold=10,w=getWeights(data$demo$dataset),k=seq(1,20,2),
##                train.func=one.vs.all.tsp,pred.func=predict.one.vs.all.tsp)


## pdf("cv.pam50.ova.chi.pdf")
## cv.plot(cv)
## dev.off()

## intrinsic <- setdiff(as.character(read.delim("~/backup/HighlySensitivePathway/src/exp-02-ml-ssp/intrinsic.gene.list.clean.txt",
##                                              sep="\t",stringsAsFactors=F,header=F)[,2]),as.character(pam50.scale$centroids.map$EntrezGene.ID))

## sel.gn <- data$probe.info$EntrezID %in% intrinsic
## D <- data$exprs[sel.gn,]
## Entrez <- data$probe.info$EntrezID[sel.gn]
## col.D <- removeDuplicatedEntrezPerPatients(D,Entrez)
## sel.rand <- 1:nrow(col.D$dataset)
## cv <- cv.hctsp(D=col.D$dataset[sel.rand,],cl=data$demo$PAM50,GeneName=col.D$EntrezID[sel.rand],k.fold=10,w=getWeights(data$demo$dataset),k=seq(1,20,4)
##                ,
##                train.func=one.vs.all.tsp,pred.func=predict.one.vs.all.tsp)
                   

## pdf("cv.intrinsic.ova.pdf")
## cv.plot(cv)
## dev.off()
## }
