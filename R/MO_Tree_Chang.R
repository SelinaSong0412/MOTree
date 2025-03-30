#' Tree-based Reinforcement Learning for estimating optimal DTR.
#'
#' a tree-based reinforcement learning (T-RL) method to directly
#' estimate optimal DTRs in a multi-stage multi-treatment setting. At
#' each stage, T-RL builds an unsupervised decision tree that directly handles
#' the problem of optimization with multiple treatment comparisons, through a
#' purity measure constructed with augmented inverse probability weighted estimators.
#'
#' @param Ys A matrix of outcome of interest, each column record one outcome.
#' @param A A vector of observed treatment options.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param pis.hat Estimated propensity score matrix.
#' @param m.method Method for calculating estimated conditional mean.
#' @param mus.reg Regression-based conditional mean outcome.
#' @param depth Maximum tree depth.
#' @param lambda.pct Minimal percent change in purity measure for split.
#' @param minsplit Minimal node size.
#' @param lookahead Whether or not to look into a further step of splitting to find the best split.
#'
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#'
#' @export
# 
# mus.hat = list(mu1.hat,mu2.hat)
# pis.hat=NULL
# mus.reg=NULL
# depth=3
# minsplit=20
# w_vec = cbind(w0,1-w0)
# weight_combine = "max"
# lambda = 0.01
# 
# mus.hat = list(Y1,Y2)
# H = matrix(x,ncol=1)
# pis.hat=NULL
# mus.reg=NULL
# depth=3
# minsplit=20
# w_vec = cbind(w0,1-w0)
# weight_combine = "max"
# lambda = 0.01
DTRtree<-function(Ys,
                  A,H,
                  pis.hat=NULL,
                  mus.reg=NULL,depth=5,
                  minsplit=20,
                  w_vec,
                  m.method = c("AIPW", "randomForest"),
                  weight_combine,lambda = 0.01){
  # initialization
  # indicator for subset data
  n<-nrow(Ys) #number of people
  n.obj<-ncol(Ys) # number of outcomes
  I.node<-rep(1,n) #indicator of nodes
  class.A<-sort(unique(A))
  output<-matrix(NA,1,5)
  colnames(output)<-c("node","X","cutoff","mEy","group")
  
  # YS update the code for calculating mus.hat fo multiple objs
  mus.hat <- lapply(1:n.obj, function(i) matrix(0, nrow = n, ncol = length(class.A)))
  for (d in c(1:n.obj)) {
    Y <- Ys[,d]
    # estimate mus.hat if not given
    if(m.method[1]=="AIPW"){
      # estimate propensity matrix if not given, using all data
      # same propensity for all subset data
      if(is.null(pis.hat)) pis.hat<-M.propen(A=A,Xs=H)
      if(is.null(mus.reg)) mus.reg<-Reg.mu(Y=Y,As=A,H=H)$mus.reg
      mus.hat.fit<-mus.AIPW(Y=Y,A=A,pis.hat=pis.hat,mus.reg=mus.reg)
    } else if(m.method[1]=="randomForest"){
      # require(randomForest)
      RF<-randomForest(Y~., data=data.frame(A,H))
      mus.hat.fit<-matrix(NA,n,length(class.A))
      for(i in 1L:length(class.A)) mus.hat.fit[,i]<-predict(RF,newdata=data.frame(A=rep(class.A[i],n),H))
    } else{
      stop("The method for estimating conditional means is not available!")
    }
    mus.hat[[d]] <- mus.hat.fit
  }
  
  for (k in 1L:depth) { # depth is the most number of split to reach one terminal node
    output <- rbind(output,matrix(NA,2^k,5) ) # 2^k??????? originally output=1*5
    output[,1] <- 1L:(2^(k+1)-1) # this does not equal to the number of rows(2^k+1)
    
    # The 1st split
    if (k==1L) {
      # apply lookahead to the first split, the most important split
      # only to first split so as to save computation time
      # use a larger minsplit for the first split
      

      # Don't consider look ahead for now
      # if (lookahead) {
      #   best.H.1<-best.H.lh(H=H,A=A,mus.hat=mus.hat,minsplit=0.15*n)
      # } else {
      #   best.H.1<-best.H(H=H,A=A,mus.hat=mus.hat,minsplit=minsplit) # choose the best covariate for splitting
      # }
      best.H.1<-best.H(H=H,A=A,mus.hat=mus.hat,
                       minsplit=minsplit,
                       w_vec=w_vec,
                       weight_combine = weight_combine) # choose the best covariate for splitting
      
      
      if (is.null(best.H.1)==F) { # meet the split criteria
        output[k,-1]<-c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1,NA)
        I.node[I.node==k & H[,best.H.1$X,drop=FALSE] <= best.H.1$X.subset]<-2*k
        output[2*k,-1]<-c(NA,NA,NA,2*k)
        I.node[I.node==k & H[,best.H.1$X,drop=FALSE] > best.H.1$X.subset]<-2*k+1
        output[2*k+1,-1]<-c(NA,NA,NA,2*k+1)
      } else { # NOT meet the split criteria (terminate node)
        break # Finish the splitting at this node
      }
      # The jth split
    } else {
      for(j in (2^(k-1)):(2^k-1)) {
        if(!is.na(output[trunc(j/2),2])){
          # best.H.j<-best.H(H=H[I.node==j,],A=A[I.node==j],
          #                  mus.hat=mus.hat[I.node==j,],minsplit=minsplit)
          tmp_mus.hat = mus.hat
          for(ii in 1:length(mus.hat)){
            tmp_mus.hat[[ii]] = mus.hat[[ii]][I.node==j,]
          }
          best.H.j<-best.H(H=H[I.node==j,,drop=FALSE],A=A[I.node==j],
                           mus.hat=tmp_mus.hat,minsplit=minsplit,
                           w_vec=w_vec,
                           weight_combine =weight_combine) # choose the best covariate for splitting
          
          if(is.null(best.H.j)==F && best.H.j$mEy.opt1 > lambda){#&& best.H.j$mEy.opt1>output[trunc(j/2),4]+lambda
            output[j,-1]<-c(best.H.j$X, best.H.j$X.subset, best.H.j$mEy.opt1,NA)
            I.node[I.node==j & H[,best.H.j$X,drop=FALSE] <= best.H.j$X.subset]<-2*j
            output[2*j,-1]<-c(NA,NA,NA,2*j)
            I.node[I.node==j & H[,best.H.j$X,drop=FALSE] > best.H.j$X.subset]<-2*j+1
            output[2*j+1,-1]<-c(NA,NA,NA,2*j+1)
          }
        }
      }
      if(sum(is.na(output[(2^(k-1)):(2^k-1),2]))==2^(k-1)) break
    }
  }
  #output<-output[!is.na(output[,2]),]
  #output<-output[!is.na(output[,2]) | !is.na(output[,5]),]
  output = as.data.frame(output)
  
  # YS update:
  output$group <- ave(output$group, is.na(output$group), 
                       FUN = function(x) ifelse(is.na(x), NA, seq_along(x)))
  
  return(output)
}



