#' Estimate propensity score by a multinomial model.
#' @param A Treatment vector.
#' @param Xs Covariate matrix.
#'
#' @importFrom nnet multinom
#' @importFrom stats predict
#' @importFrom utils capture.output
#'
#' @export
M.propen<-function(A,Xs) {
  if(ncol(as.matrix(A))!=1) stop("Cannot handle multiple stages of treatments together!")
  if(length(A)!= nrow(as.matrix(Xs))) stop("A and Xs do not match in dimension!")
  if(length(unique(A))<=1) stop("Treament options are insufficient!")
  class.A<-sort(unique(A))
  s.data <- data.frame(A,Xs)
  model<-capture.output(mlogit<-multinom(A ~., data=s.data))
  s.p  <- predict(mlogit,s.data,"probs")
  if (length(class.A)==2) {
    s.p<-cbind(1-s.p,s.p)
  }
  colnames(s.p)<-paste("pi=",class.A,sep="")
  return(s.p) # matrix of N*length(unique(A))
}

# I have tried all relevant modelling packages I can find. 
# The packages multinom, mlogit, and mclogit all can't handle random effects. 
# The packages glmmTMB and glmmADMB both can't handle multinomial models. 
# Nor can glmer, lme4, or any other similar well known packages. 
# *edit; mclogit did not work for me as the data reorganization was 
# too complex for my dataset. In reality there are 100 columns, 
# so expanding it was too taxing for my computer

# The only traditional approach that seems to meet my 
# requirements is npmlt from the mixcat package, but I 
# get the same error each time which seems to relate to my random effect(s). 
# I have tried removing all NA values, only using one random effect, 
# combining the random effects into one, and nothing works. 
# In the below code, "comb" is the column where I combined the random effect during testing. 
# I read that the data needs attaching before running for this package.



# if use gam in mgcv
# mgcv::gam(list(EEC_multinomial ~ call + duration + s(dyad, bs="re"),
          #      ~ call + duration + s(dyad, bs="re"),
          #      ~ call + duration + s(dyad, bs="re")),
          # data = dur, family = multinom(K=3))

# The s(dyad, bs="re") is the random effect term for dyad, 
# and the K value specified in multinom is the number of levels of the 
# categorical response minus 1. Because there are four levels to this response, 
# K=3 and we must repeat the formula three times, within a list. The summary() 
# output of the model will show three sections - the summary for each level of 
# the response compared to the reference level. This also allows the flexibility 
# to include different predictors for different levels of the response.


# if use npmlt in mixcat
# library(mixcat)
# model.po <- npmlt(formula = EEC_multinomial ~ call + duration, 
# formula.npo = ~ 1, 
# random = ~ 1|dyad, 
# k = 2)







