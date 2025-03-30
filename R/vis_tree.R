# install.packages("DiagrammeR")
library(data.tree)
library(dplyr)

ten2bin<-function(x,max_h){
  tmp<- x %%   2^((max_h) :1)
  tmp<-tmp %/% 2^((max_h):1 - 1)
  return(  tmp  )
}

bin2ten<-function(vec){
  n<-length(vec)
  sum( vec * 2^(n:1 - 1)  )
}

vis_tree = function(max_h,treeout,number){
  dat<-matrix(0,2^max_h,max_h)
  treeout$cutoff = round(treeout$cutoff, 3)
  for (k in 1L:(max_h)) { # depth is the most number of split to reach one terminal node
    # The 1st split
    if (k==1L) {
      tmp = rep( c(paste("X",treeout$X[k],"<=",treeout$cutoff[k],sep=""),
                   paste("X",treeout$X[k],">",treeout$cutoff[k],sep="")),
                 each = 2^(max_h-k))
      dat[,1] = tmp
    } else {
      res = vector()
      for(j in (2^(k-1)):(2^k-1)) {
        if(!is.na(treeout$X[j])){
          tmp = rep( c(paste("X",treeout$X[j],"<=",treeout$cutoff[j],sep=""),
                       paste("X",treeout$X[j],">",treeout$cutoff[j],sep="")),
                     each = 2^(max_h-k))
        }else{
          tmp = rep( c(NA,NA),
                     each = 2^(max_h-k))
        }
        res = c(res,tmp)
      }
      dat[,k] = res
    }
  }
  dat = unique(dat)
  pathString <- paste("All Data", dat[,1],sep = "/")
  for(nn in 1:nrow(dat)){
    for(ii in 2:max_h){
      if(!is.na(dat[nn,ii]))
        pathString[nn] <- paste(pathString[nn], dat[nn,ii],sep = "/")
    }  
  }
  
  pathString <- paste(pathString, number,sep = "/")
  dat = as.data.frame(dat)
  dat$pathString <- pathString
  tree0 <- as.Node(dat)  
  return(tree0)  
}



# library(igraph)
# plot(as.igraph(tree0, directed = TRUE, direction = "climb"))
# 
# library(networkD3)
# useRtreeList <- ToListExplicit(tree0, unname = TRUE)
# radialNetwork( useRtreeList)





# 
# initial_tree<-function(max_h,treeout){
#   
# 
#   dat<-matrix(0,2^max_h,max_h+5)
#   for(ii in 1:2^max_h){
#     dat[ii,1:max_h] <- ten2bin(ii-1,max_h)+1
#   }
#   pathString <- paste("dat", dat[,1],sep = "/")
#   for(ii in 2:max_h){
#     dat[,ii]<-paste(dat[,ii-1],dat[,ii],sep = "")
#     pathString <- paste(pathString, dat[,ii],sep = "/")
#   }
#   dat<-as.data.frame(dat)
#   colnames(dat)<- c(1:max_h,"judge","sign","threshold","h","n0")
#   dat$pathString <- pathString
#   
#   
# 
#   h0<- 0
#   sign0<-0
#   judge0<- 0
#   threshold0<-0
#   value0<-NA
#   
#   cur_h = 1
#   for(ii in 1:max_h){
#     for(kk in 1:(2^ii/2)){
#       tmp<-which(h0== ii-1)
#       jj<-tmp[kk]
#       new_judge<- treeout$X[cur_h]
#       thres0<- treeout$cutoff[cur_h]
#       
#       judge0<- c(judge0[1:jj],new_judge,new_judge,judge0[-(1:jj)] )
#       threshold0<-c(threshold0[1:jj],thres0,thres0,threshold0[-(1:jj)] )
#       sign0<- c(sign0[1:jj],-1,1,sign0[-(1:jj)] )
#       h0<- c(h0[1:jj],ii,ii,h0[-(1:jj)] )
#       cur_h = cur_h + 1
#     }
#   }
# 
#   tree <- as.Node(dat)  
#   tree$Set(sign=sign0,h=h0,threshold=threshold0,judge=judge0,n0=0)
#   return(tree)
# }
# 
# tree0 = initial_tree(max_h = 3,treeout)
# print(tree0,"judge","threshold","sign","h","n0")
# 
# 
# library(igraph)
# plot(as.igraph(tree0, directed = TRUE, direction = "climb"))
# 
# plot(tree0)
