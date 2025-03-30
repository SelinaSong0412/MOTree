Benefit_gain = function(w_vec,mu.hat.selected,left){
  
  m = nrow(w_vec)
  p = ncol(w_vec)
  gain_list = rep(0,m)
  
  for(ii in 1:m){
    U = 0
    for(jj in 1:p){
      U = U + mu.hat.selected[[jj]] * w_vec[ii,jj] 
    }
    
    chosen1 = left
    if(all(chosen1==1)){
      next
    }
    if(sum(chosen1)>1){
      res1 = apply(U[chosen1==1,], 2, sum)
    }else{
      res1 = U[chosen1==1,]
    }
    
    if(sum(1-chosen1)>1){
      res2 = apply(U[chosen1==0,], 2, sum)
    }else{
      res2 = U[chosen1==0,]
    }
    (res3 = apply(U, 2, sum))
    gain_list[ii] =  (max(res1) + max(res2) - max(res3)) / nrow(U)
  }
  return(gain_list)
}

