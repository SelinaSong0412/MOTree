
#,"#ae3f51"

my_vis = function(mu1.hat,mu2.hat,chosen,k,label,lwd=3){
  w0 = seq(0,1,0.02)
  res = matrix(0,length(w0),k)
  for(ii in 1:length(w0)){
    res[ii,] =  colMeans(mu1.hat[chosen,] *w0[ii] + mu2.hat[chosen,] * (1-w0[ii]  ))
  }
  for(ii in 1:nrow(res)){res[ii,] =  res[ii,] - max(res[ii,])}
  df_vis = data.frame(value = as.vector(res),
                      Treatment = rep(label,each = nrow(res)),
                      weight = rep(w0,4))
  p1 = ggplot(data = df_vis,aes(x = weight,y = value,group = Treatment,col=Treatment))+
    geom_line(lwd=lwd)+ theme_bw()+ theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),aspect.ratio = 0.5
    ) +scale_fill_manual(values = Custom.color)+
    scale_color_manual(values = Custom.color)
  
  p1
}



my_vis_pair = function(mu1.hat,mu2.hat,chosen1,chosen2,k,label,lwd=3){
  w0 = seq(0,1,0.02)
  res1 = matrix(0,length(w0),k)
  for(ii in 1:length(w0)){
    res1[ii,] =  colSums(mu1.hat[chosen1,] *w0[ii] + mu2.hat[chosen1,] * (1-w0[ii]  ))
  }
  for(ii in 1:nrow(res1)){res1[ii,] =  res1[ii,] - max(res1[ii,])}

  
  w0 = seq(0,1,0.02)
  res2 = matrix(0,length(w0),k)
  for(ii in 1:length(w0)){
    res2[ii,] =  colSums(mu1.hat[chosen2,] *w0[ii] + mu2.hat[chosen2,] * (1-w0[ii]  ))
  }
  for(ii in 1:nrow(res2)){res2[ii,] =  res2[ii,] - max(res2[ii,])}
  
  result = rep(0,length(w0))
  for(ii in 1:length(w0)){
    result[ii] = max(res1[ii,]) + max(res2[ii,]) - max(res2[ii,] + res1[ii,])
  }
  result = result / (sum(chosen1) + sum(chosen2))
  
  df_vis = data.frame(value = as.vector(result),
                      weight = rep(w0,4))
  p1 = ggplot(data = df_vis,aes(x = weight,y = value))+
    geom_line(lwd=lwd)+ theme_bw()+theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),aspect.ratio = 0.5
    )  +scale_fill_manual(values = Custom.color)+   #映射云雨图和箱线图的颜色
    scale_color_manual(values = Custom.color)
  
  p1
}




library(ggplot2)
library(ggforce)

# my_vis_point = function(chosen){
#   
#   tmp1 = mus1.reg[chosen,] 
#   tmp2 = mus2.reg[chosen,] 
#   
#   for(ii in 1:dim(tmp1)[1]){
#     tmp1[ii,] = tmp1[ii,] - tmp1[ii,4] + rnorm(1,0,0.001)
#     tmp2[ii,] = tmp2[ii,] - tmp2[ii,4] + rnorm(1,0,0.001)
#   }
#   
#   df_vis = data.frame(x = as.vector(tmp1),
#                       y = as.vector(tmp2),
#                       Treatment = rep(c("A","B","C","D"),each = nrow(tmp1)))
#   p1 = ggplot(data = df_vis,aes(x = x,y = y,group = Treatment,col=Treatment))+geom_point(alpha = 0.4) +
#     stat_ellipse(geom = "polygon",
#                  level = 1/nrow(tmp1),
#                  aes(fill = Treatment), 
#                  alpha = 0.9)
#   p1
# }


my_vis_elipse = function(mu1.hat,mu2.hat,chosen,k,ref_k,label,lwd=3,
                         textsize){
  tmp1 = mu1.hat[chosen,] 
  tmp2 = mu2.hat[chosen,] 
  
  tmp1 = tmp1 - tmp1[,ref_k]
  tmp2 = tmp2 - tmp2[,ref_k]
  

  AIPW_v1 = apply(tmp1, 2, var)
  AIPW_v2 = apply(tmp2, 2, var)
  AIPW_cov = diag(cov(tmp1,tmp2))
  AIPW_m1 = colMeans(tmp1)
  AIPW_m2 = colMeans(tmp2)
  
  
  generate_elipse = function(i,n){
    mu = c(AIPW_m1[i],AIPW_m2[i])
    sigma  = matrix(c(AIPW_v1[i],AIPW_cov[i],AIPW_cov[i],AIPW_v2[i]),2,2) 
    alpha = 0.05
    es <- eigen(sigma)
    e1 <- es$vec%*%diag( sqrt(abs(es$val))*sign(es$val) )
    r1 <- sqrt(qchisq(1-alpha,2))
    npoints = 100
    theta <- seq(0,2*pi,len=npoints)
    v1 <- cbind(r1*cos(theta),r1*sin(theta))
    pts=t(mu-(e1%*%t(v1)))
    return(pts)
  }

  
  df_elipse = vector()
  for(kk in 1:(k)){
      df_elipse =  rbind(df_elipse,
                         generate_elipse(kk,nrow(tmp1)))
  }
  df_elipse = as.data.frame(df_elipse)
  
  
  
  colnames(df_elipse) = c("x","y")
  # df_elipse$Treatment = rep(c("A","B","C"),each = nrow(df_elipse)/3)
  df_elipse$Treatment = rep(label,each = nrow(df_elipse)/length(label))
  
  # tmp1 = mus1.reg[chosen,] 
  # tmp2 = mus2.reg[chosen,] 
  tmp1 = mu1.hat[chosen,] 
  tmp2 = mu2.hat[chosen,] 
  tmp1 = tmp1 - tmp1[,ref_k]
  tmp2 = tmp2 - tmp2[,ref_k]
  
  df_vis = data.frame(x = as.vector(tmp1),
                      y = as.vector(tmp2),
                      Treatment = rep(label,each = nrow(tmp1)))
  
  p = ggplot()+geom_path(data = df_elipse,aes(x = x, y = y,group=Treatment,col = Treatment))+ 
    geom_point(data = df_vis,aes(x = x,y = y,group = Treatment,col=Treatment),alpha = 0.5)+ theme_bw() + coord_fixed()+
    theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),aspect.ratio = 0.5
    )  + xlab("Outcome 1")+ ylab("Outcome 2")+
    scale_fill_manual(values = Custom.color)+   
    scale_color_manual(values = Custom.color)
  
  p
}




my_vis_mean = function(mu1.hat,mu2.hat,chosen,k,label,lwd=3){
  
  tmp1 = mu1.hat[chosen,] 
  tmp2 = mu2.hat[chosen,] 
  
  tmp1 = tmp1 - tmp1[,k]
  tmp2 = tmp2 - tmp2[,k]
  
  
  
  df_vis = data.frame(x = colMeans(tmp1),
                      y = colMeans(tmp2),
                      Treatment = label)
  
  p = ggplot()+
    geom_point(data = df_vis,aes(x = x,y = y,group = Treatment,col=Treatment),alpha = 0.5) + theme_bw()+ 
    coord_fixed()+
    theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),aspect.ratio = 0.5
    )   + xlab("Outcome 1")+ ylab("Outcome 2")+
    scale_fill_manual(values = Custom.color)+   
    scale_color_manual(values = Custom.color)
  
  p
}



my_vis_mean_pair = function(mu1.hat,mu2.hat,chosen1,chosen2,i,j,label,k,point_size=3){

  tmp1 = mu1.hat[chosen1,] 
  tmp2 = mu2.hat[chosen1,] 
  tmp1 = tmp1 - tmp1[,k]
  tmp2 = tmp2 - tmp2[,k]
  df_vis1 = data.frame(x = colMeans(tmp1),
                      y = colMeans(tmp2),
                      Treatment = label)
  df_vis1$Group = label[i]
  
  tmp1 = mu1.hat[chosen2,] 
  tmp2 = mu2.hat[chosen2,] 
  tmp1 = tmp1 - tmp1[,k]
  tmp2 = tmp2 - tmp2[,k]
  df_vis2 = data.frame(x = colMeans(tmp1),
                      y = colMeans(tmp2),
                      Treatment = label)
  df_vis2$Group = label[j]
  
  
  df_vis = rbind(df_vis1,df_vis2)
  df_vis$Group = as.factor(df_vis$Group)
  
  
  hull <- df_vis %>% group_by(Group) %>% 
    slice(chull(x,y))
  
  # geom_polygon(aes(x = x,y = y,fill = Group, group = Group),alpha = 0.1)+
  
  p = ggplot(data = df_vis)+
    #geom_point(aes(x = x,y = y,group = Treatment,col=Treatment,shape = Group),alpha = 0.5,size=point_size) + 
    geom_text(aes(x = x,y = y,col=Group,label = Treatment))+
    geom_polygon(data = hull,aes(x=x,y=y,pattern = Group),alpha = 0.2)+
    coord_fixed()+ theme_bw()+
    theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      aspect.ratio = 0.5) + xlab("Outcome 1")+ ylab("Outcome 2")+
    scale_fill_manual(values = setNames(Custom.color,label))+   
    scale_color_manual(values = setNames(Custom.color,label))
  
  p
}





my_vis_mean_pair = function(mu1.hat,mu2.hat,chosen1,chosen2,i,j,label,k,
                            point_size=3,
                            my_magick = c("bricks","horizontalsaw","verticalsaw",
                                          "hs_horizontal","hs_vertical")){
  
  tmp1 = mu1.hat[chosen1,] 
  tmp2 = mu2.hat[chosen1,] 
  tmp1 = tmp1 - tmp1[,k]
  tmp2 = tmp2 - tmp2[,k]
  df_vis1 = data.frame(x = colMeans(tmp1),
                       y = colMeans(tmp2),
                       Treatment = label)
  df_vis1$Group = paste("Group",i)
  
  tmp1 = mu1.hat[chosen2,] 
  tmp2 = mu2.hat[chosen2,] 
  tmp1 = tmp1 - tmp1[,k]
  tmp2 = tmp2 - tmp2[,k]
  df_vis2 = data.frame(x = colMeans(tmp1),
                       y = colMeans(tmp2),
                       Treatment = label)
  df_vis2$Group = paste("Group",j)
  
  
  df_vis = rbind(df_vis1,df_vis2)
  df_vis$Group = as.factor(df_vis$Group)
  
  
  hull <- df_vis %>% group_by(Group) %>% 
    slice(chull(x,y))
  
  # geom_polygon(aes(x = x,y = y,fill = Group, group = Group),alpha = 0.1)+
  
  p = ggplot(data = df_vis)+
    #geom_point(aes(x = x,y = y,group = Treatment,col=Treatment,shape = Group),alpha = 0.5,size=point_size) + 
    geom_polygon_pattern(data = hull,
      aes(x=x,y=y, pattern_type = Group), # , pattern_fill = Group
      alpha = 0.4,
      pattern_alpha = 0.4,
      pattern_fill = "black",
      pattern                  = 'magick',
      pattern_scale            = 1,
      stat                     = "identity",
      fill                     = "#79aec1",
      pattern_aspect_ratio     = 1,
      pattern_density          = 0.1,
      show.legend = TRUE
    )+scale_pattern_type_discrete(choices = c(my_magick[i],my_magick[j])) +
    geom_text(data = df_vis,aes(x = x,y = y,label = Treatment,col=Treatment))+ theme_bw()+
    coord_fixed()+
    theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),aspect.ratio = 0.5
      )  +
    xlab("Outcome 1")+ ylab("Outcome 2")+
    scale_fill_manual(values = setNames(Custom.color,label))+   
    scale_color_manual(values = setNames(Custom.color,label))
  
  p
}






#(p = my_vis_mean_pair(mu1.hat,mu2.hat,chosen1,chosen2,i,j,label1,k,point_size=5))





my_vis_all_pair = function(mu1.hat,mu2.hat,leaf_pred,label,k,alpha = 1,
                           my_magick = c("bricks","horizontalsaw","verticalsaw",
                                         "hs_horizontal","hs_vertical")){
  
  leaf_vec = sort(unique(leaf_pred))
  m0 = length(leaf_vec)
  df_vis = data.frame()
  for(ii in 1:m0){
    chosen1 = which(leaf_pred==leaf_vec[ii])
    tmp1 = mu1.hat[chosen1,] 
    tmp2 = mu2.hat[chosen1,] 
    tmp1 = tmp1 - tmp1[,k]
    tmp2 = tmp2 - tmp2[,k]
    df_vis1 = data.frame(x = colMeans(tmp1),
                         y = colMeans(tmp2),
                         Treatment = label)
    df_vis1$Group = paste("Group",ii)
    df_vis = rbind(df_vis,df_vis1)
  }
  
  df_vis$Group = as.factor(df_vis$Group)
  
  
  hull <- df_vis %>% group_by(Group) %>% 
    slice(chull(x,y))
  
  # geom_polygon(aes(x = x,y = y,fill = Group, group = Group),alpha = 0.1)+
  
  p = ggplot(data = df_vis)+
    geom_line(aes(x = x,y = y,group = Treatment,col=Treatment)) + 
    geom_polygon_pattern(data = hull,
                         aes(x=x,y=y, pattern_type = Group), # , pattern_fill = Group
                         alpha = 0.4,
                         pattern_alpha = 0.4,
                         pattern_fill = "black",
                         pattern                  = 'magick',
                         pattern_scale            = 1,
                         stat                     = "identity",
                         fill                     = "#79aec1",
                         pattern_aspect_ratio     = 1,
                         pattern_density          = 0.1,
                         show.legend = TRUE
    )+scale_pattern_type_discrete(choices = c(my_magick[i],my_magick[j]))+
    coord_fixed()+
    theme(
      plot.title = element_text( size=textsize, face="bold.italic"),
      axis.title.x = element_text( size=textsize, face="bold"),
      axis.title.y = element_text(size=textsize, face="bold"),
      aspect.ratio = 0.5) + xlab("Outcome 1")+ ylab("Outcome 2")+
    scale_fill_manual(values = setNames(Custom.color,label))+   
    scale_color_manual(values = setNames(Custom.color,label))+ 
    theme(legend.position = "bottom")
  
  p
}


# 
# 
# 
# my_vis_all_pair = function(mu1.hat,mu2.hat,leaf_pred,label,k,alpha = 1,
#                            my_magick = c("bricks","horizontalsaw","verticalsaw",
#                                          "hs_horizontal","hs_vertical")){
#   
#   leaf_vec = sort(unique(leaf_pred))
#   m0 = length(leaf_vec)
#   df_vis = data.frame()
#   for(ii in 1:m0){
#     chosen1 = which(leaf_pred==leaf_vec[ii])
#     tmp1 = mu1.hat[chosen1,] 
#     tmp2 = mu2.hat[chosen1,] 
#     tmp1 = tmp1 - tmp1[,k]
#     tmp2 = tmp2 - tmp2[,k]
#     df_vis1 = data.frame(x = colMeans(tmp1),
#                          y = colMeans(tmp2),
#                          Treatment = label)
#     df_vis1$Group = paste("Group",ii)
#     df_vis = rbind(df_vis,df_vis1)
#   }
#   
#   df_vis$Group = as.factor(df_vis$Group)
#   
#   
#   hull <- df_vis %>% group_by(Group) %>% 
#     slice(chull(x,y))
#   
#   # geom_polygon(aes(x = x,y = y,fill = Group, group = Group),alpha = 0.1)+
#   
#   p = ggplot()+
#     geom_polygon_pattern(data = hull,
#                          aes(x=x,y=y, pattern_type = Group, pattern_fill = Group), # , pattern_fill = Group
#                          alpha = 0.4,
#                          pattern_alpha = 0.4,
#                          pattern                  = 'magick',
#                          pattern_scale            = 1,
#                          stat                     = "identity",
#                          fill                     = "#79aec1",
#                          pattern_aspect_ratio     = 1,
#                          pattern_density          = 0.1,
#                          pattern_fill = "black"
#     )+scale_pattern_type_discrete(choices = setNames(my_magick,unique(df_vis$Group)) )+
#     theme(legend.position = "bottom",legend.key.size = unit(2, 'cm'))
#   
#   p
# }
# 
# (p2 = my_vis_all_pair(mu1.hat,mu2.hat,leaf_pred,label1,k,alpha = 1))
# 
# 
# 
# df1 <- data.frame(trt = c("a", "b", "c"), outcome = c(2.3, 1.9, 3.2))
# 
# ggplot(df1, aes(trt, outcome)) +
#   geom_col_pattern(
#     aes(
#       pattern_type = trt, 
#       pattern_fill = trt
#     ),
#     pattern       = 'magick',
#     pattern_key_scale_factor = 0.7,
#     fill          = 'white',
#     colour        = 'black'
#   ) +
#   theme_bw(15) +
#   labs(
#     title    = "ggpattern::geom_col_pattern()",
#     subtitle = "pattern='magick'"
#   ) +
#   theme(legend.key.size = unit(2, 'cm')) +
#   scale_pattern_type_discrete(choices = c('bricks', 'fishscales', 'right45')) +
#   coord_fixed(ratio = 1/2)
