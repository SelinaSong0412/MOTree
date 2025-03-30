#安装ggtern包
# install.packages("ggtern")
## 加载ggtern包
library("ggtern")
library(ggplot2)
# tmp = vector()
# for(ii in seq(0,1,0.05)){
#   for(jj in seq(0,1-ii,0.05) ){
#     tmp = rbind(tmp,c(ii,jj,1-ii-jj))
#   }
# }
# tmp
# data = data.frame(x = tmp[,1],
#                   y = tmp[,2],
#                   z = tmp[,3])
# data$Treat = 0
# data$Treat[data$x + data$y > 0.8] = 1
# data$Treat[data$z + data$y > 0.5] = 2
# data$Treat = as.factor(data$Treat)
# p1<-ggtern(data=data[1:2,],aes(x=x,y=y,z=z))+
#   geom_point(aes(color=Treat),alpha=0.8)+ theme_bw()+ #For clarity
#   theme_hidegrid()
# 
# p1  + theme_showgrid_minor() + theme_showarrows()+ 
#   labs(x="X",y="Y",z="Z",title="Title")


pred_label = unique(leaf_pred)
chosen1 = which(leaf_pred == pred_label[4])

data = as.data.frame(w_vec)
colnames(data) = c("x","y","z")
data$A.opt = 0
for(ii in 1:nrow(data)){
  U = data[ii,1] * colMeans(Y_record[[1]][chosen1,]) + 
      data[ii,2] * colMeans(Y_record[[2]][chosen1,]) + 
      data[ii,3] * colMeans(Y_record[[3]][chosen1,])  
  data$Treat[ii] = which.max(U)
}
data$Treat = as.factor(data$Treat)

library(plyr)
df.tot = data
df.mu = as.data.frame(t(colMeans(data[,1:3])))
ggtern(data=df.tot,aes(x,y,z)) + 
  geom_point(aes(color=Treat),alpha=0.8) +
  geom_crosshair_tern(data=df.mu) +
  geom_point(data=df.mu,color='red') + 
  geom_label(data=100*df.mu,
             aes(
               y = y - 5,
               label=sprintf("P=%.1f, S=%.1f, T=%.1f",x,y,z)
             ),
             size=3,
             color='red'
  ) + 
  theme_bw() +
  theme_showarrows() + 
  limit_tern(1,1,1) +
  labs(
    x = 'Y1',xarrow='Y1',
    y = 'Y2',yarrow='Y2',
    z = 'Y3',zarrow='Y3'
  )

