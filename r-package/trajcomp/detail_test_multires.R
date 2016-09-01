
traj <- psplit[[24]]

dp<- DouglasPeucker(as.matrix(traj), 32)
psd<- persistence_dist(as.matrix(traj),3, 1, 2)
psm<- persistence_multires(as.matrix(traj),5, 7)
psb<- persistence_pruned(as.matrix(traj),3)

plot(traj$x,traj$y, col="yellow")
text(traj$x[1], traj$y[1], "start")
usr <- par( "usr" )
points(dp, col="black", pch=0, cex=0.7)
lines(dp, col="black")
lines(psm, col="green")
points(psm, col="green", pch=2, cex=1.5)
lines(psd, col="red")
points(psd, col="red", pch=1, cex=1.7)
lines(psb, col="aquamarine3")
points(psb, col="aquamarine3", pch=5, cex=1.0)
if(length(psm) <2)
  text(usr[1],usr[3], "psm empty")

text(usr[2],usr[4], as.character(24), adj = c( 1, 1 ))
legend('topleft',
       c("original", "DP", "PS MR", "PS D", "PS B"),
       col = c("yellow","black", "green", "red", "aquamarine3"),
       lty = c(-1, 1,1,1, 1),
       pch = c(1, 0, 2, 1, 5),
       bty="n")


for(i in 1:112){
  print(paste("ITTT:", i))
  traj <- psplit[[i]]
  
  dp<- DouglasPeucker(as.matrix(traj), 32)
  psd<- persistence_dist(as.matrix(traj),5, 32, 3)
  psm<- persistence_multires(as.matrix(traj),5, 7)
  psb<- persistence_pruned(as.matrix(traj),5)
  
  if(length(psm)> 2){
      #if(tail(psm, n=1) != tail(dp, n=1)){
    if(psm[1,] == dp[1,]){
      
      plot(traj$x,traj$y, col="yellow")
      text(traj$x[1], traj$y[1], "start")
      usr <- par( "usr" )
      points(dp, col="black", pch=0, cex=0.7)
      lines(dp, col="black")
      lines(psm, col="green")
      points(psm, col="green", pch=2, cex=1.5)
      lines(psd, col="red")
      points(psd, col="red", pch=1, cex=1.7)
      lines(psb, col="aquamarine3")
      points(psb, col="aquamarine3", pch=5, cex=1.0)
      if(length(psm) <2)
        text(usr[1],usr[3], "psm empty")
      
      text(usr[2],usr[4], as.character(i), adj = c( 1, 1 ))
      legend('topleft',
             c("original", "DP", "PS MR", "PS D", "PS B"),
             col = c("yellow","black", "green", "red", "aquamarine3"),
             lty = c(-1, 1,1,1, 1),
             pch = c(1, 0, 2, 1, 5),
             bty="n")
    }
  }
}