library(plyr)
A=read.table("clustering.csv",head=TRUE)
B=read.table("prague.XYI");
# create an empty plot of correct axis
plot(NaN,NaN, xlim=c(min(B$V1),max(B$V1)), ylim=c(min(B$V2),max(B$V2)),main="Result", xlab="X", ylab="Y");
# group on trajectory ID and add lines for the dataset in gray

ddply(B, .(V3), function(x) { lines(x$V1,x$V2, col="gray");});

# not actually elegant, but very clear for beginners:
cluster_to_color <- function(r)
{
   if (r == -2){
      return ("red");
   }
   if (r == -1){
       return("black");
   }
   if (r == 0){
       return ("yellow");
   }
   if(r == 1){
      return ("green");
   }
}
# Gruppiere nun die Line Segments (A) nach clustern und male Segmente

ddply(A,.(cluster),function(x){
   #trajectories = unique(sort(x$trajectory_index));
   points(c(x$x1,x$x2),c(x$y1,x$y2),col=cluster_to_color(x$cluster[1]));
})


