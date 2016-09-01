#library("trajcomp")

plot_sorted <-function(a,b, ylabel, leg, tit=""){
  a_sort = sort(a)
  b_sort = sort(b)
  xx <-c(a_sort, b_sort)
  xx[is.infinite(xx)] <- 0  
  ylimt <- c(0,max(xx))
  
  a_sort[is.infinite(a_sort)] <- max(xx) 
  b_sort[is.infinite(b_sort)] <- max(xx)
  print (ylimt)
  plot(a_sort, 
        ylim= ylimt,
       ylab=ylabel,
       main =ylabel,
       xlab=tit
  )
  
  legend(5,ylimt[2], leg, col = c("black", "red"), lty = c(1, 1), pch = c(1,2), bty="n")
  lines(a_sort)
  points(b_sort, col="red", pch=2)
  lines(b_sort, col="red") 
}


plot_sorted1 <-function(a, main, ylabel, tit="", extra=""){
  a_sort = sort(a)
  ylimt <- c(0,max(a_sort))
  plot(a_sort,
       main=main, 
       ylab=ylabel,
       xlab=tit)
  lines(a_sort)
  abline(h =1, untf = FALSE, lty=2, col="red")
  
  legend(5,ylimt[2], c(ylabel, "Reference Line"), col = c("black", "red"), lty = c(1, 2), pch = c(1,-1), bty="n",
         title=extra, title.col="red")
}

analyze_dataset <-function(trajectories, title, ignore, param)
{
  epsilon <-param[1]
  beta <-param[2]
  levels <-param[3]
  handle = compileSettings(
    tgroup(
      tdistance("discrete_hausdorff"),
      tdistance("dtw") #,  tdistance("discrete_frechet")
    ));
  
  ps_hau <- c()
  dp_hau <- c()
  
  ps_dtw <- c()
  dp_dtw <- c()
  
  ps_fre <- c()
  dp_fre <- c()
  
  traj_ind <- c()
  
  for(i in 1:length(trajectories)){
    
    #print(paste("ITER", i))
    traj <- trajectories[[i]]
    
    #dp<-DouglasPeucker(as.matrix(traj), epsilon)
    dp<- persistence_multires(as.matrix(traj),beta, levels)
    ps<- persistence_dist(as.matrix(traj),beta, 4, 2)
    
    if(length(ps) >= 2 && length(dp) >= 2 && !(i %in% ignore)){
      dp_dist <-  distance_vector_ddply(as.matrix(dp),as.matrix(traj),handle)
      ps_dist <-  distance_vector_ddply(as.matrix(ps),as.matrix(traj),handle)
      
      # print(paste("dp:", dp_dist)) 
      # print(paste("ps:", ps_dist))
      
      ps_hau <- c(ps_hau,ps_dist[1] )
      ps_dtw <- c(ps_dtw,ps_dist[2]/length(traj) )
      #ps_fre <- c(ps_fre,ps_dist[3] )
      
      dp_hau <- c(dp_hau,dp_dist[1] )
      dp_dtw <- c(dp_dtw,dp_dist[2]/length(traj) )
      #dp_fre <- c(dp_fre,dp_dist[3] )
      
      traj_ind <-c(traj_ind, i)
    }
  }

  dp_name <- paste("DouglasPeucker(", epsilon,")")
  ps_name <- paste("Persistence(", beta, ",", levels,")")
  pl2_name <-c(ps_name, dp_name)
  plot_sorted(ps_dtw, dp_dtw, "Dynamic Time Warping Distance", pl2_name, title)
  plot_sorted(ps_hau, dp_hau, "Discrete Hausdorff Distance", pl2_name, title)
  #plot_sorted(ps_fre, dp_fre, "Discrete Frechet Distance", pl2_name, title)

  prune_threshold <- 6
  dtw_pruned <- dp_dtw/ps_dtw
  prune_num <- length(which(dtw_pruned >prune_threshold))
  dtw_pruned <- dtw_pruned[-which(dtw_pruned >prune_threshold)]
  
  pl1_name <- paste(dp_name, "/", ps_name)
  
  #plot_sorted1(dp_dtw/ps_dtw, "Dynamic Time Warping Distance", pl1_name, title)
  # plot_sorted1(dtw_pruned, "Dynamic Time Warping Distance (Edit)", pl1_name, title,
  #              paste(prune_num,"Values >", prune_threshold,"have been removed"))
  #plot_sorted1(dp_hau/ps_hau, "Discrete Hausdorff Distance", pl1_name, title)
  #plot_sorted1(dp_fre/ps_fre, "Discrete Frechet Distance", pl1_name, title)
}

data(prague)
prague$tid = getTrajectoryIDs(prague)
prague_split = split (prague, f = prague$tid);
analyze_dataset(prague_split, "prague dataset", c(127), c(32,5,1))



data(sanfrancisco)
sanfrancisco$tid = getTrajectoryIDs(sanfrancisco)
sanfrancisco_split = split (sanfrancisco, f = sanfrancisco$tid);
analyze_dataset(sanfrancisco_split, "sanfrancisco dataset",c(), c(2,5,4))

