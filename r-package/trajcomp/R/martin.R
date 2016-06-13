#' Calculates a distance vector for a distance given by name
#'
#' @param dataset A matrix (or data frame) representing the dataset with NaNs between different trajectories
#' @param query A matrix (or data frame) representing a single query trajectory
#' @param distance_name A string defining the distance by name
#' @return A vector containing the distance of each trajectory as a column in the database (see the example)
#' @examples distance_vector(prague,tsplit(prague)[[1]],"closest_pair");
#' @family trajcomp.distances
distance_vector <- function (dataset, query, distance_name)
{
   h = compileSettings(tgroup(tdistance(distance_name)))
   return(trajectory_distance_vector(as.matrix(dataset),as.matrix(query),h));
}


# Plotting helper
cscale <- function(Q)
{
   carea = quantile(Q,c(0,1), na.rm = TRUE)
   Q = Q - carea[1];
   Q = Q / carea[2];
   Q[Q<0] = 0;
   Q[Q>1] = 1;
   return (Q*1023);
}


plot_distance_colormap <- function (Q,ref,fun,main="",...,g_xlim = NULL, g_ylim = NULL, show.query=TRUE, show.dataset = TRUE)
{

	color_table = 1024 - cscale(fun(...));

	plot(Q$x, Q$y,xlab="X", ylab="Y", pch=23,main = main, xlim=g_xlim, ylim=g_ylim,bg=color_table,col="#00000000");
	if(show.dataset){
	lines(Q$x,Q$y)
	}
	if (show.query){
	lines(ref$x,ref$y,col="red",lwd=8)
	}
}

plot_distance_colormap_defaultpalette<-function()
{
	par(bg="#cccccc");
	palette(adjustcolor(heat.colors(1024),alpha=0.2));
}

# Prague Splitter (with clean=TRUE removing NaNs for ddply, but plotting is hard / impossible)

