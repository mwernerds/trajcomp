#' Create a vector with an integer trajectory ID from a NaN separated 
#' dataset
#'
#' @param x A trajectory database separated by NaN lines
#' @return A vector containing trajectory indizes and NaN 
#' @family basic trajectory dataset tools
#' @seealso \code{\link{removeNaNs}} to create a dataset with an ID 
#' column and without NaNs (ideal for \code{\link{ddply}})
#' @examples
#' prague$tid = getTrajectoryIDs(prague);
getTrajectoryIDs <- function (Q)
{
    ret = matrix(nrow = nrow(Q));
    startPos = 1;
    k=1;
	for (i in 1:nrow(Q))
	{
	    if (is.nan(Q[i,1])){
	        ret[startPos:(i-1)] = k;
	        ret[i] = NaN;
	        k = k+1;
	        startPos = i+1;
		}
	}
	ret[startPos:nrow(Q)] = k;
    return (ret);
   
}

# DEPRECATED 
# add NaNs between the trajectories
addNaN <- function(S){
  S <- as.matrix(S);
  res <- matrix(S[1,1:3], ncol = 3)
  # colnames(res) = c("tid", "x", "y");
  for(i in 2:(length(S[,1])-1)){
  if(S[i,1] == S[i+1,1]){

          res <- rbind(res,S[i,1:3]);
    }
    else{
      res <- rbind(res,S[i,1:3]);
      res <- rbind(res, c(NaN,NaN,NaN))
    }
  }
  res <- rbind(res,S[length(S[,1]),1:3])
  return(res)
}

#' Remove separating NaN lines from a dataset
#'
#' @param x A trajectory database separated by NaN lines
#' @return The same object without NaN lines
#' @family basic trajectory dataset tools
#' @examples
#' removeNaNs(prague);
removeNaNs <- function(S){
  res = S[!is.na(S$x),] 
  return(res);
}


#' Split trajectory dataset into a list
#'
#' @param x A trajectory database separated by NaN lines
#' @return A list of trajectories
#' @family basic trajectory dataset tools
#' @seealso \code{\link{getTrajectoryIDs}} and \code{\link{ddply}} 
#' @examples
#' psplit = tsplit(prague);

tsplit<-function(x) 
{
    indizes = getTrajectoryIDs(x);
    l = split(x,f=indizes);
    l$'NaN'= NULL;
    return(l);
}


