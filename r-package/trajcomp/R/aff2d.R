#
#  A microlibrary for 2D affine transformations in R
#
#

#' Move a coordinate into its projective space 
#'
#' @param Q A vector or matrix
#' @return A vector with a column 1 added
#' @family affine
projective <- function(Q)
{
   return (cbind(Q,1));
}

#' Move a coordinate from projective space into Euclidean space
#'
#' @param Q A vector or matrix
#' @return A vector with the last column removed
#' @family affine
homogene <- function(Q)
{   
   Q[,1] = Q[,1] / Q[,3];
   Q[,2] = Q[,2] / Q[,3];
   Q = Q[,-3]
   return (Q)
}


#' Creates a translation matrix projective space
#'
#' @param s Movement in first direction
#' @param t Movement in second direction
#' @return A 3x3 matrix in projective space
#' @family affine
PM.translation<-function(s,t)
{
   Q = diag(3);
   Q[1:2,3]=c(s,t);
   return(Q);   
}


#' Creates a general affine mapping in projective space
#'
#' @param params A vector containing 6 parameters defining the affine mapping
#' @return A 3x3 matrix in projective space
#' @family affine

PM.affine<-function(params)
{
    if (length(params) != 6)
    {
       stop("An affine map in projective space needs 6 parameters")
    }
    return(rbind(matrix(params, ncol=3, byrow = TRUE),c(0,0,1)));
}



#' Applies an affine map given by params
#'
#' @param Q A vector or matrix
#' @param params A vector containing 6 parameters defining the affine mapping
#' @return A vector or matrix, after applying the transformation in projective space
#' @family affine
map_affine<-function(Q,params)
{
	# Map into projective space
	PD= projective(Q);
	# Generate an affine map
	MAP = PM.affine(params)
	# Map in projective space
	IMAGE = t(PM.affine(params) %*% t(PD))
	# And homogenify
	RESULT = homogene(IMAGE);
	return (RESULT);
}



#' Applies a translation as given by s and t
#'
#' @param Q A vector or matrix
#' @param s Parameter in first direction
#' @param t Parameter in second direction
#' @return A vector or matrix, after applying the transformation in projective space
#' @family affine
map_translate<-function(Q,s,t)
{
   return(map_affine(Q,c(1,0,s,0,1,t)));
}

#' Applies a scaling as given by a and b
#'
#' @param Q A vector or matrix
#' @param a Parameter in first direction
#' @param b Parameter in second direction
#' @return A vector or matrix, after applying the transformation in projective space
#' @family affine
map_scale<-function(Q,a,b)
{
   return(map_affine(Q,c(a,0,0,0,b,0)));
}

#' Applies a skew translation as given by a andbt
#'
#' @param Q A vector or matrix
#' @param a Parameter in first direction
#' @param b Parameter in second direction
#' @return A vector or matrix, after applying the transformation in projective space
#' @family affine

map_skew<-function(Q,a,b)
{
   return(map_affine(Q,c(1,a,0,b,1,0)));
}


#' Shows, how to use the map_affine functinos to randomly project a triangle
#'
#' @family affine
projective.example<-function()
{
	# Create a triangle
	D=matrix(c(0,0,1,0,1,1,0,0),ncol=2, byrow=TRUE);
	# Randomly affine map it around
	E = map_affine(D,runif(6)*2-1);
	# Bind it together (NaN for plot)
	DAT=rbind(D,c(NaN,NaN),E);
	# And plot
	plot(DAT[,1],DAT[,2],t="l",xlim=c(-2,2),ylim=c(-2,2))
}

