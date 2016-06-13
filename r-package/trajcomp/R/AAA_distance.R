
#' Calculates the cross product sparse matrix from two incidence matrices
#'
#' @param Matrix A: A one-sparse matrix (class matrix, not dgCMatrix, so not sparsely stored) with rownames 
#' @param Matrix B: A second one-sparse matrix (class matrix, not dgCMatrix, so not sparsely stored) with rownames 
#' @return A matrix of nrow(A) * nrow(B) rows, with rownames concatenated and combinations set to one (e.g., the cross product of the alphabets)
#' @examples matrixcross(A,B)
#' @family SVDdistance


matrixcross<-function(A,B)
{
    if (ncol(A) != ncol(B))
    {
        stop("Error: matrixcross only defined for matrices with identical numbers of columns")
    }
    C = matrix(0,nrow = nrow(A)*nrow(B), ncol=ncol(A));
    # Prepare rownames
    cn = vector("character",length=nrow(A)*nrow(B));
    for (i in 1:nrow(A))
    {
       for(j in 1:nrow(B))
       {
			cn[(i-1)*nrow(B)+j] = paste(rownames(A)[i],rownames(B)[j],sep=",")
       }
    }
    rownames(C) = cn;
    # Assuming one-sparse columns
    for (column in 1:ncol(A))
    {
        iA = which((A==1)[,column]);
        iB = which((B==1)[,column]);
        i = (iA-1)*nrow(B) + iB;
        C[i,column] = 1;
    }
    return (C);
}

#' Calculates the cross product sparse matrix from two incidence matrices (faster than the version with rownames)
#'
#' @param Matrix A: A one-sparse matrix (class matrix, not dgCMatrix, so not sparsely stored) with or without rownames 
#' @param Matrix B: A second one-sparse matrix (class matrix, not dgCMatrix, so not sparsely stored) with or without rownames 
#' @return A matrix of nrow(A) * nrow(B) rows, without rownames concatenated and combinations set to one (e.g., the cross product of the alphabets)
#' @examples matrixcross_nonames(A,B)
#' @family SVDdistance
matrixcross_nonames<-function(A,B)
{
    if (ncol(A) != ncol(B))
    {
        stop("Error: matrixcross only defined for matrices with identical numbers of columns")
    }
    C = matrix(0,nrow = nrow(A)*nrow(B), ncol=ncol(A));
    # Assuming one-sparse columns
    for (column in 1:ncol(A))
    {
        iA = which((A==1)[,column]);
        iB = which((B==1)[,column]);
        i = (iA-1)*nrow(B) + iB;
        C[i,column] = 1;
    }
    return (C);
}




# AAA distance between two trajectories, already encoded, i.e. vector of strings (same length and number of features)
# deprecated
AAA_distance <- function(t1, t2)
{
  library(corpcor)
  # number of elements
  n <- length(t1);
  # pairs of adjacents -> column headers
  nAdjPairs <- 2*n-2;
  pairs <- vector(mode="character", length=nAdjPairs);
  for(i in 1:(n-1)){
    pairs[i] <- paste(t1[i], t1[i+1], sep="");
    pairs[i+(n-1)] <- paste(t2[i], t2[i+1], sep="");
  }
  # unique adjacent pairs -> row names
  uniquePairs = unique(pairs);
  
  sparseMatrix <- data.frame(matrix(0, ncol = nAdjPairs, nrow = length(uniquePairs)));
  colnames(sparseMatrix) <- pairs;
  rownames(sparseMatrix) <- uniquePairs;

  # set one 1 in each column
  for(i in 1:nAdjPairs){
    sparseMatrix[which(colnames(sparseMatrix)[i]==rownames(sparseMatrix)),i] = 1;
  }
  # SVD on submatrices (for each trajectory)
  
  d1 = fast.svd(as.matrix(sparseMatrix[,1:(n-1)]))
  d1_long <- vector(mode="numeric", length = nrow(sparseMatrix));
  d1_long[1:length(d1$d)] <- d1$d;

  d2 = fast.svd(as.matrix(sparseMatrix[,(n:ncol(sparseMatrix))]))
  d2_long <- vector(mode="numeric", length = nrow(sparseMatrix));
  d2_long[1:length(d2$d)] <- d2$d;


  # normalization
  d1_long <- d1_long / sum(abs(d1_long));
  d2_long = d2_long / sum(abs(d2_long));
  
  # euclidean distance
  d =  sqrt(sum((d2_long-d1_long)^2));
 
  return (d);
}

singular_values_fastsvd <- function(M, r)
{
	library(corpcor)
	return ( fast.svd(as.matrix(M))$d[1:r])
}

singular_values_redsvd <- function (M,r)
{
   library(Matrix);	
   S = Matrix(M,sparse=TRUE);
   return (redSVD(S,r)$D);
}


buildMatrix <- function (t, nChar, startChar='A')
{
	library(gtools)
  t <- toString(t); 
  t <- strsplit(t,"");
  # number of elements
  n <- length(t[[1]]);
  # pairs of adjacents -> column headers
  nAdjPairs <- n-1;
  pairs <- vector(mode="character", length=nAdjPairs);
  for(i in 1:nAdjPairs){
    pairs[i] <- paste0(t[[1]][i], t[[1]][i+1]);
  }
 
  sparseMatrix <- data.frame(matrix(0, ncol = nAdjPairs, nrow = nChar*nChar));
  colnames(sparseMatrix) <- pairs;
  
  
  rowNames <- vector(mode="character", length = nChar*nChar);
  count <- 1;
  for(i in 0:(nChar-1)){
    for(j in 0:(nChar-1)){
      rowNames[count] <- paste(chr(asc(startChar)+i), chr(asc(startChar)+j), sep="");
      count <- count+1;
    }
  }
  
  rownames(sparseMatrix) <- rowNames;
  
  # set one 1 in each column
  for(i in 1:nAdjPairs){
    sparseMatrix[which(colnames(sparseMatrix)[i]==rownames(sparseMatrix)),i] = 1;
  }
  return(sparseMatrix);

}






# returns the (normalized) singular values for a given trajectory (already encoded)
# required: number of characters used for encoding and the start character (e.g. "A")
# zwei Buchstaben sind ein Buchstable != nChar

# Trajektorie => Buchstaben => Sparse Matrix => SVD => Feature Vektor

getSVs <- function(t, nChar, startChar='A', singular_value_solver=singular_values_fastsvd, numSingularValues=NA, na.rm = FALSE)
{
  
  library(gtools)
  t <- toString(t); 
  t <- strsplit(t,"");
  # number of elements
  n <- length(t[[1]]);
  # pairs of adjacents -> column headers
  nAdjPairs <- n-1;
  pairs <- vector(mode="character", length=nAdjPairs);
  for(i in 1:nAdjPairs){
    pairs[i] <- paste0(t[[1]][i], t[[1]][i+1]);
  }
 
  sparseMatrix <- data.frame(matrix(0, ncol = nAdjPairs, nrow = nChar*nChar));
  colnames(sparseMatrix) <- pairs;
  
  
  rowNames <- vector(mode="character", length = nChar*nChar);
  count <- 1;
  for(i in 0:(nChar-1)){
    for(j in 0:(nChar-1)){
      rowNames[count] <- paste(chr(asc(startChar)+i), chr(asc(startChar)+j), sep="");
      count <- count+1;
    }
  }
  
  rownames(sparseMatrix) <- rowNames;
  
  # set one 1 in each column
  for(i in 1:nAdjPairs){
    sparseMatrix[which(colnames(sparseMatrix)[i]==rownames(sparseMatrix)),i] = 1;
  }
  
  # SVD on sparse matrix
  if (is.na(numSingularValues))
  {
     numSingularValues = nrow(sparseMatrix);
  }

  d_long = singular_value_solver(sparseMatrix, numSingularValues);
  if (na.rm)
  {
       d_long[is.na(d_long)] = 0;
  }
 
   # normalization
  d_long <- d_long / sum(abs(d_long));
  return (d_long);
}
